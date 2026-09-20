#!/usr/bin/env python3
"""Account for SvABA indels and reciprocal BND records for a public SV benchmark.

Raw calls are immutable. Only a separate comparison VCF is derived. Native GT is
retained; SvABA AD/PL/GQ have caller-specific semantics and are not reinterpreted.
"""
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import re

import pysam

from build_structural_exclusions import fingerprint, save, stamp


def reverse_complement(seq):
    return seq.translate(str.maketrans('ACGTNacgtn', 'TGCANtgcan'))[::-1]


def reference_status(row, reference):
    length=reference.get_reference_length(row.contig)
    if row.start < 0 or row.start+len(row.ref) > length:
        if row.info.get('SVTYPE')=='BND' and row.pos in (0,length+1) and row.ref=='N':
            return 'virtual_boundary_no_reference_base'
        raise ValueError('native allele outside reference bounds')
    if reference.fetch(row.contig,row.start,row.start+len(row.ref)).upper()!=row.ref.upper():
        raise ValueError('native REF mismatch')
    return 'sequence_verified'


def decode_alt(alt, ref):
    m = re.fullmatch(r'([ACGTN]*)([\[\]])([^\[\]]+):(\d+)\2([ACGTN]*)', alt)
    if not m or bool(m[1]) == bool(m[5]) or (m[1] or m[5]) != ref:
        raise ValueError('unexpected SvABA breakend allele')
    return {'chrom': m[3], 'pos': int(m[4]), 'local_side': '3' if m[1] else '5',
            'remote_side': '5' if m[2] == '[' else '3'}


def paired_events(records, order):
    by_id = {row['id']: row for row in records}
    if len(by_id) != len(records): raise ValueError('duplicate breakend ID')
    seen = set(); result = []
    for row in records:
        if row['id'] in seen: continue
        other = by_id.get(row['info'].get('MATEID'))
        if other is None or other['id'] == row['id'] or other['info'].get('MATEID') != row['id']:
            raise ValueError('missing or nonreciprocal breakend mate')
        a, b = sorted((row, other), key=lambda x: (order[x['chrom']], x['pos'], x['id']))
        da, db = decode_alt(a['alt'], a['ref']), decode_alt(b['alt'], b['ref'])
        if (da['chrom'], da['pos'], da['remote_side'], da['local_side']) != (b['chrom'], b['pos'], db['local_side'], db['remote_side']):
            raise ValueError('mate coordinates or orientations disagree')
        if (db['chrom'], db['pos']) != (a['chrom'], a['pos']): raise ValueError('reverse mate coordinate mismatch')
        for key in ('SCTG', 'SPAN', 'INSERTION', 'HOMSEQ', 'IMPRECISE', 'NUMPARTS', 'EVDNC'):
            if a['info'].get(key) != b['info'].get(key): raise ValueError('mate INFO mismatch: '+key)
        if a['gt'] != b['gt'] or a['filters'] != b['filters']: raise ValueError('mate genotype or FILTER mismatch')
        if seen.intersection((a['id'], b['id'])): raise ValueError('mate reused')
        seen.update((a['id'], b['id']))
        ct = da['local_side']+'to'+db['local_side']
        kind = 'BND' if a['chrom'] != b['chrom'] else ('INV' if ct in ('3to3','5to5') else ('DUP' if ct == '5to3' else 'DEL'))
        result.append({'a': a, 'b': b, 'connection_type': ct, 'type': kind})
    if len(seen) != len(records): raise ValueError('breakend accounting incomplete')
    return result


class UnresolvedJunction(ValueError):
    """The native call is retained, but a unique linear allele is not established."""


def index_assembly_records(records):
    alignments = defaultdict(list)
    for ordinal, read in enumerate(records):
        if read.is_unmapped: continue
        item = {'ordinal': ordinal, 'start': read.reference_start, 'end': read.reference_end,
                'qstart': read.query_alignment_start, 'qend': read.query_alignment_end,
                'sequence': read.query_sequence, 'reverse': read.is_reverse,
                'secondary': read.is_secondary, 'cigar': read.cigartuples,
                'cigar_string': read.cigarstring}
        for pos, side in ((read.reference_start+1, '5'), (read.reference_end, '3')):
            alignments[(read.query_name, read.reference_name, pos, side)].append(item)
    return alignments


def assembly_alignments(path):
    with pysam.AlignmentFile(str(path)) as reader:
        return index_assembly_records(reader)


def overlap_is_ungapped(alignment, lo, hi):
    """Require a one-to-one query/reference map throughout the shared segment."""
    query = 0
    for op, length in alignment['cigar']:
        if op == 5: return False  # Hard-clipped sequence is unavailable here.
        if op in (0, 1, 4, 7, 8):
            if max(lo, query) < min(hi, query+length) and op not in (0, 7, 8): return False
            query += length
        elif op in (2, 3) and lo <= query <= hi: return False
        elif op not in (2, 3, 6): return False
    return True


def assembled_junction(event, alignments):
    """Use SAM's reference-oriented query, not potentially stale INFO fields alone.

    The pinned SvABA implementation reuses BreakPoint objects without clearing the
    opposite homology/insertion member. The active member must agree with the gap
    between the two endpoint alignments. Multiple identical reconstructions are
    harmless; different reconstructions remain unresolved, regardless of FILTER.
    """
    a, b = event['a'], event['b']; info = a['info']
    if info.get('IMPRECISE'):
        if info.get('HOMSEQ') or info.get('INSERTION'):
            raise UnresolvedJunction('imprecise junction also contains unverified sequence')
        return {'gap': 0, 'sequence': '', 'source': 'approximate_discordant_coordinates',
                'inactive_annotation_fields': [], 'matching_alignment_pairs': []}
    sa, sb = event['connection_type'].split('to'); matches = []
    for x in alignments.get((info['SCTG'], a['chrom'], a['pos'], sa), []):
        for y in alignments.get((info['SCTG'], b['chrom'], b['pos'], sb), []):
            if x['ordinal'] == y['ordinal'] or x['reverse'] != y['reverse']: continue
            if not x['sequence'] or x['sequence'] != y['sequence']: continue
            left, right = (x, y) if event['type'] == 'DEL' else (y, x)
            if left['qstart'] > right['qstart'] or left['qend'] > right['qend']: continue
            gap = right['qstart']-left['qend']; lo, hi = sorted((right['qstart'], left['qend']))
            sequence = left['sequence'][lo:hi].upper()
            active = 'INSERTION' if gap > 0 else ('HOMSEQ' if gap < 0 else None)
            native = info.get(active, '') if active else ''
            if left['reverse']: native = reverse_complement(native)
            if sequence != native.upper(): continue
            ungapped = gap >= 0 or (overlap_is_ungapped(left, lo, hi) and overlap_is_ungapped(right, lo, hi))
            matches.append({'gap': gap, 'sequence': sequence, 'ungapped_overlap': ungapped,
                            'reverse': left['reverse'], 'left_cigar': left['cigar_string'],
                            'right_cigar': right['cigar_string'], 'left_qend': left['qend'],
                            'right_qstart': right['qstart'], 'secondary': left['secondary'] or right['secondary'],
                            'active_annotation_field': active})
    identities = {(x['gap'], x['sequence'], x['ungapped_overlap']) for x in matches}
    if not identities: raise UnresolvedJunction('assembly sequence and orientation cannot be established')
    if len(identities) != 1: raise UnresolvedJunction('ambiguous assembled junction')
    first = matches[0]
    if not first['ungapped_overlap']:
        raise UnresolvedJunction('shared assembly overlap contains a gapped reference alignment')
    return {'gap': first['gap'], 'sequence': first['sequence'], 'source': 'verified_assembly_query',
            'inactive_annotation_fields': [k for k in ('HOMSEQ', 'INSERTION')
                                           if info.get(k) and k != first['active_annotation_field']],
            'matching_alignment_pairs': matches}


def linear_allele(event, reference, alignments):
    """Join retained reference sides with the actual assembled overlap/gap once."""
    a, b = event['a'], event['b']; chrom = a['chrom']; p, q = a['pos'], b['pos']
    if event['type'] not in ('DEL', 'DUP'): return None
    junction = assembled_junction(event, alignments)
    segment = junction['sequence']; hom = max(0, -junction['gap'])
    if set(segment) - set('ACGTN'): raise ValueError('non-DNA assembled sequence')
    if event['type'] == 'DEL':
        left_cut, right_cut = p-hom, q-1+hom
        assumption = 'Reference flanks outside the verified assembled junction are assumed unchanged'
    else:
        left_cut, right_cut = q-hom, p-1+hom
        assumption = 'Interpret the duplication-oriented junction using a tandem reference copy outside the assembled segment'
    n = reference.get_reference_length(chrom)
    if not 1 <= left_cut <= n or not 0 <= right_cut <= n:
        raise UnresolvedJunction('assembled junction exceeds linear reference bounds')
    anchor = left_cut
    if right_cut >= left_cut:
        ref = reference.fetch(chrom, left_cut-1, right_cut).upper(); alt = ref[0]+segment
    else:
        ref = reference.fetch(chrom, left_cut-1, left_cut).upper()
        alt = ref+segment+reference.fetch(chrom, right_cut, left_cut).upper()
    # Preserve the right-junction anchor for reference-copy duplications; trim
    # redundant sequence without choosing coordinates using benchmark outcomes.
    while len(ref)>1 and len(alt)>1 and ref[0]==alt[0]:
        anchor += 1; ref, alt = ref[1:], alt[1:]
    while len(ref)>1 and len(alt)>1 and ref[-1]==alt[-1]: ref, alt = ref[:-1], alt[:-1]
    length = len(alt)-len(ref)
    if ref == alt: raise UnresolvedJunction('assembled junction reconstructs the reference allele')
    matches_left = not hom or reference.fetch(chrom,left_cut,left_cut+hom).upper() == segment
    matches_right = not hom or reference.fetch(chrom,right_cut-hom,right_cut).upper() == segment
    return {'pos': anchor, 'ref': ref, 'alt': alt, 'length': length,
            'type': 'INS' if length > 0 else ('DEL' if length < 0 else 'REPLACEMENT'),
            'native_type': event['type'], 'homology_bases': hom, 'assembly_junction': junction,
            'overlap_matches_both_reference_flanks': matches_left and matches_right,
            'assumption': assumption}


def read_native(path, sample, reference):
    rows=[]
    with pysam.VariantFile(str(path)) as reader:
        if list(reader.header.samples) != [sample]: raise ValueError('native sample column mismatch')
        for r in reader:
            if len(r.alts or ()) != 1:raise ValueError('invalid native allele count')
            verification=reference_status(r,reference)
            rows.append({'id':r.id,'chrom':r.chrom,'pos':r.pos,'ref':r.ref,'alt':r.alts[0],
                         'qual':r.qual,'filters':list(r.filter),'info':dict(r.info),
                         'gt':list(r.samples[sample].get('GT') or ()), 'phased':r.samples[sample].phased,
                         'reference_status':verification})
    return rows


def convert(pilot, output):
    output.mkdir(parents=True,exist_ok=False)
    receipt=json.loads((pilot/'summary.json').read_text())
    if receipt['status'] != 'complete': raise ValueError('native pilot incomplete')
    reference_path=Path(receipt['inputs']['reference']['path']); sample=receipt['inputs']['alignment']['path']
    for name in ('svaba.svaba.unfiltered.indel.vcf','svaba.svaba.unfiltered.sv.vcf','svaba.contigs.bam'):
        if fingerprint(pilot/name)['sha256'] != receipt['artifacts'][name]['sha256']: raise ValueError('native artifact changed')
    alignments=assembly_alignments(pilot/'svaba.contigs.bam'); events=[]; rows=[]; exclusions=[]
    with pysam.FastaFile(str(reference_path)) as reference:
        order={name:i for i,name in enumerate(reference.references)}
        indels=read_native(pilot/'svaba.svaba.unfiltered.indel.vcf',sample,reference)
        bnds=read_native(pilot/'svaba.svaba.unfiltered.sv.vcf',sample,reference)
        pairs=paired_events(bnds,order)
        for r in indels:
            if set(r['ref'].upper()+r['alt'].upper()) - set('ACGTN'): raise ValueError('unresolved native indel')
            n=len(r['alt'])-len(r['ref']); kind='INS' if n>0 else ('DEL' if n<0 else 'REPLACEMENT')
            row={**r,'id':'SVABA_INDEL_'+r['id'],'native_ids':[r['id']], 'type':kind,'length':n,'native_type':'INDEL'}
            rows.append(row)
            events.append({'id':row['id'],'chrom':r['chrom'],'pos_1based':r['pos'],'end_chrom':r['chrom'],
                           'end_position':r['pos']+len(r['ref'])-1,'type':kind,'length':n,'connection_type':None,
                           'genotype':{'GT':r['gt']},'filter':r['filters'],'native_ids':[r['id']]})
        for e in pairs:
            a,b=e['a'],e['b']; identity='SVABA_BND_'+a['id'].split(':')[0]
            events.append({'id':identity,'chrom':a['chrom'],'pos_1based':a['pos'],'end_chrom':b['chrom'],
                           'end_position':b['pos'],'type':e['type'],'connection_type':e['connection_type'],
                           'length':a['info']['SPAN'],'genotype':{'GT':a['gt']},'filter':a['filters'],
                           'endpoint_reference_status':[a['reference_status'],b['reference_status']],
                           'native_ids':[a['id'],b['id']]})
            if e['type'] in ('INV','BND'):
                exclusions.append({'id':identity,'reason':'balanced/interchromosomal class outside DEL/INS benchmark','native_type':e['type']}); continue
            if receipt['mode']=='public' and a['chrom']!='chr20':
                exclusions.append({'id':identity,'reason':'outside the predeclared chr20 benchmark','native_type':e['type']});continue
            try:
                allele=linear_allele(e,reference,alignments)
            except UnresolvedJunction as error:
                # Retain the native pair and report representability, never silently
                # benchmark a subset of PASS linear calls.
                if a['filters']==['PASS']: raise ValueError('unresolved PASS linear junction: '+identity+': '+str(error)) from error
                exclusions.append({'id':identity,'reason':'unresolved linear representation: '+str(error),
                                   'native_type':e['type'],'filter':a['filters'],'native_ids':[a['id'],b['id']],
                                   'native_records':[a,b]});continue
            if allele is None: raise ValueError('linear event was not reconstructed')
            rows.append({**a,**allele,'id':identity,'native_ids':[a['id'],b['id']]})
        header=pysam.VariantHeader()
        for chrom,length in zip(reference.references,reference.lengths):header.contigs.add(chrom,length=length)
        for key,number,kind,description in [('SVTYPE',1,'String','Type derived from resolved allele'),('SVLEN',1,'Integer','ALT minus REF length'),
                ('NATIVE_IDS','.', 'String','Unmodified SvABA record IDs'),('NATIVE_TYPE',1,'String','Original inferred breakend or indel class')]:
            header.info.add(key,number,kind,description)
        for f in sorted({f for r in rows for f in r['filters']}):
            if f != 'PASS':header.filters.add(f,None,None,'Native SvABA filter preserved')
        header.formats.add('GT',1,'String','Native SvABA genotype, without recalibration')
        output_sample='SYNTH_SV' if receipt['mode']=='synthetic' else 'HG002';header.add_sample(output_sample)
        vcf=output/'comparison.full-context.vcf.gz'
        with pysam.VariantFile(str(vcf),'wz',header=header) as writer:
            for r in sorted(rows,key=lambda r:(order[r['chrom']],r['pos'],r['id'])):
                v=writer.new_record(contig=r['chrom'],start=r['pos']-1,id=r['id'],alleles=(r['ref'],r['alt']),qual=r['qual'])
                for f in r['filters']:v.filter.add(f)
                v.info['SVTYPE']=r['type'];v.info['SVLEN']=r['length'];v.info['NATIVE_IDS']=r['native_ids'];v.info['NATIVE_TYPE']=r['native_type']
                v.samples[output_sample]['GT']=r['gt'];v.samples[output_sample].phased=r['phased'];writer.write(v)
        pysam.tabix_index(str(vcf),preset='vcf')
    save(output/'all-native-events.json',events);save(output/'conversion-records.json',rows);save(output/'exclusions.json',exclusions)
    if len(indels)+len(bnds) != sum(len(x['native_ids']) for x in events):raise ValueError('native record accounting failed')
    summary={'completed_at':stamp(),'native_indel_records':len(indels),'native_BND_records':len(bnds),'paired_BND_events':len(pairs),
             'native_event_types':dict(Counter(e['type'] for e in events)),
             'PASS_native_event_types':dict(Counter(e['type'] for e in events if e['filter']==['PASS'])),
             'comparison_records':len(rows),'comparison_types':dict(Counter(r['type'] for r in rows)),
             'excluded_events':len(exclusions),'exclusions_by_reason':dict(Counter(x['reason'] for x in exclusions)),
             'native_nonsequence_boundary_records':sum(r['reference_status']=='virtual_boundary_no_reference_base' for r in bnds),
             'linear_junction_sources':dict(Counter(r['assembly_junction']['source'] for r in rows if 'assembly_junction' in r)),
             'overlap_differs_from_reference':sum(not r['overlap_matches_both_reference_flanks'] for r in rows if 'assembly_junction' in r),
             'inactive_annotation_fields_ignored':dict(Counter(k for r in rows if 'assembly_junction' in r for k in r['assembly_junction']['inactive_annotation_fields'])),
             'unresolved_PASS_linear_events':0,
             'unaccounted_records':0,
             'input_pilot':fingerprint(pilot/'summary.json'),'output_vcf':fingerprint(vcf),'code':fingerprint(Path(__file__)),
             'policy':['Pair by reciprocal MATEID, coordinates, orientation, genotype and FILTER; count each junction once',
                       'Verify actual query overlap/gap and orientation; include assembled overlap once even if it differs from the reference',
                       'Preserve inactive native homology/insertion annotations; do not concatenate stale fields',
                       'Retain unresolved linear pairs with exact reason; abort rather than omit any unresolved PASS linear event',
                       'Reference-copy DUP interpretation is an assumption; native calls and balanced events retained separately',
                       'Keep non-PASS and small-variant context; selection occurs in unchanged native comparator',
                       'Preserve native GT; caller-specific AD/PL/GQ are not converted to standard semantics']}
    save(output/'summary.json',summary);return summary


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('pilot','output'):p.add_argument('--'+name,required=True,type=Path)
    a=p.parse_args();print(json.dumps(convert(a.pilot,a.output)))
