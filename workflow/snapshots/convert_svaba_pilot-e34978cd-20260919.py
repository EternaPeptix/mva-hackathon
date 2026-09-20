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


def assembly_orientations(path):
    orientations = defaultdict(set)
    with pysam.AlignmentFile(str(path)) as reader:
        for read in reader:
            if read.is_unmapped or read.is_secondary: continue
            for pos, side in ((read.reference_start+1, '5'), (read.reference_end, '3')):
                orientations[(read.query_name, read.reference_name, pos, side)].add(read.is_reverse)
    return orientations


def linear_allele(event, reference, orientations):
    """Join the two retained sides, subtracting their shared contig overlap once."""
    a, b = event['a'], event['b']; chrom = a['chrom']; p, q = a['pos'], b['pos']
    if event['type'] not in ('DEL', 'DUP'): return None
    info = a['info']; insert = info.get('INSERTION', '') or ''; homology = info.get('HOMSEQ', '') or ''
    if set(insert.upper()+homology.upper()) - set('ACGTN'): raise ValueError('non-DNA insertion or homology')
    hom = len(homology); reversed_insert = None
    if insert:
        sa, sb = event['connection_type'].split('to'); sctg = info['SCTG']
        choices = orientations.get((sctg,chrom,p,sa),set()) & orientations.get((sctg,chrom,q,sb),set())
        if len(choices) != 1: raise ValueError('assembly orientation cannot be established for insertion')
        reversed_insert = choices.pop()
        if reversed_insert: insert = reverse_complement(insert)
    if event['type'] == 'DEL':
        anchor, end = p-hom, q-1
        if not 1 <= anchor <= end <= reference.get_reference_length(chrom): raise ValueError('invalid deletion/replacement bounds')
        ref = reference.fetch(chrom,anchor-1,end).upper(); alt = ref[0]+insert.upper()
        hom_equal = reference.fetch(chrom,p-hom,p).upper() == reference.fetch(chrom,q-1,q-1+hom).upper()
        assumption = 'Join reference prefix through left breakpoint minus shared homology to inserted sequence and right suffix'
    else:
        anchor = q-hom
        if not 1 <= p <= anchor <= reference.get_reference_length(chrom): raise ValueError('invalid duplication bounds')
        ref = reference.fetch(chrom,anchor-1,anchor).upper()
        alt = ref+insert.upper()+reference.fetch(chrom,p-1,anchor).upper()
        hom_equal = reference.fetch(chrom,q-hom,q).upper() == reference.fetch(chrom,p-1,p-1+hom).upper()
        assumption = 'Interpret inward duplication junction as a tandem reference copy plus any oriented assembled insertion'
    if not hom_equal: raise ValueError('breakpoint homology is inconsistent with reference flanks')
    length = len(alt)-len(ref)
    return {'pos': anchor, 'ref': ref, 'alt': alt, 'length': length,
            'type': 'INS' if length > 0 else ('DEL' if length < 0 else 'REPLACEMENT'),
            'native_type': event['type'], 'homology_bases': hom, 'insert_reverse_complemented': reversed_insert,
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
    orientations=assembly_orientations(pilot/'svaba.contigs.bam'); events=[]; rows=[]; exclusions=[]
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
            allele=linear_allele(e,reference,orientations)
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
             'unaccounted_records':0,
             'input_pilot':fingerprint(pilot/'summary.json'),'output_vcf':fingerprint(vcf),'code':fingerprint(Path(__file__)),
             'policy':['Pair by reciprocal MATEID, coordinates, orientation, genotype and FILTER; count each junction once',
                       'Subtract shared contig homology once; orient inserted sequence from matching assembly alignments',
                       'Reference-copy DUP interpretation is an assumption; native calls and balanced events retained separately',
                       'Keep non-PASS and small-variant context; selection occurs in unchanged native comparator',
                       'Preserve native GT; caller-specific AD/PL/GQ are not converted to standard semantics']}
    save(output/'summary.json',summary);return summary


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('pilot','output'):p.add_argument('--'+name,required=True,type=Path)
    a=p.parse_args();print(json.dumps(convert(a.pilot,a.output)))
