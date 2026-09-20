#!/usr/bin/env python3
"""Post hoc mitochondrial read context and reported whole-reference alternatives.

Detailed reads and sequences remain outside the source repository. Additional
read-context filters are sensitivity analyses, not validated heteroplasmy rules.
"""
import argparse
from collections import Counter, defaultdict
import datetime
import json
import os
from pathlib import Path
import subprocess
import time

import pysam

from review_cnv_context import distribution, fingerprint, require, save, sha
from review_paralog_candidate_pair import record_details, remap_details, compare_reported_mappings
from project_cluster_target_bases import reoriented_query_position, complement


def identity(read):
    require(read.has_tag('RG') and (read.is_read1 != read.is_read2), 'Read identity missing')
    return (read.get_tag('RG'), read.query_name, 1 if read.is_read1 else 2)


def base_context(read, reference_sequence, reference_start, anchor, sites):
    """CIGAR-aware raw context. Deletions, clips and missing bases stay missing."""
    aligned = {r: q for q, r in read.get_aligned_pairs() if r is not None}
    q = aligned.get(anchor - 1)
    qualities = read.query_qualities
    target = None if q is None else {
        'base': read.query_sequence[q].upper(), 'query_position_0based': q,
        'raw_BQ': None if qualities is None else qualities[q],
        'cycle_1based': read.query_length-q if read.is_reverse else q+1,
        'distance_to_nearest_read_end': min(q, read.query_length-1-q)}
    mismatches = []
    for pos, query in aligned.items():
        if query is None:
            continue
        ref = reference_sequence[pos-reference_start].upper()
        alt = read.query_sequence[query].upper()
        if ref in 'ACGT' and alt != ref:
            mismatches.append({'position_1based': pos+1, 'ref': ref, 'alt': alt,
                'raw_BQ': None if qualities is None else qualities[query]})
    extra = [m for m in mismatches if m['position_1based'] != anchor]
    indel_bases = sum(length for op, length in read.cigartuples or [] if op in (1, 2))
    site_bases = {}
    for position, ref, alt in sites:
        query = aligned.get(position-1)
        site_bases[str(position)] = None if query is None else {
            'base': read.query_sequence[query].upper(),
            'raw_BQ': None if qualities is None else qualities[query]}
    return {'anchor': target, 'sites': site_bases, 'mismatches': mismatches,
        'other_mismatches': len(extra),
        'other_Q20_mismatches': sum(m['raw_BQ'] is not None and m['raw_BQ'] >= 20 for m in extra),
        'insertion_deletion_bases': indel_bases, 'other_edit_bases': len(extra)+indel_bases,
        'read_length': read.query_length,
        'raw_quality_counts': dict(Counter(qualities or [])),
        'soft_clipped_bases': sum(n for op,n in read.cigartuples or [] if op == 4)}


def baq_observations(bam, reference, chrom, sites):
    values = defaultdict(dict)
    with pysam.FastaFile(str(reference)) as fa, pysam.AlignmentFile(str(bam), 'rb') as reader:
        for position, ref, alt in sites:
            for column in reader.pileup(chrom, position-1, position, truncate=True, stepper='samtools',
                    fastafile=fa, compute_baq=True, redo_baq=True, ignore_overlaps=False,
                    ignore_orphans=True, min_base_quality=0, min_mapping_quality=20,
                    flag_filter=0xF04, max_depth=100000):
                require(column.nsegments < 100000, 'Pileup depth cap reached')
                for item in column.pileups:
                    read = item.alignment
                    key = identity(read)
                    require(str(position) not in values[key], 'Multiple primary records per site/identity')
                    q = item.query_position
                    values[key][str(position)] = {
                        'base': None if q is None else read.query_sequence[q].upper(),
                        'BAQ_BQ': None if q is None or read.query_qualities is None else read.query_qualities[q],
                        'deletion': item.is_del, 'refskip': item.is_refskip,
                        'flag': read.flag, 'start_0based': read.reference_start, 'cigar': read.cigarstring}
    return values


def select_fragments(rows, anchor, ref, alt, policy):
    groups = defaultdict(list)
    excluded = Counter()
    for row in rows:
        reasons = []
        if row['flag'] & 0xF04: reasons.append('flag')
        if row['mapping_quality'] == 255 or row['mapping_quality'] < policy['MQ']: reasons.append('MQ')
        point = row['context']['anchor']
        if policy['quality'] == 'BAQ':
            base = row['BAQ'].get(str(anchor))
            if base is None:
                reasons.append('absent_from_BAQ_pileup')
                label = None
            elif base['refskip'] or base['deletion']:
                label = 'refskip' if base['refskip'] else 'deletion'
            else:
                label = 'REF' if base['base'] == ref else 'ALT' if base['base'] == alt else 'other'
                if base['BAQ_BQ'] is None or base['BAQ_BQ'] < policy['BQ']: reasons.append('BQ')
        else:
            label = None if point is None else 'REF' if point['base'] == ref else 'ALT' if point['base'] == alt else 'other'
            if point is None or point['raw_BQ'] is None or point['raw_BQ'] < policy['BQ']: reasons.append('missing_or_BQ')
        if policy.get('trim', 0) and (point is None or point['distance_to_nearest_read_end'] < policy['trim']):
            reasons.append('read_end')
        maximum = policy.get('max_other_edits')
        if maximum is not None and row['context']['other_edit_bases'] > maximum: reasons.append('other_edits')
        if reasons:
            excluded.update(reasons)
        else:
            groups[(row['RG'], row['query_name'])].append((label, row['id'], row['reverse']))
    fragments = []
    for (rg, name), observations in sorted(groups.items()):
        labels = {o[0] for o in observations}; strands = {o[2] for o in observations}
        fragments.append({'RG': rg, 'query_name': name,
            'label': next(iter(labels)) if len(labels) == 1 else 'mate_conflict',
            'strand': 'both' if len(strands) == 2 else 'reverse' if True in strands else 'forward',
            'observation_ids': [o[1] for o in observations]})
    counts = Counter(f['label'] for f in fragments)
    return {'policy': policy, 'counts': {label: counts[label] for label in
        ('REF','ALT','deletion','refskip','other','mate_conflict')},
        'excluded_read_reasons_nonexclusive': dict(excluded), 'fragments': fragments}


def project_anchor(old, mapped, anchor, fasta):
    """Project the actual selected base, retaining clipping and orientation changes."""
    point = old['context']['anchor']
    if point is None or mapped.is_unmapped:
        return {'status': 'anchor_absent_or_unmapped'}
    if any(op == 5 for op,n in mapped.cigartuples or []):
        return {'status': 'hard_clip_not_projected'}
    q = reoriented_query_position(point['query_position_0based'], old['context']['read_length'],
                                 old['reverse'], mapped.is_reverse)
    position = dict(mapped.get_aligned_pairs()).get(q)
    if position is None:
        return {'status': 'anchor_clipped_or_inserted'}
    base = complement(point['base']) if old['reverse'] != mapped.is_reverse else point['base']
    reference = fasta.fetch(mapped.reference_name, position, position+1).upper()
    return {'status': 'matched' if reference == base else 'mismatched',
        'chrom': mapped.reference_name, 'position_1based': position+1,
        'oriented_base': base, 'reference_base': reference}


def summarize_group(fragments, observations, mates, mappings, anchor, sites):
    rows = [observations[i] for f in fragments for i in f['observation_ids']]
    mate_rows = [mates[r['mate_id']] for r in rows if r.get('mate_id')]
    result = {'fragments': len(fragments), 'reads': len(rows),
        'by_read_group': dict(Counter(f['RG'] for f in fragments)),
        'fragment_strands': dict(Counter(f['strand'] for f in fragments)),
        'read_number_and_strand': dict(Counter(str(r['read_number'])+('_reverse' if r['reverse'] else '_forward') for r in rows)),
        'original_MQ': dict(Counter(r['mapping_quality'] for r in rows)),
        'original_NM': dict(Counter(r['tags'].get('NM') for r in rows)),
        'original_CIGAR': dict(Counter(r['cigar'] for r in rows)),
        'distinct_original_starts': len({r['start_0based'] for r in rows}),
        'raw_anchor_BQ': dict(Counter(r['context']['anchor']['raw_BQ'] for r in rows if r['context']['anchor'])),
        'BAQ_anchor_BQ': dict(Counter(r['BAQ'][str(anchor)]['BAQ_BQ'] for r in rows if str(anchor) in r['BAQ'])),
        'cycles': distribution([r['context']['anchor']['cycle_1based'] for r in rows if r['context']['anchor']]),
        'distance_to_read_end': distribution([r['context']['anchor']['distance_to_nearest_read_end'] for r in rows if r['context']['anchor']]),
        'other_edit_bases': dict(Counter(r['context']['other_edit_bases'] for r in rows)),
        'other_Q20_mismatches': dict(Counter(r['context']['other_Q20_mismatches'] for r in rows)),
        'mate_observations': len(mate_rows),
        'mate_MQ': dict(Counter(r['mapping_quality'] for r in mate_rows)),
        'mate_NM': dict(Counter(r['tags'].get('NM') for r in mate_rows)),
        'single_end_primary_contig': dict(Counter(mappings[r['id']]['primary']['chrom'] for r in rows)),
        'single_end_primary_MQ': dict(Counter(mappings[r['id']]['primary']['mapping_quality'] for r in rows)),
        'reported_mt_minus_other_AS': dict(Counter(str(mappings[r['id']]['comparison']['gene_minus_other_AS']) for r in rows)),
        'with_reported_nuclear_alternative_reads': sum(any(x['chrom'] != 'MT' and not x['flag'] & 4 for x in mappings[r['id']]['reported_alignments']) for r in rows),
        'mate_single_end_primary_contig': dict(Counter(mappings[r['id']]['primary']['chrom'] for r in mate_rows)),
    }
    result['cluster_patterns_per_read'] = {}
    for quality in ('raw','BAQ'):
        patterns = Counter()
        for row in rows:
            states=[]
            for pos,ref,alt in sites:
                base=row['context']['sites'][str(pos)] if quality == 'raw' else row['BAQ'].get(str(pos))
                q = None if base is None else base.get('raw_BQ' if quality=='raw' else 'BAQ_BQ')
                states.append('.' if q is None or q < 20 else '0' if base['base']==ref else '1' if base['base']==alt else '?')
            patterns['|'.join(states)] += 1
        result['cluster_patterns_per_read'][quality+'20'] = dict(patterns)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('root','config','output'): parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args(); root=args.root.resolve(); output=args.output.resolve()
    require(not output.exists() and not output.is_relative_to(Path(__file__).resolve().parents[2]), 'New external output required')
    config=json.loads(args.config.read_text()); output.mkdir(parents=True); os.nice(10); tick=time.monotonic()
    state={'status':'running','pid':os.getpid(),'started_at':datetime.datetime.now(datetime.timezone.utc).isoformat()}
    save(output/'status.json',state)
    try:
        run=root/config['case_directory']; reference=root/config['reference']; bam=run/'subset.bam'
        require(sha(run/'status.json')==config['case_status_sha256'],'Case receipt differs')
        case=json.loads((run/'status.json').read_text()); require(case['status']=='complete','Case incomplete')
        for name in ('subset.bam','subset.bai'): require(sha(run/name)==case['output_sha256'][name],'Subset identity differs')
        require(sha(reference)==case['config']['reference_sha256'],'Reference differs')
        prior=root/config['prior_review']; require(sha(prior)==config['prior_review_sha256'],'Prior audit differs')
        index_path=root/config['bwa_index_record']; index=json.loads(index_path.read_text())
        require(index['reference_sha256']==case['config']['reference_sha256'],'Index reference differs')
        bwa=root/'tools/bio/bin/bwa'
        require(sha(bwa)==case['source_and_tool_sha256'][str(bwa)],'BWA differs')
        paths={'case_status':run/'status.json','subset':bam,'subset_index':run/'subset.bai','reference':reference,
            'reference_index':Path(str(reference)+'.fai'),'prior_audit':prior,'bwa_index_record':index_path,'bwa':bwa,
            'config':args.config,'source':Path(__file__)}
        for name in ('review_cnv_context.py','review_paralog_candidate_pair.py','project_cluster_target_bases.py'):
            paths[name]=Path(__file__).with_name(name)
        for suffix, expected in index['index_sha256'].items():
            path=Path(str(reference)+suffix);require(sha(path)==expected,'BWA index differs');paths['index'+suffix]=path
        inputs={k:fingerprint(p) for k,p in paths.items()}
        monitored={v['path']:(Path(v['path']).stat().st_size,Path(v['path']).stat().st_mtime_ns) for v in inputs.values()}
        anchor,ref,alt=config['anchor'];sites=config['sites'];chrom=config['chrom']
        with pysam.FastaFile(str(reference)) as fa:
            for pos,r,a in sites: require(fa.fetch(chrom,pos-1,pos).upper()==r and r!=a and len(r)==len(a)==1,'Target reference differs')
        require(config['anchor'] in sites,'Anchor not in cluster')
        baq=baq_observations(bam,reference,chrom,sites)
        observations={}; originals={}; keys={}
        with pysam.FastaFile(str(reference)) as fa,pysam.AlignmentFile(str(bam),'rb') as reader:
            for read in reader.fetch(chrom,anchor-1,anchor):
                if read.is_secondary or read.is_supplementary: continue
                key=identity(read);require(key not in keys,'Repeated primary identity')
                identifier='read_'+str(len(observations));keys[key]=identifier;originals[identifier]=read
                row=record_details(read,identifier)
                row['context']=base_context(read,fa.fetch(chrom,read.reference_start,read.reference_end),read.reference_start,anchor,sites)
                row['BAQ']=baq.get(key,{})
                for value in row['BAQ'].values():
                    require((read.flag,read.reference_start,read.cigarstring)==(value['flag'],value['start_0based'],value['cigar']),'BAQ/raw identity differs')
                observations[identifier]=row
        require(0 < len(observations) <= 20000,'Unexpected observation count')
        wanted={(rg,name,3-n):identifier for (rg,name,n),identifier in keys.items()}; mates={}
        with pysam.AlignmentFile(str(bam),'rb') as reader:
            for read in reader.fetch(until_eof=True):
                if read.is_secondary or read.is_supplementary: continue
                key=identity(read)
                if key not in wanted:continue
                parent=wanted[key];identifier=parent+'_mate';require(identifier not in mates,'Repeated mate identity')
                old=originals[parent]
                require(read.reference_id==old.next_reference_id and read.reference_start==old.next_reference_start
                    and read.next_reference_id==old.reference_id and read.next_reference_start==old.reference_start,'Nonreciprocal mate coordinates')
                originals[identifier]=read;mates[identifier]=record_details(read,identifier);observations[parent]['mate_id']=identifier
        require(len(mates)==len(observations),'Primary pairs incomplete')
        queries=output/'queries.fastq'
        with queries.open('w') as f:
            for identifier,read in originals.items():
                sequence=read.get_forward_sequence();quality=read.get_forward_qualities()
                require(sequence is not None and quality is not None and len(sequence)==len(quality),'Query sequence/quality missing')
                f.write('@'+identifier+'\n'+sequence+'\n+\n'+''.join(chr(q+33) for q in quality)+'\n')
        command=[str(bwa),'mem','-a','-Y','-t','2',str(reference),str(queries)]
        save(output/'commands.json',[command]);state['current_stage']='whole_reference_single_end';save(output/'status.json',state)
        with (output/'remapped.sam').open('w') as sam,(output/'bwa.log').open('w') as log:subprocess.run(command,stdout=sam,stderr=log,check=True)
        remaps=defaultdict(list)
        with pysam.AlignmentFile(str(output/'remapped.sam'),'r') as sam,pysam.FastaFile(str(reference)) as fa:
            for read in sam:
                require(read.query_name in originals,'Unknown remapped query');row=remap_details(read)
                if read.query_name in observations:row['anchor_projection']=project_anchor(observations[read.query_name],read,anchor,fa)
                remaps[read.query_name].append(row)
        require(set(remaps)==set(originals),'Missing remapped query')
        locus={'chrom':chrom,'start_0based':0,'end_exclusive':16569};mappings={}
        for identifier,rows in remaps.items():
            comparison=compare_reported_mappings(rows,locus)
            mappings[identifier]={'primary':next(r for r in rows if not r['flag'] & (256|2048)),
                                 'comparison':comparison,'reported_alignments':rows}
        policies=[{'name':'MQ20_BAQ20','MQ':20,'BQ':20,'quality':'BAQ'},
                  {'name':'MQ60_BAQ20','MQ':60,'BQ':20,'quality':'BAQ'},
                  {'name':'MQ20_raw20','MQ':20,'BQ':20,'quality':'raw'}]
        for bq in (20,30):
            for trim in (0,10,20):
                for edits in (None,0):
                    if bq==20 and trim==0 and edits is None:continue
                    policies.append({'name':f'MQ20_BAQ{bq}_trim{trim}_otherEdits{edits}',
                        'MQ':20,'BQ':bq,'quality':'BAQ','trim':trim,'max_other_edits':edits})
        configurations=[];summaries=[]
        prior_row=next(t['evidence']['original'] for t in json.loads(prior.read_text())['targets'] if t['target']['pos']==anchor)
        for policy in policies:
            selection=select_fragments(list(observations.values()),anchor,ref,alt,policy)
            if policy['name'] in prior_row:require(selection['counts']==prior_row[policy['name']]['counts'],'Prior fragment counts differ')
            configurations.append(selection)
            summary={k:v for k,v in selection.items() if k!='fragments'}
            summary['context']={label:summarize_group([f for f in selection['fragments'] if f['label']==label],
                observations,mates,mappings,anchor,sites) for label in ('REF','ALT','other','mate_conflict')}
            summaries.append(summary)
        save(output/'read-review.json',{'observations':observations,'mates':mates,'configurations':configurations})
        save(output/'remapping.json',mappings)
        require(monitored=={p:(Path(p).stat().st_size,Path(p).stat().st_mtime_ns) for p in monitored},'Input changed')
        summary={'completed_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'configuration':config,'inputs':inputs,
            'pysam':pysam.__version__,'observations':len(observations),'mate_observations':len(mates),
            'policies':summaries,'elapsed_seconds':time.monotonic()-tick,
            'outputs':{p.name:fingerprint(p) for p in output.iterdir() if p.name not in ('summary.json','status.json')},
            'limits':['Post hoc read-context sensitivity analysis, not a calibrated artifact classifier or heteroplasmy detection limit.',
                'Same admitted primary reads and same reference. This does not recruit mitochondrial-like reads elsewhere in the genome.',
                'Single-end remapping reports BWA alternatives, not exhaustive nuclear competition or exclusion of a private NUMT.',
                'Other edits count aligned differences to the reference, including genuine linked alleles; zero-other-edit filtering can remove true variants.',
                'Read-end filters are applied symmetrically, but are exploratory. Reused reads and lanes are not independent biological replication.',
                'Raw and BAQ site patterns use explicit missingness; overlapping mates collapse and conflicts remain. No causal or treatment conclusion.']}
        save(output/'summary.json',summary);state.update(status='complete',completed_at=summary['completed_at'],summary_sha256=sha(output/'summary.json'))
        save(output/'status.json',state)
        print(json.dumps({'status':'complete','elapsed_seconds':summary['elapsed_seconds'],'observations':len(observations)}))
    except BaseException as error:
        state.update(status='failed',error_type=type(error).__name__,error=str(error));save(output/'status.json',state);raise


if __name__=='__main__':main()
