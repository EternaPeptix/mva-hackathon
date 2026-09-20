#!/usr/bin/env python3
"""Describe mate anchoring without treating single-end uncertainty as a false call."""
import argparse
from collections import Counter
import datetime
import json
from pathlib import Path

import pysam
from review_cnv_context import fingerprint, require, save, sha, verify


def summarize_mates(review, remapping, primary, locus):
    observed = {r['id']:r for r in review['observations']}
    mates = {r['id']:r for r in review['observed_mates']}
    phase = review['phase_thresholds']['20']
    rows = []
    for fragment in phase['fragments']:
        kept = {}
        missing = []
        for identifier in fragment['observation_ids']:
            observation = observed[identifier]
            if 'mate_id' not in observation:
                missing.append(observation.get('mate_status', 'unavailable'))
                continue
            mate = mates[observation['mate_id']]
            require((mate['RG'],mate['query_name']) == (fragment['RG'],fragment['query_name']), 'Mate fragment differs')
            reciprocal = (mate['mate_chrom'] == observation['chrom']
                          and mate['mate_start_0based'] == observation['start_0based']
                          and observation['mate_chrom'] == mate['chrom']
                          and observation['mate_start_0based'] == mate['start_0based'])
            p = primary[mate['id']]
            comparison = remapping['comparisons'][mate['id']]
            key = (mate['RG'],mate['query_name'],mate['read_number'])
            kept[key] = {'mate_id':mate['id'], 'reciprocal_coordinates':reciprocal,
                         'original_flag':mate['flag'], 'original_MQ':mate['mapping_quality'],
                         'original_NM':mate['tags'].get('NM'),
                         'original_eligible':not mate['flag'] & 0xF04 and mate['mapping_quality'] != 255 and mate['mapping_quality'] >= 20,
                         'single_end_primary':p, 'score_comparison':comparison}
        rows.append({'RG':fragment['RG'], 'query_name':fragment['query_name'], 'label':fragment['label'],
                     'mates':list(kept.values()), 'unavailable_observation_mates':missing})
    groups = {}
    for label in sorted({r['label'] for r in rows}):
        fragments = [r for r in rows if r['label'] == label]
        all_mates = [m for f in fragments for m in f['mates']]
        group = {'fragments':len(fragments), 'retrieved_mate_observations':len(all_mates),
                 'fragments_with_nonreciprocal_mate':sum(any(not m['reciprocal_coordinates'] for m in f['mates']) for f in fragments),
                 'original_mate_MQ':dict(Counter(m['original_MQ'] for m in all_mates)),
                 'single_end_mate_MQ':dict(Counter(m['single_end_primary']['MQ'] for m in all_mates)),
                 'gene_minus_other_AS':dict(Counter('unreported_comparison' if m['score_comparison']['gene_minus_other_AS'] is None else str(m['score_comparison']['gene_minus_other_AS']) for m in all_mates)),
                 'thresholds':{}}
        for threshold in (20,40,60):
            def qualifies(m, inside):
                p = m['single_end_primary']
                in_gene = (p['chrom'] == locus['chrom'] and p['end_exclusive'] is not None
                           and p['start_0based'] < locus['end_exclusive'] and p['end_exclusive'] > locus['start_0based'])
                return (m['original_eligible'] and m['reciprocal_coordinates'] and not p['flag'] & 4
                        and p['MQ'] != 255 and p['MQ'] >= threshold and in_gene == inside)
            group['thresholds'][str(threshold)] = {
                'fragments_with_any_eligible_mate_primary_overlapping_gene':sum(any(qualifies(m,True) for m in f['mates']) for f in fragments),
                'fragments_with_any_eligible_mate_primary_outside_gene':sum(any(qualifies(m,False) for m in f['mates']) for f in fragments),
                'fragments_with_gene_mate_and_at_least_90_percent_available_query_aligned':sum(any(qualifies(m,True) and m['single_end_primary']['query_length'] and m['single_end_primary']['aligned_query_length'] >= 0.9*m['single_end_primary']['query_length'] for m in f['mates']) for f in fragments)}
        groups[label] = group
    return {'groups':groups, 'fragments':rows}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--analysis',type=Path,required=True)
    parser.add_argument('--expected-summary-sha256',required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    require(not args.output.exists() and not args.output.resolve().is_relative_to(Path(__file__).resolve().parents[2]), 'Use new output outside source tree')
    require(sha(args.analysis/'summary.json') == args.expected_summary_sha256, 'Parent summary differs')
    parent=json.loads((args.analysis/'summary.json').read_text())
    for name in ('read-review.json','single-end-remap.json','single-end-remap.sam'):
        verify(args.analysis/name,parent['outputs'][name])
    review=json.loads((args.analysis/'read-review.json').read_text())
    remapping=json.loads((args.analysis/'single-end-remap.json').read_text())
    primary={}
    with pysam.AlignmentFile(str(args.analysis/'single-end-remap.sam'),'r') as sam:
        for r in sam:
            if r.is_secondary or r.is_supplementary:continue
            require(r.query_name not in primary,'Duplicate primary query')
            primary[r.query_name]={'chrom':r.reference_name,'start_0based':r.reference_start,
                'end_exclusive':r.reference_end,'flag':r.flag,'MQ':r.mapping_quality,'cigar':r.cigarstring,
                'query_length':r.query_length,'aligned_query_length':r.query_alignment_length,'tags':dict(r.tags)}
    result=summarize_mates(review,remapping,primary,parent['configuration']['gene_locus'])
    args.output.mkdir(parents=True)
    save(args.output/'mate-evidence.json',result)
    summary={'completed_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'parent_summary':fingerprint(args.analysis/'summary.json'),'script':fingerprint(Path(__file__)),
        'groups':result['groups'],'mate_evidence':fingerprint(args.output/'mate-evidence.json'),
        'scope':'Companion sensitivity from existing read/remap outputs; no realignment repeated, new specimen or independent validation.',
        'limitations':['Both mate and spanning-read information may come from the same fragment. Counts across categories must not be added.',
            'A favorable mate mapping is supporting context, not demonstrated unique locus origin. Reference incompleteness, mismatches, clipping and library artifacts remain possible.',
            'The 90 percent threshold uses available sequence length, not missing hard-clipped bases. MQ/score changes are not a clinical classification.']}
    save(args.output/'summary.json',summary)
    print(json.dumps(summary))


if __name__=='__main__':main()
