#!/usr/bin/env python3
"""Audit multiallelic CNV denominators and emit a compact internal review table."""
import argparse
from collections import Counter
import datetime
import gzip
import json
from pathlib import Path
from review_cnv_context import fingerprint, require, sha, verify, save


def copy_state_distribution(record):
    require(record['SVTYPE'] == 'CNV', 'Copy-state distribution requires a multiallelic CNV')
    info = record['INFO']
    n, counts, frequencies = (info.get(k) for k in ('CN_NUMBER', 'CN_COUNT', 'CN_FREQ'))
    require(type(n) is int and n > 0 and counts and frequencies and len(counts) == len(frequencies), 'Missing copy-state distribution')
    require(all(type(c) is int and c >= 0 for c in counts) and sum(counts) == n, 'Copy-state count total differs')
    require(all(f is not None and abs(c/n - f) < 1e-6 for c, f in zip(counts, frequencies)), 'Copy-state frequencies differ beyond rounding')
    return {'sample_denominator': n, 'copy_states_starting_at_zero': list(range(len(counts))),
            'sample_counts': counts, 'reported_sample_frequencies': frequencies,
            'count_divided_by_sample_denominator': [c/n for c in counts],
            'scope': 'Population sample copy states for this locus, not allele frequency or measured case copy state.'}


def population_example(relation, record):
    return {k: relation[k] for k in ('id', 'minimum_reciprocal_overlap', 'fraction_of_call', 'fraction_of_population_span',
            'start_distance_bp', 'end_distance_bp', 'both_endpoints_within_one_bin', 'population_PASS', 'same_simple_type')} | {
        'chrom': record['chrom'], 'VCF_POS': record['POS'], 'VCF_END': record['END'], 'SVTYPE': record['SVTYPE'],
        'FILTER': record['FILTER'], 'biallelic_frequency': record['biallelic_frequency'],
        'population_copy_states': copy_state_distribution(record) if record['SVTYPE'] == 'CNV' else None}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--context', type=Path, required=True)
    p.add_argument('--expected-summary-sha256', required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    require(not args.output.exists(), 'Existing compact review preserved')
    require(sha(args.context / 'summary.json') == args.expected_summary_sha256, 'Context summary differs')
    source = json.loads((args.context / 'summary.json').read_text())
    for name, record in source['outputs'].items():
        verify(args.context / name, record)
    with gzip.open(args.context / 'all-call-context.json.gz', 'rt') as f:
        calls = json.load(f)['calls']
    with gzip.open(args.context / 'population-overlapping-records.json.gz', 'rt') as f:
        records = json.load(f)['records']
    require(len(calls) == source['call_census']['all']['calls'], 'Call count differs')
    require(len(records) == source['population_census']['totals']['overlapping_records'], 'Population count differs')
    mcnv = [r for r in records.values() if r['SVTYPE'] == 'CNV']
    for record in mcnv:
        copy_state_distribution(record)
    primary, dosage_alerts, mcnv_counts = [], [], Counter()
    for call in calls:
        if call['triage_failures']['primary']:
            continue
        item = {k: call[k] for k in ('call_id', 'chrom', 'start_1based', 'end_inclusive', 'type', 'bin_size',
                'relative_depth', 'MQ_zero_fraction', 'reference_gap_fraction', 'triage_failures', 'gene_overlaps', 'other_resolution_matches')}
        item['marker_context'] = call['markers_inside_nominal_call']
        item['flanks'] = call['flanks']
        item['dosage_context'] = [{k: d[k] for k in ('kind', 'id', 'name', 'chrom', 'start_1based', 'end_inclusive',
            'overlap_bp', 'fraction_of_curated_span', 'fraction_of_call', 'relevant_score_field', 'relevant_score', 'relevant_description')} | {
                'last_evaluated': d['curation']['Date Last Evaluated'],
                'haploinsufficiency_score': d['curation']['Haploinsufficiency Score'],
                'triplosensitivity_score': d['curation']['Triplosensitivity Score']} for d in call['dosage_context']]
        for d in item['dosage_context']:
            if d['relevant_score'] in ('2', '3'):
                dosage_alerts.append({'call_id': call['call_id'], 'context': d,
                    'limit': 'A small partial overlap does not establish disruption of the dosage-sensitive gene or critical syndrome region.'})
        pop = call['population_context']
        item['simple_population_summary'] = {k: pop[k] for k in ('any_nominal_span_overlap', 'same_simple_type_PASS',
            'same_simple_type_PASS_reciprocal_overlap', 'same_simple_type_PASS_both_endpoints_within_one_bin', 'other_SVTYPE_or_multiallelic', 'non_PASS')}
        item['simple_population_examples'] = [population_example(m, m['population_record']) for m in pop['best_overlap_examples']]
        multiallelic = sorted([m for m in call['all_population_relations'] if records[m['id']]['SVTYPE'] == 'CNV'],
                             key=lambda m: (-m['population_PASS'], -m['minimum_reciprocal_overlap'], m['id']))
        item['multiallelic_population_examples'] = [population_example(m, records[m['id']]) for m in multiallelic[:5]]
        for threshold in (.5, .8):
            count = sum(m['population_PASS'] and m['minimum_reciprocal_overlap'] >= threshold for m in multiallelic)
            item['PASS_multiallelic_RO' + str(threshold)] = count
            mcnv_counts[str(threshold)] += count > 0
        primary.append(item)
    require(len(primary) == source['call_census']['primary']['calls'], 'Primary count differs')
    result = {'completed_at': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'scope': 'Compact internal context for every primary-review call; all raw calls and population relations remain in the verified upstream artifacts.',
        'source_summary': fingerprint(args.context / 'summary.json'), 'script': fingerprint(Path(__file__)),
        'verified_source_outputs': source['outputs'], 'source_inputs': source['inputs'],
        'population_census': source['population_census'], 'marker_census': source['marker_census'], 'call_census': source['call_census'],
        'dosage_source_census': {k: (len(v) if isinstance(v, list) else v) for k, v in source['dosage_source_census'].items()},
        'multiallelic_population_audit': {'overlapping_records': len(mcnv), 'all_count_totals_and_frequencies_reconciled': True,
            'FILTER_counts': dict(Counter(','.join(r['FILTER']) for r in mcnv)),
            'primary_calls_with_PASS_multiallelic_match_by_reciprocal_overlap': dict(mcnv_counts)},
        'primary_relevant_dosage_score_2_or_3_overlaps': dosage_alerts, 'primary_call_context': primary,
        'limitations': source['limitations'] + ['Primary subset and overlap cutoffs are descriptive. No case CNV has been assigned ACMG class, EPCR or a drug-response claim by this summary.']}
    save(args.output, result)
    print(json.dumps({'status': 'complete', 'primary_calls':len(primary), 'mCNV_audit':result['multiallelic_population_audit'],
                      'output':fingerprint(args.output)}))


if __name__ == '__main__':
    main()
