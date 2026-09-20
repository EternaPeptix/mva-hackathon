#!/usr/bin/env python3
"""Reconcile completed structural outputs and describe CNV/breakpoint agreement.

Different signatures from the same reads do not constitute biological replication.
Native PASS, genotype scores and breakpoint precision are not clinical probability.
"""
import argparse
from collections import Counter
import datetime
import json
import os
from pathlib import Path
import time
from review_cnv_context import fingerprint, require, sha, verify, save, overlap


def usable_support(row):
    gt = row['genotype'].get('GT')
    return (row['filter'] == ['PASS'] and row['genotype'].get('FT') == 'PASS'
            and gt is not None and len(gt) == 2 and all(x in (0, 1) for x in gt) and 1 in gt
            and all(row['endpoint_reference_bases_in_query_intervals'])
            and not row['native_diploid_sex_genotype_uninterpreted'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('analysis', 'cnv-review', 'gff', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    for name in ('summary', 'cnv-review', 'gff'):
        parser.add_argument('--expected-' + name + '-sha256', required=True)
    args = parser.parse_args(); tick = time.monotonic(); os.nice(10)
    from run_structural_variants import load_mask, review_calls
    from review_four_lane_copy_number import load_genes
    require(not args.output.exists(), 'Existing structural review preserved')
    state = json.loads((args.analysis / 'status.json').read_text())
    require(state['status'] == 'complete' and state['summary_sha256'] == args.expected_summary_sha256
            == sha(args.analysis / 'summary.json'), 'Structural terminal summary differs')
    native = json.loads((args.analysis / 'summary.json').read_text())
    require(native['analysis_status'] == 'complete_with_calls', 'Unexpected native completion scope')
    for name, record in native['outputs'].items():
        verify(args.analysis / name, record)
    verify(args.analysis / 'provenance.json', native['provenance'])
    proof = json.loads((args.analysis / 'provenance.json').read_text())
    for command in native['commands']:
        require(command['exit_code'] == 0, 'Failed native command')
        verify(Path(command['log']['path']), command['log'])
    index = proof['inputs']['reference_index']; verify(Path(index['path']), index)
    dictionary = {}
    for line in Path(index['path']).read_text().splitlines():
        fields = line.split('\t'); require(fields[0] not in dictionary, 'Duplicate reference contig')
        dictionary[fields[0]] = int(fields[1])
    mask = Path(proof['exclusion_summary']['path']).parent
    _, intervals, whole = load_mask(mask, dictionary, native['config'], proof['inputs']['reference'], index)
    require(sha(args.cnv_review / 'summary.json') == args.expected_cnv_review_sha256, 'CNV review differs')
    cnv_review = json.loads((args.cnv_review / 'summary.json').read_text())
    verify(args.cnv_review / 'annotated-calls.json', cnv_review['output'])
    cnv_native_path = Path(cnv_review['inputs']['native_summary']['path'])
    verify(cnv_native_path, cnv_review['inputs']['native_summary'])
    cnv_native = json.loads(cnv_native_path.read_text())
    require(cnv_native['inputs']['alignment']['sha256'] == proof['inputs']['original_alignment']['sha256'],
            'CNV and structural reads differ')
    require(cnv_native['inputs']['reference']['sha256'] == proof['inputs']['reference']['sha256'], 'Reference differs')
    require(sha(args.gff) == args.expected_gff_sha256 == cnv_review['inputs']['gff']['sha256'], 'Gene annotation differs')
    annotation = load_genes(args.gff, set(native['config']['chromosomes']))
    args.output.mkdir(parents=True)
    # Reuse the previously tested native decoder/serializer and verify every byte;
    # explicitly not an independent decoder or caller validation.
    counts = review_calls(args.analysis / 'calls.bcf', args.output / 'reconciled-calls.jsonl', native['config']['sample'],
                         dictionary, native['config']['chromosomes'], intervals, whole)
    require(counts == native['native_discovery'], 'Native count/genotype census differs')
    require(sha(args.output / 'reconciled-calls.jsonl') == native['outputs']['calls-review.jsonl']['sha256'],
            'Native fields differ from original serialized review')
    rows = [json.loads(line) for line in (args.output / 'reconciled-calls.jsonl').read_text().splitlines()]
    cnvs = json.loads((args.cnv_review / 'annotated-calls.json').read_text())['calls']
    comparisons = []
    for call in cnvs:
        if call['triage_failures']['primary']:
            continue
        matches = []
        for row in rows:
            if row['chrom'] != call['chrom'] or row['end_chrom'] != call['chrom'] or row['type'] != {'deletion':'DEL','duplication':'DUP'}[call['type']]:
                continue
            a, b = row['pos_1based'], row['end_position']
            n = overlap(a, b, call['start_1based'], call['end_inclusive'])
            reciprocal = min(n / call['size_bp'], n / max(1, b-a+1))
            if reciprocal >= .5:
                matches.append({'record': row, 'minimum_reciprocal_overlap': reciprocal,
                    'start_distance_bp': abs(a-call['start_1based']), 'end_distance_bp': abs(b-call['end_inclusive']),
                    'both_endpoints_within_one_bin': max(abs(a-call['start_1based']),abs(b-call['end_inclusive'])) <= call['bin_size'],
                    'PASS_nonreference_autosomal_in_query': usable_support(row)})
        comparisons.append({k:call[k] for k in ('call_id','type','chrom','start_1based','end_inclusive','bin_size','relative_depth','gene_overlaps')} | {
            'matches_at_RO50': sorted(matches, key=lambda m:-m['minimum_reciprocal_overlap'])})
    genes = {}
    for name, gene in cnv_review['selected_gene_profiles'].items():
        g = annotation['gene:' + gene['id']]
        require((g['chrom'],g['start'],g['end']) == (gene['chrom'],gene['start'],gene['end']), 'Selected gene annotation differs')
        endpoints, spans = [], []
        for row in rows:
            if ((row['chrom'] == gene['chrom'] and gene['start'] <= row['pos_1based'] <= gene['end']) or
                    (row['end_chrom'] == gene['chrom'] and gene['start'] <= row['end_position'] <= gene['end'])):
                item = dict(row)
                item['gene_CDS_context'] = {
                    'union_CDS_length_bp':sum(b-a+1 for a,b in g['cds']),
                    'endpoint_native_reference_positions_in_CDS':[
                        row['chrom']==g['chrom'] and any(a<=row['pos_1based']<=b for a,b in g['cds']),
                        row['end_chrom']==g['chrom'] and any(a<=row['end_position']<=b for a,b in g['cds'])],
                    'DEL_DUP_INV_nominal_span_union_CDS_overlap_bp':sum(overlap(row['pos_1based'],row['end_position'],a,b) for a,b in g['cds'])
                        if row['type'] in ('DEL','DUP','INV') and row['chrom']==row['end_chrom']==g['chrom'] else None,
                    'limit':'Union across all annotated transcripts; no MANE/splicing, exact affected sequence or transcript consequence is inferred.'}
                endpoints.append(item)
            if row['chrom'] == row['end_chrom'] == gene['chrom'] and row['type'] in ('DEL', 'DUP', 'INV'):
                if overlap(row['pos_1based'], row['end_position'], gene['start'], gene['end']):
                    spans.append(row['id'])
        genes[name] = {'gene_span':{k:gene[k] for k in ('id','chrom','start','end')}, 'native_endpoints_inside':endpoints,
                      'same_chromosome_DEL_DUP_INV_nominal_span_overlap_ids':spans,
                      'limit':'A long imprecise span is not proof of gene disruption; lack of endpoints does not exclude a structural/regulatory mechanism.'}
    save(args.output / 'primary-cnv-structural-context.json', {'calls':comparisons, 'selected_gene_context':genes})
    summary = {'completed_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'elapsed_seconds':time.monotonic()-tick,
        'scope':'Read-only full native output reconciliation and descriptive agreement with the pre-existing primary CNV review subset.',
        'inputs':{'native_summary':fingerprint(args.analysis/'summary.json'),'cnv_review':fingerprint(args.cnv_review/'summary.json'),
                  'gff':fingerprint(args.gff), 'script':fingerprint(Path(__file__))},'native_census':counts,
        'PASS_nonreference_autosomal_in_query_records':sum(usable_support(r) for r in rows),
        'PASS_records_with_native_homozygous_reference_GT':sum(r['filter']==['PASS'] and r['genotype'].get('GT')==[0,0] for r in rows),
        'primary_CNV_call_records':len(comparisons),
        'primary_CNV_records_with_any_RO50_match':sum(bool(c['matches_at_RO50']) for c in comparisons),
        'primary_CNV_records_with_PASS_nonreference_autosomal_in_query_RO50_match':sum(any(m['PASS_nonreference_autosomal_in_query'] for m in c['matches_at_RO50']) for c in comparisons),
        'primary_CNV_records_with_PASS_nonreference_autosomal_in_query_RO80_match':sum(any(m['PASS_nonreference_autosomal_in_query'] and m['minimum_reciprocal_overlap']>=.8 for m in c['matches_at_RO50']) for c in comparisons),
        'primary_CDS_matches':[c for c in comparisons if c['matches_at_RO50'] and any(g['union_CDS_overlap_bp'] for g in c['gene_overlaps'])],
        'selected_gene_endpoint_counts':{g:len(d['native_endpoints_inside']) for g,d in genes.items()},
        'selected_gene_endpoint_CDS_context':{g:d['native_endpoints_inside'] for g,d in genes.items() if d['native_endpoints_inside']},
        'outputs':{name:fingerprint(args.output/name) for name in ('reconciled-calls.jsonl','primary-cnv-structural-context.json')},
        'limitations':['The serializer and decoder are reused; byte equality verifies output consistency, not independent biological accuracy.',
            'Only same-type intrachromosomal DEL/DUP nominal spans enter CNV reciprocal comparisons. Complex/BND components are not merged into inferred dosage alleles.',
            'Site PASS, sample FT=PASS, nonreference diploid GT and retained autosomal endpoints define a descriptive subset; saturated GQ or native PRECISE is not clinical probability.',
            'All raw structural records, homozygous-reference genotypes, failed/sex/end-boundary records and unmatched CNVs remain available.',
            'Public DELLY/SvABA controls showed substantial unresolved SV sensitivity. No match is not a clinical negative, and neither match nor missing match establishes molecular phase, parental origin or causal probability.']}
    save(args.output/'summary.json',summary)
    print(json.dumps({k:summary[k] for k in ('completed_at','elapsed_seconds','PASS_nonreference_autosomal_in_query_records',
        'PASS_records_with_native_homozygous_reference_GT','primary_CNV_records_with_PASS_nonreference_autosomal_in_query_RO50_match',
        'primary_CNV_records_with_PASS_nonreference_autosomal_in_query_RO80_match','selected_gene_endpoint_counts')}))


if __name__ == '__main__':
    main()
