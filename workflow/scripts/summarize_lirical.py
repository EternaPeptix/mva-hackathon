#!/usr/bin/env python3
"""Validate LIRICAL receipts and export compact comparisons, not genome-wide variants."""
import argparse
from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import time

TARGETS = ('OMIM:257300', 'OMIM:241200', 'OMIM:601678', 'OMIM:218040',
           'OMIM:615471', 'OMIM:617402', 'OMIM:277440')


def digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def warning_summary(path):
    counts, unknown = Counter(), Counter()
    failed_positions = set()
    inheritance_ids = set()
    for line in path.read_text().splitlines():
        if 'WARN' not in line:
            continue
        message = line.split(' - ', 1)[-1]
        if message.startswith('Unknown contig '):
            counts['unknown_contig'] += 1
            unknown[message.split()[-1]] += 1
        elif message.startswith('Error annotating variant '):
            counts['functional_annotation_error'] += 1
            match = re.search(r':([^: ]+):(\d+):[^:]+>', message)
            if match:
                failed_positions.add((match[1], int(match[2])))
        elif message.startswith('No inheritance mode annotation found'):
            counts['missing_inheritance_default'] += 1
            inheritance_ids.update(re.findall(r'NCBIGene:\d+', message))
        else:
            counts['other'] += 1
    return {
        'counts': dict(counts), 'unknown_contig_distinct_count': len(unknown),
        'unknown_primary_contig_records': {c: n for c, n in unknown.items()
            if re.fullmatch(r'(?:chr)?(?:[1-9]|1\d|2[0-2]|X|Y|M|MT)', c)},
        'missing_inheritance_gene_ids': sorted(inheritance_ids),
        'unique_failed_annotation_positions': len(failed_positions),
        'bub1b_candidate_annotation_failures': sorted(p for c, p in failed_positions
            if c in ('15', 'chr15') and p in (40209701, 40220612)),
        'scope': 'Warning counts, not proof that other annotations are correct.'}


def summarize_run(directory):
    state_path = directory / 'run-status.json'
    state = json.loads(state_path.read_text())
    if state['status'] != 'complete':
        return {'status': state['status'], 'state_sha256': digest(state_path),
                'error_type': state.get('error_type'), 'observed': state['observed'],
                'log_sha256': digest(directory / 'lirical.log')}
    for name in ('lirical.tsv', 'lirical.json', 'lirical.html'):
        if digest(directory / name) != state['output_hashes'][name]:
            raise ValueError(f'{directory.name}: output checksum mismatch: {name}')
    data = json.loads((directory / 'lirical.json').read_text())
    observed, excluded = state['observed'], state['negated']
    subject = data['analysisData']
    if (subject['observedPhenotypicFeatures'] != observed or
            subject['excludedPhenotypicFeatures'] != excluded):
        raise ValueError('output HPO features differ from recorded input')
    genotype_aware = state['mode'] == 'genotype_aware'
    command = state['command']
    if genotype_aware and subject['sampleId'] != command[command.index('--sample-id') + 1]:
        raise ValueError('output sample differs from VCF sample')
    if data['analysisMetadata']['isGlobalAnalysisMode'] != ('--global' in command):
        raise ValueError('global analysis metadata disagrees with command')
    results = data['analysisResults']
    if not results:
        raise ValueError('empty disease analysis')
    by_id = {}
    for result in results:
        ident = result['diseaseId']
        if ident in by_id:
            raise ValueError('duplicate disease ID requires explicit handling')
        by_id[ident] = result
        if ([f['query'] for f in result['observedPhenotypicFeatures']] != observed or
                [f['query'] for f in result['excludedPhenotypicFeatures']] != excluded):
            raise ValueError('per-disease HPO features differ from input')
        prior, lr, post = (result[k] for k in
                          ('pretestProbability', 'compositeLR', 'posttestProbability'))
        if not (all(math.isfinite(x) for x in (prior, lr, post)) and
                0 < prior < 1 and lr >= 0 and 0 <= post <= 1):
            raise ValueError('invalid probability or likelihood ratio')
        expected = lr * prior / (1 - prior + lr * prior)
        if not math.isclose(post, expected, rel_tol=1e-10, abs_tol=1e-15):
            raise ValueError('post-test probability does not match prior and LR')
        individual = [f['lr'] for f in result['observedPhenotypicFeatures'] +
                      result['excludedPhenotypicFeatures']]
        if result['genotypeLR'] is not None:
            individual.append(result['genotypeLR']['lr'])
        if not math.isclose(math.prod(individual), lr, rel_tol=1e-10, abs_tol=1e-100):
            raise ValueError('composite LR does not equal recorded feature/genotype LRs')
    csv.field_size_limit(64 * 1024 * 1024)
    with (directory / 'lirical.tsv').open() as stream:
        rows = list(csv.DictReader((l for l in stream if not l.startswith('!')), delimiter='\t'))
    tsv_by_id, previous_post = {}, float('inf')
    for i, row in enumerate(rows, 1):
        ident = row['diseaseCurie']
        if ident in tsv_by_id or ident not in by_id:
            raise ValueError('duplicate TSV ID or disease missing from JSON')
        tsv_by_id[ident] = row
        result = by_id[ident]
        post, lr = result['posttestProbability'], result['compositeLR']
        if int(row['rank'].replace(',', '')) != i or post > previous_post + 1e-14:
            raise ValueError('TSV ranks/probabilities are not sorted as expected')
        previous_post = post
        if abs(float(row['posttestprob'].rstrip('%')) / 100 - post) > 0.00005001:
            raise ValueError('TSV rounded probability disagrees with JSON')
        if lr > 0 and abs(float(row['compositeLR']) - math.log10(lr)) > 0.00050001:
            raise ValueError('TSV log10 LR disagrees with JSON')
    ordered = sorted(results, key=lambda x: (-x['posttestProbability'], x['diseaseId']))
    ranks = {}
    last_probability, competition_rank = None, None
    for index, result in enumerate(ordered, 1):
        if result['posttestProbability'] != last_probability:
            competition_rank = index
            last_probability = result['posttestProbability']
        ranks[result['diseaseId']] = competition_rank

    def brief(result, details=False):
        ident = result['diseaseId']
        row = tsv_by_id.get(ident)
        value = {k: result[k] for k in ('diseaseId', 'pretestProbability',
                                       'posttestProbability', 'compositeLR', 'genotypeLR')}
        value.update(json_competition_rank=ranks[ident],
                     native_tsv_rank=int(row['rank'].replace(',', '')) if row else None,
                     disease_name=row['diseaseName'] if row else None)
        if details:
            value['phenotype_contributions'] = result['observedPhenotypicFeatures']
            value['displayed_nonzero_score_variants'] = []
            if row and row.get('variants'):
                for variant in row['variants'].split('; '):
                    match = re.search(r'pathogenicity:([0-9.eE+-]+)', variant)
                    if match and float(match[1]) > 0:
                        value['displayed_nonzero_score_variants'].append(variant)
            value['variant_display_caveat'] = ('TSV scores are rounded model annotations; '
                'nonzero does not imply a clinically pathogenic or contributing allele.')
        return value

    return {'status': 'complete', 'mode': state['mode'], 'metadata': data['analysisMetadata'],
            'started_at': state['started_at'], 'completed_at': state['completed_at'],
            'elapsed_seconds': state['elapsed_seconds'], 'command': command,
            'observed': observed, 'excluded': excluded,
            'input_hashes': {k: state[k] for k in ('config_sha256', 'runner_sha256',
                'jar_sha256', 'resource_hashes', 'input_vcf_sha256', 'exomiser_files') if k in state},
            'state_sha256': digest(state_path), 'output_hashes': state['output_hashes'],
            'log_sha256': digest(directory / 'lirical.log'),
            'json_disease_count': len(results), 'tsv_disease_count': len(rows),
            'json_diseases_absent_from_tsv': len(results) - len(rows),
            'sum_of_model_posttest_probabilities': sum(x['posttestProbability'] for x in results),
            'sample_id_note': ('VCF sample identity validated.' if genotype_aware else
                'Upstream phenotype-only output uses subject; HPO inputs independently checked.'),
            'top_json_results': [brief(x) for x in ordered[:15]],
            'target_results': [brief(by_id[t], details=True) for t in TARGETS if t in by_id],
            'warnings': warning_summary(directory / 'lirical.log')}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runs-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    summaries = {p.parent.name: summarize_run(p.parent)
                 for p in sorted(args.runs_root.glob('*/run-status.json'))}
    if not summaries:
        raise ValueError('no recorded runs')
    result = {'schema_version': 1, 'created_at': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
              'script_sha256': digest(Path(__file__)),
              'scope': 'Known-candidate algorithm comparison; model probabilities are not EPCR, '
                       'mutually exclusive diagnosis probabilities, phase evidence or clinical confirmation.',
              'rank_definition': 'JSON rank is competition rank with exact-probability ties; '
                                 'native TSV rank is retained separately and may filter out diseases.',
              'runs': summaries}
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({name: {k: value[k] for k in ('status', 'json_disease_count',
          'tsv_disease_count', 'elapsed_seconds') if k in value} for name, value in summaries.items()}))


if __name__ == '__main__':
    main()
