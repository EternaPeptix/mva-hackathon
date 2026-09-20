#!/usr/bin/env python3
"""Export bounded baseline evidence and exact public-allele follow-up queries.

Full coding/model tables remain in the compute output directory. This is a
post hoc review receipt, not an automated clinical classification.
"""
import argparse
import datetime
import gzip
import json
from pathlib import Path
import sys

import pysam

from baseline_transcript_sensitivity import checked_file
from screen_secondary_findings import fingerprint


def completed_run(directory):
    state = json.loads((directory / 'run-status.json').read_text())
    if state['status'] != 'complete':
        raise ValueError('analysis did not complete')
    checked_file(directory / 'summary.json', state['output_hashes']['summary.json'])
    return state, json.loads((directory / 'summary.json').read_text())


def short_effects(effects):
    result = {}
    for transcript, values in effects.items():
        fields = [entry.split('|') for entry in values['annotations']]
        amino = sorted({row[5] for row in fields if len(row) > 5 and row[5]})
        result[transcript] = {'group': values['group'],
            'effect_names': sorted({effect for row in fields for effect in row[0].split('&')}),
            'amino_acid_changes': [value[:200] for value in amino],
            'amino_acid_text_truncated': any(len(value) > 200 for value in amino)}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--sensitivity', type=Path, required=True)
    parser.add_argument('--clinvar', type=Path, required=True)
    parser.add_argument('--example-gene', action='append', default=[])
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('existing receipt preserved; choose a new output')
    state, summary = completed_run(args.baseline)
    sensitivity_state, sensitivity_summary = completed_run(args.sensitivity)
    checked_file(args.baseline / 'run-status.json', sensitivity_state['inputs']['baseline_status']['sha256'])
    variants_path = args.baseline / 'coding-variants.json.gz'
    checked_file(variants_path, state['output_hashes'][variants_path.name])
    with gzip.open(variants_path, 'rt') as stream:
        variants = {variant['key']: variant for variant in json.load(stream)}
    limit = str(state['config']['frequency']['primary_maximum_AF_fraction'])
    models_path = args.baseline / ('models-AF%s.json.gz' % limit)
    checked_file(models_path, state['output_hashes'][models_path.name])
    with gzip.open(models_path, 'rt') as stream:
        models = json.load(stream)
    selected = {key for row in models if row['frequency_group'] == 'observed_rare'
                or row['gene'] in args.example_gene for key in row['allele_keys']}
    if len(selected) > 50:
        raise ValueError('review exceeds bounded 50-allele export limit')
    public = fingerprint(args.clinvar)
    public_index = fingerprint(Path(str(args.clinvar) + '.tbi'))
    records = []
    with pysam.VariantFile(args.clinvar) as reader:
        for key in sorted(selected):
            variant = variants[key]
            allele = variant['allele']
            matches = []
            for row in reader.fetch(allele['chrom'], allele['pos'] - 1, allele['pos']):
                if row.pos == allele['pos'] and row.ref == allele['ref'] and allele['alt'] in (row.alts or ()):
                    matches.append({'record_id': row.id, 'chrom': row.contig, 'pos': row.pos,
                        'ref': row.ref, 'alts': list(row.alts),
                        'info': {name: row.info[name] for name in ('ALLELEID', 'CLNSIG', 'CLNREVSTAT',
                            'CLNDN', 'CLNDISDB', 'CLNSIGINCL', 'CLNDNINCL', 'CLNDISDBINCL') if name in row.info}})
            call = variant['source_call']
            records.append({'key': key, 'source_call': {name: call[name] for name in (
                'GT_values', 'GT_phased', 'DP', 'GQ', 'PGT', 'PID', 'PS', 'alt_copies', 'source_alt_index')},
                'frequency': variant['frequency'],
                'genes': {gene: short_effects(effects) for gene, effects in variant['genes'].items()},
                'exact_ClinVar_VCF_records': matches})
    for record in (public, public_index):
        stat = Path(record['path']).stat()
        if (stat.st_size, stat.st_mtime_ns) != (record['bytes'], record['mtime_ns']):
            raise ValueError('public input changed during query')
    result = {'created_at': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'baseline': {'run': state, 'summary': summary},
        'transcript_sensitivity': {'run': sensitivity_state, 'summary': sensitivity_summary},
        'candidate_review': {'selection': 'All primary-threshold observed-rare model alleles, plus explicitly named frequency-unresolved example genes; post hoc, not blind.',
            'example_genes': args.example_gene, 'primary_AF_fraction': float(limit),
            'public_ClinVar': public, 'public_index': public_index, 'records': records,
            'limits': ['Exact coordinates/REF/ALT in this dated public VCF only; no equivalent-indel or XML assertion search here.',
                       'No exact record is not evidence of benignity.', 'No new diagnosis, pathogenicity, phase or functional evidence is inferred.']},
        'exporter': fingerprint(Path(__file__)), 'command': [sys.executable] + sys.argv}
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps({'output': str(args.output), 'candidate_alleles': len(records),
                      'sha256': fingerprint(args.output)['sha256']}))


if __name__ == '__main__':
    main()
