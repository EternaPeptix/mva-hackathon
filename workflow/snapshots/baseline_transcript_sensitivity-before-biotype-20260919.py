#!/usr/bin/env python3
"""Post-baseline sensitivity to transcript tags in the same retained Ensembl GFF.

This does not replace the all-transcript primary analysis. Full model tables stay
in the compute output directory; summary.json is a bounded review receipt.
"""
import argparse
from collections import Counter
import gzip
import json
import os
from pathlib import Path
import time
from urllib.parse import unquote

from rare_biallelic_baseline import models_at_frequency, summarize_models
from screen_secondary_findings import fingerprint


TAGS = ('MANE_Select', 'Ensembl_canonical')


def tagged_transcripts(lines):
    transcripts = {}
    for line in lines:
        if line.startswith('##FASTA'):
            break
        if line.startswith('#') or not line.strip():
            continue
        fields = line.rstrip('\n').split('\t')
        if len(fields) != 9:
            raise ValueError('invalid GFF row')
        if fields[2] != 'mRNA':
            continue
        attributes = {}
        for entry in fields[8].split(';'):
            if not entry:
                continue
            key, value = entry.split('=', 1)
            if key in attributes:
                raise ValueError('duplicate GFF attribute')
            attributes[key] = unquote(value)
        if attributes.get('biotype') != 'protein_coding':
            continue
        transcript = attributes['transcript_id']
        if attributes['ID'] != 'transcript:' + transcript or transcript in transcripts:
            raise ValueError('ambiguous GFF transcript identity')
        transcripts[transcript] = {
            'chrom': fields[0], 'gene_id': attributes['Parent'].removeprefix('gene:'),
            'name': attributes.get('Name'), 'version': attributes.get('version'),
            'tags': attributes.get('tag', '').split(',') if attributes.get('tag') else []}
    selected = {tag: {tr for tr, row in transcripts.items() if tag in row['tags']} for tag in TAGS}
    if not all(selected.values()):
        raise ValueError('expected transcript tags missing from GFF')
    return transcripts, selected


def restrict_transcripts(variants, allowed):
    result = []
    for variant in variants:
        genes = {gene: {tr: effect for tr, effect in effects.items() if tr in allowed}
                 for gene, effects in variant['genes'].items()}
        genes = {gene: effects for gene, effects in genes.items() if effects}
        if genes:
            result.append(dict(variant, genes=genes))
    return result


def model_identity(model):
    return model['gene'], tuple(sorted(model['allele_keys']))


def verify_subset(original, restricted):
    by_id = {model_identity(row): row for row in original}
    if len(by_id) != len(original) or len({model_identity(row) for row in restricted}) != len(restricted):
        raise ValueError('duplicate model identity')
    tier_changes = Counter()
    for row in restricted:
        before = by_id.get(model_identity(row))
        if before is None:
            raise ValueError('transcript restriction created a new model')
        if any(before[key] != row[key] for key in ('phase', 'frequency_group')):
            raise ValueError('transcript restriction changed phase or frequency status')
        if row['tier'] < before['tier']:
            raise ValueError('transcript restriction improved a consequence tier')
        tier_changes[str(before['tier']) + '->' + str(row['tier'])] += 1
    return {'all_retained_models_present_in_primary': True,
            'unchanged_phase_and_frequency_groups': True,
            'tier_transitions': dict(tier_changes),
            'removed_models': len(original) - len(restricted)}


def checked_file(path, expected_sha):
    result = fingerprint(path)
    if result['sha256'] != expected_sha:
        raise ValueError('checksum mismatch: ' + str(path))
    return result


def run(args):
    os.nice(10)
    started = time.monotonic()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    baseline = args.baseline.resolve()
    state = {'status': 'running', 'started_at': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
             'pid': os.getpid(), 'scope': 'Post hoc transcript-tag sensitivity; all-transcript baseline remains primary.'}
    try:
        status = json.loads((baseline / 'run-status.json').read_text())
        if status['status'] != 'complete':
            raise ValueError('baseline did not complete')
        paths = {'baseline_status': fingerprint(baseline / 'run-status.json'),
                 'script': fingerprint(Path(__file__))}
        paths['gff'] = checked_file(args.gff, status['inputs']['gff']['sha256'])
        for name in ('script', 'normalization_helpers'):
            source = status['inputs'][name]
            paths['baseline_' + name] = checked_file(Path(__file__).with_name(Path(source['path']).name), source['sha256'])
        for name in ('coding-variants.json.gz', 'summary.json'):
            paths[name] = checked_file(baseline / name, status['output_hashes'][name])
        state['inputs'] = paths
        with gzip.open(baseline / 'coding-variants.json.gz', 'rt') as stream:
            variants = json.load(stream)
        primary_summary = json.loads((baseline / 'summary.json').read_text())
        with gzip.open(args.gff, 'rt') as stream:
            transcripts, selected = tagged_transcripts(stream)
        observed_transcripts = {tr for variant in variants for effects in variant['genes'].values() for tr in effects}
        if observed_transcripts - set(transcripts):
            raise ValueError('baseline has coding transcripts absent from this GFF parser')
        with gzip.open(output / 'transcript-tags.json.gz', 'wt') as stream:
            json.dump(transcripts, stream)
        counts = {'protein_coding_mRNA_transcripts': len(transcripts),
                  'protein_coding_genes': len({row['gene_id'] for row in transcripts.values()}),
                  'quality_passing_coding_alleles': len(variants),
                  'affected_coding_transcripts': len(observed_transcripts),
                  'selected': {tag: {'transcripts': len(ids),
                      'genes': len({transcripts[tr]['gene_id'] for tr in ids}),
                      'autosomal_transcripts': sum(transcripts[tr]['chrom'] in set(map(str, range(1, 23))) for tr in ids)}
                      for tag, ids in selected.items()}}
        frequencies = status['config']['frequency']
        limits = [frequencies['primary_maximum_AF_fraction']] + frequencies['sensitivity_maximum_AF_fractions']
        primary_models = {}
        for limit in limits:
            filename = 'models-AF%s.json.gz' % limit
            paths[filename] = checked_file(baseline / filename, status['output_hashes'][filename])
            models, audit = models_at_frequency(variants, limit)
            with gzip.open(baseline / filename, 'rt') as stream:
                if models != json.load(stream):
                    raise ValueError('baseline model recomputation mismatch')
            # JSON round-trip normalizes integer keys in tier counters.
            if json.loads(json.dumps(summarize_models(models, audit))) != primary_summary['profiles'][str(limit)]:
                raise ValueError('baseline summary recomputation mismatch')
            primary_models[str(limit)] = models
        profiles = {}
        for tag, allowed in selected.items():
            restricted = restrict_transcripts(variants, allowed)
            tagged_profiles = {}
            for limit in limits:
                models, audit = models_at_frequency(restricted, limit)
                check = verify_subset(primary_models[str(limit)], models)
                tagged_profiles[str(limit)] = dict(summarize_models(models, audit), subset_validation=check)
                with gzip.open(output / ('models-%s-AF%s.json.gz' % (tag, limit)), 'wt') as stream:
                    json.dump(models, stream)
            profiles[tag] = {'retained_coding_alleles': len(restricted), 'profiles': tagged_profiles}
        focus = {}
        primary_limit = str(frequencies['primary_maximum_AF_fraction'])
        for gene in args.review_gene:
            implicated = {key for row in primary_models[primary_limit] if row['gene'] == gene for key in row['allele_keys']}
            rows = []
            for variant in variants:
                if variant['key'] not in implicated:
                    continue
                effects = variant['genes'][gene]
                rows.append({'allele_key': variant['key'], 'frequency_status': variant['frequency']['status'],
                             'max_AF_fraction': variant['frequency']['max_AF_fraction'],
                             'all_affected_transcripts': sorted(effects),
                             'selected_affected_transcripts': {tag: {tr: dict(transcripts[tr], effect=effects[tr])
                                 for tr in sorted(set(effects) & ids)} for tag, ids in selected.items()}})
            focus[gene] = rows
        summary = {'scope': state['scope'], 'counts': counts, 'profiles': profiles,
                   'baseline_models_and_summaries_reproduced_exactly': True,
                   'review_genes': list(args.review_gene), 'review_selection': 'Post hoc; primary observed-rare genes plus CRYAA and EIF2S2 as transcript-dependent unresolved examples, not an exhaustive review.',
                   'review_gene_alleles': focus,
                   'limitations': ['Uses MANE_Select and Ensembl_canonical tags in Ensembl 116, not a separately updated MANE release.',
                       'Restricting transcripts can miss clinically relevant alternative transcripts or genes without these primary-assembly tags.',
                       'No new frequency, read, phase, expression or functional evidence is created.',
                       'Disappearance under a transcript restriction is not evidence of benignity.',
                       'No phenotype, ClinVar or pathogenicity predictor enters model construction.']}
        (output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
        for record in paths.values():
            current = Path(record['path']).stat()
            if (current.st_size, current.st_mtime_ns) != (record['bytes'], record['mtime_ns']):
                raise ValueError('input changed during sensitivity run')
        state['output_hashes'] = {p.name: fingerprint(p)['sha256'] for p in output.iterdir() if p.is_file()}
        state.update(status='complete', completed_at=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
                     elapsed_seconds=time.monotonic() - started)
        print(json.dumps({'status': state['status'], 'elapsed_seconds': state['elapsed_seconds'], 'counts': counts,
                          'profiles': {tag: {'retained_coding_alleles': data['retained_coding_alleles'],
                              'profiles': {limit: {group: {'models': values['models'], 'genes': values['genes'],
                                  'genes_by_best_tier': values['genes_by_best_tier']} for group, values in profile['groups'].items()}
                                  for limit, profile in data['profiles'].items()}} for tag, data in profiles.items()}}))
    except BaseException as error:
        state.update(status='failed', error=str(error), elapsed_seconds=time.monotonic() - started)
        raise
    finally:
        (output / 'run-status.json').write_text(json.dumps(state, indent=2) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--gff', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--review-gene', action='append', default=[])
    run(parser.parse_args())


if __name__ == '__main__':
    main()
