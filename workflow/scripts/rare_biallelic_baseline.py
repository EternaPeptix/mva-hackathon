#!/usr/bin/env python3
"""Phenotype-free autosomal coding baseline with explicit missingness and phase.

Full variant/model tables stay in the selected output directory on the compute
host. Export only the bounded summary for review. Tiers are not probabilities.
"""
import argparse
from collections import Counter, defaultdict
import gzip
import itertools
import json
import os
from pathlib import Path
import re
import subprocess
import time

import pysam

from screen_secondary_findings import fingerprint, original_call, quality_failures, concrete_alleles


def effects_by_gene(entries, config):
    result = defaultdict(dict)
    groups = config['consequence_groups']
    disruptive = set(groups['predicted_disruptive'])
    altering = set(groups['other_protein_altering'])
    for entry in entries:
        if entry.startswith('@'):
            raise ValueError('haplotype pointer unexpected in local-csq output')
        fields = entry.split('|')
        if len(fields) < 4 or fields[3] != 'protein_coding':
            continue
        effects = set(fields[0].split('&'))
        group = 2 if effects & disruptive else 1 if effects & altering else 0
        if not group:
            continue
        gene, transcript = fields[1:3]
        if not gene or not transcript:
            raise ValueError('coding consequence lacks gene/transcript identity')
        previous = result[gene].get(transcript, {'group': 0, 'annotations': []})
        result[gene][transcript] = {'group': max(group, previous['group']),
                                  'annotations': previous['annotations'] + [entry]}
    return dict(result)


def frequency_state(record):
    observed, unavailable = [], []
    for entry in record['frequencies']:
        ac, an = entry['ac'], entry['an']
        if not isinstance(ac, int) or not isinstance(an, int) or ac < 0 or an < 0 or ac > an:
            raise ValueError('invalid public allele counts')
        if an == 0:
            unavailable.append(entry['source'])
        else:
            observed.append(dict(entry, AF_fraction=ac / an))
    maximum = max((row['AF_fraction'] for row in observed), default=None)
    return {'status': 'observed' if observed else 'unresolved', 'max_AF_fraction': maximum,
            'counts': observed, 'unavailable_count_sources': unavailable,
            'allele_found_in_database': record['allele_found']}


def phase_token(call):
    gt = call['GT_values']
    index = call['source_alt_index']
    if len(gt) != 2 or None in gt or gt.count(index) != 1:
        return None
    pgt, pid = call['PGT'], call['PID']
    if isinstance(pgt, str) and re.fullmatch(r'\d+\|\d+', pgt) and pid not in (None, '.', ''):
        phased = list(map(int, pgt.split('|')))
        if sorted(phased) == sorted(gt):
            return (call['source_allele']['chrom'], 'PID', str(pid), phased.index(index))
    ps = call['PS']
    if call['GT_phased'] and ps not in (None, '.', ''):
        return (call['source_allele']['chrom'], 'PS', str(ps), gt.index(index))
    return None


def source_record_identity(call):
    allele = call['source_allele']
    return (allele['chrom'], allele['pos'], allele['ref'], tuple(call['source_all_alts']))


def pair_phase(a, b):
    first, second = a['source_call'], b['source_call']
    if source_record_identity(first) == source_record_identity(second):
        if (first['GT_values'] == second['GT_values'] and first['called_ploidy'] == 2
                and first['source_alt_index'] != second['source_alt_index']
                and first['alt_copies'] == second['alt_copies'] == 1):
            return 'trans_from_diploid_multiallelic_GT'
        raise ValueError('inconsistent alternate alleles from one source record')
    x, y = phase_token(first), phase_token(second)
    if x is not None and y is not None and x[:3] == y[:3]:
        return 'cis' if x[3] == y[3] else 'trans_from_explicit_source_block'
    return 'unknown'


def overlapping(a, b):
    x, y = a['allele'], b['allele']
    return x['chrom'] == y['chrom'] and max(x['pos'], y['pos']) < min(
        x['pos'] + len(x['ref']), y['pos'] + len(y['ref']))


def models_at_frequency(variants, maximum):
    by_gene = defaultdict(list)
    counts = Counter()
    for variant in variants:
        frequency = variant['frequency']['max_AF_fraction']
        if frequency is not None and frequency > maximum:
            counts['variant_over_frequency_threshold'] += 1
            continue
        counts['retained_' + variant['frequency']['status'] + '_variants'] += 1
        for gene in variant['genes']:
            by_gene[gene].append(variant)
    models = []
    for gene, rows in by_gene.items():
        homo = [x for x in rows if x['source_call']['alt_copies'] == 2]
        het = [x for x in rows if x['source_call']['alt_copies'] == 1]
        proposed = [(x, x, 'homozygous_alt') for x in homo]
        for a, b in itertools.combinations(het, 2):
            shared = set(a['genes'][gene]) & set(b['genes'][gene])
            if not shared:
                counts['heterozygous_pair_without_shared_transcript'] += 1
                continue
            relation = pair_phase(a, b)
            if relation == 'cis':
                counts['explicit_cis_pairs_excluded'] += 1
                continue
            if overlapping(a, b) and relation != 'trans_from_diploid_multiallelic_GT':
                counts['overlapping_pair_requiring_review'] += 1
                continue
            proposed.append((a, b, relation))
        for a, b, relation in proposed:
            common = sorted(set(a['genes'][gene]) & set(b['genes'][gene]))
            if not common:
                raise ValueError('model without a shared affected transcript')
            sums = {tr: a['genes'][gene][tr]['group'] + b['genes'][gene][tr]['group'] for tr in common}
            best = max(sums.values())
            tier = 5 - best  # 1: two predicted-disruptive; 2: one; 3: neither.
            bucket = 'observed_rare' if all(x['frequency']['status'] == 'observed' for x in (a, b)) else 'frequency_unresolved'
            models.append({'gene': gene, 'tier': tier, 'frequency_group': bucket,
                'phase': relation, 'allele_keys': [a['key']] if a is b else [a['key'], b['key']],
                'shared_transcripts_at_best_tier': [tr for tr in common if sums[tr] == best]})
    models.sort(key=lambda x: (x['frequency_group'], x['tier'], x['gene'], x['allele_keys']))
    return models, dict(counts)


def summarize_models(models, counters):
    groups = {}
    for group in ('observed_rare', 'frequency_unresolved'):
        rows = [x for x in models if x['frequency_group'] == group]
        best_by_gene = {}
        for row in rows:
            best_by_gene[row['gene']] = min(row['tier'], best_by_gene.get(row['gene'], 99))
        groups[group] = {'models': len(rows), 'genes': len(best_by_gene),
            'best_gene_tier_counts': dict(Counter(best_by_gene.values())),
            'phase_counts': dict(Counter(x['phase'] for x in rows)),
            'genes_by_best_tier': {str(tier): sorted(g for g, t in best_by_gene.items() if t == tier) for tier in (1, 2, 3)},
            'BUB1B_models': [x for x in rows if x['gene'] == 'BUB1B']}
    return {'groups': groups, 'screen_counts': counters,
            'interpretation': 'Tied descriptive consequence tiers, not clinical probabilities or validated recessive diagnoses.'}


def run_pipeline(commands, directory):
    processes, streams = [], []
    try:
        for number, command in enumerate(commands):
            log = (directory / ('pipeline-%d.log' % number)).open('wb')
            streams.append(log)
            previous = processes[-1].stdout if processes else None
            process = subprocess.Popen(command, stdin=previous, stdout=subprocess.PIPE, stderr=log)
            if previous is not None:
                previous.close()
            processes.append(process)
        # The final command writes to --output; consume any informational stdout.
        out, _ = processes[-1].communicate()
        if out:
            (directory / 'pipeline-stdout.txt').write_bytes(out)
        codes = [p.wait() for p in processes]
        if any(codes):
            raise RuntimeError('annotation pipeline failed: ' + str(codes))
    finally:
        for p in processes:
            if p.poll() is None:
                p.terminate()
                try:
                    p.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    p.kill()
                    p.wait()
        for stream in streams:
            stream.close()


def run(args):
    os.nice(10)
    started = time.monotonic()
    root, output = args.root.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    paths = {
        'vcf': root / 'data/WGS_EX2312012_HGWCNDSX7.vcf.gz',
        'reference': root / 'ref/realignment/GRCh38_nochr.fa',
        'gff': root / 'ref/baseline-20260919/Homo_sapiens.GRCh38.116.chr.gff3.gz',
        'frequency_db': root / 'exomiser/2512_hg38/2512_hg38_variants.mv.db',
        'lirical_jar': root / 'lirical/lirical-cli-2.4.1/lirical-cli-2.4.1.jar',
        'frequency_probe': root / 'repo/workflow/scripts/AlleleFrequencyProbe.java',
        'config': args.config.resolve(), 'script': Path(__file__).resolve(),
        'normalization_helpers': Path(__file__).with_name('screen_secondary_findings.py'),
        'bcftools': root / 'tools/bio/bin/bcftools',
        'java': root / 'tools/bio/lib/jvm/bin/java'}
    paths['vcf_index'] = Path(str(paths['vcf']) + '.tbi')
    paths['reference_index'] = Path(str(paths['reference']) + '.fai')
    state = {'status': 'running', 'started_at': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
             'pid': os.getpid(), 'scope': 'Autosomal coding baseline without phenotype, ClinVar or pathogenicity-predictor scoring.'}
    def save():
        (output / 'run-status.json').write_text(json.dumps(state, indent=2) + '\n')
    save()
    try:
        config = json.loads(paths['config'].read_text())
        state['config'] = config
        state['inputs'] = {name: fingerprint(path) for name, path in paths.items()}
        state['stage'] = 'normalize and annotate autosomal variants'
        save()
        annotated = output / 'autosomal-local-csq.bcf'
        commands = [
            [str(paths['bcftools']), 'view', '-r', ','.join(map(str, range(1, 23))), '-Ou', str(paths['vcf'])],
            [str(paths['bcftools']), 'norm', '-f', str(paths['reference']), '-c', 'e', '-m', '-any', '--old-rec-tag', 'SF_ORIG', '-Ou'],
            [str(paths['bcftools']), 'csq', '--local-csq', '--samples', '-', '-f', str(paths['reference']),
             '-g', str(paths['gff']), '-Ob', '-o', str(annotated)]]
        state['commands'] = commands
        run_pipeline(commands, output)
        counts, excluded, seen = Counter(), Counter(), {}
        duplicates = set()
        quality = config['quality']
        with pysam.VariantFile(paths['vcf']) as source, pysam.VariantFile(annotated) as reader:
            if list(source.header.samples) != [config['sample']]:
                raise ValueError('unexpected sample')
            if 'Local consequence' not in reader.header.info['BCSQ'].description:
                raise ValueError('expected local, not inferred-haplotype, consequences')
            for record in reader:
                counts['autosomal_normalized_allele_records'] += 1
                if not concrete_alleles(record):
                    excluded['non_concrete'] += 1
                    continue
                if len(record.alts) != 1:
                    raise ValueError('expected split biallelic record')
                if 'BCSQ' not in record.info:
                    counts['records_without_BCSQ'] += 1
                    continue
                genes = effects_by_gene(record.info['BCSQ'], config)
                if not genes:
                    continue
                counts['protein_altering_allele_records'] += 1
                call = original_call(record, source, config['sample'])
                reasons = quality_failures(call, quality['minimum_DP'], quality['minimum_GQ'])
                if call['called_ploidy'] != 2:
                    reasons.append('not_diploid')
                if reasons:
                    excluded.update(reasons)
                    continue
                key = '%s:%d:%s:%s' % (record.contig, record.pos, record.ref, record.alts[0])
                if key in seen:
                    duplicates.add(key)
                    continue
                seen[key] = {'key': key, 'allele': {'chrom': record.contig, 'pos': record.pos,
                    'ref': record.ref, 'alt': record.alts[0]}, 'source_call': call, 'genes': genes}
        for key in duplicates:
            del seen[key]
        counts['duplicate_normalized_keys_quarantined'] = len(duplicates)
        counts['quality_passing_unique_coding_alleles'] = len(seen)
        (output / 'duplicate-keys.json').write_text(json.dumps(sorted(duplicates), indent=2) + '\n')
        queries = output / 'frequency-queries.txt'
        queries.write_text(''.join(key.replace(':', '\t') + '\n' for key in seen))
        state['stage'] = 'query retained public frequency counts'
        save()
        command = [str(paths['java']), '-Xmx4g', '-XX:ActiveProcessorCount=2', '-cp', str(paths['lirical_jar']),
                   str(paths['frequency_probe']), str(paths['frequency_db'])]
        state['frequency_command'] = command
        with queries.open('rb') as inp, (output / 'frequencies.jsonl').open('wb') as out, (output / 'frequency.log').open('wb') as log:
            subprocess.run(command, stdin=inp, stdout=out, stderr=log, check=True)
        answered = set()
        with (output / 'frequencies.jsonl').open() as stream:
            for line in stream:
                row = json.loads(line)
                key = row['key']
                if key not in seen or key in answered:
                    raise ValueError('unexpected or duplicated frequency result')
                answered.add(key)
                seen[key]['frequency'] = frequency_state(row)
        if answered != set(seen):
            raise ValueError('missing frequency responses')
        variants = list(seen.values())
        counts.update({'frequency_' + k: v for k, v in Counter(x['frequency']['status'] for x in variants).items()})
        with gzip.open(output / 'coding-variants.json.gz', 'wt') as stream:
            json.dump(variants, stream)
        profiles = {}
        limits = [config['frequency']['primary_maximum_AF_fraction']] + config['frequency']['sensitivity_maximum_AF_fractions']
        for limit in limits:
            models, audit = models_at_frequency(variants, limit)
            profiles[str(limit)] = summarize_models(models, audit)
            with gzip.open(output / ('models-AF%s.json.gz' % limit), 'wt') as stream:
                json.dump(models, stream)
        summary = {'scope': state['scope'], 'sample': config['sample'], 'counts': dict(counts),
            'quality_exclusion_reasons': dict(excluded), 'profiles': profiles,
            'BUB1B_variants': [v for v in variants if 'BUB1B' in v['genes']],
            'limitations': config['limitations'], 'config': config}
        (output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
        for name, record in state['inputs'].items():
            current = Path(record['path']).stat()
            if (current.st_size, current.st_mtime_ns) != (record['bytes'], record['mtime_ns']):
                raise ValueError('input changed during run: ' + name)
        state['output_hashes'] = {p.name: fingerprint(p)['sha256'] for p in output.iterdir()
                                 if p.is_file() and p.name != 'run-status.json'}
        state.update(status='complete', stage='completed', completed_at=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
                     elapsed_seconds=time.monotonic() - started)
        print(json.dumps({'status': 'complete', 'elapsed_seconds': state['elapsed_seconds'], 'counts': dict(counts),
                          'profiles': {k: {g: {'models': v['models'], 'genes': v['genes'], 'tier_counts': v['best_gene_tier_counts'],
                                              'BUB1B_models': v['BUB1B_models']} for g, v in value['groups'].items()}
                                       for k, value in profiles.items()}}), flush=True)
    except BaseException as error:
        state.update(status='failed', error=str(error), elapsed_seconds=time.monotonic() - started)
        raise
    finally:
        save()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', required=True, type=Path)
    p.add_argument('--output', required=True, type=Path)
    p.add_argument('--config', required=True, type=Path)
    run(p.parse_args())


if __name__ == '__main__':
    main()
