#!/usr/bin/env python3
"""Bounded SNP audit with explicit reference genotypes and fragment-level evidence.

Missing variant records never become reference calls. Raw reads and regional VCFs
stay outside the repository; only target-level aggregates are exported.
"""
import argparse
from collections import Counter, defaultdict
import datetime
import json
import os
from pathlib import Path
import statistics
import subprocess
import sys
import time

import pysam

from audit_baseline_reads import collapse_fragments
from screen_secondary_findings import fingerprint


def timestamp():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def validate_targets(targets, reference):
    seen = set()
    with pysam.FastaFile(str(reference)) as fasta:
        for row in targets:
            chrom, pos, ref, alt = (row[k] for k in ('chrom', 'pos', 'ref', 'alt'))
            if not isinstance(pos, int) or pos < 1 or ref not in 'ACGT' or alt not in 'ACGT' or len(ref) != 1 or len(alt) != 1 or ref == alt:
                raise ValueError('require distinct concrete single-base alleles')
            if (chrom, pos) in seen:
                raise ValueError('duplicate target position')
            seen.add((chrom, pos))
            if row['key'] != f'{chrom}:{pos}:{ref}:{alt}':
                raise ValueError('target key mismatch')
            if fasta.fetch(chrom, pos - 1, pos).upper() != ref:
                raise ValueError('target reference mismatch')


def fragment_evidence(bam, reference, target, filters):
    chrom, pos, ref, alt = (target[k] for k in ('chrom', 'pos', 'ref', 'alt'))
    groups, excluded, features = defaultdict(list), Counter(), defaultdict(list)
    depth = 0
    with pysam.FastaFile(str(reference)) as fasta, pysam.AlignmentFile(str(bam), 'rb') as reader:
        for column in reader.pileup(chrom, pos - 1, pos, truncate=True, stepper='samtools',
                fastafile=fasta, compute_baq=True, redo_baq=True, ignore_overlaps=False,
                ignore_orphans=True, min_base_quality=0,
                min_mapping_quality=filters['minimum_mapping_quality'],
                flag_filter=int(filters['exclude_flags'], 0), max_depth=100000):
            if column.reference_pos != pos - 1:
                continue
            depth = column.nsegments
            if depth >= 100000:
                raise ValueError('pileup depth cap reached; evidence may be truncated')
            for item in column.pileups:
                read = item.alignment
                group = read.get_tag('RG') if read.has_tag('RG') else ''
                key = (group, read.query_name)
                if item.is_refskip or item.is_del or item.query_position is None:
                    label = 'refskip' if item.is_refskip else 'deletion'
                    groups[key].append((label, read.is_reverse))
                    excluded['no_base_' + label] += 1
                    continue
                query = item.query_position
                if read.query_qualities is None or read.query_qualities[query] < filters['minimum_base_quality']:
                    excluded['low_or_missing_BAQ_base_quality'] += 1
                    continue
                base = read.query_sequence[query].upper()
                label = 'REF' if base == ref else 'ALT' if base == alt else 'other'
                groups[key].append((label, read.is_reverse))
                features[label].append((read.mapping_quality, read.query_qualities[query], read.reference_start))
    result = collapse_fragments(groups, ['REF', 'ALT', 'deletion', 'refskip'])
    result.update(key=target['key'], pileup_segments_before_BQ=depth, excluded_read_counts=dict(excluded),
                  counting='One observation per (read group, query name); disagreeing eligible mates remain conflicts.',
                  BAQ=True, overlapping_mates_premerged_by_pysam=False, filters=filters,
                  orphan_pairs_excluded=True)
    result['read_quality_by_allele'] = {
        label: {'eligible_reads': len(rows), 'median_MQ': statistics.median(x[0] for x in rows),
                'median_BAQ_BQ': statistics.median(x[1] for x in rows),
                'distinct_aligned_starts': len({x[2] for x in rows})}
        for label, rows in features.items()}
    return result


def record_evidence(record, target, sample):
    call = record.samples[sample]
    gt = list(call.get('GT') or [])
    exact = record.pos == target['pos'] and record.ref == target['ref']
    alts = list(record.alts or [])
    target_index = alts.index(target['alt']) + 1 if exact and target['alt'] in alts else None
    if len(gt) != 2 or any(x is None for x in gt):
        state = 'no_call_or_non_diploid'
    elif gt == [0, 0]:
        state = 'hom_reference'
    elif target_index is not None and target_index in gt:
        state = 'target_alt_present'
    else:
        state = 'other_alt_or_representation'
    return {'pos': record.pos, 'end': record.stop, 'ref': record.ref, 'alts': alts,
            'exact_site_ref': exact, 'GT': gt, 'phased': call.phased, 'state': state,
            'DP': call.get('DP'), 'GQ': call.get('GQ'), 'AD': call.get('AD'),
            'PL': call.get('PL'), 'MIN_DP': call.get('MIN_DP'),
            'FILTER': list(record.filter), 'QUAL': record.qual,
            'target_ALT_index': target_index,
            'target_ALT_copies': None if target_index is None or len(gt) != 2 or any(x is None for x in gt) else gt.count(target_index)}


def target_calls(path, targets, sample):
    result = {}
    with pysam.VariantFile(str(path)) as reader:
        if list(reader.header.samples) != [sample]:
            raise ValueError('unexpected VCF sample')
        for row in targets:
            records = [record_evidence(r, row, sample) for r in reader.fetch(row['chrom'], row['pos'] - 1, row['pos'])]
            result[row['key']] = {'records': records, 'missing_record_is_reference': False,
                'exact_site_records': [r for r in records if r['exact_site_ref']],
                'spanning_or_different_ref_records': [r for r in records if not r['exact_site_ref']]}
    return result


def reference_support(reads, caller_rows, policy):
    reasons = []
    counts = reads['counts']
    if counts['REF'] < policy['minimum_reference_fragments']:
        reasons.append('insufficient_reference_fragments')
    if any(n for label, n in counts.items() if label != 'REF'):
        reasons.append('nonreference_or_conflicting_fragments')
    directions = reads['strand_counts'].get('REF', {})
    if policy['require_both_strands'] and not (directions.get('both', 0) or (directions.get('forward', 0) and directions.get('reverse', 0))):
        reasons.append('reference_not_observed_on_both_strands')
    for name, row in caller_rows.items():
        exact = row['exact_site_records']
        if len(exact) != 1:
            reasons.append(name + '_missing_or_multiple_exact_records')
        else:
            call = exact[0]
            if call['state'] != 'hom_reference':
                reasons.append(name + '_not_explicit_homref')
            if call['DP'] is None or call['DP'] < policy['minimum_caller_DP']:
                reasons.append(name + '_low_or_unknown_DP')
            if name in policy['GQ_required_modes'] and (call['GQ'] is None or call['GQ'] < policy['minimum_reference_GQ']):
                reasons.append(name + '_low_or_unknown_GQ')
            if name in policy['target_PL_gap_required_modes']:
                gap = reference_PL_gap(call)
                if gap is None or gap < policy['minimum_target_reference_PL_gap']:
                    reasons.append(name + '_low_or_unknown_reference_PL_gap')
            if any(f != 'PASS' for f in call['FILTER']):
                reasons.append(name + '_filtered')
        if row['spanning_or_different_ref_records']:
            reasons.append(name + '_overlapping_representation_requires_review')
    for name in policy['required_modes']:
        if name not in caller_rows:
            reasons.append(name + '_missing_mode')
    return {'status': 'reference_supported_in_retained_lanes' if not reasons else 'reference_support_not_established',
            'reasons': reasons, 'meaning': policy['meaning']}


def reference_PL_gap(call):
    """Likelihood separation of 0/0 from 0/1 and 1/1, explicitly not GQ."""
    pl = call.get('PL')
    if (call.get('GT') != [0, 0] or call.get('target_ALT_index') != 1 or len(call.get('alts', [])) != 1
            or pl is None or len(pl) != 3 or any(x is None or x < 0 for x in pl)):
        return None
    return min(pl[1:]) - pl[0]


def run(args):
    os.nice(10)
    started = time.monotonic()
    root, output = args.root.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    config = json.loads(args.config.read_text())
    manifest = json.loads((args.slices / 'slice-provenance.json').read_text())
    state = {'status': 'running', 'started_at': timestamp(), 'pid': os.getpid(), 'commands': [],
             'scope': config['scope'], 'invocation': [sys.executable] + sys.argv}
    def save():
        (output / 'run-status.json').write_text(json.dumps(state, indent=2) + '\n')
    save()
    try:
        paths = {'config': args.config, 'slice_manifest': args.slices / 'slice-provenance.json',
                 'reference': root / 'ref/vendor/GRCh38_nochr.fa',
                 'reference_index': root / 'ref/vendor/GRCh38_nochr.fa.fai',
                 'reference_dictionary': root / 'ref/vendor/GRCh38_nochr.dict',
                 'script': Path(__file__), 'read_helpers': Path(__file__).with_name('audit_baseline_reads.py'),
                 'provenance_helpers': Path(__file__).with_name('screen_secondary_findings.py'),
                 'samtools': root / 'tools/bio/bin/samtools', 'bcftools': root / 'tools/bio/bin/bcftools',
                 'java': root / 'tools/bio/lib/jvm/bin/java',
                 'gatk_jar': root / 'gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar'}
        inputs = {name: fingerprint(path) for name, path in paths.items()}
        state['inputs'] = inputs
        if manifest['status'] != 'complete' or inputs['config']['sha256'] != manifest['config_sha256']:
            raise ValueError('slice/config mismatch')
        if inputs['reference']['sha256'] != config['reference_sha256']:
            raise ValueError('reference checksum mismatch')
        reference = paths['reference']
        validate_targets(config['targets'], reference)
        with pysam.FastaFile(str(reference)) as fasta:
            dictionary = dict(zip(fasta.references, fasta.lengths))
        slices = []
        for row in manifest['files']:
            path = args.slices / row['slice_name']
            inputs[row['slice_name']] = fingerprint(path)
            if inputs[row['slice_name']]['sha256'] != row['slice_sha256']:
                raise ValueError('slice checksum mismatch')
            with pysam.AlignmentFile(str(path), 'rb') as reader:
                if dict(zip(reader.references, reader.lengths)) != dictionary:
                    raise ValueError('reference dictionaries differ')
                if {r['SM'] for r in reader.header.to_dict().get('RG', [])} != {config['sample']}:
                    raise ValueError('slice sample differs')
            slices.append(str(path))
        env = dict(os.environ, PATH=str(root / 'tools/bio/bin') + os.pathsep + os.environ.get('PATH', ''))
        def execute(command, label):
            command = list(map(str, command))
            state.update(current_stage=label)
            state['commands'].append(command)
            save()
            with (output / (label + '.log')).open('w') as log:
                subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
        execute(['samtools', 'merge', '-@', '1', '-o', output / 'merged.bam'] + slices, 'merge')
        for command, label in [(['samtools', 'collate', '-@', '1', '-o', output / 'collated.bam', output / 'merged.bam'], 'collate'),
            (['samtools', 'fixmate', '-@', '1', '-m', output / 'collated.bam', output / 'fixmate.bam'], 'fixmate'),
            (['samtools', 'sort', '-@', '1', '-m', '128M', '-o', output / 'sorted.bam', output / 'fixmate.bam'], 'sort'),
            (['samtools', 'markdup', '-@', '1', '-f', output / 'duplicates.txt', output / 'sorted.bam', output / 'markdup.bam'], 'markdup'),
            (['samtools', 'quickcheck', '-v', output / 'markdup.bam'], 'quickcheck'),
            (['samtools', 'index', output / 'markdup.bam'], 'index')]:
            execute(command, label)
        bam = output / 'markdup.bam'
        evidence = {'config': config, 'source_slices': manifest,
                    'reads': {row['key']: fragment_evidence(bam, reference, row, config['read_filters']) for row in config['targets']}}
        regions = [f'{chrom}:{start + 1}-{end}' for chrom, start, end in manifest['regions_BED_0based_half_open']]
        intervals = [x for region in regions for x in ('-L', region)]
        java = [paths['java'], '-XX:ActiveProcessorCount=4', '-Xmx8g', '-jar', paths['gatk_jar']]
        execute(java + ['HaplotypeCaller', '-R', reference, '-I', bam, '-O', output / 'gatk.g.vcf.gz',
            '--pair-hmm-implementation', 'LOGLESS_CACHING', '-ERC', 'BP_RESOLUTION'] + intervals, 'gatk-reference-confidence')
        execute(java + ['GenotypeGVCFs', '-R', reference, '-V', output / 'gatk.g.vcf.gz',
            '-O', output / 'gatk.allsites.vcf.gz', '--include-non-variant-sites', 'true'] + intervals, 'gatk-genotype-allsites')
        execute(['bcftools', 'mpileup', '-Ob', '-o', output / 'pileup.bcf', '-f', reference, '-r', ','.join(regions),
                 '-a', 'AD,DP', '-q', '20', '-Q', '20', '-d', '100000', '--ff', '0xF04', bam], 'bcftools-pileup')
        execute(['bcftools', 'call', '-m', '-A', '-a', 'GQ', '-Oz', '-o', output / 'bcftools.allsites.vcf.gz', output / 'pileup.bcf'], 'bcftools-allsites')
        targets_file = output / 'alleles.tsv'
        targets_file.write_text(''.join(f"{r['chrom']}\t{r['pos']}\t{r['ref']},{r['alt']}\n" for r in config['targets']))
        pysam.tabix_compress(str(targets_file), str(targets_file) + '.gz', force=False)
        pysam.tabix_index(str(targets_file) + '.gz', seq_col=0, start_col=1, end_col=1, zerobased=False)
        execute(['bcftools', 'call', '-m', '-A', '-C', 'alleles', '-T', str(targets_file) + '.gz', '-a', 'GQ',
                 '-Oz', '-o', output / 'bcftools.targeted.vcf.gz', output / 'pileup.bcf'], 'bcftools-targeted')
        evidence['calls'] = {}
        for name, filename in [('GATK_BP_gVCF', 'gatk.g.vcf.gz'), ('GATK_genotyped', 'gatk.allsites.vcf.gz'),
                               ('bcftools_allsites', 'bcftools.allsites.vcf.gz'), ('bcftools_targeted', 'bcftools.targeted.vcf.gz')]:
            path = output / filename
            execute(['bcftools', 'index', '-t', '-f', path], name + '-index')
            evidence['calls'][name] = target_calls(path, config['targets'], config['sample'])
        evidence['reference_support'] = {row['key']: reference_support(evidence['reads'][row['key']],
            {name: data[row['key']] for name, data in evidence['calls'].items()},
            config['reference_support_policy']) for row in config['targets']}
        evidence['bcftools_target_reference_PL_gaps'] = {row['key']: [reference_PL_gap(call)
            for call in evidence['calls']['bcftools_targeted'][row['key']]['exact_site_records']] for row in config['targets']}
        evidence['limitations'] = [
            'Same retained biological reads support both callers and all modes; they are not independent biological replication.',
            'Historical two-lane reads, native diploid calls and no BQSR; full four-lane alignment remains separate.',
            'Reference support concerns the selected SNPs, not exclusion of mosaicism, copy-number changes, or other noncoding alleles.',
            'Targeted bcftools constrains alleles; unconstrained calls and read-level other-allele counts are retained to expose alternatives.',
            'GATK BP-resolution gVCF supplies the intermediate reference-confidence GQ, alongside explicit GenotypeGVCFs all-sites GT and DP.',
            'Reference-only GQ is often omitted from final callers; missing values stay null. The constrained bcftools PL gap is a likelihood separation, not GQ or a posterior probability.',
            'Caller DP and fragment counts apply different filters and need not agree; GQ is not a clinical causal probability.']
        (output / 'evidence.json').write_text(json.dumps(evidence, indent=2) + '\n')
        for record in inputs.values():
            stat = Path(record['path']).stat()
            if (stat.st_size, stat.st_mtime_ns) != (record['bytes'], record['mtime_ns']):
                raise ValueError('input changed during audit')
        state.update(status='complete', completed_at=timestamp(), elapsed_seconds=time.monotonic() - started,
                     output_hashes={p.name: fingerprint(p)['sha256'] for p in output.iterdir() if p.is_file() and p.name != 'run-status.json'})
        print(json.dumps({'status': state['status'], 'elapsed_seconds': state['elapsed_seconds'], 'targets': len(config['targets'])}))
    except BaseException as error:
        state.update(status='failed', error=str(error), elapsed_seconds=time.monotonic() - started)
        raise
    finally:
        save()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--slices', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args())
