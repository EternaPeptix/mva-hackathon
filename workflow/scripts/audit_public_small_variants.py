#!/usr/bin/env python3
"""Audit completed genotype-aware public comparisons without changing calls or matching rules.

Reports separate-caller truth overlap, fixed descriptive quality bins and bounded
coordinate-selected discordance examples. It does not construct a merged callset.
"""
import argparse
from collections import Counter, defaultdict
import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import statistics

import pysam


def fingerprint(path):
    path = Path(path)
    before = path.stat()
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(block)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError('Input changed during audit')
    return {'path': str(path), 'bytes': after.st_size, 'sha256': h.hexdigest()}


def key(record):
    sample = next(iter(record.samples))
    return (record.contig, record.pos, record.ref, tuple(record.alts or ()),
            tuple(record.samples[sample].get('GT') or ()))


def kind(record):
    alts = record.alts or ()
    if not alts or any(not alt or any(base.upper() not in 'ACGTN' for base in alt) for alt in alts):
        return 'non_sequence'
    if len(record.ref) == 1 and all(len(alt) == 1 for alt in alts):
        return 'snv'
    if all(len(alt) == len(record.ref) for alt in alts):
        return 'mnp'
    if all(len(alt) != len(record.ref) for alt in alts):
        return 'indel'
    return 'mixed'


def scalar(value):
    return value if isinstance(value, (int, float)) and math.isfinite(value) else None


def quantiles(values):
    values = sorted(value for value in values if value is not None)
    if not values:
        return {'available': 0, 'min': None, 'median': None, 'max': None}
    return {'available': len(values), 'min': values[0], 'median': statistics.median(values), 'max': values[-1]}


def quality_bin(quality):
    if quality is None:
        return 'missing_or_nonfinite'
    if quality < 20:
        return 'QUAL_lt_20'
    if quality < 100:
        return 'QUAL_20_to_lt_100'
    if quality < 500:
        return 'QUAL_100_to_lt_500'
    return 'QUAL_ge_500'


def inspect(path):
    records = Counter()
    classes = Counter()
    quality_bins = Counter()
    genotype_types = Counter()
    values = {name: [] for name in ('QUAL', 'DP', 'GQ')}
    with pysam.VariantFile(str(path)) as vcf:
        if list(vcf.header.samples) != ['HG002']:
            raise ValueError('Public HG002 sample required')
        for record in vcf:
            records[key(record)] += 1
            classes[kind(record)] += 1
            sample = record.samples['HG002']
            quality = scalar(record.qual)
            quality_bins[quality_bin(quality)] += 1
            values['QUAL'].append(quality)
            for field in ('DP', 'GQ'):
                values[field].append(scalar(sample.get(field)))
            gt = sample.get('GT')
            genotype_types[str(tuple(gt) if gt is not None else None)] += 1
    return records, {'records': sum(records.values()), 'record_classes': dict(classes),
                     'fixed_QUAL_bins': dict(quality_bins), 'genotypes': dict(genotype_types),
                     'format_summaries': {name: quantiles(v) for name, v in values.items()}}


def compact(record):
    sample = record.samples['HG002']
    return {'chrom': record.contig, 'pos': record.pos, 'ref': record.ref,
            'alt': list(record.alts or ()), 'GT': sample.get('GT'), 'phased': sample.phased,
            'QUAL': scalar(record.qual), 'DP': scalar(sample.get('DP')),
            'GQ': scalar(sample.get('GQ')), 'AD': sample.get('AD'), 'record_class': kind(record)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--comparison', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--sites-per-class', type=int, default=2)
    args = parser.parse_args()
    if args.output.exists() or not 0 <= args.sites_per_class <= 5:
        parser.error('Use a new output directory and zero to five examples per record class')
    summary_path = args.comparison / 'summary.json'
    summary = json.loads(summary_path.read_text())
    if summary['status'] != 'complete' or not summary['identical_scored_truth_universe']:
        parser.error('Completed comparisons over the identical truth universe required')
    os.nice(10)
    inputs = {}
    for method, result in summary['callers'].items():
        for label, pin in result['files'].items():
            checked = fingerprint(pin['path'])
            if checked['sha256'] != pin['sha256']:
                raise ValueError('Completed native output changed')
            inputs[method + '/' + label] = checked
    for label in ('truth', 'gatk', 'bcftools'):
        pin = summary['inputs'][label]
        checked = fingerprint(pin['path'])
        if checked['sha256'] != pin['sha256']:
            raise ValueError('Original comparison input changed')
        inputs[label] = checked
    args.output.mkdir(parents=True)
    counters = {}
    diagnostics = {}
    class_by_key = {}
    details = {}
    raw_alleles_to_gt = defaultdict(set)
    with pysam.VariantFile(summary['inputs']['truth']['path']) as baseline:
        for record in baseline:
            identity = key(record)
            class_by_key[identity] = kind(record)
            details[identity] = compact(record)
            raw_alleles_to_gt[identity[:4]].add(tuple(sorted(identity[4])))
    for method, result in summary['callers'].items():
        counters[method] = {}
        diagnostics[method] = {}
        for label, pin in result['files'].items():
            counts, stats = inspect(pin['path'])
            if sum(counts.values()) != result['counts'][label]:
                raise ValueError('Independent output census disagrees with completed comparison')
            if stats['record_classes'] != result['record_class_counts'][label]:
                raise ValueError('Record-class census disagrees')
            counters[method][label] = counts
            diagnostics[method][label] = stats
        same_alleles_wrong_gt = 0
        for identity, copies in counters[method]['fp'].items():
            gt_options = raw_alleles_to_gt.get(identity[:4])
            if gt_options and tuple(sorted(identity[4])) not in gt_options:
                same_alleles_wrong_gt += copies
        diagnostics[method]['fp_exact_raw_alleles_with_discordant_GT'] = same_alleles_wrong_gt
    truth_a = counters['gatk']['tp_baseline'] + counters['gatk']['fn']
    truth_b = counters['bcftools']['tp_baseline'] + counters['bcftools']['fn']
    if truth_a != truth_b:
        raise ValueError('Caller truth universes differ')
    match_a = counters['gatk']['tp_baseline']
    match_b = counters['bcftools']['tp_baseline']
    partitions = {'matched_both': match_a & match_b, 'gatk_only': match_a - match_b,
                  'bcftools_only': match_b - match_a,
                  'missed_both': counters['gatk']['fn'] & counters['bcftools']['fn']}
    if sum(partitions.values(), Counter()) != truth_a:
        raise ValueError('Truth overlap does not partition the complete scored universe')
    overlap = {name: {'records': sum(rows.values()),
                       'record_classes': dict(sum((Counter({class_by_key[k]: n}) for k, n in rows.items()), Counter()))}
               for name, rows in partitions.items()}
    combined = (match_a | match_b)
    combined_count = sum(combined.values())
    truth_count = sum(truth_a.values())
    examples = []
    call_files = {name: pysam.VariantFile(summary['inputs'][name]['path']) for name in ('gatk', 'bcftools')}
    try:
        for partition in ('gatk_only', 'bcftools_only', 'missed_both'):
            for cls in ('snv', 'indel'):
                chosen = sorted((k for k in partitions[partition] if class_by_key[k] == cls),
                                key=lambda k: (k[0], k[1], k[2], k[3]))[:args.sites_per_class]
                for identity in chosen:
                    chrom, pos, ref, alts, gt = identity
                    entry = {'truth_match_partition': partition, 'truth': details[identity], 'nearby_original_calls': {}}
                    for name, calls in call_files.items():
                        nearby = [compact(r) for r in calls.fetch(chrom, max(0, pos - 26), pos + max(25, len(ref)))]
                        entry['nearby_original_calls'][name] = {'total_records_in_window': len(nearby),
                                                               'displayed_records': nearby[:20],
                                                               'truncated': len(nearby) > 20}
                    examples.append(entry)
    finally:
        for vcf in call_files.values():
            vcf.close()
    result = {'completed_at': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'scope': summary['scope'], 'comparison': fingerprint(summary_path),
              'code': fingerprint(Path(__file__)), 'verified_inputs': inputs,
              'scored_truth_records': truth_count, 'truth_partitions': overlap,
              'matched_by_at_least_one_caller': {'records': combined_count,
                  'fraction_of_scored_truth': combined_count / truth_count if truth_count else None,
                  'is_validated_merged_callset': False, 'merged_precision': None},
              'output_diagnostics': diagnostics,
              'examples': examples,
              'example_selection': {'sites_per_record_class_per_partition': args.sites_per_class,
                                    'classes': ['snv', 'indel'], 'order': 'ascending chromosome/position/allele',
                                    'window': '25 bases around the truth REF span',
                                    'purpose': 'bounded illustrative raw-call context; not a representative sample or read-level validation'},
              'limits': ['All matching labels come from the unchanged genotype-aware public comparison.',
                         'Overlap of separate truth matches is not a merged callset or merged precision estimate.',
                         'Record classes can differ between equivalent truth and call representations.',
                         'QUAL bins are fixed descriptive bins, not chosen filtering thresholds or comparable probability scales across callers.',
                         'Exact raw-allele GT disagreement is a limited diagnostic; nearby representation differences are not resolved by that count.',
                         'Public benchmark errors are not diagnoses, parental phase or causal probabilities for the child.']}
    (args.output / 'summary.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'scored_truth': truth_count, 'truth_partitions': overlap, 'examples': len(examples)}))


if __name__ == '__main__':
    main()
