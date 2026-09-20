#!/usr/bin/env python3
"""Reconcile frozen and exploratory public results, including scoring scope.

Reads completed artifacts only. BED geometry is a descriptive cross-tab against
native labels, never a substitute for haplotype-aware boundary matching.
"""
import argparse
from bisect import bisect_right
from collections import Counter
import json
import math
import os
from pathlib import Path

import pysam

from audit_public_small_variants import compact, fingerprint, key, kind
from rtg_small_variant_helpers import save, stamp


def read_regions(path):
    intervals = []
    for line in Path(path).read_text().splitlines():
        if not line or line.startswith('#'):
            continue
        chrom, left, right = line.split()[:3]
        left, right = int(left), int(right)
        if chrom != 'chr20' or left < 0 or right <= left:
            raise ValueError('Unexpected evaluation interval')
        if intervals and left < intervals[-1][1]:
            raise ValueError('Intervals must be sorted and nonoverlapping')
        intervals.append((left, right))
    return intervals, [left for left, right in intervals]


def overlaps(intervals, starts, left, right):
    index = bisect_right(starts, right - 1) - 1
    return index >= 0 and intervals[index][1] > left


def checked(pin):
    current = fingerprint(pin['path'])
    if current['sha256'] != pin['sha256']:
        raise ValueError('Completed input/output changed: ' + pin['path'])
    return current


def counts(path):
    with pysam.VariantFile(str(path)) as variants:
        return Counter(key(record) for record in variants)


def scope(path, scored, regions):
    intervals, starts = regions
    remaining = scored.copy()
    categories, totals = Counter(), Counter()
    examples = []
    with pysam.VariantFile(str(path)) as variants:
        for record in variants:
            identity = key(record)
            label = 'scored' if remaining[identity] else 'unscored'
            if remaining[identity]:
                remaining[identity] -= 1
            left, right = record.start, record.start + len(record.ref)
            ref_overlap = overlaps(intervals, starts, left, right)
            anchor_overlap = overlaps(intervals, starts, left, left + 1)
            sample = record.samples['HG002']
            gt, ft = sample.get('GT'), sample.get('FT')
            reasons = []
            if not ref_overlap:
                reasons.append('no_REF_span_overlap')
            if not anchor_overlap:
                reasons.append('anchor_outside')
            if len(record.ref) > 1000 or any(len(alt) > 1000 for alt in (record.alts or ())):
                reasons.append('allele_length_gt_1000')
            if kind(record) == 'non_sequence':
                reasons.append('nonsequence_ALT')
            if ft not in (None, '.', 'PASS', (), ('.',), ('PASS',)):
                reasons.append('sample_FT_not_PASS_or_dot')
            if set(record.filter.keys()) - {'PASS'}:
                reasons.append('record_FILTER_not_PASS_or_dot')
            if gt is None or len(gt) != 2 or any(allele is None for allele in gt):
                reasons.append('not_fully_called_diploid')
            elif not any(gt):
                reasons.append('reference_GT')
            if record.qual is None or not math.isfinite(record.qual) or record.qual < 0:
                reasons.append('QUAL_missing_nonfinite_or_negative')
            categories[(label, tuple(reasons))] += 1
            totals['input_records'] += 1
            totals[label] += 1
            if label == 'unscored':
                totals['unscored_with_REF_overlap' if ref_overlap else 'unscored_without_REF_overlap'] += 1
                if ref_overlap and kind(record) == 'non_sequence':
                    totals['unscored_REF_overlapping_nonsequence'] += 1
                if ref_overlap and not reasons:
                    totals['unscored_REF_overlapping_without_listed_flag'] += 1
            if (label == 'unscored' and ref_overlap) or (label == 'scored' and not anchor_overlap):
                if len(examples) < 6:
                    examples.append(dict(compact(record), native_status=label, descriptive_flags=reasons))
    if any(remaining.values()):
        raise ValueError('Scored record not found in original input')
    return {'totals': dict(totals), 'categories': [
        {'native_status': label, 'descriptive_flags': list(flags), 'records': number}
        for (label, flags), number in sorted(categories.items())], 'bounded_examples': examples}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ('baseline', 'sensitivity', 'output'):
        parser.add_argument('--' + arg, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Existing study preserved; use a new output directory')
    os.nice(10)
    comparisons = {name: json.loads((folder / 'summary.json').read_text())
                   for name, folder in [('frozen_split', args.baseline), ('retained_multiallelic', args.sensitivity)]}
    original, retained = comparisons.values()
    for summary in comparisons.values():
        if summary['status'] != 'complete' or not summary['identical_scored_truth_universe']:
            raise ValueError('Completed, reconciled comparisons required')
    for label in ('config', 'code', 'helpers', 'runtime'):
        if original[label] != retained[label]:
            raise ValueError('Comparator settings or implementation changed: ' + label)
    for label in ('truth', 'truth_index', 'bed', 'reference'):
        if original['inputs'][label] != retained['inputs'][label]:
            raise ValueError('Benchmark inputs changed')
    verified = {'config': checked(original['config']),
                'truth': checked(original['inputs']['truth']), 'bed': checked(original['inputs']['bed'])}
    regions = read_regions(original['inputs']['bed']['path'])
    classes, details = {}, {}
    with pysam.VariantFile(original['inputs']['truth']['path']) as variants:
        for record in variants:
            classes[key(record)] = kind(record)
            details[key(record)] = compact(record)
    all_counters, all_scopes, all_counts = {}, {}, {}
    universe = None
    for arm, summary in comparisons.items():
        all_counters[arm], all_scopes[arm], all_counts[arm] = {}, {}, {}
        for method, result in summary['callers'].items():
            counters = {}
            for label, pin in result['files'].items():
                verified[arm + '/' + method + '/' + label] = checked(pin)
                counters[label] = counts(pin['path'])
                if sum(counters[label].values()) != result['counts'][label]:
                    raise ValueError('Native output recount disagrees')
            truth = counters['tp_baseline'] + counters['fn']
            if universe is None:
                universe = truth
            if universe != truth:
                raise ValueError('Actual truth identities differ across study arms')
            all_counters[arm][method] = counters
            all_counts[arm][method] = result['counts']
            verified[arm + '/' + method + '/input'] = checked(summary['inputs'][method])
            all_scopes[arm][method] = scope(summary['inputs'][method]['path'],
                                           counters['tp_calls'] + counters['fp'], regions)
        all_scopes[arm]['truth'] = scope(summary['inputs']['truth']['path'], universe, regions)
    transitions = {}
    for method in ('gatk', 'bcftools'):
        old = all_counters['frozen_split'][method]['tp_baseline']
        new = all_counters['retained_multiallelic'][method]['tp_baseline']
        parts = {'matched_both_arms': old & new, 'gained_truth_matches': new - old,
                 'lost_truth_matches': old - new, 'missed_both_arms': universe - (old | new)}
        transitions[method] = {label: {'records': sum(rows.values()),
            'record_classes': dict(sum((Counter({classes[k]: n}) for k, n in rows.items()), Counter()))}
            for label, rows in parts.items()}
        transitions[method]['first_five_gained_truth_records'] = [details[k] for k in sorted(new - old)[:5]]
    focal = []
    for identity, detail in details.items():
        if identity[1] == 685696:
            focal.append({'truth': detail, 'matches': {arm: {method: bool(counters['tp_baseline'][identity])
                for method, counters in methods.items()} for arm, methods in all_counters.items()}})
    summary = {'status': 'complete', 'completed_at': stamp(), 'code': fingerprint(Path(__file__)),
        'source_comparisons': {name: fingerprint(folder / 'summary.json')
            for name, folder in [('frozen_split', args.baseline), ('retained_multiallelic', args.sensitivity)]},
        'verified_inputs': verified, 'identical_truth_universe_all_four_comparisons': True,
        'scored_truth_records': sum(universe.values()), 'counts': all_counts,
        'truth_match_transitions': transitions, 'scope_census': all_scopes, 'motivating_site': focal,
        'limits': ['Retained-multiallelic arm was prompted by observed discordances and is exploratory, not prespecified.',
            'No caller was rerun or tuned. Reference, comparator, genotype rules, evaluation region and score threshold are identical.',
            'Call-record denominators differ when multiallelic records are split; FP differences are not one-to-one corrected biological errors.',
            'Scoring-scope flags describe geometry/metadata of native labels, not an independent reimplementation of vcfeval exclusions.',
            'Native matching can cross BED boundaries; preclipping calls by these descriptive flags would be incorrect.',
            'Neither comparison is a patient diagnosis, clinical sensitivity, family phase result or causal probability.']}
    args.output.mkdir(parents=True)
    save(args.output / 'summary.json', summary)
    print(json.dumps({'scored_truth': summary['scored_truth_records'], 'transitions': transitions,
                      'scopes': {a: {m: v['totals'] for m, v in scopes.items()} for a, scopes in all_scopes.items()}}))


if __name__ == '__main__':
    main()
