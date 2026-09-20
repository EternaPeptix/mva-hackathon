#!/usr/bin/env python3
"""Reconcile completed marker/ROH artifacts and write a bounded research summary.

This reads existing outputs; it does not rerun a caller or establish accuracy,
ancestry, parental origin, copy neutrality, clinical status, or causal probability.
The descriptive length bins are post-run review choices, not diagnostic cutoffs.
"""
import argparse
from collections import Counter
import datetime
import hashlib
from itertools import zip_longest
import json
from pathlib import Path
import re

import pysam


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def check_file(path, expected):
    require(path.is_file() and path.stat().st_size == expected['bytes'], 'file size differs: ' + str(path))
    require(sha256(path) == expected['sha256'], 'file hash differs: ' + str(path))


def signature(record, sample):
    return (record.contig, record.pos, record.alleles, dict(record.samples[sample]),
            record.samples[sample].phased)


def category(record, sample):
    gt = record.samples[sample].get('GT')
    if not gt or any(a is None for a in gt):
        return 'missing'
    require(len(gt) == 2 and all(a in (0, 1) for a in gt), 'non-biallelic diploid genotype')
    return ('homref', 'heterozygous', 'homalt')[sum(gt)]


def review(root, expected_sha, review_regions):
    state = json.loads((root / 'status.json').read_text())
    require(state['status'] == 'complete', 'analysis must be complete')
    actual_sha = sha256(root / 'summary.json')
    require(actual_sha == expected_sha == state['summary_sha256'], 'terminal summary identity differs')
    source = json.loads((root / 'summary.json').read_text())
    require(source['analysis_status'] == 'completed_requested_models', 'models did not complete')
    require(source['config']['scope'] == 'verified_four_lane_autosomal_marker_genotyping', 'unexpected input scope')
    require(source['completed_at'] == state['completed_at'], 'completion timestamps differ')
    sample = source['config']['sample']
    for name, record in source['outputs'].items():
        require(Path(record['path']).resolve() == (root / name).resolve(), 'output escapes selected analysis')
        check_file(root / name, record)
    for item in source['commands']:
        require(item['exit_code'] == 0, 'native stage did not succeed')
        check_file(Path(item['log']['path']), item['log'])

    requested, qualified, excluded, missing_gq, raw_categories = Counter(), Counter(), Counter(), Counter(), Counter()
    qualified_categories, qualified_missing_gq, reasons = Counter(), Counter(), Counter()
    markers, positions = {}, {}
    with pysam.VariantFile(str(root / 'calls.bcf')) as raw, \
            pysam.VariantFile(str(root / 'annotated-all-markers.vcf.gz')) as annotated, \
            pysam.VariantFile(str(root / 'qualified-markers.vcf.gz')) as usable:
        for reader in (raw, annotated, usable):
            require(list(reader.header.samples) == [sample], 'sample identity differs')
        usable_iter = iter(usable)
        for call, row in zip_longest(raw, annotated):
            require(call is not None and row is not None, 'raw/annotated marker count differs')
            require(signature(call, sample) == signature(row, sample), 'annotation changed native genotype evidence')
            label = category(row, sample)
            requested[row.contig] += 1; raw_categories[label] += 1
            if row.samples[sample].get('GQ') is None:
                missing_gq[label] += 1
            if row.info.get('MVA_USABLE', False):
                require('MVA_EXCLUDE' not in row.info, 'usable marker also excluded')
                retained = next(usable_iter, None)
                require(retained is not None and signature(row, sample) == signature(retained, sample),
                        'qualified evidence differs or marker missing')
                require(dict(retained.info) == dict(row.info), 'qualified annotation differs')
                require(label != 'missing', 'qualified marker has missing genotype')
                qualified[row.contig] += 1; qualified_categories[label] += 1
                if row.samples[sample].get('GQ') is None:
                    qualified_missing_gq[label] += 1
                key = (row.contig, row.pos)
                require(key not in markers, 'duplicate qualified position')
                markers[key] = label
                positions.setdefault(row.contig, []).append(row.pos)
            else:
                excluded[row.contig] += 1
                require(bool(row.info.get('MVA_EXCLUDE')), 'excluded marker has no reason')
                reasons.update(row.info['MVA_EXCLUDE'])
        require(next(usable_iter, None) is None, 'extra qualified marker')
    counts = source['marker_qualification']['counts']
    require(sum(requested.values()) == counts['requested_and_returned_markers'], 'requested count differs')
    require(sum(qualified.values()) == counts['qualified_markers'], 'qualified count differs')
    require(sum(excluded.values()) == counts['excluded_markers'], 'excluded count differs')
    require(dict(reasons) == source['marker_qualification']['exclusion_reason_counts_nonexclusive'], 'exclusion counts differ')
    require(sum(missing_gq.values()) == counts['raw_GQ_missing'], 'missing GQ count differs')
    for label, n in raw_categories.items():
        require(n == counts['raw_genotype_' + label], 'raw genotype census differs')
    for chrom in map(str, range(1, 23)):
        saved = source['marker_qualification']['per_chromosome'][chrom]
        require((requested[chrom], qualified[chrom], excluded[chrom]) == (
            saved['requested_and_returned_markers'], saved['qualified_markers'], saved.get('excluded_markers', 0)),
            'chromosome marker counts differ')
        require(positions[chrom] == sorted(set(positions[chrom])), 'qualified positions unordered')
    extraction = source['read_extraction']
    native = extraction['native_counts']
    require(native['records_processed'] == native['records_filter_accepted'] + native['records_filter_rejected'],
            'read extraction does not conserve records')
    require(native['records_filter_accepted'] == sum(extraction['accepted_by_chromosome'].values())
            == sum(extraction['accepted_by_read_group'].values()), 'read-group counts differ')

    models, native_files, observed_models = [], 0, set()
    for model in source['homozygosity_models']:
        tag, gap = model['population_af_tag'], model['maximum_gap_bp']
        require((tag, gap) not in observed_models, 'duplicate model scenario')
        observed_models.add((tag, gap))
        folder = root / (tag + '-gap-' + str(gap))
        check_file(folder / 'summary.json', model['summary'])
        detail = json.loads((folder / 'summary.json').read_text())
        require(detail['input_sha256'] == source['outputs']['qualified-markers.vcf.gz']['sha256'], 'model input differs')
        require(detail['counts'] == model['counts'], 'model counts differ')
        require(detail['population_af_tag'] == tag and detail['max_marker_gap_bp'] == gap, 'model scenario differs')
        expected_gaps = [{'chrom': c, 'left_marker_1based': a, 'right_marker_1based': b, 'marker_distance_bp': b-a}
                         for c, ps in positions.items() for a, b in zip(ps, ps[1:]) if b-a > gap]
        require(detail['gap_boundaries'] == expected_gaps and len(expected_gaps) == model['gap_boundaries'],
                'large marker gap accounting differs')
        covered, regions, roh_labels = set(), [], Counter()
        for i, segment in enumerate(detail['segments']):
            segment_dir = folder / ('segment-%04d' % i)
            for name, key in [('input.vcf', 'input_sha256'), ('roh.txt', 'output_sha256'), ('stderr.txt', 'stderr_sha256')]:
                require(sha256(segment_dir / name) == segment[key], 'native model artifact changed')
                native_files += 1
            with pysam.VariantFile(str(segment_dir / 'input.vcf')) as reader:
                expected_sites = [(r.contig, r.pos) for r in reader]
            require(len(expected_sites) == segment['markers'], 'segment input count differs')
            actual_sites, states, segment_regions = [], Counter(), []
            for line in (segment_dir / 'roh.txt').read_text().splitlines():
                fields = line.split('\t')
                if fields[0] == 'ST':
                    require(len(fields) == 6 and fields[1] == sample and fields[4] in ('0', '1'), 'invalid native site')
                    key = (fields[2], int(fields[3]))
                    require(key in markers and key not in covered, 'unqualified or repeated native marker')
                    actual_sites.append(key); covered.add(key); states[fields[4]] += 1
                    if fields[4] == '1':
                        roh_labels[markers[key]] += 1
                elif fields[0] == 'RG':
                    require(len(fields) == 8 and fields[1] == sample, 'invalid native region')
                    item = dict(chrom=fields[2], start_1based=int(fields[3]), end_inclusive=int(fields[4]),
                                length_bp=int(fields[5]), markers=int(fields[6]), mean_model_phred=float(fields[7]))
                    require(item['length_bp'] == item['end_inclusive'] - item['start_1based'] + 1, 'invalid region span')
                    require(segment['start_1based'] <= item['start_1based'] <= item['end_inclusive'] <= segment['end_inclusive'],
                            'region crosses independent segment boundary')
                    segment_regions.append(item)
            require(actual_sites == expected_sites and dict(states) == segment['site_states'], 'native site census differs')
            require(segment_regions == segment['regions'], 'native regions differ from saved summary')
            regions.extend(segment_regions)
        skipped = set()
        for segment in detail['skipped_segments']:
            selected = [(segment['chrom'], p) for p in positions[segment['chrom']]
                        if segment['start_1based'] <= p <= segment['end_inclusive']]
            require(len(selected) == segment['markers'] < detail['minimum_segment_markers'], 'invalid skipped segment')
            require(not set(selected) & (covered | skipped), 'overlapping skipped segment')
            skipped.update(selected)
        require(covered | skipped == set(markers), 'model does not account for all qualified markers')
        require(len(covered) == detail['counts']['analyzed_markers'] and len(skipped) == detail['counts']['skipped_markers'],
                'analyzed/skipped census differs')
        require(regions == model['regions'], 'top-level region list differs')
        require(len(detail['segments']) == model['segments'], 'segment total differs')
        models.append(dict(population_af_tag=tag, maximum_gap_bp=gap, counts=detail['counts'],
            independent_segments=len(detail['segments']), gap_boundaries=len(expected_gaps),
            raw_region_count=len(regions), raw_region_span_bp=sum(r['length_bp'] for r in regions),
            maximum_region_bp=max((r['length_bp'] for r in regions), default=0),
            posthoc_descriptive_length_bins=[dict(minimum_bp=n, count=sum(r['length_bp'] >= n for r in regions),
                span_bp=sum(r['length_bp'] for r in regions if r['length_bp'] >= n)) for n in (1000000, 3000000, 5000000, 10000000)],
            genotype_labels_in_native_autozygous_state=dict(roh_labels),
            longest_five_regions=sorted(regions, key=lambda r: r['length_bp'], reverse=True)[:5],
            selected_region_overlaps={name: [r for r in regions if r['chrom'] == c and r['start_1based'] <= b and r['end_inclusive'] >= a]
                                     for name, (c, a, b) in review_regions.items()}))
    expected_models = {(tag, source['config']['primary_maximum_gap_bp']) for tag in source['config']['population_profiles']}
    expected_models.update(('POP_AF', gap) for gap in source['config']['pooled_gap_sensitivity_bp'])
    require(observed_models == expected_models, 'declared model scenarios missing or extra')
    return dict(reviewed_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        analysis_completed_at=source['completed_at'], analysis_root=str(root), native_summary_sha256=actual_sha,
        review_script_sha256=sha256(Path(__file__)),
        verification=dict(top_level_outputs=len(source['outputs']), native_segment_artifacts=native_files,
            raw_annotated_markers=sum(requested.values()), qualified_markers=sum(qualified.values()), scenarios=len(models),
            scope='Output identity and internal evidence consistency; not an independent model, decoder or biological validation.'),
        read_extraction=extraction, marker_counts=counts, raw_genotypes=dict(raw_categories),
        raw_missing_GQ_by_genotype=dict(missing_gq), qualified_genotypes=dict(qualified_categories),
        qualified_missing_GQ_by_genotype=dict(qualified_missing_gq), exclusion_reasons_nonexclusive=dict(reasons),
        markers_by_chromosome={c: dict(requested=requested[c], qualified=qualified[c], excluded=excluded[c],
             first_qualified_1based=ps[0], last_qualified_1based=ps[-1]) for c, ps in positions.items()},
        reviewed_intervals_1based_inclusive=review_regions, models=models,
        limitations=source['limitations'] + [
            'All native regions retained; no additional region-quality threshold. A minimum of 20 markers governs input segments, not reported regions.',
            'Length bins and selected-region overlap are post-run descriptive reviews, not preregistered diagnostic thresholds.',
            'Native interval spans include sequence between markers. They are not measured base-level autozygosity or an inbreeding coefficient.',
            'A native autozygous state may include heterozygous maximum-likelihood calls; model scores are not clinical probabilities.'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--analysis', type=Path, required=True)
    parser.add_argument('--expected-summary-sha256', required=True)
    parser.add_argument('--review-region', action='append', default=[], help='NAME=CHROM:START-END, 1-based inclusive')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require(not args.output.exists(), 'output already exists')
    regions = {}
    for value in args.review_region:
        match = re.fullmatch(r'([A-Za-z0-9_]+)=([0-9]+):([0-9]+)-([0-9]+)', value)
        require(match is not None, 'invalid review region')
        name, chrom, first, last = match.groups()
        require(name not in regions and 1 <= int(chrom) <= 22 and 1 <= int(first) <= int(last), 'invalid region bounds')
        regions[name] = (chrom, int(first), int(last))
    result = review(args.analysis.resolve(), args.expected_summary_sha256, regions)
    with args.output.open('x') as stream:
        stream.write(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps(result['verification']))


if __name__ == '__main__':
    main()
