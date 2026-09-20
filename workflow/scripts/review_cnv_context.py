#!/usr/bin/env python3
"""Add reproducible population, dosage and common-marker context to every CNV.

Full call/marker/population details remain outside the source repository. Read
balance, overlap and dosage curation are descriptive evidence, not a diagnosis.
"""
import argparse
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
import csv
import datetime
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import statistics
import time


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def fingerprint(path):
    return {'path': str(Path(path).resolve()), 'bytes': Path(path).stat().st_size, 'sha256': sha(path)}


def verify(path, record):
    require(Path(path).stat().st_size == record['bytes'] and sha(path) == record['sha256'], 'Input identity differs: ' + str(path))


def save(path, value):
    tmp = Path(str(path) + '.tmp')
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    os.replace(tmp, path)


def overlap(a, b, c, d):
    return max(0, min(b, d) - max(a, c) + 1)


class IntervalIndex:
    """Inclusive intervals; prefix maximum ends retain nested/long overlaps."""
    def __init__(self, rows):
        self.rows = sorted(rows, key=lambda r: r['start_1based'])
        self.starts = [r['start_1based'] for r in self.rows]
        self.maxends = []
        last = 0
        for row in self.rows:
            last = max(last, row['end_inclusive'])
            self.maxends.append(last)

    def fetch(self, start, end):
        if end < start:
            return []
        return [r for r in self.rows[bisect_left(self.maxends, start):bisect_right(self.starts, end)]
                if r['end_inclusive'] >= start]


def distribution(values):
    values = sorted(values)
    if not values:
        return {'n': 0, 'minimum': None, 'median': None, 'maximum': None, 'mean': None}
    return {'n': len(values), 'minimum': values[0], 'median': statistics.median(values),
            'maximum': values[-1], 'mean': statistics.mean(values)}


def gt_category(gt):
    if not gt or any(x is None for x in gt):
        return 'missing'
    require(len(gt) == 2 and all(x in (0, 1) for x in gt), 'Unexpected panel genotype')
    return ('homref', 'heterozygous', 'homalt')[sum(gt)]


def allele_fraction(ad):
    if ad is None or len(ad) != 2 or any(type(x) is not int or x < 0 for x in ad) or not sum(ad):
        return None
    return ad[1] / sum(ad)


def marker_summary(rows, minimum_ad):
    result = {'markers': len(rows), 'usable_for_previous_ROH_model': sum(r['usable'] for r in rows),
        'native_diploid_genotypes': dict(Counter(r['category'] for r in rows)),
        'exclusion_reasons_nonexclusive': dict(Counter(x for r in rows for x in r['excluded'])),
        'DP': distribution([r['DP'] for r in rows if r['DP'] is not None]),
        'missing_AD': sum(r['AD'] is None for r in rows),
        'missing_PL': sum(r['PL'] is None for r in rows),
        'missing_GQ': sum(r['GQ'] is None for r in rows),
        'observed_marker_span_bp': rows[-1]['pos'] - rows[0]['pos'] + 1 if rows else None}
    for label, subset in [('all_markers', rows), ('ROH_usable', [r for r in rows if r['usable']]),
            ('native_heterozygous', [r for r in rows if r['category'] == 'heterozygous']),
            ('native_heterozygous_ROH_usable', [r for r in rows if r['category'] == 'heterozygous' and r['usable']])]:
        fractions = [allele_fraction(r['AD']) for r in subset
                     if allele_fraction(r['AD']) is not None and sum(r['AD']) >= minimum_ad]
        hist = [0] * 10
        for f in fractions:
            hist[min(int(f * 10), 9)] += 1
        result[label + '_AD_at_least_' + str(minimum_ad)] = {
            'ALT_fraction': distribution(fractions), 'absolute_distance_from_half': distribution([abs(f - .5) for f in fractions]),
            'ALT_fraction_tenth_bins_left_closed_last_includes_one': hist}
    return result


def read_markers(root, expected_sha):
    import pysam
    require(sha(root / 'summary.json') == expected_sha, 'ROH summary differs')
    native = json.loads((root / 'summary.json').read_text())
    state = json.loads((root / 'status.json').read_text())
    require(state['status'] == 'complete' and state['summary_sha256'] == expected_sha, 'ROH terminal status differs')
    require(native['config']['scope'] == 'verified_four_lane_autosomal_marker_genotyping'
            and native['analysis_status'] == 'completed_requested_models', 'Unexpected marker analysis scope')
    for name in ('annotated-all-markers.vcf.gz', 'annotated-all-markers.vcf.gz.tbi'):
        verify(root / name, native['outputs'][name])
    sample = native['config']['sample']
    markers, lengths = defaultdict(list), {}
    total, qualified, excluded, gq_missing, genotypes, reasons = 0, 0, 0, 0, Counter(), Counter()
    with pysam.VariantFile(str(root / 'annotated-all-markers.vcf.gz')) as reader:
        require(list(reader.header.samples) == [sample], 'Panel sample differs')
        lengths = {c: reader.header.contigs[c].length for c in native['config']['chromosomes']}
        for record in reader:
            require(record.contig in lengths and len(record.alleles) == 2 and all(len(a) == 1 for a in record.alleles), 'Unexpected panel allele')
            s = record.samples[sample]; cat = gt_category(s.get('GT'))
            row = {'chrom': record.contig, 'pos': record.pos, 'REF': record.ref, 'ALT': record.alts[0],
                   'GT': s.get('GT'), 'AD': s.get('AD'), 'DP': s.get('DP'), 'PL': s.get('PL'), 'GQ': s.get('GQ'),
                   'phased': s.phased, 'category': cat, 'usable': bool(record.info.get('MVA_USABLE', False)),
                   'excluded': record.info.get('MVA_EXCLUDE', ())}
            require(row['usable'] != bool(row['excluded']), 'Marker usability/exclusion mismatch')
            require(not markers[record.contig] or row['pos'] > markers[record.contig][-1]['pos'], 'Repeated/unordered marker')
            require(1 <= row['pos'] <= lengths[record.contig], 'Marker outside chromosome')
            markers[record.contig].append(row)
            total += 1; qualified += row['usable']; excluded += not row['usable']
            gq_missing += row['GQ'] is None; genotypes[cat] += 1; reasons.update(row['excluded'])
    counts = native['marker_qualification']['counts']
    require((total, qualified, excluded, gq_missing) == (counts['requested_and_returned_markers'], counts['qualified_markers'],
            counts['excluded_markers'], counts['raw_GQ_missing']), 'Panel census differs')
    require(all(n == counts['raw_genotype_' + c] for c, n in genotypes.items()), 'Genotype census differs')
    require(dict(reasons) == native['marker_qualification']['exclusion_reason_counts_nonexclusive'], 'Marker exclusion census differs')
    for chrom, rows in markers.items():
        require(len(rows) == native['marker_qualification']['per_chromosome'][chrom]['requested_and_returned_markers'], 'Chromosome count differs')
    return markers, lengths, native, {'markers': total, 'usable': qualified, 'excluded': excluded,
        'genotypes': dict(genotypes), 'exclusion_reasons_nonexclusive': dict(reasons), 'missing_GQ': gq_missing}


def read_clingen(path, kind):
    lines = path.read_text().splitlines()
    header = next(i for i, line in enumerate(lines) if line.startswith('#') and '\tGenomic Location\t' in line)
    require(any('GRCh38' in line for line in lines[:header]), 'ClinGen assembly differs')
    rows, unmapped = [], []
    id_field = 'Gene ID' if kind == 'gene' else 'ISCA ID'
    name_field = 'Gene Symbol' if kind == 'gene' else 'ISCA Region Name'
    seen = set()
    for raw in csv.DictReader([lines[header][1:]] + lines[header + 1:], delimiter='\t'):
        require(None not in raw and all(v is not None for v in raw.values()), 'ClinGen column count differs')
        # Sex chromosome pseudoautosomal genes may have more than one locus.
        key = (raw[id_field], raw['Genomic Location'])
        require(key not in seen, 'Repeated ClinGen identity/location'); seen.add(key)
        base = {'kind': kind, 'id': raw[id_field], 'name': raw[name_field], 'curation': raw}
        match = re.fullmatch(r'chr([0-9XY]+):(\d+)-(\d+)', raw['Genomic Location'])
        if match is None:
            unmapped.append(base); continue
        chrom, start, end = match.groups()
        require(1 <= int(start) <= int(end), 'ClinGen interval reversed')
        rows.append(base | {'chrom': chrom, 'start_1based': int(start), 'end_inclusive': int(end)})
    return rows, unmapped


def annotate_dosage(call, matches):
    field = 'Haploinsufficiency' if call['type'] == 'deletion' else 'Triplosensitivity'
    result = []
    for row in matches:
        n = overlap(call['start_1based'], call['end_inclusive'], row['start_1based'], row['end_inclusive'])
        item = dict(row, overlap_bp=n, fraction_of_curated_span=n / (row['end_inclusive'] - row['start_1based'] + 1),
                    fraction_of_call=n / call['size_bp'], relevant_score_field=field + ' Score',
                    relevant_score=row['curation'][field + ' Score'], relevant_description=row['curation'][field + ' Description'])
        result.append(item)
    return result


def population_relation(call, record):
    """Nominal span follows official gnomAD BED: VCF POS..END, inclusive.

    It is not an assertion of exact affected sequence or breakpoint resolution.
    """
    a, b = record['affected_start_1based'], record['affected_end_inclusive']
    n = overlap(call['start_1based'], call['end_inclusive'], a, b)
    require(n > 0, 'Population relation requested without span overlap')
    same = record['SVTYPE'] == ('DEL' if call['type'] == 'deletion' else 'DUP')
    simple = same and len(record['ALT']) == 1 and not record['INFO'].get('MULTIALLELIC', False) and record['SVTYPE'] in ('DEL', 'DUP')
    passed = record['FILTER'] == ['PASS']
    min_reciprocal = min(n / call['size_bp'], n / (b - a + 1))
    return {'id': record['ID'], 'overlap_bp': n, 'fraction_of_call': n / call['size_bp'],
        'fraction_of_population_span': n / (b - a + 1), 'minimum_reciprocal_overlap': min_reciprocal,
        'same_simple_type': simple, 'population_PASS': passed,
        'start_distance_bp': abs(call['start_1based'] - a), 'end_distance_bp': abs(call['end_inclusive'] - b),
        'both_endpoints_within_one_bin': max(abs(call['start_1based'] - a), abs(call['end_inclusive'] - b)) <= call['bin_size']}


def biallelic_frequency(record):
    info = record['INFO']
    if record['SVTYPE'] not in ('DEL', 'DUP') or len(record['ALT']) != 1 or info.get('MULTIALLELIC', False):
        return None
    ac, an, af = info.get('AC'), info.get('AN'), info.get('AF')
    if ac is None or an is None or af is None or len(ac) != 1 or len(af) != 1 or ac[0] is None or af[0] is None or an == 0:
        return None
    require(type(ac[0]) is int and type(an) is int and 0 <= ac[0] <= an, 'Invalid population allele counts')
    exact = ac[0] / an
    require(math.isfinite(af[0]) and abs(exact - af[0]) < 1e-6, 'Population AF disagrees with AC/AN beyond rounding')
    return {'AC': ac[0], 'AN': an, 'reported_AF': af[0], 'AC_divided_by_AN': exact,
            'scope': 'This population record only; matching its interval does not assign this AF to the case allele.'}


def population_context(vcf, calls, out):
    import pysam
    indexes = {c: IntervalIndex([r for r in calls if r['chrom'] == c]) for c in map(str, range(1, 23))}
    by_call, records = defaultdict(list), {}
    totals, types, filters, omitted = Counter(), Counter(), Counter(), Counter()
    keep_prefixes = ('AC', 'AN', 'AF', 'CN_NUMBER', 'CN_COUNT', 'CN_FREQ', 'N_HOMREF', 'N_HET', 'N_HOMALT')
    keep = {'SVTYPE', 'SVLEN', 'CHR2', 'END2', 'CPX_TYPE', 'CPX_INTERVALS', 'MULTIALLELIC', 'ALGORITHMS', 'EVIDENCE', 'NCR'}
    with pysam.VariantFile(str(vcf)) as reader:
        require(len(reader.header.samples) == 0, 'Expected public sites-only VCF')
        require(all(k in reader.header.info for k in ('SVTYPE', 'END', 'AC', 'AN', 'AF', 'CN_NUMBER', 'CN_COUNT', 'CN_FREQ')), 'Population header differs')
        (out / 'population-header.txt').write_text(str(reader.header))
        for row in reader:
            totals['records'] += 1; svtype = row.info.get('SVTYPE'); types[svtype] += 1
            filters.update(list(row.filter) or ['.'])
            chrom = row.contig.removeprefix('chr')
            if chrom not in indexes:
                omitted['nonautosomal'] += 1; continue
            # A nominal CPX/INS/INV/BND span is never treated as a same-type CNV.
            if row.stop < row.pos or row.info.get('CHR2', row.contig) != row.contig:
                omitted['nonpositive_or_interchromosomal_nominal_span'] += 1; continue
            a, b = row.pos, row.stop
            candidates = indexes[chrom].fetch(a, b)
            if not candidates:
                omitted['no_raw_call_span_overlap'] += 1; continue
            require(row.id and row.id not in records, 'Missing or duplicate population variant ID')
            info = {k: v for k, v in row.info.items() if k in keep or any(k == p or k.startswith(p + '_') for p in keep_prefixes)
                    or reader.header.info[k].type == 'Flag'}
            record = {'ID': row.id, 'chrom': chrom, 'POS': row.pos, 'END': row.stop, 'REF': row.ref,
                'ALT': list(row.alts or ()), 'QUAL': row.qual, 'FILTER': list(row.filter), 'SVTYPE': svtype,
                'affected_start_1based': a, 'affected_end_inclusive': b, 'INFO': info}
            record['biallelic_frequency'] = biallelic_frequency(record)
            records[row.id] = record
            for call in candidates:
                by_call[call['call_id']].append(population_relation(call, record))
            if totals['records'] % 100000 == 0:
                print(json.dumps({'population_records_seen': totals['records'], 'overlapping_records': len(records)}), flush=True)
    totals['overlapping_records'] = len(records)
    totals['call_population_relations'] = sum(map(len, by_call.values()))
    require(totals['records'] == sum(omitted.values()) + len(records), 'Population row accounting differs')
    with gzip.open(out / 'population-overlapping-records.json.gz', 'wt') as f:
        json.dump({'records': records}, f, allow_nan=False)
    return by_call, records, {'totals': dict(totals), 'SVTYPE_counts': dict(types),
            'FILTER_counts_nonexclusive': dict(filters), 'omitted_from_nominal_span_context': dict(omitted)}


def summarize_population(matches, records, config):
    same_pass = [m for m in matches if m['same_simple_type'] and m['population_PASS']]
    bins = {str(threshold): [m for m in same_pass if m['minimum_reciprocal_overlap'] >= threshold]
            for threshold in config['population_reciprocal_overlap_thresholds']}
    selected = sorted(same_pass, key=lambda m: (-m['minimum_reciprocal_overlap'], m['start_distance_bp'] + m['end_distance_bp'], m['id']))
    example_ids = [m['id'] for m in selected[:config['maximum_population_examples_per_call']]]
    def largest_frequency(group):
        values = [records[m['id']]['biallelic_frequency']['AC_divided_by_AN'] for m in group if records[m['id']]['biallelic_frequency'] is not None]
        return max(values, default=None)
    return {'any_nominal_span_overlap': len(matches), 'same_simple_type_PASS': len(same_pass),
        'same_simple_type_PASS_reciprocal_overlap': {k: {'count': len(v), 'maximum_population_record_AF': largest_frequency(v)} for k, v in bins.items()},
        'same_simple_type_PASS_both_endpoints_within_one_bin': sum(m['both_endpoints_within_one_bin'] for m in same_pass),
        'best_overlap_examples': [m | {'population_record': records[m['id']]} for m in selected[:config['maximum_population_examples_per_call']]],
        'other_SVTYPE_or_multiallelic': sum(not m['same_simple_type'] for m in matches),
        'non_PASS': sum(not m['population_PASS'] for m in matches),
        'example_ids': example_ids}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('cnv-review', 'roh-analysis', 'resources', 'config', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    for name in ('cnv-review', 'roh-summary', 'resources'):
        parser.add_argument('--expected-' + name + '-sha256', required=True)
    args = parser.parse_args(); tick = time.monotonic(); os.nice(10)
    require(not args.output.exists(), 'Existing review preserved')
    require(sha(args.cnv_review / 'summary.json') == args.expected_cnv_review_sha256, 'CNV review differs')
    review = json.loads((args.cnv_review / 'summary.json').read_text())
    verify(args.cnv_review / 'annotated-calls.json', review['output'])
    calls = json.loads((args.cnv_review / 'annotated-calls.json').read_text())['calls']
    require(len(calls) == review['native_call_records_equal'] == len({c['call_id'] for c in calls}), 'CNV call census differs')
    require(sha(args.resources / 'summary.json') == args.expected_resources_sha256, 'Public resource receipt differs')
    resources = json.loads((args.resources / 'summary.json').read_text())
    require(resources['status'] == 'complete', 'Public sources incomplete')
    for record in resources['files']:
        verify(args.resources / record['filename'], record)
    config = json.loads(args.config.read_text())
    markers, lengths, native_roh, marker_census = read_markers(args.roh_analysis, args.expected_roh_summary_sha256)
    cnv_native_path = Path(review['inputs']['native_summary']['path'])
    verify(cnv_native_path, review['inputs']['native_summary'])
    cnv_native = json.loads(cnv_native_path.read_text())
    for name in ('alignment', 'alignment_index', 'reference', 'reference_index'):
        require(cnv_native['inputs'][name]['sha256'] == native_roh['inputs'][name]['sha256'],
                'CNV and marker input identities differ: ' + name)
    genes, unmapped_genes = read_clingen(args.resources / 'ClinGen_gene_curation_list_GRCh38.tsv', 'gene')
    regions, unmapped_regions = read_clingen(args.resources / 'ClinGen_region_curation_list_GRCh38.tsv', 'region')
    dosage_index = {c: IntervalIndex([r for r in genes + regions if r['chrom'] == c]) for c in lengths}
    primary_index = {c: IntervalIndex([r for r in calls if r['chrom'] == c and not r['triage_failures']['primary']]) for c in lengths}
    positions = {c: [r['pos'] for r in rows] for c, rows in markers.items()}
    def select(c, a, b):
        return markers[c][bisect_left(positions[c], a):bisect_right(positions[c], b)] if a <= b else []
    args.output.mkdir(parents=True)
    population, records, population_census = population_context(args.resources / 'gnomad.v4.1.sv.sites.vcf.gz', calls, args.output)
    all_context, cds_context, cross_context = [], [], []
    for call in calls:
        c, a, b = call['chrom'], call['start_1based'], call['end_inclusive']
        item = {k: call[k] for k in ('call_id', 'chrom', 'start_1based', 'end_inclusive', 'size_bp', 'type', 'bin_size',
                'relative_depth', 'MQ_zero_fraction', 'reference_gap_fraction', 'within_reference_bounds', 'triage_failures', 'gene_overlaps', 'other_resolution_matches')}
        item['dosage_context'] = annotate_dosage(call, dosage_index[c].fetch(a, b))
        item['population_context'] = summarize_population(population[call['call_id']], records, config)
        item['all_population_relations'] = population[call['call_id']]
        interior = select(c, a, b)
        item['markers_inside_nominal_call'] = marker_summary(interior, config['marker_minimum_allelic_depth'])
        item['markers_outside_reference_portion'] = 'No markers exist outside reference; out-of-bounds raw calls remain flagged.' if not call['within_reference_bounds'] else None
        width = max(config['flank_minimum_bp'], min(config['flank_maximum_bp'], call['size_bp']))
        item['flanks'] = {}
        for label, start, end in [('left', max(1, a-width), min(lengths[c], a-1)), ('right', max(1, b+1), min(lengths[c], b+width))]:
            rows = select(c, start, end)
            clean = [r for r in rows if not primary_index[c].fetch(r['pos'], r['pos'])]
            item['flanks'][label] = {'start_1based': start, 'end_inclusive': end, 'empty_interval': end < start,
                'all': marker_summary(rows, config['marker_minimum_allelic_depth']),
                'outside_all_primary_review_CNVs': marker_summary(clean, config['marker_minimum_allelic_depth'])}
        all_context.append(item)
        # Bounded receipt keeps every primary coding overlap and every primary
        # cross-resolution call, including zero-marker and zero-match cases.
        bounded = {k: v for k, v in item.items() if k != 'all_population_relations'}
        if not call['triage_failures']['primary'] and any(g['union_CDS_overlap_bp'] for g in call['gene_overlaps']):
            cds_context.append(bounded)
        if not call['triage_failures']['primary'] and any(m['passes_primary'] for m in call['other_resolution_matches']):
            cross_context.append(bounded)
    with gzip.open(args.output / 'all-call-context.json.gz', 'wt') as f:
        json.dump({'config': config, 'calls': all_context}, f, allow_nan=False)
    # Native records are retained for every marker, including missing/excluded;
    # do not turn an excluded marker into a homozygous-reference call.
    with gzip.open(args.output / 'all-marker-context.jsonl.gz', 'wt') as f:
        for c in native_roh['config']['chromosomes']:
            for row in markers[c]:
                f.write(json.dumps(row, allow_nan=False) + '\n')
    def census(rows):
        return {'calls': len(rows), 'no_marker_inside': sum(r['markers_inside_nominal_call']['markers'] == 0 for r in rows),
            'no_native_heterozygous_marker_inside': sum(r['markers_inside_nominal_call']['native_diploid_genotypes'].get('heterozygous', 0) == 0 for r in rows),
            'any_ClinGen_dosage_overlap': sum(bool(r['dosage_context']) for r in rows),
            'any_same_type_PASS_population_RO50': sum(r['population_context']['same_simple_type_PASS_reciprocal_overlap']['0.5']['count'] > 0 for r in rows),
            'any_same_type_PASS_population_RO80': sum(r['population_context']['same_simple_type_PASS_reciprocal_overlap']['0.8']['count'] > 0 for r in rows)}
    primary = [r for r in all_context if not r['triage_failures']['primary']]
    summary = {'completed_at': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'elapsed_seconds': time.monotonic()-tick,
        'scope': config['scope'], 'config': config, 'inputs': {
            'cnv_review_summary': fingerprint(args.cnv_review / 'summary.json'), 'annotated_cnv_calls': review['output'],
            'roh_summary': fingerprint(args.roh_analysis / 'summary.json'), 'annotated_all_markers': native_roh['outputs']['annotated-all-markers.vcf.gz'],
            'resource_summary': fingerprint(args.resources / 'summary.json'), 'configuration': fingerprint(args.config), 'script': fingerprint(Path(__file__))},
        'marker_census': marker_census, 'population_census': population_census,
        'dosage_source_census': {'mapped_genes': len(genes), 'mapped_regions': len(regions), 'unmapped_genes': unmapped_genes, 'unmapped_regions': unmapped_regions},
        'call_census': {'all': census(all_context), 'primary': census(primary), 'primary_CDS': census(cds_context), 'primary_cross_resolution': census(cross_context)},
        'primary_CDS_call_context': cds_context, 'primary_cross_resolution_call_context': cross_context,
        'outputs': {name: fingerprint(args.output / name) for name in ('population-header.txt', 'population-overlapping-records.json.gz', 'all-call-context.json.gz', 'all-marker-context.jsonl.gz')},
        'limitations': [f'All {len(calls)} raw call records and all {marker_census["markers"]} panel markers retained; interval overlap across resolutions reuses the same reads.',
            'Population nominal span follows the official gnomAD BED convention: VCF POS..END, inclusive. This is not an exact affected-sequence or breakpoint assertion. Other nominal SV spans are context only; complex component alleles are not reconstructed.',
            'PASS and reciprocal overlap are descriptive. Even both endpoints within one native bin does not establish the same allele, topology, genotype, benignity or a clinical population frequency for the case.',
            'Multiallelic CNV CN_COUNT/CN_FREQ are sample copy-state distributions, not biallelic AC/AN/AF. Missing frequencies are not zero.',
            'ClinGen scores remain raw strings: no evidence, recessive, unlikely, not evaluated and missing are distinct. A partial overlap does not inherit the full-gene or recurrent-region classification.',
            'Common SNP ascertainment, sparse spacing, mapping bias, diploid genotype priors and depth filtering limit read-balance interpretation. No calibrated copy-number or mosaic estimate is reported.',
            'Flanking markers outside primary calls are not proven copy-neutral. Absence of markers or heterozygosity does not prove deletion or autozygosity.',
            'This does not independently validate breakpoints, phasing, parental origin, diagnosis, clinical sensitivity or drug response.']}
    save(args.output / 'summary.json', summary)
    print(json.dumps({'status': 'complete', 'call_census': summary['call_census'], 'population_census': population_census,
                      'elapsed_seconds': summary['elapsed_seconds']}), flush=True)


if __name__ == '__main__':
    main()
