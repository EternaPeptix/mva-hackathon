#!/usr/bin/env python3
"""Describe unchanged saved public truth decisions; never call, align or rematch."""
import argparse
from collections import Counter, defaultdict
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import statistics
import sys

OUTCOMES = ('matched_by_both', 'DELLY_only', 'SvABA_only', 'missed_by_both')
CLASSES = ('DEL', 'INS')
BINS = ('50-99', '100-499', '500-999', '1000-9999', '10000-50000')
ANNOTATIONS = {
    'TRF': 'Flag', 'TRFdiff': 'Float', 'TRFrepeat': 'String', 'TRFovl': 'Float',
    'TRFstart': 'Integer', 'TRFend': 'Integer', 'TRFperiod': 'Integer',
    'TRFcopies': 'Float', 'TRFscore': 'Integer', 'TRFentropy': 'Float', 'TRFsim': 'Float',
    'RM_score': 'Integer', 'RM_repeat': 'String', 'RM_clsfam': 'String',
    'LCR': 'Float', 'REMAP': 'String',
}
ROLES = ('delly_final_truth', 'svaba_final_truth', 'original_truth_chr20_no_star',
         'official_chr20_benchmark', 'delly_summary', 'delly_status', 'svaba_summary',
         'svaba_status', 'delly_variant_summary', 'svaba_variant_summary',
         'comparison_config', 'prior_final_audit')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def no_duplicates(pairs):
    out = {}
    for key, value in pairs:
        require(key not in out, 'duplicate JSON key')
        out[key] = value
    return out


def parse_json(text):
    return json.loads(text, object_pairs_hook=no_duplicates,
                      parse_constant=lambda x: (_ for _ in ()).throw(ValueError('nonfinite JSON')))


def fingerprint(path):
    path = Path(path)
    h = hashlib.sha256()
    with path.open('rb') as reader:
        for block in iter(lambda: reader.read(1024 * 1024), b''):
            h.update(block)
    return {'path': str(path.resolve()), 'bytes': path.stat().st_size, 'sha256': h.hexdigest()}


def bind(path, expected):
    actual = fingerprint(path)
    require(type(expected.get('bytes')) is int and expected['bytes'] == actual['bytes'], 'input byte count differs')
    require(expected.get('sha256') == actual['sha256'], 'input SHA differs')
    require(str(Path(expected['path']).resolve()) == actual['path'], 'input path differs')
    return actual


def integer(value):
    require(isinstance(value, str) and re.fullmatch(r'-?[0-9]+', value) is not None, 'malformed integer')
    return int(value)


def info_values(text):
    out = {}
    if text == '.':
        return out
    for token in text.split(';'):
        key, sep, value = token.partition('=')
        require(key and key not in out, 'empty or duplicate INFO key')
        out[key] = value if sep else None
    return out


def header_definition(line):
    values = {}
    for key in ('ID', 'Number', 'Type'):
        found = re.findall(r'(?:<|,)' + key + r'=([^,>]+)', line)
        require(len(found) == 1, 'malformed VCF field declaration')
        values[key] = found[0]
    return values


def stream_vcf(path):
    """Yield a header then records with original one-based data-row ordinals."""
    headers, infos, formats = [], {}, {}
    with gzip.open(path, 'rt') as reader:
        for line in reader:
            require(len(line) <= 2 * 1024 * 1024, 'VCF line exceeds finite parser bound')
            line = line.rstrip('\r\n')
            require(line.startswith('#'), 'VCF data before column header')
            headers.append(line)
            for prefix, definitions in (('##INFO=', infos), ('##FORMAT=', formats)):
                if line.startswith(prefix):
                    definition = header_definition(line)
                    key = definition['ID']
                    require(key not in definitions, 'duplicate VCF field declaration')
                    definitions[key] = definition
            if line.startswith('#CHROM\t'):
                cols = line.split('\t')
                require(cols == ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'HG002'],
                        'unexpected public VCF sample/columns')
                break
        else:
            raise ValueError('VCF column header absent')
        yield {'headers': headers, 'info_definitions': infos, 'format_definitions': formats}
        for ordinal, raw in enumerate(reader, 1):
            require(ordinal <= 250000 and len(raw) <= 2 * 1024 * 1024, 'VCF finite bound exceeded')
            raw = raw.rstrip('\r\n')
            fields = raw.split('\t')
            require(len(fields) == 10 and not raw.startswith('#'), 'malformed VCF row')
            pos = integer(fields[1])
            require(pos >= 1 and fields[0] == 'chr20' and fields[3] and fields[4], 'unexpected chromosome/allele/position')
            fmt, sample = fields[8].split(':'), fields[9].split(':')
            require(len(fmt) == len(set(fmt)) and len(sample) <= len(fmt), 'malformed FORMAT/sample values')
            yield {'ordinal': ordinal, 'chrom': fields[0], 'pos': pos, 'id': fields[2],
                   'ref': fields[3], 'alts': fields[4].split(','), 'qual': fields[5], 'filter': fields[6],
                   'info': info_values(fields[7]), 'sample': dict(zip(fmt, sample)), 'raw_record': raw}


def key(row):
    return row['chrom'], row['pos'], row['id'], row['ref'], tuple(row['alts'])


def genotype(row):
    raw = row['sample'].get('GT')
    require(raw is not None and re.fullmatch(r'(?:\.|[0-9]+)(?:[|/](?:\.|[0-9]+))*', raw), 'malformed/absent GT')
    require(not ('|' in raw and '/' in raw), 'mixed GT separators')
    alleles = tuple(None if x == '.' else int(x) for x in re.split(r'[|/]', raw))
    require(all(x is None or x <= len(row['alts']) for x in alleles), 'GT allele outside original ALT list')
    # Match the old audit's allele-content policy; raw order/phase stays in the ledger.
    return tuple(sorted(alleles, key=lambda x: -1 if x is None else x))


def annotations(row, header):
    out = {}
    for name, kind in ANNOTATIONS.items():
        definition = header['info_definitions'].get(name)
        require(definition and definition['Type'] == kind and definition['Number'] == ('0' if kind == 'Flag' else '1'),
                'missing or changed annotation schema: ' + name)
        present = name in row['info']
        raw = row['info'].get(name)
        value = None
        if kind == 'Flag':
            require(not present or raw is None, 'flag must use bare VCF encoding')
            value = True if present else None
        elif present and raw != '.':
            require(raw is not None and raw != '' and ',' not in raw, 'invalid scalar annotation')
            if kind == 'Integer':
                value = integer(raw)
            elif kind == 'Float':
                value = float(raw)
                require(math.isfinite(value), 'nonfinite annotation')
            else:
                value = raw
        out[name] = {'present': present, 'value': value, 'raw': raw}
    return out


def size_bin(size):
    require(type(size) is int and 50 <= size <= 50000, 'truth outside frozen size range')
    for edge, label in zip((100, 500, 1000, 10000, 50001), BINS):
        if size < edge:
            return label
    raise ValueError('unreachable size')


def truth_shape(row):
    require(len(row['alts']) == 1 and all(re.fullmatch('[ACGTNacgtn]+', x) for x in [row['ref'], *row['alts']]),
            'unsupported selected truth allele shape')
    kind = row['info'].get('SVTYPE')
    require(kind in CLASSES, 'unsupported selected truth class')
    size = abs(integer(row['info'].get('SVLEN')))
    return kind, size, size_bin(size)


def read_final(path):
    stream = stream_vcf(path)
    header = next(stream)
    bd = header['format_definitions'].get('BD')
    require(bd and bd['Number'] == '1' and bd['Type'] == 'String', 'missing/changed BD schema')
    rows = {}
    for row in stream:
        require(key(row) not in rows, 'duplicate final allele identity; do not collapse multiplicity')
        require(row['sample'].get('BD') in ('TP', 'FN'), 'unsupported final truth decision')
        genotype(row)
        truth_shape(row)
        annotations(row, header)
        rows[key(row)] = row
        require(len(rows) <= 10000, 'final record bound exceeded')
    return header, rows


def outcome(a, b):
    require(a in ('TP', 'FN') and b in ('TP', 'FN'), 'invalid decision')
    return {('TP', 'TP'): 'matched_by_both', ('TP', 'FN'): 'DELLY_only',
            ('FN', 'TP'): 'SvABA_only', ('FN', 'FN'): 'missed_by_both'}[(a, b)]


def margin(rows):
    counts = Counter(row['outcome'] for row in rows)
    require(set(counts) <= set(OUTCOMES), 'unknown outcome')
    total = len(rows)
    return {'truth_records': total, 'outcomes': {k: counts[k] for k in OUTCOMES},
            'shared_miss_fraction': counts['missed_by_both'] / total if total else None,
            'union_matched_fraction': (total - counts['missed_by_both']) / total if total else None}


def value_state(field):
    if not field['present']:
        return ('not_encoded', None)
    if field['value'] is None:
        return ('explicit_missing', None)
    return ('value', field['value'])


def repeat_state(row):
    ann = row['annotations']
    trf = 'flag_present' if ann['TRF']['present'] else 'flag_not_encoded'
    rm = [ann[k] for k in ('RM_score', 'RM_repeat', 'RM_clsfam')]
    state = 'value_present' if any(x['value'] is not None for x in rm) else (
        'fields_explicit_missing' if any(x['present'] for x in rm) else 'not_encoded')
    return trf, state


def summarize(rows):
    class_size, repeats, families, remaps, lcr = [], [], [], [], []
    for kind in CLASSES:
        for bin_name in BINS:
            subset = [r for r in rows if r['type'] == kind and r['size_bin'] == bin_name]
            base = {'type': kind, 'size_bin': bin_name}
            class_size.append({**base, **margin(subset)})
            for trf in ('flag_present', 'flag_not_encoded'):
                for rm in ('value_present', 'fields_explicit_missing', 'not_encoded'):
                    selected = [r for r in subset if repeat_state(r) == (trf, rm)]
                    repeats.append({**base, 'TRF_annotation_state': trf, 'RM_annotation_state': rm, **margin(selected)})
            for name, target in [('RM_clsfam', families), ('REMAP', remaps)]:
                groups = defaultdict(list)
                # Missing/absent strata exist explicitly even when their count is zero.
                groups[('not_encoded', None)] = []
                groups[('explicit_missing', None)] = []
                for row in subset:
                    groups[value_state(row['annotations'][name])].append(row)
                for (state, value), members in sorted(groups.items(), key=lambda item: str(item[0])):
                    target.append({**base, 'state': state, 'literal_value': value, **margin(members)})
            for status in OUTCOMES:
                selected = [r for r in subset if r['outcome'] == status]
                states = Counter(value_state(r['annotations']['LCR'])[0] for r in selected)
                vals = [r['annotations']['LCR']['value'] for r in selected if r['annotations']['LCR']['value'] is not None]
                lcr.append({**base, 'outcome': status, 'truth_records': len(selected),
                            'not_encoded': states['not_encoded'], 'explicit_missing': states['explicit_missing'],
                            'finite_values': len(vals), 'minimum': min(vals) if vals else None,
                            'median': statistics.median(vals) if vals else None, 'maximum': max(vals) if vals else None})
    for table in (class_size, repeats, families, remaps):
        require(sum(r['truth_records'] for r in table) == len(rows), 'table denominator mismatch')
        for label in OUTCOMES:
            require(sum(r['outcomes'][label] for r in table) == sum(r['outcome'] == label for r in rows),
                    'table outcome mismatch')
    require(sum(r['truth_records'] for r in lcr) == len(rows), 'LCR denominator mismatch')
    return {'total': margin(rows), 'class_size': class_size, 'TRF_RM_by_class_size': repeats,
            'RM_class_family_by_class_size': families, 'REMAP_by_class_size': remaps,
            'raw_LCR_by_class_size_and_outcome': lcr}


def load_bed(path):
    intervals = []
    for line in Path(path).read_text().splitlines():
        fields = line.split('\t')
        require(len(fields) == 3 and fields[0] == 'chr20', 'unexpected benchmark BED')
        start, end = integer(fields[1]), integer(fields[2])
        require(0 <= start < end and (not intervals or intervals[-1][1] <= start), 'overlapping/unsorted BED')
        intervals.append((start, end))
    return intervals


def build_ledger(paths, expected):
    def check_count_types(value):
        if isinstance(value, dict):
            for item in value.values():
                check_count_types(item)
        else:
            require(type(value) is int and value >= 0, 'expected count must be nonnegative integer')
    check_count_types(expected)
    dh, delly = read_final(paths['delly_final_truth'])
    sh, svaba = read_final(paths['svaba_final_truth'])
    require(set(delly) == set(svaba) and len(delly) == expected['truth_records'], 'final truth universes differ')
    stream = stream_vcf(paths['original_truth_chr20_no_star'])
    oh = next(stream)
    for name in ANNOTATIONS:
        require(oh['info_definitions'].get(name) == dh['info_definitions'].get(name) == sh['info_definitions'].get(name),
                'annotation header lineage differs')
    intervals = load_bed(paths['official_chr20_benchmark'])
    require(len(intervals) == expected['benchmark_intervals'] and sum(b-a for a, b in intervals) == expected['benchmark_bases'],
            'benchmark BED scope differs')
    rows, seen, original_count = [], set(), 0
    for original in stream:
        original_count += 1
        identity = key(original)
        if identity not in delly:
            continue
        require(identity not in seen, 'duplicate original selected identity; do not collapse multiplicity')
        seen.add(identity)
        dr, sr = delly[identity], svaba[identity]
        kind, size, bin_name = truth_shape(original)
        require(truth_shape(dr) == truth_shape(sr) == (kind, size, bin_name), 'type/size reclassification')
        require(genotype(original) == genotype(dr) == genotype(sr), 'original/final genotype differs')
        ann = annotations(original, oh)
        # Numeric encoding may be reformatted, but typed value and presence must agree.
        typed = lambda a: {k: (v['present'], v['value']) for k, v in a.items()}
        require(typed(ann) == typed(annotations(dr, dh)) == typed(annotations(sr, sh)), 'annotation value/presence changed')
        start, end = original['pos'] - 1, original['pos'] - 1 + len(original['ref'])
        contained = any(a <= start and end <= b for a, b in intervals)
        rows.append({'original_ordinal': original['ordinal'], 'chrom': original['chrom'], 'pos_1based': original['pos'],
                     'id': original['id'], 'ref': original['ref'], 'alts': original['alts'],
                     'allele_sha256': hashlib.sha256(json.dumps([original['ref'], original['alts']], separators=(',', ':')).encode()).hexdigest(),
                     'ref_span_0based_half_open': [start, end], 'type': kind, 'size': size, 'size_bin': bin_name,
                     'REF_span_contained_in_one_benchmark_interval': contained,
                     'original_raw_record': original['raw_record'], 'original_GT': original['sample']['GT'],
                     'GT_encoding_matches_original': {name: r['sample']['GT'] == original['sample']['GT']
                                                      for name, r in [('DELLY', dr), ('SvABA', sr)]},
                     'annotations': ann, 'outcome': outcome(dr['sample']['BD'], sr['sample']['BD']),
                     'callers': {name: {'refined_baseline_ordinal': r['ordinal'], 'decision': r['sample']['BD'],
                                       'GT': r['sample']['GT'], 'raw_info': r['info'], 'raw_sample': r['sample']}
                                 for name, r in [('DELLY', dr), ('SvABA', sr)]}})
    require(original_count == expected['original_full_context_records'], 'original full-context count differs')
    require(seen == set(delly) and len(rows) == expected['truth_records'], 'selected original truth missing')
    require(dict(Counter(r['type'] for r in rows)) == expected['classes'], 'class margin differs')
    summary = summarize(rows)
    require(summary['total']['outcomes'] == expected['outcomes'], 'four-outcome margin differs')
    for kind in CLASSES:
        require(margin([r for r in rows if r['type'] == kind])['outcomes'] == expected['outcomes_by_class'][kind],
                'class outcome margin differs')
    for caller in ('DELLY', 'SvABA'):
        observed = Counter(('INS_like' if r['type'] == 'INS' else 'DEL') + ':' + r['size_bin'] + ':' + r['callers'][caller]['decision'] for r in rows)
        require(dict(observed) == expected['per_caller_by_size'][caller], 'published class/size margins differ')
    return rows, summary, {'original': oh, 'DELLY': dh, 'SvABA': sh}


def check_receipts(paths, pins, rows):
    read = lambda role: parse_json(Path(paths[role]).read_text())
    for role in ('delly_status', 'svaba_status'):
        require(read(role)['status'] == 'complete', 'historical comparison incomplete')
    ds, ss, audit = read('delly_summary'), read('svaba_summary'), read('prior_final_audit')
    require(audit['comparison_status'] == 'complete', 'historical final audit incomplete')
    for flag in ('truth_alleles_equal_between_callers', 'truth_alleles_unique', 'truth_genotypes_unchanged_between_callers',
                 'every_final_truth_allele_occurs_exactly_in_original_input', 'counts_reconcile', 'overlap_totals_reconcile'):
        require(audit[flag] is True, 'historical allele audit did not pass')
    for name, role in [('benchmark/refine.base.vcf.gz', 'svaba_final_truth'), ('summary.json', 'svaba_summary')]:
        require(all(audit['artifacts'][name][k] == pins[role][k] for k in ('bytes', 'sha256')), 'historical audit artifact pin differs')
    for name, role in [('baseline', 'original_truth_chr20_no_star'), ('benchmark_bed', 'official_chr20_benchmark')]:
        require(all(ds['inputs'][name][k] == pins[role][k] for k in ('path', 'bytes', 'sha256')), 'original comparison input binding differs')
    for receipt, field in [(ds, 'config'), (ss, 'configuration')]:
        require(all(receipt[field][k] == pins['comparison_config'][k] for k in ('path', 'bytes', 'sha256')), 'configuration binding differs')
    config = read('comparison_config')
    require(config['truvari_version'] == '5.5.0', 'frozen comparator version differs')
    for caller, summary, role in [('DELLY', ds, 'delly_variant_summary'), ('SvABA', ss, 'svaba_variant_summary')]:
        counts = Counter(r['callers'][caller]['decision'] for r in rows)
        refined = read(role)
        require(counts['TP'] == refined['TP-base'] and counts['FN'] == refined['FN'], 'refined count receipt differs')
        nested = summary['comparison']['refined']
        require(nested['TP-base'] == counts['TP'] and nested['FN'] == counts['FN'], 'summary count receipt differs')
    observed = Counter(('INS_like' if r['type'] == 'INS' else 'DEL') + ':' + r['size_bin'] + ':' + r['callers']['SvABA']['decision'] for r in rows)
    require(dict(observed) == ss['by_size']['baseline'], 'published SvABA class/size table differs')


def run(config_path, expected_config_sha256, output):
    config_path, output = Path(config_path), Path(output).resolve()
    require(not output.exists(), 'existing output preserved')
    config_pin = fingerprint(config_path)
    require(config_pin['sha256'] == expected_config_sha256, 'config SHA differs')
    config = parse_json(config_path.read_text())
    require(config.get('ready') is True, 'real input analysis not approved')
    require(set(config['inputs']) == set(ROLES), 'input role membership differs')
    pins = {role: bind(item['path'], item) for role, item in config['inputs'].items()}
    require(sum(x['bytes'] for x in pins.values()) <= 8 * 1024**2, 'compressed input bound exceeded')
    paths = {k: v['path'] for k, v in pins.items()}
    rows, summary, headers = build_ledger(paths, config['expected'])
    check_receipts(paths, pins, rows)
    for role, pin in pins.items():
        bind(paths[role], pin)
    summary.update(scope=config['scope'], inputs=pins, config=config_pin, source=fingerprint(__file__),
                   annotation_policy='Literal existing public annotations only. Absence is not assessed absence of repeats; LCR remains continuous.',
                   genotype_policy='Allele-content equality ignores GT order/phase under the old audit policy. Exact raw GT and encoding differences are retained; unchanged phase is not asserted.',
                   coordinate_policy='VCF POS is one-based; REF span is zero-based half-open and includes the anchor. BED containment is descriptive and never a new eligibility filter.',
                   limitations=['No new caller, comparator, alignment, read support, normalization or nearby-equivalence analysis.',
                                'Conditional counts describe a fixed draft public chromosome benchmark; association is not cause.',
                                'Variants in repeat regions are not independent statistical observations; no confidence intervals.',
                                'Union truth detection is not merged-call precision, clinical sensitivity or a negative patient finding.'])
    payloads = {'ledger.jsonl': ''.join(json.dumps(r, sort_keys=True, allow_nan=False) + '\n' for r in rows),
                'summary.json': json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + '\n',
                'headers.json': json.dumps(headers, indent=2, sort_keys=True, allow_nan=False) + '\n'}
    require(sum(len(x.encode()) for x in payloads.values()) <= 64 * 1024**2, 'output bound exceeded')
    output.mkdir(parents=True, exist_ok=False)
    for name, data in payloads.items():
        (output / name).write_text(data)
    status = {'status': 'complete', 'truth_records': len(rows), 'outputs': {name: fingerprint(output/name) for name in payloads}}
    (output/'status.json').write_text(json.dumps(status, indent=2, sort_keys=True) + '\n')
    return status


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--expected-config-sha256', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    try:
        result = run(args.config, args.expected_config_sha256, args.output)
    except (ValueError, KeyError, OSError, TypeError) as error:
        print(str(error), file=sys.stderr)
        return 2
    print(json.dumps({'status': result['status'], 'truth_records': result['truth_records']}))
    return 0


if __name__ == '__main__':
    sys.exit(main())
