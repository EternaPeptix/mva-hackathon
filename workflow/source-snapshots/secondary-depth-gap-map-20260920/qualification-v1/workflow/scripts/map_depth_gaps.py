#!/usr/bin/env python3
"""Map saved aligned-read depth to the unchanged all-transcript CDS geometry."""
import argparse
from collections import Counter, defaultdict
import gzip
import json
import math
from pathlib import Path
import re
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent / 'vendor'))
from coding_depth_targets import prepare_targets
from guard_helpers import fingerprint, parse, require, save, stamp

MODES = ('MQ20_BQ20', 'MQ60_BQ20')
MODE_VALUES = ({'name': MODES[0], 'minimum_MQ': 20, 'minimum_BQ': 20},
               {'name': MODES[1], 'minimum_MQ': 60, 'minimum_BQ': 20})
ROLES = {'gff', 'fai', 'public_design', 'public_status', 'public_summary',
         'selection', 'summary', 'status', 'verification'} | {
             f'{m}_{kind}' for m in MODES for kind in ('depth', 'summary')}
SOURCE_FILES = ('map_depth_gaps.py', 'vendor/coding_depth_targets.py', 'vendor/guard_helpers.py')
LIMITATIONS = [
    'Saved aligned-read depth only; no alignments opened and no depth/caller rerun.',
    'All annotated CDS transcripts under the frozen Ensembl 116 policy; no selected clinical transcript policy, translation or reportability inference.',
    'Depth <20 and zero are descriptive bins, not clinical callability or a negative screen.',
    'Missing positions are errors, never zero or reference genotypes. X/Y depth is not diploid-rescaled.',
    'Both modes have BQ20; they reuse the same reads and differ in the MQ threshold. No additional independent evidence or locus-origin validation.',
    'Transcript totals overlap. Gene totals may overlap. Only the panel union counts each genomic base once.',
    'No promoter, UTR, intronic splice flank, regulatory, repeat-expansion, structural or low-fraction sensitivity claim.']


def typed_equal(a, b):
    return json.dumps(a, sort_keys=True, allow_nan=False) == json.dumps(b, sort_keys=True, allow_nan=False)


def require_equal(a, b, label):
    require(typed_equal(a, b), label + ' differs')


def read_json(path):
    value = parse(Path(path).read_text())
    def finite(item):
        if isinstance(item,float):require(math.isfinite(item),'nonfinite JSON number')
        elif isinstance(item,dict):
            for v in item.values():finite(v)
        elif isinstance(item,list):
            for v in item:finite(v)
    finite(value)
    return value


def require_pin(actual, saved, label):
    for key in ('path', 'bytes', 'sha256'):
        require(key in saved and typed_equal(actual[key], saved[key]), label + ' pin ' + key + ' differs')


def bind(config_path, expected_sha):
    config_pin = fingerprint(config_path, {'sha256': expected_sha})
    config = read_json(config_path)
    require(config.get('schema_version') == 1 and type(config['schema_version']) is int, 'schema')
    require(config.get('ready') is True and config.get('assembly') == 'GRCh38', 'not approved')
    require(config.get('scope') == '84-gene-all-transcript-CDS-saved-depth', 'scope')
    require_equal(config.get('depth_modes'), list(MODE_VALUES), 'depth modes')
    require(set(config['inputs']) == ROLES, 'input membership')
    require(set(config['source_sha256']) == set(SOURCE_FILES), 'source membership')
    pins = {'config': config_pin}
    for rel, digest in config['source_sha256'].items():
        pins['source/' + rel] = fingerprint(Path(__file__).resolve().parent / rel, {'sha256': digest})
    # This finite role allowlist deliberately contains no alignment/index/runtime executable.
    for role in sorted(ROLES):
        pins[role] = fingerprint(config['inputs'][role]['path'], config['inputs'][role])
    doc = {role: read_json(pins[role]['path']) for role in ROLES if role not in {'gff', 'fai'} and not role.endswith('_depth')}
    for role in ('status', 'summary', 'public_status', 'public_summary', 'verification'):
        require(doc[role].get('status') == 'complete', role + ' incomplete')
    require(doc['verification'].get('all_verified_inputs_stable') is True
            and doc['verification'].get('public_target_reconciled') is True, 'independent verification incomplete')
    require_pin(pins['status'], doc['verification']['producer_status'], 'independent producer status')
    require_pin(pins['summary'], doc['verification']['producer_summary'], 'independent producer summary')
    require(doc['status']['summary_sha256'] == pins['summary']['sha256'], 'terminal summary hash')
    for role in ('selection', 'summary') + tuple(f'{m}_{k}' for m in MODES for k in ('depth', 'summary')):
        basename = Path(pins[role]['path']).name
        require_pin(pins[role], doc['status']['outputs'][basename], 'terminal ' + role)
    for role in ('public_design', 'public_summary'):
        require_pin(pins[role], doc['public_status']['outputs'][Path(pins[role]['path']).name], role)
    for role in ('status', 'summary', 'selection', 'public_design', 'public_status', 'gff', 'fai') + tuple(
            f'{m}_{k}' for m in MODES for k in ('depth', 'summary')):
        require_pin(pins[role], doc['verification']['verified_files'][pins[role]['path']], 'independent ' + role)
    for source_role, saved_role in (('gff', 'gff'), ('fai', 'reference_index')):
        for role in ('public_design', 'public_status', 'summary'):
            require_pin(pins[source_role], doc[role]['inputs'][saved_role], role + ' annotation input')
    require(doc['public_design']['reads_inspected'] is False and doc['public_status']['reads_inspected'] is False,
            'public design read scope')
    return config, pins, doc


def reference_lengths(path):
    result = {}
    for line in Path(path).read_text().splitlines():
        f = line.split('\t')
        require(len(f) >= 2 and re.fullmatch(r'[1-9][0-9]*', f[1]) is not None and f[0] not in result, 'invalid FAI')
        result[f[0]] = int(f[1])
    require(result, 'empty FAI')
    return result


def geometry(config, pins, doc):
    public = doc['public_design']
    require(public['coordinates'] == '0-based half-open', 'public coordinate convention')
    names = [g['name'] for g in public['genes']]
    rebuilt = prepare_targets(pins['gff']['path'], names, reference_lengths(pins['fai']['path']))
    require(len(rebuilt['genes']) == len(public['genes']), 'public gene count')
    for expected, actual in zip(rebuilt['genes'], public['genes']):
        require_equal(expected, {k: actual[k] for k in expected}, 'public gene geometry')
    require_equal(rebuilt['genes'], doc['selection']['genes'], 'selected gene geometry')
    for field in ('panel_union_bases', 'gene_union_bases_sum'):
        require_equal(rebuilt[field], public[field], 'public ' + field)
        require_equal(rebuilt[field], doc['selection'][field], 'selection ' + field)
    stripped = [{k: w[k] for k in ('chrom', 'start', 'end', 'gene_ids')} for w in public['windows']]
    require_equal(rebuilt['windows'], stripped, 'public windows')
    require_equal(rebuilt['windows'], doc['selection']['windows'], 'selected windows')
    require(doc['selection']['coordinates'] == '0-based half-open; all CDS phases retained, not translated', 'selection coordinates')
    # This map has no per-base FASTA dependency. Require the already qualified all-ACGT scope.
    require(all(type(w['reference_non_ACGT_bases']) is int and w['reference_non_ACGT_bases'] == 0
                and type(w['reference_ACGT_bases']) is int and w['reference_ACGT_bases'] == w['end'] - w['start']
                for w in public['windows']), 'non-ACGT scope unsupported')
    require(doc['public_summary']['panel_non_ACGT_bases'] == 0
            and doc['public_summary']['panel_ACGT_bases'] == rebuilt['panel_union_bases'], 'public ACGT totals')
    scope = dict(genes=len(rebuilt['genes']), transcripts=len({t for g in rebuilt['genes'] for t in g['transcript_ids']}),
                 windows=len(rebuilt['windows']), panel_union_bases=rebuilt['panel_union_bases'],
                 gene_union_bases_sum=rebuilt['gene_union_bases_sum'])
    require_equal(scope, config['expected_scope'], 'expected scope')
    require_equal([{k: g[k] for k in ('id', 'name', 'chrom', 'strand', 'CDS_union_bases', 'transcript_ids')}
                   for g in rebuilt['genes']], doc['summary']['genes'], 'summary genes')
    for key in ('panel_union_bases', 'gene_union_bases_sum'):
        require_equal(rebuilt[key], doc['summary'][key], 'summary ' + key)
    return rebuilt


def atoms(design):
    """CDS source-row membership boundaries; genomic order on either strand."""
    events = defaultdict(lambda: defaultdict(Counter))
    for g in design['genes']:
        for row in g['cds_rows']:
            require(type(row['ordinal']) is int and row['ordinal'] > 0, 'CDS source ordinal')
            key = (g['id'], row['transcript_id'], row['ordinal'])
            events[row['chrom']][row['start']][key] += 1
            events[row['chrom']][row['end']][key] -= 1
    result = []
    order = list(dict.fromkeys(w['chrom'] for w in design['windows']))
    for chrom in order:
        active, previous = Counter(), None
        for pos, delta in sorted(events[chrom].items()):
            members = sorted(k for k, n in active.items() if n > 0)
            if members and previous is not None and previous < pos:
                result.append(dict(chrom=chrom, start=previous, end=pos, members=members))
            active.update(delta)
            require(all(n >= 0 for n in active.values()), 'negative CDS membership')
            previous = pos
        require(not any(active.values()), 'unclosed CDS membership')
    require(sum(a['end'] - a['start'] for a in result) == design['panel_union_bases'], 'atomic denominator')
    return result


def depth_values(path, windows):
    with gzip.open(path, 'rt') as stream:
        for w in windows:
            for pos in range(w['start'], w['end']):
                row = next(stream, '').rstrip('\n').split('\t')
                require(len(row) == 3 and row[0] == w['chrom']
                        and re.fullmatch(r'[1-9][0-9]*', row[1]) is not None
                        and int(row[1]) == pos + 1 and re.fullmatch(r'(0|[1-9][0-9]*)', row[2]) is not None,
                        'missing, malformed, extra, duplicate or out-of-order depth position')
                yield w['chrom'], pos, int(row[2])
        require(next(stream, None) is None, 'extra depth after target')


def category(depth):
    return 'zero' if depth == 0 else '1-19' if depth < 20 else 'at-least-20'


def metric(hist):
    return dict(CDS_union_bases=sum(hist.values()), zero_depth_bases=hist[0],
                below20_bases=sum(n for d, n in hist.items() if d < 20),
                sum_depth=sum(d*n for d, n in hist.items()),
                histogram={str(d): n for d, n in sorted(hist.items()) if n})


def normalized_histogram(value):
    require(isinstance(value, dict), 'histogram object')
    require(all(isinstance(k, str) and re.fullmatch(r'(0|[1-9][0-9]*)', k)
                and type(v) is int and v >= 0 for k,v in value.items()), 'invalid histogram')
    return {k:v for k,v in value.items() if v}


def calculate(design, depth_paths, emit):
    streams = [iter(depth_values(p, design['windows'])) for p in depth_paths]
    hist = {m: defaultdict(Counter) for m in MODES}
    windows = {m: [Counter() for _ in design['windows']] for m in MODES}
    transitions, increases = Counter(), 0
    window_index, interval_count = 0, 0
    for atom in atoms(design):
        genes = sorted({x[0] for x in atom['members']})
        txs = sorted({x[1] for x in atom['members']})
        keys = ['panel_union'] + genes + txs
        run = None
        for pos in range(atom['start'], atom['end']):
            values = [next(s, None) for s in streams]
            require(all(v is not None and v[:2] == (atom['chrom'], pos) for v in values), 'atomic/depth identity')
            depths = [v[2] for v in values]
            cats = tuple(category(d) for d in depths)
            transitions[cats] += 1
            increases += depths[1] > depths[0]
            while not (design['windows'][window_index]['chrom'] == atom['chrom'] and
                       design['windows'][window_index]['start'] <= pos < design['windows'][window_index]['end']):
                window_index += 1
                require(window_index < len(design['windows']), 'position outside windows')
            require_equal(genes, design['windows'][window_index]['gene_ids'], 'atom gene membership')
            for m, d in zip(MODES, depths):
                windows[m][window_index][d] += 1
                for key in keys:
                    hist[m][key][d] += 1
            is_gap = any(d < 20 for d in depths)
            if run is not None and (not is_gap or run['categories'] != list(cats)):
                emit(run); interval_count += 1; run = None
            if is_gap:
                if run is None:
                    run = dict(chrom=atom['chrom'], start=pos, end=pos, coordinates='0-based half-open',
                               gene_ids=genes, transcript_ids=txs,
                               source_CDS_members=[dict(gene_id=g, transcript_id=t, gff_line_ordinal=o) for g,t,o in atom['members']],
                               categories=list(cats), depth_modes=list(MODES),
                               depth_statistics=[dict(minimum=d, maximum=d, sum_depth=0) for d in depths])
                run['end'] = pos + 1
                for s, d in zip(run['depth_statistics'], depths):
                    s['minimum'] = min(s['minimum'], d); s['maximum'] = max(s['maximum'], d); s['sum_depth'] += d
        if run is not None:
            emit(run); interval_count += 1
    require(all(next(s, None) is None for s in streams), 'extra streamed positions')
    tx_metadata = {}
    for g in design['genes']:
        for tx in g['transcript_ids']:
            require(tx not in tx_metadata, 'transcript attached to multiple genes')
            tx_metadata[tx] = dict(gene_id=g['id'], gene_name=g['name'], chrom=g['chrom'], strand=g['strand'],
                                  cds_rows=[r for r in g['cds_rows'] if r['transcript_id'] == tx])
    result = dict(panel_union={m: metric(hist[m]['panel_union']) for m in MODES},
                  genes=[dict(gene_id=g['id'], name=g['name'], chrom=g['chrom'], strand=g['strand'],
                              modes={m: metric(hist[m][g['id']]) for m in MODES}) for g in design['genes']],
                  transcripts=[dict(transcript_id=t, **meta, modes={m: metric(hist[m][t]) for m in MODES})
                               for t, meta in sorted(tx_metadata.items())],
                  paired_categories=[dict(MQ20_BQ20=a, MQ60_BQ20=b, bases=transitions[(a,b)])
                                     for a in ('zero','1-19','at-least-20') for b in ('zero','1-19','at-least-20')],
                  MQ60_depth_greater_than_MQ20_bases=increases, gap_intervals=interval_count,
                  limitations=LIMITATIONS)
    return result, hist, windows


def reconcile(design, doc, hist, windows):
    for m, mode in zip(MODES, MODE_VALUES):
        saved = doc[m + '_summary']
        require_equal(saved['mode'], mode, m + ' mode')
        require_equal(doc['summary']['depth_modes'][m], {k: saved[k] for k in ('mode','position_count','groups')}, m + ' saved summary')
        require_equal(saved['position_count'], design['panel_union_bases'], m + ' position count')
        require(len(saved['windows']) == len(design['windows']), 'saved window count')
        for expected, actual, h in zip(design['windows'], saved['windows'], windows[m]):
            require_equal({k: actual[k] for k in expected}, expected, 'window identity')
            for name in ('span_depth', 'ACGT_depth'):
                measured = metric(h)
                require_equal(normalized_histogram(actual[name]['histogram']), measured['histogram'], m + ' window histogram')
                for field, value in dict(bases=measured['CDS_union_bases'], sum_depth=measured['sum_depth'],
                        mean=measured['sum_depth']/measured['CDS_union_bases'],
                        bases_at_least={str(t):sum(n for d,n in h.items() if d>=t) for t in (1,10,20,30)}).items():
                    require_equal(actual[name][field],value,m+' window '+field)
        require(set(saved['groups']) == {'panel_union'} | {g['id'] for g in design['genes']}, 'group membership')
        for key, group in saved['groups'].items():
            h = hist[m][key]; measured = metric(h)
            for name in ('span', 'ACGT'):
                for field, value in dict(bases=measured['CDS_union_bases'], zero_depth_bases=measured['zero_depth_bases'],
                                         sum_depth=measured['sum_depth'], mean=measured['sum_depth']/measured['CDS_union_bases']).items():
                    require_equal(group[name][field], value, m + ' group ' + field)
                require_equal(normalized_histogram(group[name]['histogram']),measured['histogram'],m+' group histogram')
                require_equal(group[name]['bases_at_least'], {str(t): sum(n for d,n in h.items() if d >= t) for t in (1,10,20,30)}, 'depth bins')
            require_equal(group['non_ACGT_reference_bases'],0,'group ACGT scope')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config', type=Path, required=True)
    ap.add_argument('--expected-config-sha256', required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    require(not args.output.exists(), 'existing output refused')
    config, pins, doc = bind(args.config, args.expected_config_sha256)
    design = geometry(config, pins, doc)
    args.output.mkdir(parents=True, exist_ok=False)
    state = dict(status='running', started_at=stamp(), inputs=pins, reads_opened=False)
    save(args.output/'status.json', state)
    try:
        with gzip.open(args.output/'gap-intervals.jsonl.gz', 'wt') as stream:
            result, hist, windows = calculate(design, [pins[m+'_depth']['path'] for m in MODES],
                lambda row: stream.write(json.dumps(row, sort_keys=True, allow_nan=False)+'\n'))
        reconcile(design, doc, hist, windows)
        for pin in pins.values():
            fingerprint(pin['path'], pin)
        result.update(status='complete', completed_at=stamp(), expected_scope=config['expected_scope'], inputs=pins)
        save(args.output/'summary.json', result)
        state.update(status='complete', completed_at=stamp(),
                     outputs={p.name: fingerprint(p) for p in args.output.iterdir() if p.is_file() and p.name!='status.json'})
        save(args.output/'status.json', state)
    except BaseException as exc:
        state.update(status='failed', failed_at=stamp(), error=repr(exc));save(args.output/'status.json',state)
        raise


if __name__ == '__main__':
    main()
