#!/usr/bin/env python3
"""Reconcile native CNV outputs, retain every call, and annotate exploratory triage.

Full annotated calls stay outside Git. Overlap across bin sizes reuses the same
reads/method; gene/CDS overlap and descriptive thresholds do not classify a CNV.
"""
import argparse
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import time
from urllib.parse import unquote


def require(value, message):
    if not value:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for data in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(data)
    return h.hexdigest()


def overlap(a, b, c, d):
    return max(0, min(b, d) - max(a, c) + 1)


def union(intervals):
    out = []
    for a, b in sorted(set(intervals)):
        require(1 <= a <= b, 'invalid annotation interval')
        if out and a <= out[-1][1] + 1:
            out[-1][1] = max(b, out[-1][1])
        else:
            out.append([a, b])
    return out


def finite(v):
    return isinstance(v, (float, int)) and math.isfinite(v)


def reasons(row, policy):
    failures = []
    if not row['within_reference_bounds']:
        failures.append('outside_reference_bounds')
    q = row['MQ_zero_fraction']
    if not finite(q) or not 0 <= q <= 1:
        failures.append('q0_missing_or_invalid')
    elif q > policy['maximum_q0']:
        failures.append('q0_above_limit')
    g = row['reference_gap_fraction']
    if not finite(g) or not 0 <= g <= 1:
        failures.append('gap_metric_missing_or_invalid')
    elif g > policy['maximum_gap_fraction']:
        failures.append('gap_fraction_above_limit')
    e = row['e1']
    if not finite(e) or e < 0 or e >= policy['maximum_e1_exclusive']:
        failures.append('e1_outside_limit')
    d = row['relative_depth']
    if not finite(d) or d < 0 or (d > policy['maximum_deletion_ratio'] if row['type'] == 'deletion'
                                      else d < policy['minimum_duplication_ratio']):
        failures.append('relative_depth_effect_below_limit')
    return failures


def load_genes(path, chromosomes):
    genes, parents, pending = {}, {}, []
    cds = defaultdict(set)
    with gzip.open(path, 'rt') as stream:
        for line in stream:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            require(len(f) == 9, 'GFF row schema differs')
            if f[0] not in chromosomes or f[2] in ('exon', 'five_prime_UTR', 'three_prime_UTR'):
                continue
            a = {k: unquote(v) for item in f[8].split(';') if '=' in item for k, v in [item.split('=', 1)]}
            ident = a.get('ID', '')
            if ident.startswith('gene:'):
                require(ident not in genes, 'duplicate gene identity')
                genes[ident] = {'id': ident.removeprefix('gene:'), 'name': a.get('Name', ident),
                    'biotype': a.get('biotype'), 'chrom': f[0], 'start': int(f[3]), 'end': int(f[4]), 'strand': f[6]}
            elif ident.startswith('transcript:'):
                parent = a.get('Parent', '').split(',')
                require(len(parent) == 1 and parent[0].startswith('gene:'), 'ambiguous transcript parent')
                require(ident not in parents or parents[ident] == parent[0], 'transcript parent differs')
                parents[ident] = parent[0]
            elif f[2] == 'CDS':
                for parent in a.get('Parent', '').split(','):
                    require(parent.startswith('transcript:'), 'CDS parent is not a transcript')
                    if parent in parents:
                        cds[parents[parent]].add((int(f[3]), int(f[4])))
                    else:
                        pending.append((parent, int(f[3]), int(f[4])))
    for parent, a, b in pending:
        require(parent in parents, 'unresolved CDS transcript')
        cds[parents[parent]].add((a, b))
    require(set(cds) <= genes.keys(), 'CDS has unrecognized gene')
    for ident, g in genes.items():
        g['cds'] = union(cds[ident])
        require(all(g['start'] <= a <= b <= g['end'] for a, b in g['cds']), 'CDS outside gene span')
    return genes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('analysis', 'reference', 'gff', 'config', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--expected-summary-sha256', required=True)
    parser.add_argument('--expected-gff-sha256', required=True)
    args = parser.parse_args(); os.nice(10)
    import datetime
    import numpy as np
    import pysam
    from cnvpytor.io import IO, FLAG_AUTO, FLAG_GC_CORR
    from run_copy_number import fingerprint, save, clean
    tick = time.monotonic()
    root, out = args.analysis.resolve(), args.output.resolve()
    require(not out.exists(), 'existing review preserved')
    state = json.loads((root / 'status.json').read_text())
    require(state['status'] == 'complete', 'CNV analysis is incomplete')
    require(sha(root / 'summary.json') == args.expected_summary_sha256 == state['summary_sha256'], 'summary identity differs')
    native = json.loads((root / 'summary.json').read_text())
    require(native['config']['scope'] == 'verified_four_lane_autosomal_discovery', 'unexpected analysis scope')
    for name, record in native['outputs'].items():
        require(Path(record['path']).resolve() == root / name and sha(root / name) == record['sha256'], 'native output differs')
    require(sha(args.reference) == native['inputs']['reference']['sha256'], 'reference differs')
    require(sha(args.gff) == args.expected_gff_sha256, 'annotation differs')
    config = json.loads(args.config.read_text())
    raw = json.loads((root / 'calls.json').read_text())['calls']
    chroms = native['config']['chromosomes']
    require(chroms == list(map(str, range(1, 23))), 'autosomal scope differs')
    genes = load_genes(args.gff, set(chroms))
    gindex = {}
    for c in chroms:
        rows = sorted((g for g in genes.values() if g['chrom'] == c), key=lambda g: g['start'])
        ends = []; last = 0
        for g in rows:
            last = max(last, g['end']); ends.append(last)
        gindex[c] = (rows, [g['start'] for g in rows], ends)
    io = IO(str(root / 'depth.pytor'), ro=True)
    fasta = pysam.FastaFile(str(args.reference))
    lengths = dict(io.get_chromosome_lengths())
    allrows, totals, gene_profiles = [], Counter(), {}
    maximum_gap_difference = 0.0
    try:
        for c in chroms:
            require(lengths[c] == fasta.get_reference_length(c), 'reference length differs')
            p, u = io.read_rd(c); observed = native['read_midpoint_import'][c]
            require(int(p.sum()) == observed['included_read_midpoints'], 'midpoint total differs')
            require(int(u.sum()) == observed['MQ_positive'] and int((p-u).sum()) == observed['MQ_zero'], 'MQ channel totals differ')
            require(np.all(u <= p), 'positive MQ exceeds all reads')
            require(observed['observed_records'] == observed['included_read_midpoints'] + observed.get('excluded_flags', 0)
                    + observed.get('excluded_unknown_MQ_255', 0) + observed.get('excluded_reference_skip', 0), 'import census does not conserve records')
            totals.update(observed)
        for b, rows in raw.items():
            size = int(b)
            require(Counter(r['type'] for r in rows) == native['calls_by_bin_size'][b], 'call counts differ')
            stat = io.get_signal(None, size, 'RD stat', FLAG_AUTO | FLAG_GC_CORR)
            require(float(stat[4]) == native['normalization'][b]['fitted_mean_read_midpoints_per_bin'], 'normalization differs')
            native_calls = {c: io.read_calls(c, size, 'calls', FLAG_GC_CORR) for c in chroms}
            offsets = Counter()
            for idx, row in enumerate(rows):
                c = row['chrom']; a, z = row['start_1based'], row['end_inclusive']
                require(type(a) is int and type(z) is int and row['size_bp'] == z-a+1, 'call coordinates differ')
                require(row['within_reference_bounds'] == (1 <= a <= z <= lengths[c]), 'reference-bound flag differs')
                item = native_calls[c][offsets[c]]; offsets[c] += 1
                require(item['type'] == (-1 if row['type'] == 'deletion' else 1), 'native call type differs')
                for src, dst in [('start','start_1based'),('end','end_inclusive'),('size','size_bp'),('cnv','relative_depth'),
                        ('p_val','e1'),('p_val_2','e2'),('p_val_3','e3'),('p_val_4','e4'),('Q0','MQ_zero_fraction'),('pN','reference_gap_fraction'),('dG','distance_to_gap_bp')]:
                    require(clean(item[src]) == row[dst], 'native call field differs: ' + dst)
                annotated = dict(row, call_id=b + ':' + str(idx+1), bin_size=size)
                annotated['triage_failures'] = {label: reasons(row, config[label]) for label in ('primary','stricter')}
                if row['within_reference_bounds']:
                    seq = fasta.fetch(c, a-1, z).upper()
                    gap = 1 - sum(seq.count(base) for base in 'ACGT') / len(seq)
                    annotated['direct_reference_non_ACGT_fraction'] = gap
                    maximum_gap_difference = max(maximum_gap_difference, abs(gap-row['reference_gap_fraction']))
                grows, starts, maxends = gindex[c]
                annotated['gene_overlaps'] = []
                for g in grows[bisect_left(maxends, a):bisect_right(starts, z)]:
                    n = overlap(a,z,g['start'],g['end'])
                    if n:
                        annotated['gene_overlaps'].append({k:g[k] for k in ('id','name','biotype')} | {
                            'gene_span_overlap_bp':n,'gene_span_length_bp':g['end']-g['start']+1,
                            'union_CDS_overlap_bp':sum(overlap(a,z,x,y) for x,y in g['cds']),
                            'union_CDS_length_bp':sum(y-x+1 for x,y in g['cds'])})
                allrows.append(annotated)
            require(all(offsets[c] == len(native_calls[c]) for c in chroms), 'unaccounted native calls')
            require(sum(not r['within_reference_bounds'] for r in rows) == native['out_of_reference_calls_by_bin_size'][b], 'reference-bound census differs')
            for name in config['review_genes']:
                matches = [g for g in genes.values() if g['name'] == name]
                require(len(matches) == 1, 'review gene missing or ambiguous: ' + name)
                g = matches[0]
                profile = gene_profiles.setdefault(name, {k:g[k] for k in ('id','chrom','start','end','biotype')} | {'by_bin_size':{}})
                signal = np.asarray(io.get_signal(g['chrom'],size,'RD',FLAG_GC_CORR),dtype=float)
                values = signal[(g['start']-1)//size:(g['end']-1)//size+1] / stat[4]
                profile['by_bin_size'][b] = {'intersecting_bins':len(values), 'finite_bins':int(np.isfinite(values).sum()),
                    'relative_density_quantiles':clean(np.nanquantile(values,[0,0.25,0.5,0.75,1]))}
        require(maximum_gap_difference < 1e-9, 'native gap metrics differ from exact reference sequence')
    finally:
        io.file.close(); fasta.close()
    for row in allrows:
        row['other_resolution_matches'] = []
        for other in allrows:
            if row['bin_size'] == other['bin_size'] or row['chrom'] != other['chrom'] or row['type'] != other['type']:
                continue
            n = overlap(row['start_1based'],row['end_inclusive'],other['start_1based'],other['end_inclusive'])
            if min(n/row['size_bp'],n/other['size_bp']) >= config['reciprocal_overlap_minimum']:
                row['other_resolution_matches'].append({'call_id':other['call_id'], 'bin_size':other['bin_size'],
                    'overlap_fraction_self':n/row['size_bp'],'overlap_fraction_other':n/other['size_bp'],
                    'passes_primary':not other['triage_failures']['primary'], 'passes_stricter':not other['triage_failures']['stricter']})
    count_summary = {}
    for b, rows in raw.items():
        items = [r for r in allrows if r['bin_size'] == int(b)]
        count_summary[b] = {'raw':len(items),'raw_types':dict(Counter(r['type'] for r in items)),
            'outside_reference':sum(not r['within_reference_bounds'] for r in items)}
        for label in ('primary','stricter'):
            passed = [r for r in items if not r['triage_failures'][label]]
            count_summary[b][label] = {'passed':len(passed),'types':dict(Counter(r['type'] for r in passed)),
                'failure_reasons_nonexclusive':dict(Counter(x for r in items for x in r['triage_failures'][label])),
                'with_CDS_overlap':sum(any(g['union_CDS_overlap_bp'] for g in r['gene_overlaps']) for r in passed),
                'with_same_type_reciprocal_match_also_passing':sum(any(m['passes_'+label] for m in r['other_resolution_matches']) for r in passed)}
    for name, profile in gene_profiles.items():
        profile['raw_gene_overlapping_calls'] = [{k:r[k] for k in ('call_id','type','bin_size','start_1based','end_inclusive','relative_depth','MQ_zero_fraction','reference_gap_fraction','triage_failures')}
            for r in allrows if any(g['name'] == name for g in r['gene_overlaps'])]
    out.mkdir(parents=True)
    save(out/'annotated-calls.json', {'configuration':config,'calls':allrows})
    selected = [r for r in allrows if not r['triage_failures']['primary'] and any(m['passes_primary'] for m in r['other_resolution_matches'])]
    result = {'completed_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'elapsed_seconds':time.monotonic()-tick,
        'scope':'Read-only output reconciliation and exploratory post-run CNV/gene triage; all raw calls retained outside Git.',
        'inputs':{'native_summary':fingerprint(root/'summary.json'),'gff':fingerprint(args.gff),'reference':fingerprint(args.reference),
                  'configuration':fingerprint(args.config),'script':fingerprint(Path(__file__))},
        'native_outputs_verified':native['outputs'],'native_call_records_equal':len(allrows),'read_midpoint_census':dict(totals),
        'normalization':native['normalization'],'exact_reference_gap_max_difference':maximum_gap_difference,
        'annotation_genes':len(genes),'counts':count_summary,'selected_gene_profiles':gene_profiles,
        'same_type_cross_resolution_primary_calls':selected,'output':fingerprint(out/'annotated-calls.json'),
        'limits':['Post-run thresholds describe a review subset; rejected intervals remain available and are not proven false.',
            'All resolutions reuse the same reads and method. Reciprocal overlap is not independent breakpoint validation.',
            'Gene overlap uses Ensembl gene spans and union of all annotated CDS; it does not establish a disrupted transcript, dosage sensitivity, pathogenicity or clinical probability.',
            'Profiles summarize intersecting whole bins, not exon-level coverage, absolute copy number, molecular phase or a clinical negative result.',
            'No BAF, parental origin, contamination assessment or public real-read CNV sensitivity is added here.']}
    save(out/'summary.json',result)
    print(json.dumps({'status':'complete','native_calls':len(allrows),'counts':count_summary,'cross_resolution_primary_calls':len(selected),'elapsed_seconds':result['elapsed_seconds']}))


if __name__ == '__main__':
    main()
