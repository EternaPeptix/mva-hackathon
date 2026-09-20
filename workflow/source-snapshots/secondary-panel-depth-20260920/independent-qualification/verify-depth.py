#!/usr/bin/env python3
"""Independent target/payload/bookkeeping audit plus an explicit same-engine recount.

No producer module is imported. Only a new audit directory is written.
"""
import argparse
from array import array
from bisect import bisect_right
from collections import Counter, defaultdict
import datetime
import gzip
import hashlib
import json
from pathlib import Path
import platform
import subprocess
from urllib.parse import unquote

CONFIG_SHA = '224763eda065ed0f658c4809aeacc7f4568ee190341086fcda7a0aa1d50fbcd8'
CONFIG_REL = 'results/secondary-panel-cds-depth-20260920/review-v2-launch/producer-config.json'
RUN_REL = 'results/secondary-panel-cds-depth-20260920/review-v2'
PUBLIC_REL = 'methods/secondary-panel-depth-preparation-20260920/public-target-v1'
PUBLIC_STATUS_SHA = '3673bafb827a6b1d6f6ac90b08b837c7053373fa5c3c21e4645f113b2ecf5fc2'
EXPECTED_GENES = 84
SOURCE_NAMES = {'measure_coding_depth.py', 'coding_depth_targets.py', 'survey_chromosome_depth.py', 'prepare_gatk_compatible_cram.py'}
QUALIFICATION_SHA = '74398ab873453551e1b2b3aafdf42f032949fcdbb8cc9ca5a9357a1f63f43fa6'
MODES = [dict(name='MQ0_BQ1', minimum_MQ=0, minimum_BQ=1),
         dict(name='MQ20_BQ20', minimum_MQ=20, minimum_BQ=20),
         dict(name='MQ60_BQ20', minimum_MQ=60, minimum_BQ=20)]


def require(ok, message):
    if not ok:
        raise ValueError(message)


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def same(actual, expected, message):
    require(canonical(actual) == canonical(expected), message)


def parse(text):
    def pairs(items):
        result = {}
        for k, v in items:
            require(k not in result, 'duplicate JSON key')
            result[k] = v
        return result
    def invalid(value):
        raise ValueError('nonfinite JSON value ' + value)
    return json.loads(text, object_pairs_hook=pairs, parse_constant=invalid)


class Pins:
    def __init__(self):
        self.files = {}

    def add(self, path, expected=None):
        path = Path(path).resolve()
        key = str(path)
        if expected and 'path' in expected:
            same(str(Path(expected['path']).resolve()), key, 'bound path differs')
        if key not in self.files:
            before = path.stat()
            h = hashlib.sha256()
            with path.open('rb') as stream:
                for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
                    h.update(block)
            after = path.stat()
            same([before.st_size, before.st_mtime_ns], [after.st_size, after.st_mtime_ns], 'file changed during hash')
            self.files[key] = dict(path=key, bytes=after.st_size, mtime_ns=after.st_mtime_ns, sha256=h.hexdigest())
        got = self.files[key]
        if expected:
            for field in ('sha256', 'bytes', 'mtime_ns'):
                if field in expected:
                    same(got[field], expected[field], 'bound ' + field + ' differs: ' + key)
        return got

    def read(self, path, expected=None):
        self.add(path, expected)
        return parse(Path(path).read_text())

    def stable(self):
        for pin in self.files.values():
            st = Path(pin['path']).stat()
            same([st.st_size, st.st_mtime_ns], [pin['bytes'], pin['mtime_ns']], 'input changed during audit')


def independent_targets(gff, names, reference_lengths):
    """Use per-base membership sets, independently of production interval sweep."""
    require(type(names) is list and names and len(names) == len(set(names)), 'invalid gene selection')
    genes, transcripts, cds = {}, {}, []
    with gzip.open(gff, 'rt') as stream:
        for ordinal, line in enumerate(stream, 1):
            if not line.strip() or line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            require(len(f) == 9, 'invalid annotation columns')
            if f[0] not in set(map(str, range(1, 23))) | {'X', 'Y'}:
                continue
            attrs = {}
            for pair in f[8].split(';'):
                if '=' in pair:
                    k, v = pair.split('=', 1)
                    require(k not in attrs, 'duplicate annotation attribute')
                    attrs[k] = unquote(v)
            ident = attrs.get('ID', '')
            if not (ident.startswith(('gene:', 'transcript:')) or f[2] == 'CDS'):
                continue
            a, b = int(f[3]) - 1, int(f[4])
            require(f[0] in reference_lengths and 0 <= a < b <= reference_lengths[f[0]], 'annotation outside reference')
            row = dict(chrom=f[0], start=a, end=b, strand=f[6], ordinal=ordinal)
            if ident.startswith('gene:'):
                require(ident not in genes, 'duplicate gene')
                genes[ident] = dict(row, id=ident, name=attrs.get('Name'), biotype=attrs.get('biotype'))
            elif ident.startswith('transcript:'):
                require(ident not in transcripts, 'duplicate transcript')
                transcripts[ident] = dict(row, parent=attrs.get('Parent'))
            else:
                require(f[7] in ('0', '1', '2'), 'invalid CDS phase')
                cds.append(dict(row, parents=attrs.get('Parent', '').split(','), phase=int(f[7])))
    selected = []
    positions = defaultdict(set)
    for name in names:
        matches = [g for g in genes.values() if g['name'] == name]
        require(len(matches) == 1, 'missing or ambiguous gene')
        g = dict(matches[0])
        require(g['biotype'] == 'protein_coding', 'gene not protein coding')
        tids = {tid for tid, t in transcripts.items() if t['parent'] == g['id']}
        used, rows, bases = set(), [], set()
        for c in cds:
            for tid in c['parents']:
                if tid not in tids:
                    continue
                t = transcripts[tid]
                require(c['chrom'] == t['chrom'] == g['chrom'] and c['strand'] == t['strand'] == g['strand'], 'CDS parent identity differs')
                require(g['start'] <= t['start'] <= c['start'] < c['end'] <= t['end'] <= g['end'], 'CDS parent geometry differs')
                used.add(tid)
                rows.append(dict(c, transcript_id=tid))
                bases.update(range(c['start'], c['end']))
        require(bases, 'gene without CDS')
        intervals = []
        for p in sorted(bases):
            positions[(g['chrom'], p)].add(g['id'])
            if intervals and intervals[-1][1] == p:
                intervals[-1][1] += 1
            else:
                intervals.append([p, p + 1])
        g.update(transcript_ids=sorted(used), cds_rows=rows, cds_union=intervals, CDS_union_bases=len(bases))
        selected.append(g)
    rank = {c: i for i, c in enumerate(reference_lengths)}
    windows = []
    for (chrom, p), ids in sorted(positions.items(), key=lambda row: (rank[row[0][0]], row[0][1])):
        ids = sorted(ids)
        if windows and windows[-1]['chrom'] == chrom and windows[-1]['end'] == p and windows[-1]['gene_ids'] == ids:
            windows[-1]['end'] = p + 1
        else:
            windows.append(dict(chrom=chrom, start=p, end=p + 1, gene_ids=ids))
    return dict(genes=selected, windows=windows, panel_union_bases=len(positions),
                gene_union_bases_sum=sum(g['CDS_union_bases'] for g in selected))


def metrics(hist, include_zero=False):
    total = sum(hist.values())
    result = dict(bases=total, sum_depth=sum(d * n for d, n in hist.items()),
        mean=sum(d * n for d, n in hist.items()) / total if total else None,
        bases_at_least={str(t): sum(n for d, n in hist.items() if d >= t) for t in (1, 10, 20, 30)},
        histogram={str(d): n for d, n in sorted(hist.items())})
    if include_zero:
        result['zero_depth_bases'] = hist.get(0, 0)
    return result


def positions_from_lines(lines, windows):
    iterator = iter(lines)
    depths = {}
    for w in windows:
        for p in range(w['start'], w['end']):
            fields = next(iterator, '').split()
            require(len(fields) == 3, 'missing or malformed depth position')
            require(fields[0] == w['chrom'] and int(fields[1]) == p + 1, 'depth position/order differs')
            d = int(fields[2])
            require(d >= 0, 'negative depth')
            key = (w['chrom'], p)
            require(key not in depths, 'duplicate target position')
            depths[key] = d
    require(next(iterator, None) is None, 'extra depth position')
    return depths


def aggregate(depths, design, sequences):
    rows = []
    group_hist = {name: [Counter(), Counter()] for name in ['panel_union'] + [g['id'] for g in design['genes']]}
    for w, sequence in zip(design['windows'], sequences):
        require(len(sequence) == w['end'] - w['start'], 'reference sequence length differs')
        hs = [Counter(), Counter()]
        for offset, base in enumerate(sequence):
            d = depths[(w['chrom'], w['start'] + offset)]
            hs[0][d] += 1
            if base in 'ACGT':
                hs[1][d] += 1
        rows.append(dict(w, span_depth=metrics(hs[0]), ACGT_depth=metrics(hs[1])))
        for name in ['panel_union'] + w['gene_ids']:
            for index in (0, 1):
                group_hist[name][index].update(hs[index])
    groups = {}
    for name, hs in group_hist.items():
        groups[name] = dict(span=metrics(hs[0], True), ACGT=metrics(hs[1], True),
            non_ACGT_reference_bases=sum(hs[0].values()) - sum(hs[1].values()))
    return rows, groups


def quality_available(read):
    if read.flag & 0xF04:
        return
    require(read.mapping_quality != 255, 'unknown MQ')
    q = read.query_qualities
    require(q is not None and len(q) == read.query_length and all(type(x) is int and 0 <= x < 255 for x in q), 'unknown BQ')


def record_payload(read):
    """Retain semantic tag types, arrays, all core fields and missing qualities.

    BAM integer storage widths c/C/s/S/i/I represent the same SAM integer type.
    Array subtype, noninteger types and every tag value remain exact.
    """
    tags = []
    for name, value, dtype in read.get_tags(with_value_type=True):
        if isinstance(value, array):
            value = dict(subtype=value.typecode, values=list(value))
        tags.append([name, 'i' if dtype in 'cCsSiI' else dtype, value])
    tags.sort(key=canonical)
    return dict(qname=read.query_name, flag=read.flag, ref=read.reference_name,
        start=read.reference_start, MQ=read.mapping_quality, cigar=read.cigartuples,
        mate_ref=read.next_reference_name, mate_start=read.next_reference_start,
        template_length=read.template_length, sequence=read.query_sequence,
        qualities=list(read.query_qualities) if read.query_qualities is not None else None, tags=tags)


def payload_key(read):
    return hashlib.sha256(canonical(record_payload(read)).encode()).hexdigest()


def first_overlapping_window(start, end, windows):
    index = bisect_right([w['end'] for w in windows], start)
    return index if index < len(windows) and windows[index]['start'] < end else None


def compare_extraction(original, selected, reference, index, design, config):
    import pysam
    expected, actual = Counter(), Counter()
    by_chrom = defaultdict(list)
    for w in design['windows']:
        by_chrom[w['chrom']].append(w)
    with pysam.AlignmentFile(str(original), reference_filename=str(reference), index_filename=str(index)) as source, \
            pysam.AlignmentFile(str(selected)) as extracted:
        same(source.header.to_dict(), extracted.header.to_dict(), 'extraction header differs')
        rg = extracted.header.to_dict().get('RG', [])
        require(len(rg) == config['expected_read_groups'] and len({x['ID'] for x in rg}) == len(rg), 'RG count or identity differs')
        require({x.get('SM') for x in rg} == {config['sample']}, 'sample differs')
        allowed = {x['ID'] for x in rg}
        groups, flags, mq = Counter(), Counter(), Counter()
        for read in extracted:
            require(read.has_tag('RG') and read.get_tag('RG') in allowed, 'unregistered RG')
            actual[payload_key(read)] += 1
            quality_available(read)
            groups[read.get_tag('RG')] += 1
            flags[read.flag] += 1
            if not read.flag & 0xF04:
                mq[read.mapping_quality] += 1
        for chrom, windows in by_chrom.items():
            for number, w in enumerate(windows):
                for read in source.fetch(chrom, w['start'], w['end']):
                    # Select a stored record once even if it spans multiple target
                    # windows. Identical source records retain their multiplicity.
                    end = read.reference_end or read.reference_start + 1
                    if first_overlapping_window(read.reference_start, end, windows) == number:
                        expected[payload_key(read)] += 1
        same(dict(actual), dict(expected), 'original-to-selected complete record payload or multiplicity differs')
    return dict(records=sum(actual.values()), unique_payloads=len(actual),
        canonical_multiset_sha256=hashlib.sha256(canonical(dict(sorted(actual.items()))).encode()).hexdigest(),
        groups={g: groups[g] for g in sorted(allowed)}, flags={str(k): v for k, v in sorted(flags.items())},
        eligible_MQ={str(k): v for k, v in sorted(mq.items())},
        header_equal=True, complete_core_and_tag_payload_equal=True, identical_record_multiplicity_retained=True,
        scope='Independent indexed original regional selection and complete semantic record multiset equality; tag order and BAM scalar-integer storage width are immaterial.')


def public_bindings(root, pins, config):
    """Bind the fixed public-only 84-gene selection before any read audit."""
    directory = root / PUBLIC_REL
    state = pins.read(directory / 'status.json', {'sha256': PUBLIC_STATUS_SHA})
    same(state['status'], 'complete', 'public preparation incomplete')
    same(state['reads_inspected'], False, 'public preparation inspected reads')
    same(sorted(state['outputs']), ['reference-design.json', 'regions.bed', 'target-summary.json'], 'public output membership differs')
    for name, pin in state['outputs'].items():
        pins.add(directory / name, pin)
    require(not ({'alignment', 'alignment_index'} & set(state['inputs'])), 'public preparation contains read inputs')
    for name, pin in state['inputs'].items():
        pins.add(pin['path'], pin)
        if name in ('reference', 'reference_index', 'gff'):
            pins.add(pin['path'], config['inputs'][name])
    panel = pins.read(state['inputs']['panel']['path'], state['inputs']['panel'])
    same(type(panel['total']).__name__, 'int', 'panel count type differs')
    same(panel['total'], EXPECTED_GENES, 'panel count differs')
    symbols = [row['symbol'] for row in panel['rows']]
    require(len(symbols) == len(set(symbols)) == EXPECTED_GENES, 'panel symbols missing or repeated')
    same(sorted(symbols), config['genes'], 'configuration panel membership/order differs')
    if 'gene_list_source' in config:
        pins.add(state['inputs']['panel']['path'], config['gene_list_source'])
    design = pins.read(directory / 'reference-design.json')
    summary = pins.read(directory / 'target-summary.json')
    same(design['status'], 'prepared', 'public design incomplete')
    same(design['reads_inspected'], False, 'public design inspected reads')
    same(design['inputs'], state['inputs'], 'public design input membership differs')
    same(design['panel_rows'], panel['rows'], 'public literal panel rows differ')
    same(summary['status'], 'complete', 'public summary incomplete')
    same(summary['reads_inspected'], False, 'public summary inspected reads')
    same(summary['genes'], EXPECTED_GENES, 'public summary gene count differs')
    bed = ''.join(f"{w['chrom']}\t{w['start']}\t{w['end']}\n" for w in design['windows'])
    same((directory / 'regions.bed').read_text(), bed, 'public BED differs')
    return design, summary


def reconcile_public_design(public, summary, independent, sequences):
    """Reconstruct every target base and reference sequence, not just totals."""
    for key in ('panel_union_bases', 'gene_union_bases_sum'):
        same(public[key], independent[key], 'public target denominator differs')
        same(summary[key], independent[key], 'public summary denominator differs')
    same(len(public['genes']), len(independent['genes']), 'public gene count differs')
    same(len(public['windows']), len(independent['windows']), 'public window count differs')
    same(summary['windows'], len(independent['windows']), 'public summary window count differs')
    gene_sequences = defaultdict(list)
    total_acgt = 0
    for actual, expected, sequence in zip(public['windows'], independent['windows'], sequences):
        same({k: actual[k] for k in expected}, expected, 'public window membership differs')
        acgt = sum(base in 'ACGT' for base in sequence)
        same(actual['reference_ACGT_bases'], acgt, 'public window ACGT differs')
        same(actual['reference_non_ACGT_bases'], len(sequence) - acgt, 'public window non-ACGT differs')
        same(actual['reference_sequence_sha256'], hashlib.sha256(sequence.encode()).hexdigest(), 'public window reference sequence differs')
        total_acgt += acgt
        for ident in expected['gene_ids']:
            gene_sequences[ident].append(sequence)
    gene_summary = []
    for actual, expected in zip(public['genes'], independent['genes']):
        same({k: actual[k] for k in expected}, expected, 'public complete gene/transcript/CDS payload differs')
        sequence = ''.join(gene_sequences[expected['id']])
        acgt = sum(base in 'ACGT' for base in sequence)
        same(actual['concatenated_CDS_union_sha256'], hashlib.sha256(sequence.encode()).hexdigest(), 'public gene reference sequence differs')
        same(actual['reference_ACGT_bases'], acgt, 'public gene ACGT differs')
        same(actual['reference_non_ACGT_bases'], len(sequence)-acgt, 'public gene non-ACGT differs')
        gene_summary.append({k:actual[k] for k in ('name','id','chrom','CDS_union_bases','reference_ACGT_bases','reference_non_ACGT_bases')})
    same(summary['gene_summary'], gene_summary, 'public gene summary differs')
    same(summary['panel_ACGT_bases'], total_acgt, 'public panel ACGT differs')
    same(summary['panel_non_ACGT_bases'], independent['panel_union_bases']-total_acgt, 'public panel non-ACGT differs')


def bind_launch(root, pins, config, expected_sha):
    require(type(expected_sha) is str and len(expected_sha) == 64 and all(c in '0123456789abcdef' for c in expected_sha), 'explicit completed launch SHA required')
    run = root / RUN_REL
    launch = run.with_name(run.name + '-launch')
    state = pins.read(launch / 'status.json', {'sha256': expected_sha})
    same(state['status'], 'complete', 'launch incomplete')
    same(state['exit_code'], 0, 'launch failed')
    binding = pins.read(launch / 'worker-binding.json', state['worker_binding'])
    same(binding['status'], 'complete', 'worker binding incomplete')
    same(binding['before_case_access'], True, 'worker binding was not before access')
    pins.add(root / CONFIG_REL, binding['producer_config'])
    same(binding['producer_config']['sha256'], CONFIG_SHA, 'worker inner configuration differs')
    for pin in binding['inputs'].values():
        pins.add(pin['path'], pin)
    for name in ('config', 'source', 'qualification'):
        pins.add(state[name]['path'], state[name])
    for name, key in [('config','outer_config'), ('source','wrapper_source')]:
        same(state[name]['sha256'], binding['inputs'][key]['sha256'], 'launch/worker '+name+' differs')
        same(state[name]['path'], binding['inputs'][key]['path'], 'launch/worker '+name+' path differs')
    controls = pins.read(state['qualification']['path'])
    same(controls['status'], 'complete', 'guard qualification incomplete')
    same(controls['tests_run'], 21, 'guard qualification count differs')
    same(controls['skipped'], [], 'guard qualification skipped')
    for field in ('errors','failures','exit_code'):
        same(controls[field], 0, 'guard qualification '+field+' differs')
    expected_sources = ['scripts/launch_secondary_depth_v2.py','scripts/qualify_secondary_depth_guard_v2.py','scripts/vendor/guard_helpers.py','tests/test_secondary_depth_guard_v2.py']
    same(sorted(controls['source_files']), sorted(expected_sources), 'guard qualification membership differs')
    workflow = Path(state['source']['path']).resolve().parents[1]
    for name, pin in controls['source_files'].items():
        pins.add(workflow / name, pin)
    same(controls['source_files']['scripts/launch_secondary_depth_v2.py']['sha256'], state['source']['sha256'], 'guard executed source qualification differs')
    outer = pins.read(state['config']['path'])
    same(outer['producer_configuration'], config, 'outer/inner configuration differs')
    same(outer['producer_config_sha256'], CONFIG_SHA, 'outer config SHA differs')
    pins.add(run / 'status.json', state['producer_status'])
    for name in ('status', 'summary'):
        pins.add(run / (name+'.json'), binding['producer_results'][name])
    return state


def audit(root, out, source, expected_launch_sha):
    import pysam
    root, out = root.resolve(), out.resolve()
    method = root / 'methods/coding-depth-20260920'
    run = root / RUN_REL
    pins = Pins()
    pins.add(source)
    config = pins.read(root / CONFIG_REL, {'sha256': CONFIG_SHA})
    same(config['depth_modes'], MODES, 'depth modes differ')
    bind_launch(root, pins, config, expected_launch_sha)
    qualification = pins.read(method / 'qualification-v2/controls.json', {'sha256': QUALIFICATION_SHA})
    require(qualification['status'] == 'complete' and type(qualification['exit_code']) is int and qualification['exit_code'] == 0, 'producer qualification not complete')
    for name, digest in qualification['source_sha256'].items():
        pins.add(method / 'qualification-v2' / name, {'sha256': digest})
    for name, digest in config['source_sha256'].items():
        same(digest, qualification['source_sha256']['workflow/scripts/' + name], 'config/source qualification differs')
    same(sorted(config['source_sha256']), sorted(SOURCE_NAMES), 'producer source membership differs')
    public, public_summary = public_bindings(root, pins, config)
    state = pins.read(run / 'status.json')
    same(state['status'], 'complete', 'producer not complete')
    summary = pins.read(run / 'summary.json', {'sha256': state['summary_sha256']})
    same(summary['status'], 'complete', 'summary not complete')
    same(state['inputs'], summary['inputs'], 'status/summary input bindings differ')
    for name, pin in summary['inputs'].items():
        pins.add(pin['path'], pin)
        if name in config['inputs']:
            pins.add(pin['path'], config['inputs'][name])
        elif name == 'config':
            pins.add(root / CONFIG_REL, pin)
        else:
            require(name.startswith('source/'), 'unexpected scientific input')
            short = name.split('/', 1)[1]
            pins.add(method / 'qualification-v2/workflow/scripts' / short, pin)
            same(pin['sha256'], config['source_sha256'][short], 'executed source pin differs')
    same(sorted(summary['inputs']), sorted(list(config['inputs']) + ['config'] + ['source/' + n for n in config['source_sha256']]), 'scientific input membership differs')
    expected_outputs = {'selection.json', 'regions.bed', 'selected.bam', 'selected.bam.bai', 'commands.json',
        'extraction-counts.json', 'extract-preserving-records.log', 'quickcheck.log', 'index.log', 'summary.json'}
    for mode in MODES:
        expected_outputs.update(mode['name'] + suffix for suffix in ('.json', '.depth.tsv.gz', '.log'))
    same(sorted(state['outputs']), sorted(expected_outputs), 'producer output membership differs')
    for name, pin in state['outputs'].items():
        pins.add(run / name, pin)
    paths = {n: Path(config['inputs'][n]['path']).resolve() for n in config['inputs']}
    require(str(paths['reference_index']) == str(paths['reference']) + '.fai', 'reference index binding differs')
    require(str(paths['alignment_index']) == str(paths['alignment']) + ('.crai' if paths['alignment'].suffix == '.cram' else '.bai'), 'alignment index binding differs')
    with pysam.FastaFile(str(paths['reference'])) as fasta:
        lengths = dict(zip(fasta.references, fasta.lengths))
        design = independent_targets(paths['gff'], config['genes'], lengths)
        sequences = [fasta.fetch(w['chrom'], w['start'], w['end']).upper() for w in design['windows']]
        reconcile_public_design(public, public_summary, design, sequences)
    with pysam.AlignmentFile(str(run / 'selected.bam')) as alignment:
        same(list(zip(alignment.references, alignment.lengths)), list(lengths.items()), 'alignment/reference dictionaries differ')
    selection = pins.read(run / 'selection.json')
    same(selection['config'], config, 'selection configuration differs')
    for key, value in design.items():
        same(selection[key], value, 'independent target ' + key + ' differs')
    for key in ('panel_union_bases', 'gene_union_bases_sum'):
        same(summary[key], design[key], 'summary denominator differs')
    gene_fields = ('id', 'name', 'chrom', 'strand', 'CDS_union_bases', 'transcript_ids')
    same(summary['genes'], [{k: g[k] for k in gene_fields} for g in design['genes']], 'summary gene identities differ')
    bed_text = ''.join(f"{w['chrom']}\t{w['start']}\t{w['end']}\n" for w in design['windows'])
    same((run / 'regions.bed').read_text(), bed_text, 'selected BED differs')
    own_bed = out / 'independent-regions.bed'
    own_bed.write_text(bed_text)
    extraction = compare_extraction(paths['alignment'], run / 'selected.bam', paths['reference'], paths['alignment_index'], design, config)
    same(summary['read_records_by_group'], extraction['groups'], 'RG read counts differ')
    same(summary['flags'], extraction['flags'], 'flag counts differ')
    same(summary['eligible_flag_mapping_qualities'], extraction['eligible_MQ'], 'eligible MQ counts differ')
    commands = pins.read(run / 'commands.json')
    same([r['label'] for r in commands], ['extract-preserving-records', 'quickcheck', 'index'] + [m['name'] for m in MODES], 'command membership differs')
    require(all(type(r['exit_code']) is int and r['exit_code'] == 0 for r in commands), 'producer command failure')
    modes, recounts = {}, []
    for mode in MODES:
        name = mode['name']
        expected_command = [str(paths['samtools']), 'depth', '-aa', '-b', str(run / 'regions.bed'), '-q', str(mode['minimum_BQ']), '-Q', str(mode['minimum_MQ']), '-G', '0xF04', '-s', str(run / 'selected.bam')]
        same(next(r['command'] for r in commands if r['label'] == name), expected_command, 'producer depth command differs')
        with gzip.open(run / (name + '.depth.tsv.gz'), 'rt') as stream:
            depths = positions_from_lines(stream, design['windows'])
        rows, groups = aggregate(depths, design, sequences)
        detail = pins.read(run / (name + '.json'))
        expected = dict(mode=mode, position_count=len(depths), groups=groups, windows=rows)
        same(detail, expected, 'complete per-window detail differs: ' + name)
        same(summary['depth_modes'][name], {k: expected[k] for k in ('mode', 'position_count', 'groups')}, 'complete summary groups differ: ' + name)
        command = expected_command.copy()
        command[4] = str(own_bed)
        log_path = out / (name + '.recount.stderr.txt')
        with log_path.open('x') as log:
            process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=log, text=True)
            try:
                recounted = positions_from_lines(process.stdout, design['windows'])
            except BaseException:
                process.terminate()
                process.wait()
                raise
            finally:
                process.stdout.close()
            code = process.wait()
        require(code == 0, 'same-engine technical recount failed')
        require(recounted == depths, 'same-engine per-position recount differs')
        recounts.append(dict(command=command, exit_code=code, all_positions_equal=True, positions=len(depths)))
        modes[name] = dict(position_count=len(depths), groups=groups, all_window_details_reconciled=True)
    same(sorted(summary['depth_modes']), sorted(m['name'] for m in MODES), 'summary mode membership differs')
    pins.stable()
    return dict(status='complete', completed_at=datetime.datetime.now(datetime.timezone.utc).isoformat(), host=platform.node(),
        producer_status=pins.add(run / 'status.json'), producer_summary=pins.add(run / 'summary.json'),
        reference_only_target_reconstruction=dict(genes=len(design['genes']), atomic_segments=len(design['windows']),
            panel_union_bases=design['panel_union_bases'], gene_union_bases_sum=design['gene_union_bases_sum']),
        extraction=extraction, depth_modes=modes, same_engine_technical_recounts=recounts,
        verified_files=pins.files, all_verified_inputs_stable=True,
        public_target_reconciled=True, producer_configuration_sha256=CONFIG_SHA,
        limitations=['Independent GFF/base membership, reference denominator, decoded-record extraction, and descriptive bookkeeping checks.',
            'Depth recount uses the identical pinned samtools engine and mate-overlap suppression policy; it is not an independent depth algorithm or biological validation.',
            'No genotype, locus-origin, phase, clinical callability, sensitivity, negative genetic finding, or independent-molecule inference.'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--expected-launch-status-sha256', required=True)
    args = parser.parse_args()
    root, out = args.root.resolve(), args.output.resolve()
    allowed = root / 'results/secondary-panel-cds-depth-20260920/independent-verification'
    require(out.is_relative_to(allowed) and out != allowed and not out.exists(), 'new owned audit subdirectory required')
    out.mkdir(parents=True)
    state = dict(status='running', started_at=datetime.datetime.now(datetime.timezone.utc).isoformat())
    def save(name, value):
        (out / name).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')
    save('status.json', state)
    try:
        report = audit(root, out, Path(__file__), args.expected_launch_status_sha256)
        save('verification.json', report)
        state.update(status='complete', completed_at=report['completed_at'],
            verification_sha256=hashlib.sha256((out / 'verification.json').read_bytes()).hexdigest())
        save('status.json', state)
    except BaseException as error:
        state.update(status='failed', error_type=type(error).__name__, error=str(error))
        save('status.json', state)
        raise
    print(json.dumps(state, sort_keys=True))


if __name__ == '__main__':
    main()
