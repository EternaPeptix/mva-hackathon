#!/usr/bin/env python3
"""Prepare explicit primary-nuclear DELLY query exclusions for the exact reference.

The input template already contains identical chr-prefixed and unprefixed primary
entries. Retain exact matching entries and exclude every other reference contig.
Query coverage is not hard read filtering, final breakpoint filtering or callability.
"""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import subprocess
import time

import pysam

PRIMARY = list(map(str, range(1, 23))) + ['X', 'Y']
REFERENCE_SHA = '187917c00fa7373a054f9768009d19e4f6e5fc87d3a6c125ce131164ca69fcd6'
TEMPLATE_SHA = 'caef8593f82f513694ae41b429680eab8abe39c057a119dae7760715aef006ee'
PARSER_UTIL_SHA = 'ae492aafea144c7746da230f60141aaca6a5a9ff78c27264aa5fb27686a77554'


def stamp(): return datetime.datetime.now(datetime.timezone.utc).isoformat()


def save(path, value): path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def fingerprint(path):
    path = path.resolve(); before = path.stat(); digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''): digest.update(block)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError('input changed while hashing')
    return {'path': str(path), 'bytes': after.st_size, 'mtime_ns': after.st_mtime_ns, 'sha256': digest.hexdigest()}


def merge(intervals):
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]: merged[-1][1] = max(end, merged[-1][1])
        else: merged.append([start, end])
    return merged


def native_query_intervals(length, excluded, whole=False):
    """Mirror the pinned, independently executed parser's actual iterator bounds."""
    if whole: return []
    result = []; start = 0
    for left, right in merge(excluded):
        if start + 1 < left: result.append([start, left - 1])
        start = right
    if start + 1 < length: result.append([start, length])
    return result


def parse_exclusions(path, dictionary):
    intervals = {chrom: [] for chrom in dictionary}; whole = set()
    for number, line in enumerate(path.read_text().splitlines(), 1):
        fields = line.split('\t')
        if len(fields) not in (1, 3) or fields[0] not in dictionary:
            raise ValueError('unknown contig or malformed exclusion row at line %d' % number)
        chrom = fields[0]
        if len(fields) == 1: whole.add(chrom); continue
        try: start, end = map(int, fields[1:])
        except ValueError: raise ValueError('nonnumeric exclusion boundary') from None
        if not 0 <= start < end <= dictionary[chrom]: raise ValueError('exclusion boundary outside reference')
        intervals[chrom].append([start, end])
    return {chrom: merge(rows) for chrom, rows in intervals.items()}, whole


def query_rows(dictionary, intervals, whole):
    return [[chrom, start, end] for chrom, length in dictionary.items()
            for start, end in native_query_intervals(length, intervals[chrom], chrom in whole)]


def write_bed(path, rows):
    path.write_text(''.join('%s\t%d\t%d\n' % tuple(row) for row in rows))


def build(reference, template, parser, output):
    output.mkdir(parents=True, exist_ok=False); started = time.monotonic()
    state = {'status': 'running', 'started_at': stamp()}; save(output / 'status.json', state)
    try:
        inputs = {'reference': fingerprint(reference), 'reference_index': fingerprint(Path(str(reference) + '.fai')),
                  'template': fingerprint(template), 'parser': fingerprint(parser),
                  'parser_build': fingerprint(parser.parent / 'build.json'), 'builder': fingerprint(Path(__file__))}
        if inputs['reference']['sha256'] != REFERENCE_SHA or inputs['template']['sha256'] != TEMPLATE_SHA:
            raise ValueError('reference or publisher template differs from reviewed version')
        receipt = json.loads((parser.parent / 'build.json').read_text())
        if (receipt['source_util_sha256'] != PARSER_UTIL_SHA or receipt['exit_code'] != 0
                or receipt['binary_sha256'] != inputs['parser']['sha256']):
            raise ValueError('native parser build provenance mismatch')
        with pysam.FastaFile(str(reference)) as fasta: dictionary = dict(zip(fasta.references, fasta.lengths))
        if any(chrom not in dictionary for chrom in PRIMARY): raise ValueError('required primary nuclear contig absent')
        raw = [line.split() for line in template.read_text().splitlines() if line.strip()]
        selected = []; naming = []
        for chrom in PRIMARY + ['M']:
            bare = sorted(row[1:] for row in raw if row[0] == chrom)
            prefixed = sorted(row[1:] for row in raw if row[0] == 'chr' + chrom)
            if not bare or bare != prefixed: raise ValueError('publisher primary naming variants differ')
            naming.append({'contig': chrom, 'exact_name_rows': len(bare), 'chr_name_rows': len(prefixed),
                           'entries_identical': True})
            if chrom in PRIMARY:
                for fields in bare:
                    if len(fields) < 2: raise ValueError('whole primary chromosome excluded by template')
                    start, end = map(int, fields[:2])
                    if not 0 <= start < end <= dictionary[chrom]: raise ValueError('template primary boundary outside reference')
                    selected.append([chrom, start, end, ' '.join(fields[2:])])
        whole = {chrom for chrom in dictionary if chrom not in PRIMARY}
        native_path = output / 'native-exclusions.tsv'
        native_path.write_text(''.join('%s\t%d\t%d\n' % tuple(row[:3]) for row in selected)
                               + ''.join(chrom + '\n' for chrom in dictionary if chrom in whole))
        intervals, parsed_whole = parse_exclusions(native_path, dictionary)
        if whole != parsed_whole: raise ValueError('whole-contig exclusion accounting mismatch')
        expected = query_rows(dictionary, intervals, whole)
        args = [str(parser), str(native_path), str(reference) + '.fai']
        native = subprocess.run(args, capture_output=True, text=True)
        (output / 'native-parser.stdout').write_text(native.stdout)
        (output / 'native-parser.stderr').write_text(native.stderr)
        observed = [[fields[0], int(fields[1]), int(fields[2])] for line in native.stdout.splitlines()
                    if (fields := line.split('\t'))]
        if native.returncode or observed != expected: raise ValueError('native query intervals differ from explicit accounting')
        write_bed(output / 'query-intervals.bed', observed)
        intended = [[chrom, a, b] for chrom in dictionary if chrom in PRIMARY for a, b in intervals[chrom]]
        write_bed(output / 'primary-template-intervals.bed', intended)
        stats = []
        for chrom in PRIMARY:
            valid = native_query_intervals(dictionary[chrom], intervals[chrom])
            original_bases = sum(b - a for a, b in intervals[chrom])
            omitted = dictionary[chrom] - sum(b - a for a, b in valid)
            stats.append({'chrom': chrom, 'length': dictionary[chrom], 'merged_template_intervals': len(intervals[chrom]),
                          'template_excluded_bases': original_bases, 'native_query_intervals': len(valid),
                          'native_omitted_query_bases': omitted, 'extra_omitted_query_bases': omitted - original_bases})
        result = {'completed_at': stamp(), 'elapsed_seconds': time.monotonic() - started,
                  'scope': 'verified_reference_primary_nuclear_query_exclusions', 'chromosomes': PRIMARY,
                  'inputs': inputs, 'reference_dictionary': dictionary,
                  'template_naming_audit': naming, 'selected_template_rows': selected,
                  'original_template_row_count': len(raw),
                  'literal_template_matches': sum(row[0] in dictionary for row in raw),
                  'unmatched_template_rows': sum(row[0] not in dictionary for row in raw),
                  'whole_contigs_excluded': sorted(whole), 'whole_contig_count': len(whole),
                  'whole_contig_bases_excluded': sum(dictionary[chrom] for chrom in whole),
                  'primary_statistics': stats, 'native_parser_command': args,
                  'native_parser_exit_code': native.returncode, 'native_parser_matches_all_expected_intervals': True,
                  'outputs': {name: fingerprint(output / name) for name in
                              ('native-exclusions.tsv', 'query-intervals.bed', 'primary-template-intervals.bed')},
                  'limitations': ['Exact-name primary rows already match the reference; earlier aggregate-name interpretation was too broad.',
                                 'M and every nonprimary reference contig are excluded; primary-to-excluded-contig rearrangements are outside this scan.',
                                 'Template includes Y pseudoautosomal exclusions; raw X/Y diploid discovery genotypes require separate ploidy interpretation.',
                                 'Query intervals can retrieve reads overlapping their edges; they are not hard breakpoint/read filters or measured callability.',
                                 'No patient reads or clinical mask sensitivity is assessed.']}
        save(output / 'summary.json', result); state.update(status='complete', completed_at=result['completed_at']); save(output / 'status.json', state)
        return result
    except BaseException as error:
        state.update(status='failed', failed_at=stamp(), error_type=type(error).__name__); save(output / 'status.json', state); raise


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('reference', 'template', 'parser', 'output'): p.add_argument('--' + name, type=Path, required=True)
    a = p.parse_args(); r = build(a.reference, a.template, a.parser, a.output)
    print(json.dumps({'completed_at': r['completed_at'], 'elapsed_seconds': r['elapsed_seconds'],
                      'selected_template_rows': len(r['selected_template_rows']), 'whole_contig_count': r['whole_contig_count'],
                      'extra_omitted_query_bases': sum(row['extra_omitted_query_bases'] for row in r['primary_statistics'])}, indent=2))


if __name__ == '__main__': main()
