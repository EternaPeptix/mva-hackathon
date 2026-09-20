#!/usr/bin/env python3
"""Guarded native structural discovery on a completed alignment.

Retain the unmodified alignment, supplementary evidence, complete native BCF and
per-record review outside Git. Native diploid discovery is not clinical genotype,
phase, mosaicism, causal probability or comprehensive structural exclusion.
"""
import argparse
from collections import Counter
import json
import os
from pathlib import Path
import shlex
import subprocess
import time

import pysam

from build_structural_exclusions import (PRIMARY, fingerprint, native_query_intervals,
                                        parse_exclusions, query_rows, save, stamp)

PATIENT_SCOPE = 'verified_four_lane_primary_nuclear_structural_discovery'
NATIVE_OPTIONS = ['-t', 'ALL', '-q', '20', '-u', '20', '-r', '20', '-z', '3',
                  '-s', '9', '-c', '25', '-m', '25', '-n', '40', '-p', '20', '-a', '250']


def validate_config(config):
    if config['scope'] not in ('synthetic', PATIENT_SCOPE): raise ValueError('unreviewed structural scope')
    chroms = config['chromosomes']
    if not chroms or len(chroms) != len(set(chroms)) or any(chrom not in PRIMARY for chrom in chroms):
        raise ValueError('explicit distinct primary nuclear chromosomes required')
    if config['scope'] == PATIENT_SCOPE and chroms != PRIMARY:
        raise ValueError('patient discovery requires all 24 primary nuclear chromosomes')
    if config['scope'] == 'synthetic' and config['sample'] != 'SYNTH_SV':
        raise ValueError('synthetic mode requires the reserved synthetic sample')
    if config['read_policy'] != 'native_delly_sr_unmodified_alignment_reject_unknown_mapped_MQ':
        raise ValueError('unreviewed structural read policy')
    if config['ploidy_policy'] != 'native_diploid_discovery_sex_genotypes_uninterpreted':
        raise ValueError('unreviewed discovery ploidy policy')
    if type(config['threads']) is not int or not 1 <= config['threads'] <= 4:
        raise ValueError('invalid structural thread count')
    if config['native_options'] != NATIVE_OPTIONS: raise ValueError('unreviewed native options')


def verify_completion(config, alignment, identity, path):
    if config['scope'] == 'synthetic': return {'scope': 'synthetic; no patient completion claim'}
    if path is None: raise ValueError('verified final four-lane completion receipt required')
    receipt = fingerprint(path); state = json.loads(path.read_text())
    if (state.get('status') != 'complete' or state.get('expected_lanes') != 4 or alignment.suffix != '.cram'
            or Path(state.get('cram_path', '')).resolve() != alignment.resolve()
            or state.get('cram_bytes') != identity['bytes'] or state.get('cram_sha256') != identity['sha256']):
        raise ValueError('final four-lane completion receipt does not match input')
    return receipt


def bed_rows(path):
    rows = []
    for line in path.read_text().splitlines():
        fields = line.split('\t')
        if len(fields) != 3: raise ValueError('malformed query-interval sidecar')
        rows.append([fields[0], int(fields[1]), int(fields[2])])
    return rows


def load_mask(folder, dictionary, config, reference, reference_index):
    receipt = fingerprint(folder / 'summary.json')
    if receipt['sha256'] != config['exclusion_summary_sha256']: raise ValueError('exclusion-summary checksum mismatch')
    summary = json.loads((folder / 'summary.json').read_text())
    expected_scope = 'synthetic_query_exclusions' if config['scope'] == 'synthetic' else 'verified_reference_primary_nuclear_query_exclusions'
    if summary['scope'] != expected_scope or summary['chromosomes'] != config['chromosomes']:
        raise ValueError('exclusion scope differs from discovery scope')
    if summary['reference_dictionary'] != dictionary: raise ValueError('exclusion reference dictionary mismatch')
    if (summary['inputs']['reference']['sha256'] != reference['sha256']
            or summary['inputs']['reference_index']['sha256'] != reference_index['sha256']):
        raise ValueError('exclusions use a different reference/index')
    for name in ('native-exclusions.tsv', 'query-intervals.bed', 'primary-template-intervals.bed'):
        actual = fingerprint(folder / name); expected = summary['outputs'][name]
        if actual['sha256'] != expected['sha256'] or actual['bytes'] != expected['bytes']:
            raise ValueError('exclusion sidecar changed')
    intervals, whole = parse_exclusions(folder / 'native-exclusions.tsv', dictionary)
    if whole != set(dictionary) - set(config['chromosomes']): raise ValueError('whole-contig exclusions do not match scope')
    if query_rows(dictionary, intervals, whole) != bed_rows(folder / 'query-intervals.bed'):
        raise ValueError('query intervals disagree with native parser semantics')
    expected = [[chrom, a, b] for chrom in dictionary if chrom in config['chromosomes'] for a, b in intervals[chrom]]
    if expected != bed_rows(folder / 'primary-template-intervals.bed'):
        raise ValueError('template intervals disagree with native exclusions')
    if config['scope'] == PATIENT_SCOPE and summary.get('native_parser_matches_all_expected_intervals') is not True:
        raise ValueError('patient exclusions lack native boundary verification')
    return receipt, intervals, whole


def alignment_census(reader, groups):
    counts = Counter(); by_rg = Counter(); mq = Counter(); last = (-1, -2); total_contigs = len(reader.references)
    for read in reader:
        counts['alignment_records'] += 1
        key = (read.reference_id if read.reference_id >= 0 else total_contigs, read.reference_start)
        if key < last: raise ValueError('actual alignment order is not coordinate sorted')
        last = key
        if not read.has_tag('RG') or read.get_tag('RG') not in groups: raise ValueError('read has missing or undeclared read group')
        by_rg[read.get_tag('RG')] += 1
        if read.is_unmapped: counts['unmapped_records'] += 1
        else:
            counts['mapped_records'] += 1; mq[str(read.mapping_quality)] += 1
            if read.mapping_quality == 255: counts['mapped_unknown_MQ_255'] += 1
        if read.is_duplicate: counts['duplicate_records'] += 1
        if read.is_qcfail: counts['QC_failed_records'] += 1
        if read.is_secondary: counts['secondary_records'] += 1
        if read.is_supplementary: counts['supplementary_records'] += 1
        if not read.flag & 0x900: counts['primary_read_records'] += 1
    return {'counts': dict(counts), 'records_by_read_group': dict(by_rg), 'mapped_MQ_histogram': dict(mq)}


def in_interval(position, intervals):
    return any(start <= position < end for start, end in intervals)


def review_calls(path, output, sample, dictionary, chroms, intervals, whole):
    counts = Counter(); types = Counter(); filters = Counter(); genotype_labels = Counter(); ids = set()
    with pysam.VariantFile(str(path)) as reader, output.open('w') as writer:
        samples = list(reader.header.samples)
        if samples not in ([], [sample]): raise ValueError('native caller sample identity mismatch')
        for record in reader:
            if not samples: raise ValueError('nonempty native callset has no sample column')
            info = dict(record.info)
            kind = info.get('SVTYPE'); partner = info.get('CHR2', record.contig)
            end = info.get('POS2') if kind == 'BND' else record.stop
            if kind not in ('DEL', 'DUP', 'INV', 'INS', 'BND'): raise ValueError('unexpected structural record type')
            if record.contig not in chroms or partner not in chroms: raise ValueError('native breakpoint outside declared chromosome scope')
            if end is None or not 0 <= record.pos <= dictionary[record.contig] + 1 or not 0 <= end <= dictionary[partner] + 1:
                raise ValueError('native endpoint missing or outside reference/telomere positions')
            if not record.id or record.id in ids: raise ValueError('missing or duplicate structural record ID')
            ids.add(record.id); counts['records'] += 1; types[kind] += 1
            filter_label = ';'.join(record.filter) or '.'; filters[filter_label] += 1
            gt = record.samples[sample].get('GT'); label = '/'.join('.' if x is None else str(x) for x in gt) if gt else '.'
            genotype_labels[label] += 1
            queries_a = native_query_intervals(dictionary[record.contig], intervals[record.contig], record.contig in whole)
            queries_b = native_query_intervals(dictionary[partner], intervals[partner], partner in whole)
            row = {'id': record.id, 'chrom': record.contig, 'pos_1based': record.pos, 'ref': record.ref,
                   'alt': list(record.alts or ()), 'type': kind, 'end_chrom': partner, 'end_position': end,
                   'connection_type': info.get('CT'), 'length': info.get('SVLEN'),
                   'filter': list(record.filter), 'qual': record.qual, 'precise': info.get('PRECISE', False),
                   'evidence': {key: info.get(key) for key in ('PE','SR','MAPQ','SRMAPQ','CIPOS','CIEND')},
                   'genotype': dict(record.samples[sample]), 'genotype_phased': record.samples[sample].phased,
                   'native_diploid_sex_genotype_uninterpreted': record.contig in ('X','Y') or partner in ('X','Y'),
                   'endpoint_reference_bases_in_query_intervals': [in_interval(record.pos - 1, queries_a), in_interval(end - 1, queries_b)],
                   'coordinate_note': 'Native positions retained. Endpoint-adjacent reference-base flags use POS-1/POS2-1, not normalized interbase breakpoints or hard exclusion.'}
            if row['native_diploid_sex_genotype_uninterpreted']: counts['records_involving_sex_chromosomes'] += 1
            if not all(row['endpoint_reference_bases_in_query_intervals']): counts['records_with_endpoint_base_outside_query_intervals'] += 1
            writer.write(json.dumps(row, allow_nan=False) + '\n')
    return {'counts': dict(counts), 'types': dict(types), 'filters': dict(filters), 'raw_diploid_genotypes': dict(genotype_labels),
            'native_header_samples': samples, 'native_sample_column_status': 'present' if samples else 'absent_empty_callset'}


def run_analysis(alignment, reference, exclusions, config, delly, samtools, output, completion_state=None):
    validate_config(config)
    if config['scope'] == PATIENT_SCOPE and output.resolve().is_relative_to(Path(__file__).resolve().parents[2]):
        raise ValueError('patient structural output must remain outside Git')
    output.mkdir(parents=True, exist_ok=False); start = time.monotonic(); commands = []
    state = {'status': 'running', 'pid': os.getpid(), 'started_at': stamp(), 'stage': 'verify completed input'}
    save(output / 'status.json', state)
    try:
        if config['scope'] == PATIENT_SCOPE:
            if completion_state is None: raise ValueError('verified final four-lane completion receipt required')
            completion = json.loads(completion_state.read_text())
            if completion.get('status') != 'complete' or completion.get('expected_lanes') != 4:
                raise ValueError('patient alignment is not complete')
        inputs = {'alignment': fingerprint(alignment), 'reference': fingerprint(reference),
                  'reference_index': fingerprint(Path(str(reference) + '.fai')),
                  'alignment_index': fingerprint(Path(str(alignment) + ('.crai' if alignment.suffix == '.cram' else '.bai'))),
                  'delly': fingerprint(delly), 'samtools': fingerprint(samtools), 'runner': fingerprint(Path(__file__)),
                  'exclusion_code': fingerprint(Path(__file__).with_name('build_structural_exclusions.py'))}
        for key in ('reference', 'reference_index', 'delly', 'samtools'):
            if inputs[key]['sha256'] != config[key + '_sha256']: raise ValueError(key + ' checksum mismatch')
        completion_receipt = verify_completion(config, alignment, inputs['alignment'], completion_state)
        def command(stage, args):
            state.update(stage=stage); save(output / 'status.json', state)
            log = output / ('step-%02d.log' % len(commands)); tick = time.monotonic()
            with log.open('w') as stream: result = subprocess.run(list(map(str, args)), stdout=stream, stderr=subprocess.STDOUT)
            commands.append({'stage':stage, 'command':list(map(str,args)), 'exit_code':result.returncode,
                             'elapsed_seconds':time.monotonic()-tick, 'log_sha256':fingerprint(log)['sha256']})
            save(output / 'commands.json', commands)
            if result.returncode: raise RuntimeError('native structural stage failed: ' + stage)
        command('alignment format check', [samtools, 'quickcheck', '-v', alignment])
        with pysam.FastaFile(str(reference)) as fasta, pysam.AlignmentFile(str(alignment), reference_filename=str(reference)) as reader:
            dictionary = dict(zip(fasta.references, fasta.lengths))
            if list(dictionary.items()) != list(zip(reader.references, reader.lengths)):
                raise ValueError('alignment/reference dictionary or order differs')
            header = reader.header.to_dict(); groups = header.get('RG', [])
            if not reader.has_index() or header.get('HD', {}).get('SO') != 'coordinate': raise ValueError('indexed coordinate-sorted alignment required')
            if (not groups or any(not g.get('ID') or g.get('SM') != config['sample'] for g in groups)
                    or len({g['ID'] for g in groups}) != len(groups)):
                raise ValueError('alignment sample/read-group metadata mismatch')
            if config['scope'] == PATIENT_SCOPE and len(groups) != 4: raise ValueError('four read groups required for patient input')
            if not any('markdup' in shlex.split(pg.get('CL', '')) for pg in header.get('PG', [])):
                raise ValueError('duplicate-marking provenance absent from alignment header')
            mask_receipt, intervals, whole = load_mask(exclusions, dictionary, config, inputs['reference'], inputs['reference_index'])
            state.update(stage='complete read-metadata census; alignment unchanged'); save(output / 'status.json', state)
            census = alignment_census(reader, {g['ID'] for g in groups}); save(output / 'alignment-census.json', census)
        if census['counts'].get('mapped_unknown_MQ_255'): raise ValueError('mapped reads have unknown MQ 255; explicit review required')
        if not census['counts'].get('alignment_records'): raise ValueError('alignment has no read records')
        if set(census['records_by_read_group']) != {g['ID'] for g in groups}: raise ValueError('declared read group has no observed records')
        command('native structural discovery', [delly, 'sr', '-h', str(config['threads']), *config['native_options'],
                                               '-x', exclusions / 'native-exclusions.tsv', '-g', reference,
                                               '-o', output / 'calls.bcf', alignment])
        if not (output / 'calls.bcf.csi').is_file(): raise ValueError('native structural index missing')
        reviewed = review_calls(output / 'calls.bcf', output / 'calls-review.jsonl', config['sample'],
                                dictionary, config['chromosomes'], intervals, whole)
        for key, identity in inputs.items():
            after = Path(identity['path']).stat()
            if (after.st_size, after.st_mtime_ns) != (identity['bytes'], identity['mtime_ns']): raise ValueError('input changed during structural analysis: ' + key)
        if fingerprint(exclusions / 'summary.json')['sha256'] != mask_receipt['sha256']: raise ValueError('exclusion receipt changed during analysis')
        load_mask(exclusions, dictionary, config, inputs['reference'], inputs['reference_index'])
        result = {'completed_at': stamp(), 'elapsed_seconds': time.monotonic()-start, 'scope': config['scope'],
                  'config': config, 'inputs': inputs, 'alignment_completion': completion_receipt,
                  'exclusion_summary': mask_receipt, 'commands': commands, 'alignment_census': census,
                  'native_discovery': reviewed, 'analysis_status': 'complete_with_calls' if reviewed['counts'].get('records') else 'complete_without_calls',
                  'outputs': {name:fingerprint(output/name) for name in ('calls.bcf','calls.bcf.csi','calls-review.jsonl','alignment-census.json')},
                  'limitations': ['Exploratory native diploid discovery; no genotype/phase/clinical/EPCR assignment.',
                                 'X/Y genotypes remain uninterpreted; M, nonprimary contigs and their cross-contig rearrangements are outside this scope.',
                                 'Query exclusions do not prove final breakpoint exclusion or measured callability; raw records and endpoint context are retained.',
                                 'Synthetic insertion miss and real-read sensitivity remain material limits; no-call output is not a negative clinical test.',
                                 'Native flag/MQ policies apply to the unmodified alignment, including supplementary evidence; input QC and downstream interpretation remain separate.']}
        save(output/'summary.json', result);state.update(status='complete',stage='complete',completed_at=result['completed_at']);save(output/'status.json',state)
        return result
    except BaseException as error:
        state.update(status='failed',failed_at=stamp(),error_type=type(error).__name__);save(output/'status.json',state);raise


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('alignment','reference','exclusions','config','delly','samtools','output'): p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--completion-state',type=Path)
    a=p.parse_args();os.nice(10);r=run_analysis(a.alignment,a.reference,a.exclusions,json.loads(a.config.read_text()),a.delly,a.samtools,a.output,a.completion_state)
    print(json.dumps({'completed_at':r['completed_at'],'elapsed_seconds':r['elapsed_seconds'],'analysis_status':r['analysis_status'],
                      'counts':r['native_discovery']['counts'],'scope':r['scope']},indent=2))


if __name__ == '__main__': main()
