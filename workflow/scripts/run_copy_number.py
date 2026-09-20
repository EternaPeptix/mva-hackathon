#!/usr/bin/env python3
"""Autosomal, exploratory CNVpytor depth analysis with explicit input safeguards.

The input adapter counts mapped read midpoints, not covered bases or independent
DNA molecules. Full depth profiles and calls must stay outside the repository.
"""
import argparse
from collections import Counter
import datetime
import hashlib
import importlib.metadata
import json
import logging
import math
import os
from pathlib import Path
import time

# These settings affect this analysis process, not the active alignment/runtime.
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('MPLBACKEND', 'Agg')

import numpy as np
import pysam


EXCLUDE_FLAGS = 0xF04


def stamp():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def fingerprint(path):
    path = Path(path).resolve()
    before = path.stat()
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(block)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError('input changed during hashing')
    return {'path': str(path), 'bytes': after.st_size, 'mtime_ns': after.st_mtime_ns,
            'sha256': h.hexdigest()}


def clean(value):
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [clean(v) for v in value]
    if isinstance(value, np.generic):
        return clean(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return {'nonfinite': str(value)}
    return value


def save(path, value):
    path.write_text(json.dumps(clean(value), indent=2, allow_nan=False) + '\n')


def validate_definition(config, reference_dictionary):
    chroms = config['chromosomes']
    if len(chroms) < 2 or len(chroms) != len(set(chroms)):
        raise ValueError('at least two distinct autosomes required for normalization')
    if any(c not in set(map(str, range(1, 23))) or c not in reference_dictionary for c in chroms):
        raise ValueError('explicit no-chr autosomal selection required')
    bins = config['bin_sizes']
    if not bins or bins != sorted(set(bins)) or any(type(b) is not int or b < 100 or b % 100 for b in bins):
        raise ValueError('bin sizes must be distinct increasing multiples of 100')
    if config['scope'] not in ('synthetic', 'verified_four_lane_autosomal_discovery'):
        raise ValueError('unsupported input scope')
    if config['scope'] != 'synthetic' and chroms != list(map(str, range(1, 23))):
        raise ValueError('patient analysis requires all 22 autosomes')
    if config['mask_policy'] != 'unmasked_discovery_with_q0_and_reference_gap_metrics':
        raise ValueError('unreviewed mask policy')


def verify_completion(config, alignment, alignment_fingerprint, completion_state):
    if config['scope'] == 'synthetic':
        return {'scope': 'synthetic; no patient completion claim'}
    if completion_state is None:
        raise ValueError('verified final alignment completion receipt required')
    state = json.loads(completion_state.read_text())
    if (state.get('status') != 'complete' or state.get('expected_lanes') != 4
            or Path(state.get('cram_path', '')).resolve() != alignment.resolve()
            or state.get('cram_bytes') != alignment_fingerprint['bytes']
            or state.get('cram_sha256') != alignment_fingerprint['sha256']):
        raise ValueError('final four-lane completion receipt does not match input')
    return fingerprint(completion_state)


def count_chromosome(reader, chrom, length, read_groups):
    """Match CNVpytor midpoint bins, with unconditional flags and start-zero handling."""
    all_mapped = np.zeros(length // 100 + 1, dtype=np.uint64)
    positive_mq = np.zeros_like(all_mapped)
    counts = Counter()
    for read in reader.fetch(chrom):
        counts['observed_records'] += 1
        if read.flag & EXCLUDE_FLAGS:
            counts['excluded_flags'] += 1
            continue
        if not read.has_tag('RG') or read.get_tag('RG') not in read_groups:
            raise ValueError('eligible read has missing or unrecognized read group')
        start, end = read.reference_start, read.reference_end
        if start is None or end is None or not 0 <= start < end <= length:
            raise ValueError('eligible read has invalid aligned coordinates')
        if read.mapping_quality == 255:
            counts['excluded_unknown_MQ_255'] += 1
            continue
        if any(op == 3 for op, _ in read.cigartuples or []):
            counts['excluded_reference_skip'] += 1
            continue
        index = (start + end) // 200
        all_mapped[index] += 1
        if read.mapping_quality > 0:
            positive_mq[index] += 1
        counts['included_read_midpoints'] += 1
        counts['MQ_zero' if read.mapping_quality == 0 else 'MQ_positive'] += 1
        if read.mapping_quality < 20:
            counts['included_MQ_below_20'] += 1
        if start == 0:
            counts['included_reference_start_zero'] += 1
    if not counts['included_read_midpoints']:
        raise ValueError('selected chromosome has no eligible reads')
    if np.max(all_mapped) > np.iinfo(np.uint32).max:
        raise ValueError('count exceeds CNVpytor storage range')
    if int(all_mapped.sum()) != counts['included_read_midpoints'] or np.any(positive_mq > all_mapped):
        raise ValueError('read count conservation failed')
    return all_mapped, positive_mq, dict(counts)


def run_analysis(alignment, reference, config, output, completion_state=None):
    from cnvpytor import __version__
    from cnvpytor.root import Root
    from cnvpytor.genome import Genome
    from cnvpytor.io import FLAG_AUTO, FLAG_GC_CORR
    import cnvpytor

    if __version__ != '1.3.2':
        raise ValueError('unreviewed CNVpytor version')
    output.mkdir(parents=True, exist_ok=False)
    state = {'status': 'running', 'started_at': stamp(), 'pid': os.getpid(), 'stage': 'validate inputs'}
    save(output / 'status.json', state)
    started = time.monotonic()
    app = None
    old_reference_genomes, old_detected = Genome.reference_genomes, Genome.detected_genome
    try:
        inputs = {'alignment': fingerprint(alignment), 'reference': fingerprint(reference),
                  'reference_index': fingerprint(Path(str(reference) + '.fai')),
                  'runner': fingerprint(Path(__file__))}
        index = Path(str(alignment) + ('.crai' if alignment.suffix == '.cram' else '.bai'))
        inputs['alignment_index'] = fingerprint(index)
        if inputs['reference']['sha256'] != config['reference_sha256']:
            raise ValueError('reference checksum differs from configuration')
        completion = verify_completion(config, alignment, inputs['alignment'], completion_state)
        with pysam.FastaFile(str(reference)) as fasta, pysam.AlignmentFile(str(alignment),
                reference_filename=str(reference)) as reader:
            dictionary = dict(zip(fasta.references, fasta.lengths))
            validate_definition(config, dictionary)
            if dictionary != dict(zip(reader.references, reader.lengths)):
                raise ValueError('alignment and reference dictionaries differ')
            if not reader.has_index() or reader.header.to_dict().get('HD', {}).get('SO') != 'coordinate':
                raise ValueError('indexed coordinate-sorted input required')
            groups = reader.header.to_dict().get('RG', [])
            if (not groups or any(g.get('SM') != config['sample'] for g in groups)
                    or len({g['ID'] for g in groups}) != len(groups)):
                raise ValueError('sample metadata mismatch or duplicate read-group IDs')
            read_groups = {g['ID'] for g in groups}
            chroms = config['chromosomes']
            lengths = {c: dictionary[c] for c in chroms}
            # Explicit custom reference: use exact input FASTA GC, never bundled hg38 resources.
            Genome.reference_genomes = {'mva_exact_reference': {
                'name': 'mva_exact_reference', 'species': 'analysis scope only',
                'chromosomes': {c: (lengths[c], 'A') for c in chroms}}}
            Genome.detected_genome = 'mva_exact_reference'
            app = Root(str(output / 'depth.pytor'), create=True, max_cores=1)
            app.io.create_signal(None, None, 'reference genome', np.array([np.bytes_('mva_exact_reference')]))
            app.io.create_signal(None, None, 'use reference', np.array([0, 0], dtype=np.uint8))
            import_counts = {}
            for chrom in chroms:
                state.update(stage='read midpoint import', current_chromosome=chrom)
                save(output / 'status.json', state)
                parity, positive, counts = count_chromosome(reader, chrom, lengths[chrom], read_groups)
                app.io.save_rd(chrom, parity, positive, chromosome_length=lengths[chrom])
                actual_p, actual_u = app.io.read_rd(chrom)
                if not np.array_equal(parity, actual_p) or not np.array_equal(positive, actual_u):
                    raise ValueError('depth storage round-trip changed counts')
                import_counts[chrom] = counts
        state.update(stage='exact-reference GC and histograms'); save(output / 'status.json', state)
        app.gc(str(reference), chroms=chroms)
        bins = config['bin_sizes']
        app.calculate_histograms(bins, chroms=chroms, gc_auto=False, gc_corr_rm_ol=False)
        normalization = {}
        for size in bins:
            stat = np.asarray(app.io.get_signal(None, size, 'RD stat', FLAG_AUTO | FLAG_GC_CORR))
            if len(stat) != 6 or not np.all(np.isfinite(stat)) or min(stat[4:6]) <= 0:
                raise ValueError('missing or invalid autosomal depth normalization')
            normalization[size] = {'fitted_mean_read_midpoints_per_bin': stat[4], 'fitted_sd': stat[5]}
            for chrom in chroms:
                signal = app.io.get_signal(chrom, size, 'RD', FLAG_GC_CORR)
                if len(signal) < math.ceil(lengths[chrom] / size):
                    raise ValueError('corrected depth profile does not cover selected chromosome')
        state.update(stage='partition'); save(output / 'status.json', state)
        genome_size = sum(lengths.values())
        app.partition(bins, chroms=chroms, use_gc_corr=True, use_mask=False, genome_size=genome_size)
        state.update(stage='call relative depth changes'); save(output / 'status.json', state)
        calls = app.call(bins, chroms=chroms, print_calls=False, use_gc_corr=True,
                         use_mask=False, genome_size=genome_size, genome_cnv_fraction=0.01)
        if set(calls) != set(bins):
            raise ValueError('caller omitted a requested bin size')
        for size in bins:
            for chrom in chroms:
                if not app.io.signal_exists(chrom, size, 'calls', FLAG_GC_CORR):
                    raise ValueError('missing per-chromosome call completion signal')
        columns = ['type', 'chrom', 'start_1based', 'end_inclusive', 'size_bp', 'relative_depth',
                   'e1', 'e2', 'e3', 'e4', 'MQ_zero_fraction', 'reference_gap_fraction', 'distance_to_gap_bp']
        call_rows = {}
        for size, values in calls.items():
            call_rows[size] = []
            for values_row in values:
                if len(values_row) != len(columns):
                    raise ValueError('unrecognized call schema')
                row = dict(zip(columns, values_row))
                row['within_reference_bounds'] = 1 <= row['start_1based'] <= row['end_inclusive'] <= lengths[row['chrom']]
                call_rows[size].append(row)
        save(output / 'calls.json', {'scope': config['scope'], 'coordinates': '1-based inclusive, bin-level estimates',
             'interpretation': 'Unfiltered relative-depth discovery calls; not diagnoses or calibrated causal probabilities.',
             'calls': call_rows})
        app.io.file.flush(); app.io.file.close(); app = None
        for row in inputs.values():
            stat = Path(row['path']).stat()
            if (stat.st_size, stat.st_mtime_ns) != (row['bytes'], row['mtime_ns']):
                raise ValueError('input changed during analysis')
        package = Path(cnvpytor.__file__).parent
        result = {'completed_at': stamp(), 'elapsed_seconds': time.monotonic() - started,
            'config': config, 'inputs': inputs, 'alignment_completion': completion,
            'runtime': {p: importlib.metadata.version(p) for p in ('CNVpytor','numpy','scipy','pysam','h5py')},
            'package_source_hashes': {p.name: fingerprint(p)['sha256'] for p in sorted(package.glob('*.py'))},
            'read_midpoint_import': import_counts, 'normalization': normalization,
            'counts_round_trip_verified': True,
            'calls_by_bin_size': {size: dict(Counter(r['type'] for r in rows)) for size, rows in call_rows.items()},
            'out_of_reference_calls_by_bin_size': {size: sum(not r['within_reference_bounds'] for r in rows) for size, rows in call_rows.items()},
            'outputs': {name: fingerprint(output / name) for name in ('depth.pytor','calls.json')},
            'limitations': ['Read midpoint density is not base-level callability or independent molecule depth.',
                'Autosomal discovery only; baseline normalization assumes a dominant copy-neutral population.',
                'No strict accessibility mask, independent breakpoint evidence, BAF, contamination or clinical classification is supplied.',
                'MQ-zero reads retained for q0; the positive-MQ channel does not guarantee unique mapping.',
                'No balanced-SV, UPD, low-level mosaicism, sex/mitochondrial, repeat or methylation exclusion.',
                'No patient-negative or clinical sensitivity claim follows from empty calls or synthetic tests.']}
        save(output / 'summary.json', result)
        state.update(status='complete', stage='complete', completed_at=result['completed_at'],
                     summary_sha256=fingerprint(output / 'summary.json')['sha256'])
        save(output / 'status.json', state)
        return result
    except BaseException as error:
        state.update(status='failed', failed_at=stamp(), error_type=type(error).__name__)
        save(output / 'status.json', state)
        raise
    finally:
        if app is not None:
            app.io.file.close()
        Genome.reference_genomes, Genome.detected_genome = old_reference_genomes, old_detected


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--alignment', type=Path, required=True)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--completion-state', type=Path)
    args = parser.parse_args()
    os.nice(10)
    config = json.loads(args.config.read_text())
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
    result = run_analysis(args.alignment.resolve(), args.reference.resolve(), config,
                          args.output.resolve(), args.completion_state)
    print(json.dumps({'scope': config['scope'], 'elapsed_seconds': result['elapsed_seconds'],
                      'calls_by_bin_size': result['calls_by_bin_size']}))


if __name__ == '__main__':
    main()
