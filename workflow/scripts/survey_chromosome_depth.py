#!/usr/bin/env python3
"""Fixed-region depth survey for copy-number assumptions, not a clinical karyotype.

Reference-only windows are fixed before inspecting reads. Missing depth output is
an error; explicit zero depth is retained. No variant calls or source edits occur.
"""
import argparse
from collections import Counter, defaultdict
import gzip
import json
import os
from pathlib import Path
import statistics
import subprocess
import time

import pysam

from prepare_gatk_compatible_cram import digest, save, stamp


def grid_windows(chrom, start, end, count, width, category):
    if count < 1 or width < 1 or start < 0 or end - start < count * width:
        raise ValueError('invalid or overlapping grid design')
    rows = []
    for i in range(count):
        center = start + (2 * i + 1) * (end - start) // (2 * count)
        left = center - width // 2
        rows.append({'chrom': chrom, 'start': left, 'end': left + width, 'category': category})
    return rows


def prepare_windows(fasta, config):
    windows = []
    for chrom in map(str, range(1, 23)):
        windows.extend(grid_windows(chrom, 0, fasta.get_reference_length(chrom),
            config['autosome_windows_per_chromosome'], config['window_bp'], 'autosome'))
    par_proof = []
    for chrom in ('X', 'Y'):
        intervals = sorted(config['PAR_1based_inclusive'][chrom], key=lambda x: x['start'])
        if len(intervals) != 2 or intervals[0]['end'] >= intervals[1]['start']:
            raise ValueError('require ordered nonoverlapping PAR1 and PAR2')
        windows.extend(grid_windows(chrom, intervals[0]['end'], intervals[1]['start'] - 1,
            config['nonPAR_windows_per_sex_chromosome'], config['window_bp'], 'nonPAR_' + chrom))
        for row in intervals:
            windows.extend(grid_windows(chrom, row['start'] - 1, row['end'],
                config['windows_per_PAR_copy'], config['window_bp'], row['name'] + '_' + chrom))
            sequence = fasta.fetch(chrom, row['start'] - 1, row['end']).upper()
            import hashlib
            par_proof.append({'chrom': chrom, 'name': row['name'], 'start_1based': row['start'],
                'end_1based': row['end'], 'length': len(sequence),
                'ACGT_bases': sum(sequence.count(base) for base in 'ACGT'),
                'sequence_sha256': hashlib.sha256(sequence.encode()).hexdigest()})
    windows.append({'chrom': 'MT', 'start': 0, 'end': fasta.get_reference_length('MT'), 'category': 'MT'})
    order = {chrom: i for i, chrom in enumerate(fasta.references)}
    windows.sort(key=lambda w: (order[w['chrom']], w['start']))
    for i, row in enumerate(windows):
        if i and windows[i - 1]['chrom'] == row['chrom'] and windows[i - 1]['end'] > row['start']:
            raise ValueError('overlapping target windows')
        sequence = fasta.fetch(row['chrom'], row['start'], row['end']).upper()
        n = sum(sequence.count(base) for base in 'ACGT')
        gc = (sequence.count('G') + sequence.count('C')) / n if n else None
        row.update(index=i, span=len(sequence), ACGT_bases=n, GC_fraction=gc,
                   reference_eligible=n / len(sequence) >= config['minimum_reference_ACGT_fraction'] and
                   gc is not None and config['reference_GC_bounds'][0] <= gc <= config['reference_GC_bounds'][1])
    return windows, par_proof


def collect_depth(lines, windows, sequences, raw_output=None):
    if len(windows) != len(sequences) or any(len(s) != w['end'] - w['start'] for w, s in zip(windows, sequences)):
        raise ValueError('sequence/window lengths differ')
    rows, line_count = [], 0
    stream = iter(lines)
    for window, sequence in zip(windows, sequences):
        all_depth, acgt_depth = Counter(), Counter()
        for position in range(window['start'], window['end']):
            try:
                raw = next(stream)
            except StopIteration as error:
                raise ValueError('missing explicit depth positions') from error
            fields = raw.split()
            if len(fields) != 3 or fields[0] != window['chrom'] or int(fields[1]) != position + 1:
                raise ValueError('depth order/position/column mismatch')
            depth = int(fields[2])
            if depth < 0:
                raise ValueError('negative depth')
            if raw_output is not None:
                raw_output.write(raw)
            all_depth[depth] += 1
            if sequence[position - window['start']] in 'ACGT':
                acgt_depth[depth] += 1
            line_count += 1
        row = dict(window)
        for label, histogram in [('span', all_depth), ('ACGT', acgt_depth)]:
            total = sum(histogram.values())
            row[label + '_depth'] = {'bases': total, 'sum_depth': sum(d * n for d, n in histogram.items()),
                'mean': sum(d * n for d, n in histogram.items()) / total if total else None,
                'bases_at_least': {str(t): sum(n for d, n in histogram.items() if d >= t) for t in (1, 10, 20, 30)},
                'histogram': dict(sorted(histogram.items()))}
        rows.append(row)
    if next(stream, None) is not None:
        raise ValueError('unexpected depth outside selected windows')
    return rows, line_count


def summarize(rows):
    groups = defaultdict(list)
    for row in rows:
        groups[row['category']].append(row)
    result = {}
    autosomal = [r['ACGT_depth']['mean'] for r in groups['autosome'] if r['reference_eligible']]
    denominator = statistics.median(autosomal) if autosomal else None
    for category, values in groups.items():
        eligible = [r for r in values if r['reference_eligible']]
        means = [r['ACGT_depth']['mean'] for r in eligible]
        median = statistics.median(means) if means else None
        result[category] = {'selected_windows': len(values), 'reference_eligible_windows': len(eligible),
            'reference_ineligible_windows': len(values) - len(eligible),
            'span_bases': sum(r['span'] for r in values), 'ACGT_bases': sum(r['ACGT_bases'] for r in values),
            'median_reference_eligible_window_mean': median,
            'median_over_autosomal_median': median / denominator if median is not None and denominator else None,
            'reference_eligible_window_mean_min': min(means) if means else None,
            'reference_eligible_window_mean_max': max(means) if means else None}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('root', 'alignment-state', 'config', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    os.nice(10)
    root, out = args.root.resolve(), args.output.resolve()
    if out.exists():
        raise ValueError('existing survey preserved')
    started = time.monotonic()
    alignment_state = json.loads(args.alignment_state.read_text())
    if alignment_state['status'] != 'complete' or alignment_state['expected_lanes'] != 4:
        raise ValueError('completed four-lane alignment required')
    source = Path(alignment_state['cram_path'])
    if source.stat().st_size != alignment_state['cram_bytes'] or digest(source) != alignment_state['cram_sha256']:
        raise ValueError('alignment differs from its completed identity')
    config = json.loads(args.config.read_text())
    reference = root / 'ref/realignment/GRCh38_nochr.fa'
    if digest(reference) != config['reference_sha256']:
        raise ValueError('reference differs')
    with pysam.FastaFile(str(reference)) as fasta:
        windows, par_proof = prepare_windows(fasta, config)
        sequences = [fasta.fetch(w['chrom'], w['start'], w['end']).upper() for w in windows]
    out.mkdir(parents=True, exist_ok=False)
    samtools = root / 'tools/bio/bin/samtools'
    state = {'status': 'running', 'pid': os.getpid(), 'started_at': stamp(),
             'source_cram_sha256': alignment_state['cram_sha256'], 'reference_sha256': config['reference_sha256'],
             'source_and_tool_sha256': {str(p): digest(p) for p in (Path(__file__), args.config, args.alignment_state,
                 samtools, Path(str(source) + '.crai'), Path(str(reference) + '.fai'),
                 Path(__file__).with_name('prepare_gatk_compatible_cram.py'))}}
    commands = []

    def command_record(command, label):
        state.update(current_stage=label)
        save(out / 'status.json', state)
        record = {'label': label, 'command': list(map(str, command)), 'started_at': stamp()}
        commands.append(record)
        save(out / 'commands.json', commands)
        return record

    def run(command, label):
        record = command_record(command, label)
        with (out / (label + '.log')).open('w') as log:
            p = subprocess.run(record['command'], stdout=log, stderr=subprocess.STDOUT)
        record.update(exit_code=p.returncode, completed_at=stamp())
        save(out / 'commands.json', commands)
        if p.returncode:
            raise RuntimeError('stage failed: ' + label)

    try:
        save(out / 'selection.json', {'config': config, 'windows': windows, 'PAR_reference_copies': par_proof})
        bed = out / 'regions.bed'
        bed.write_text(''.join(f"{w['chrom']}\t{w['start']}\t{w['end']}\n" for w in windows))
        bam = out / 'selected.bam'
        run([samtools, 'view', '-@', '1', '--no-PG', '-T', reference, '-b', '-M', '-L', bed,
             '--save-counts', out / 'extraction-counts.json', '-o', bam, source], 'extract-preserving-records')
        run([samtools, 'quickcheck', '-v', bam], 'quickcheck')
        run([samtools, 'index', bam], 'index')
        counts, mq_counts = Counter(), Counter()
        with pysam.AlignmentFile(str(bam)) as reads:
            groups = reads.header.to_dict().get('RG', [])
            if len(groups) != 4 or {g.get('SM') for g in groups} != {config['sample']}:
                raise ValueError('unexpected sample/lane metadata')
            for read in reads:
                counts[read.get_tag('RG')] += 1
                if not read.flag & 0xF04:
                    mq_counts[read.mapping_quality] += 1
        if set(counts) != {g['ID'] for g in groups} or mq_counts[255]:
            raise ValueError('missing lane or unknown mapping qualities require separate handling')
        results = {}
        for mode in config['depth_modes']:
            label = mode['name']
            command = [samtools, 'depth', '-aa', '-b', bed, '-q', str(mode['minimum_BQ']),
                       '-Q', str(mode['minimum_MQ']), '-G', '0xF04', '-s', bam]
            record = command_record(command, label)
            with (out / (label + '.log')).open('w') as log, gzip.open(out / (label + '.depth.tsv.gz'), 'wt') as raw:
                process = subprocess.Popen(record['command'], stdout=subprocess.PIPE, stderr=log, text=True)
                try:
                    rows, n = collect_depth(process.stdout, windows, sequences, raw)
                except BaseException:
                    process.terminate()
                    process.wait()
                    raise
                finally:
                    process.stdout.close()
                code = process.wait()
            record.update(exit_code=code, completed_at=stamp())
            save(out / 'commands.json', commands)
            if code:
                raise RuntimeError('depth process failed')
            results[label] = {'mode': mode, 'position_count': n, 'groups': summarize(rows), 'windows': rows}
            save(out / (label + '.json'), results[label])
        summary = {'completed_at': stamp(), 'config': config, 'window_count': len(windows),
                   'selected_bases': sum(w['span'] for w in windows), 'read_records_by_lane': dict(counts),
                   'eligible_flag_read_records_by_MQ': dict(mq_counts), 'PAR_reference_copies': par_proof,
                   'depth_modes': {k: {'mode': v['mode'], 'position_count': v['position_count'], 'groups': v['groups']} for k, v in results.items()},
                   'limitations': ['Reference-only systematic sample, not whole-genome or gene-level callability.',
                       'Ratios use medians of equally sized reference-eligible windows, not estimated chromosome copies or a clinical karyotype.',
                       'Reference eligibility uses sequence composition, never observed depth; zero and low-mapping windows are preserved.',
                       'Native depth excludes deletions/reference skips, applies no BAQ, and suppresses overlapping mates using samtools rules.',
                       'No GC correction, external mappability mask, contamination estimate, cell-level mosaicism test or parental inference.',
                       'PAR copies remain separately aligned; mitochondrial coverage does not validate heteroplasmy, circular-boundary or NUMT handling.']}
        save(out / 'summary.json', summary)
        state.update(status='complete', completed_at=stamp(), elapsed_seconds=time.monotonic() - started,
                     output_sha256={p.name: digest(p) for p in out.iterdir() if p.is_file() and p.name != 'status.json'})
        save(out / 'status.json', state)
        print(json.dumps({'status': 'complete', 'windows': len(windows), 'bases': summary['selected_bases']}))
    except BaseException as error:
        state.update(status='failed', failed_at=stamp(), error_type=type(error).__name__)
        save(out / 'status.json', state)
        raise


if __name__ == '__main__':
    main()
