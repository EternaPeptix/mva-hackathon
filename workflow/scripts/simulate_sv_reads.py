#!/usr/bin/env python3
"""Construct diploid synthetic SV haplotypes and unaligned paired FASTQs.

Uniform random sequence and simple substitutions are deliberately limited controls,
not a model of human repeats, GC bias, contamination or mosaicism.
"""
import argparse
from collections import Counter
import gzip
import hashlib
import io
import json
from pathlib import Path
import time

import numpy as np


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''): h.update(block)
    return h.hexdigest()


def write_fasta(path, sequences):
    with path.open('w') as stream:
        for chrom, seq in sequences.items():
            stream.write('>' + chrom + '\n')
            for start in range(0, len(seq), 80): stream.write(seq[start:start+80] + '\n')


def simulate(output, coverage=30, seed=20260919, error_rate=0.001):
    if not 5 <= coverage <= 100 or not 0 <= error_rate <= 0.02:
        raise ValueError('simulation coverage/error rate outside declared bounds')
    output.mkdir(parents=True, exist_ok=False); started = time.monotonic()
    rng = np.random.default_rng(seed)
    def random_dna(length):
        return np.array(list('ACGT'), dtype='S1')[rng.integers(0, 4, length)].tobytes().decode('ascii')
    reference = {c: random_dna(1000000) for c in ('1', '2', '3')}
    r1, r2 = reference['1'], reference['2']; inserted = random_dna(100)
    alternate = {
        '1': r1[:150000] + r1[155000:356000] + r1[350000:356000]
             + r1[356000:550000] + reverse_complement(r1[550000:558000])
             + r1[558000:750000] + inserted + r1[750000:900000] + r2[500000:],
        '2': r2[:500000] + r1[900000:], '3': reference['3']}
    if sum(map(len, alternate.values())) != 3000000 - 5000 + 6000 + 100:
        raise AssertionError('alternate haplotype length accounting failed')
    write_fasta(output / 'reference.fa', reference); write_fasta(output / 'alternate-haplotype.fa', alternate)
    events = [
        {'id': 'deletion', 'type': 'DEL', 'chrom': '1', 'left_cut_0based': 150000, 'right_cut_0based': 155000,
         'length_bp': 5000, 'expected_genotype': [0, 1]},
        {'id': 'tandem_duplication', 'type': 'DUP', 'chrom': '1', 'left_cut_0based': 350000, 'right_cut_0based': 356000,
         'length_bp': 6000, 'expected_genotype': [0, 1]},
        {'id': 'inversion', 'type': 'INV', 'chrom': '1', 'left_cut_0based': 550000, 'right_cut_0based': 558000,
         'length_bp': 8000, 'expected_genotype': [0, 1], 'expected_connection_types': ['3to3', '5to5']},
        {'id': 'insertion', 'type': 'INS', 'chrom': '1', 'left_cut_0based': 750000, 'right_cut_0based': 750000,
         'length_bp': 100, 'inserted_sequence': inserted, 'expected_genotype': [0, 1]},
        {'id': 'reciprocal_translocation', 'type': 'BND', 'chrom': '1', 'left_cut_0based': 900000,
         'partner_chrom': '2', 'partner_cut_0based': 500000, 'expected_genotype': [0, 1],
         'expected_connection_types': ['3to5', '5to3']}]
    read_length = 150; mean_fragment = 350; fragment_sd = 35
    quality = '?' * read_length  # Fixed Q30, corresponding to the default substitution probability.
    counts = Counter(); fragment_counts = Counter(); pair_id = 0
    def add_errors(sequence):
        n = int(rng.binomial(read_length, error_rate))
        if not n: return sequence
        bases = list(sequence)
        for pos in rng.choice(read_length, n, replace=False):
            options = 'ACGT'.replace(bases[pos], '')
            bases[pos] = options[int(rng.integers(0, 3))]
        counts['introduced_substitution_errors'] += n
        return ''.join(bases)
    p1 = output / 'SYNTH_SV_R1.fastq.gz'; p2 = output / 'SYNTH_SV_R2.fastq.gz'
    with p1.open('wb') as f1, p2.open('wb') as f2, \
            gzip.GzipFile(fileobj=f1, mode='wb', filename='', mtime=0, compresslevel=3) as g1, \
            gzip.GzipFile(fileobj=f2, mode='wb', filename='', mtime=0, compresslevel=3) as g2, \
            io.TextIOWrapper(g1, encoding='ascii') as w1, io.TextIOWrapper(g2, encoding='ascii') as w2:
        for haplotype, sequences in [('reference', reference), ('alternate', alternate)]:
            for chrom, sequence in sequences.items():
                pairs = round((coverage / 2) * len(sequence) / (2 * read_length))
                for _ in range(pairs):
                    fragment = int(np.clip(round(rng.normal(mean_fragment, fragment_sd)), read_length+1, 700))
                    start = int(rng.integers(0, len(sequence)-fragment+1))
                    left = sequence[start:start+read_length]
                    right = reverse_complement(sequence[start+fragment-read_length:start+fragment])
                    if rng.integers(0, 2): left, right = right, left
                    left, right = add_errors(left), add_errors(right)
                    name = f'SYNTH_SV.{pair_id}'
                    w1.write(f'@{name}/1\n{left}\n+\n{quality}\n')
                    w2.write(f'@{name}/2\n{right}\n+\n{quality}\n')
                    counts[f'{haplotype}_chr{chrom}_pairs'] += 1; fragment_counts[str(fragment)] += 1; pair_id += 1
    truth = {'scope': 'Synthetic diploid simple-SV reference; no patient sequence or observed biological effect.',
        'seed': seed, 'random_generator': 'numpy.default_rng/PCG64', 'numpy_version': np.__version__,
        'nominal_diploid_depth': coverage, 'read_length': read_length, 'fragment_mean_before_truncation': mean_fragment,
        'fragment_sd_before_truncation': fragment_sd, 'fragment_limits': [read_length+1, 700],
        'independent_substitution_probability': error_rate, 'base_quality': 30,
        'reference_lengths': {c: len(s) for c, s in reference.items()},
        'alternate_lengths': {c: len(s) for c, s in alternate.items()}, 'reference_haplotype_unchanged': True,
        'unaltered_control_chromosome': '3', 'events_on_one_haplotype': events,
        'read_pairs': pair_id, 'read_records': pair_id*2, 'read_bases': pair_id*2*read_length,
        'counts': dict(counts), 'fragment_length_histogram': dict(sorted(fragment_counts.items(), key=lambda r: int(r[0]))),
        'elapsed_seconds_including_output_hashing': None,
        'files': {p.name: {'bytes': p.stat().st_size, 'sha256': sha(p)} for p in (output/'reference.fa', output/'alternate-haplotype.fa', p1, p2)},
        'script_sha256': sha(Path(__file__)),
        'limitations': ['Uniform random reference, simple events, one seed and no human repeat/GC/mappability model.',
            'Only independent substitutions are simulated, with fixed Q30; alternative error rates do not alter quality calibration.',
            'Fragments are sampled uniformly from each derivative chromosome; no reads are injected specifically at junctions.',
            'No simulated PCR duplicates; coordinate collisions can still be marked by the native duplicate algorithm.',
            'No clinical sensitivity, diagnostic probability or mosaic detection limit follows.']}
    truth['elapsed_seconds_including_output_hashing'] = time.monotonic() - started
    (output/'truth.json').write_text(json.dumps(truth, indent=2) + '\n')
    return truth


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True)
    p.add_argument('--coverage', type=float, default=30); p.add_argument('--seed', type=int, default=20260919)
    p.add_argument('--error-rate', type=float, default=0.001); args = p.parse_args()
    result = simulate(args.output, args.coverage, args.seed, args.error_rate)
    print(json.dumps({k: result[k] for k in ('read_pairs', 'counts', 'elapsed_seconds_including_output_hashing')}, indent=2))


if __name__ == '__main__': main()
