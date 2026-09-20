#!/usr/bin/env python3
"""Explain selected exact-window losses without relaxing the evidence filters."""
import argparse
from collections import Counter
import json
from pathlib import Path

import pysam

from audit_final_alignment_indels import window_hypotheses
from audit_candidate_phase import spanning_sequence
from prepare_gatk_compatible_cram import digest, save, stamp


def window_state(read, start, end):
    if read.flag & 0xF04:
        return 'excluded_flags'
    if read.mapping_quality < 20 or read.mapping_quality == 255:
        return 'low_or_unknown_mapping_quality'
    if any(op == 3 for op, _ in read.cigartuples or []):
        return 'reference_skip'
    pairs = {r: q for q, r in read.get_aligned_pairs() if r in (start, end - 1)}
    left, right = pairs.get(start), pairs.get(end - 1)
    if left is None or right is None:
        return 'missing_left' if left is None and right is not None else 'missing_right' if right is None and left is not None else 'missing_both'
    if read.query_qualities is None:
        return 'missing_qualities'
    if right < left:
        return 'invalid_boundary_order'
    if min(read.query_qualities[left:right + 1], default=0) < 20:
        return 'low_window_quality'
    return 'eligible'


def target_cycle(read, position_0based):
    query = next((q for q, r in read.get_aligned_pairs() if r == position_0based), None)
    if query is None:
        return None
    return read.query_length - query if read.is_reverse else query + 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('run', 'reference', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError('existing diagnostic preserved')
    status = json.loads((args.run / 'status.json').read_text())
    if status['status'] != 'complete':
        raise ValueError('completed indel review required')
    for name in ('targets.bam', 'targets.bam.bai', 'evidence.json'):
        if digest(args.run / name) != status['output_sha256'][name]:
            raise ValueError('completed input differs: ' + name)
    if digest(args.reference) != status['reference_sha256']:
        raise ValueError('reference differs')
    evidence = json.loads((args.run / 'evidence.json').read_text())
    targets = []
    with pysam.FastaFile(str(args.reference)) as fasta, pysam.AlignmentFile(str(args.run / 'targets.bam')) as reader:
        for target in evidence['config']['targets']:
            row = {'key': target['key'], 'reference_101bp': fasta.fetch(target['chrom'], target['pos'] - 51, target['pos'] + 50), 'windows': []}
            for left, right in evidence['config']['window_flanks_bp']:
                start, end, hypotheses = window_hypotheses(fasta, target, left, right)
                lookup = {sequence: label for label, sequence in hypotheses.items()}
                counts, geometry = Counter(), Counter()
                for read in reader.fetch(target['chrom'], start, end):
                    state = window_state(read, start, end)
                    counts[state] += 1
                    if state in ('excluded_flags', 'low_or_unknown_mapping_quality'):
                        continue
                    if (left, right) == (1, 2):
                        sequence = spanning_sequence(read, start, end, 20)
                        label = lookup.get(sequence, 'other') if sequence is not None else 'not_eligible'
                        geometry[(label, read.is_reverse, read.reference_start, read.cigarstring,
                                  read.mapping_quality, read.get_tag('NM') if read.has_tag('NM') else None,
                                  target_cycle(read, target['pos'] - 1), read.query_length, read.get_tag('RG'))] += 1
                row['windows'].append({'flanks_bp': [left, right], 'per_read_reasons': dict(counts),
                    'geometries': [{'label': k[0], 'reverse': k[1], 'start_0based': k[2], 'cigar': k[3],
                        'MQ': k[4], 'NM': k[5], 'anchor_sequencing_cycle': k[6], 'read_length': k[7],
                        'RG': k[8], 'read_count': n} for k, n in geometry.items()]})
            targets.append(row)
    result = {'completed_at': stamp(), 'status_sha256': digest(args.run / 'status.json'),
              'reference_sha256': status['reference_sha256'], 'targets': targets,
              'source_sha256': {p.name: digest(p) for p in [Path(__file__), Path(__file__).with_name('audit_final_alignment_indels.py'),
                                  Path(__file__).with_name('audit_candidate_phase.py')]},
              'scope': 'Descriptive per-read window and alignment geometry, not a genotype, artifact diagnosis or independent fragment count.',
              'limitations': ['Different window intervals include different overlapping records; their totals are not a shared denominator.',
                  'Every state is exclusive in the documented filter order; a record can have more than one underlying problem.',
                  'Sequencing cycle is based on stored query orientation and length; hard-clipped/unrecorded cycles are not reconstructed.',
                  'No base-quality threshold relaxation, read realignment or clinical classification.']}
    save(args.output, result)
    print(json.dumps({'output': str(args.output), 'sha256': digest(args.output), 'targets': len(targets)}))


if __name__ == '__main__':
    main()
