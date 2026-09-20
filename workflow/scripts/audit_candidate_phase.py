#!/usr/bin/env python3
"""Audit one nearby candidate pair from retained reads; emit aggregate counts only."""
import argparse
from collections import Counter, defaultdict
import hashlib
import itertools
import json
import os
from pathlib import Path
import subprocess
import time

import pysam


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def haplotypes(reference, start, sites):
    """Apply non-overlapping VCF alleles, retaining all matching hypotheses."""
    sequences = defaultdict(list)
    for first, second in itertools.product((0, 1), repeat=2):
        sequence = reference
        for (pos, ref, alt), use_alt in sorted(zip(sites, (first, second)), reverse=True):
            offset = pos - 1 - start
            if offset < 0 or reference[offset:offset + len(ref)] != ref:
                raise ValueError('candidate reference mismatch')
            if use_alt:
                sequence = sequence[:offset] + alt + sequence[offset + len(ref):]
        sequences[sequence].append(f'{first}|{second}')
    return sequences


def spanning_sequence(read, start, end, minimum_base_quality=20):
    """Query substring between aligned reference boundaries; reject skips and poor bases."""
    if any(operation == 3 for operation, _ in read.cigartuples or []):
        return None
    boundary = {r: q for q, r in read.get_aligned_pairs() if r in (start, end - 1)}
    left, right = boundary.get(start), boundary.get(end - 1)
    if left is None or right is None or right < left or read.query_qualities is None:
        return None
    if min(read.query_qualities[left:right + 1], default=0) < minimum_base_quality:
        return None
    return read.query_sequence[left:right + 1].upper()


def count_haplotypes(bam, reference, contig, sites):
    if len(sites) != 2 or sites[0][0] + len(sites[0][1]) > sites[1][0]:
        raise ValueError('exactly two ordered, non-overlapping sites required')
    # One extra flank before the first VCF anchor and two after the last REF allele.
    start = sites[0][0] - 2
    end = sites[-1][0] - 1 + len(sites[-1][1]) + 2
    with pysam.FastaFile(str(reference)) as fasta:
        expected = haplotypes(fasta.fetch(contig, start, end).upper(), start, sites)
    fragments = defaultdict(list)
    rejected = Counter()
    with pysam.AlignmentFile(str(bam), 'rb') as alignment:
        for read in alignment.fetch(contig, start, end):
            if read.flag & 0xF04 or read.mapping_quality < 20:
                rejected['flags_or_mapping_quality'] += 1
                continue
            sequence = spanning_sequence(read, start, end)
            if sequence is None:
                rejected['not_spanning_or_low_quality_or_refskip'] += 1
                continue
            labels = expected.get(sequence, [])
            label = labels[0] if len(labels) == 1 else ('ambiguous' if labels else 'other')
            rg = read.get_tag('RG') if read.has_tag('RG') else ''
            fragments[(rg, read.query_name)].append((label, read.is_reverse))
    counts, strand_counts = Counter(), defaultdict(Counter)
    for observations in fragments.values():
        labels = {label for label, _ in observations}
        label = next(iter(labels)) if len(labels) == 1 else 'mate_conflict'
        counts[label] += 1
        strands = {reverse for _, reverse in observations}
        strand = 'both' if len(strands) > 1 else ('reverse' if True in strands else 'forward')
        strand_counts[label][strand] += 1
    return {'window': {'contig': contig, 'start_1based': start + 1, 'end_1based': end},
            'counts': {key: counts[key] for key in ('0|0', '1|0', '0|1', '1|1', 'other', 'ambiguous', 'mate_conflict')},
            'label_definition': 'Each label describes alleles at site 1 and site 2 on one observed fragment; 0=REF, 1=ALT.',
            'strand_counts': dict(strand_counts), 'eligible_fragments': len(fragments),
            'excluded_read_counts': dict(rejected),
            'filters': {'minimum_mapping_quality': 20, 'minimum_base_quality': 20,
                        'exclude_flags': '0xF04', 'one_count_per_read_group_and_template': True,
                        'require_one_read_spanning_both_sites': True, 'baq_applied': False}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root, output = args.root.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    os.nice(10)
    started = time.monotonic()
    source = root / 'work/fancd2-read-validation-20260919'
    manifest = json.loads((source / 'slice-provenance.json').read_text())
    env = dict(os.environ, PATH=str(root/'tools/bio/bin') + os.pathsep + os.environ.get('PATH', ''))
    reference = root / 'ref/vendor/GRCh38_nochr.fa'
    commands = []
    def run(command):
        commands.append(command)
        with (output / 'processing.log').open('a') as log:
            subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    with pysam.FastaFile(str(reference)) as fasta:
        dictionary = dict(zip(fasta.references, fasta.lengths))
    inputs = []
    for item in manifest['files']:
        path = source / item['slice_name']
        if digest(path) != item['slice_sha256']:
            raise ValueError('slice checksum mismatch')
        with pysam.AlignmentFile(str(path), 'rb') as bam:
            if dict(zip(bam.references, bam.lengths)) != dictionary:
                raise ValueError('reference dictionary mismatch')
        inputs.append(str(path))
    run(['samtools', 'merge', '-o', str(output/'merged.bam')] + inputs)
    run(['samtools', 'collate', '-o', str(output/'collated.bam'), str(output/'merged.bam')])
    run(['samtools', 'fixmate', '-m', str(output/'collated.bam'), str(output/'fixmate.bam')])
    run(['samtools', 'sort', '-m', '128M', '-o', str(output/'sorted.bam'), str(output/'fixmate.bam')])
    bam = output / 'markdup.bam'
    run(['samtools', 'markdup', str(output/'sorted.bam'), str(bam)])
    run(['samtools', 'quickcheck', str(bam)])
    run(['samtools', 'index', str(bam)])
    sites = [(10046723, 'AG', 'A'), (10046725, 'TAAG', 'T')]
    annotations = []
    with pysam.VariantFile(str(root/'data/WGS_EX2312012_HGWCNDSX7.vcf.gz')) as vcf:
        for rec in vcf.fetch('3', sites[0][0]-1, sites[-1][0]):
            if (rec.pos, rec.ref, rec.alts[0]) in sites:
                sample = rec.samples['WGS_EX2312012']
                annotations.append({'contig': rec.chrom, 'position': rec.pos, 'ref': rec.ref,
                    'alt': rec.alts[0], 'GT': sample['GT'], 'GT_phased': sample.phased,
                    'PGT': sample.get('PGT'), 'PID': sample.get('PID')})
    result = {'scope': manifest['scope'], 'source_slices': manifest,
              'script_sha256': digest(Path(__file__)), 'vendor_annotations': annotations,
              'read_evidence': count_haplotypes(bam, reference, '3', sites), 'commands': commands,
              'completed_at': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
              'elapsed_seconds': round(time.monotonic()-started, 2),
              'limitations': ['Two historical lane alignments only; not the new four-lane run.',
                             'Exact local sequence matching; unmatched or ambiguous fragments remain separate.',
                             'Local phase does not establish clinical benignity or exclude other FANCD2 variants.']}
    (output/'phase-evidence.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'vendor_annotations': annotations, 'read_evidence': result['read_evidence']}))


if __name__ == '__main__':
    main()
