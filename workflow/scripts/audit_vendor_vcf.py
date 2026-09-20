#!/usr/bin/env python3
"""Compute aggregate VCF/reference checks on owned compute; emit no genome-wide rows."""
import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import random
import time

import pysam


def inspect(vcf_path, reference_path, sample_size=3000, seed=20260918):
    rng = random.Random(seed)
    counts = Counter()
    filters = Counter()
    reservoir = []
    blocks = {}
    candidate_phase = []
    with pysam.VariantFile(str(vcf_path)) as vcf:
        if len(vcf.header.samples) != 1:
            raise ValueError('expected a single-sample VCF')
        for record in vcf:
            counts['records'] += 1
            filter_names = tuple(record.filter.keys())
            filters[';'.join(filter_names) or '.'] += 1
            if filter_names != ('PASS',):
                continue
            counts['pass_records'] += 1
            sample = record.samples[0]
            counts['pass_phased_genotypes'] += bool(sample.phased)
            pgt, pid = sample.get('PGT'), sample.get('PID')
            if pgt not in (None, '.', '') and pid not in (None, '.', ''):
                counts['pass_records_with_pgt_and_pid'] += 1
                key = (record.contig, str(pid))
                if key not in blocks:
                    blocks[key] = [record.pos, record.pos]
                else:
                    blocks[key][0] = min(blocks[key][0], record.pos)
                    blocks[key][1] = max(blocks[key][1], record.pos)
            if record.contig.removeprefix('chr') == '15' and record.pos in (40209701, 40220612):
                candidate_phase.append({'chromosome': '15', 'position': record.pos,
                    'genotype_phased': bool(sample.phased), 'pgt_present': pgt not in (None, '.', ''),
                    'pid_present': pid not in (None, '.', ''), 'ps_present': sample.get('PS') not in (None, '.', '')})
            item = (record.contig, record.pos, record.ref)
            if len(reservoir) < sample_size:
                reservoir.append(item)
            else:
                index = rng.randrange(counts['pass_records'])
                if index < sample_size:
                    reservoir[index] = item
    comparison = Counter()
    with pysam.FastaFile(str(reference_path)) as reference:
        available = set(reference.references)
        for contig, position, expected in reservoir:
            if contig not in available:
                comparison['missing_contig'] += 1
                continue
            observed = reference.fetch(contig, position - 1, position - 1 + len(expected))
            comparison['matched' if observed.upper() == expected.upper() else 'mismatched'] += 1
    return {'record_counts': dict(counts), 'filter_counts': dict(filters),
            'reference_check': {'sample_size': len(reservoir), 'sampling': 'reservoir over all PASS records',
                                'seed': seed, **dict(comparison)},
            'vendor_phase': {'pass_phase_blocks': len(blocks),
                            'maximum_observed_variant_span_bp': max((end - start for start, end in blocks.values()), default=0),
                            'candidate_phase_fields': candidate_phase,
                            'interpretation': 'Vendor annotations only; absence of an observed phase block does not prove cis or trans.'}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--vcf', required=True, type=Path)
    parser.add_argument('--reference', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--sample-size', type=int, default=3000)
    parser.add_argument('--seed', type=int, default=20260918)
    args = parser.parse_args()
    if args.sample_size < 1:
        parser.error('sample-size must be positive')
    if args.output.exists():
        parser.error('existing audit preserved; choose a new output path')
    os.nice(10)
    started = time.monotonic()
    result = inspect(args.vcf, args.reference, args.sample_size, args.seed)
    result.update(completed_at=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
                  elapsed_seconds=round(time.monotonic() - started, 2),
                  vcf_file=args.vcf.name, reference_file=args.reference.name,
                  audit_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    check = result['reference_check']
    result['status'] = 'passed' if check.get('matched', 0) == check['sample_size'] and check['sample_size'] > 0 else 'failed'
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix('.json.tmp')
    temporary.write_text(json.dumps(result, indent=2) + '\n')
    os.replace(temporary, args.output)
    print(json.dumps(result, indent=2))
    if result['status'] != 'passed':
        raise SystemExit(1)


if __name__ == '__main__':
    main()
