#!/usr/bin/env python3
"""Export bounded SNP evidence plus vendor calls and raw read-quality context."""
import argparse
from collections import Counter
import datetime
import json
from pathlib import Path

import pysam

from audit_targeted_snvs import target_calls
from screen_secondary_findings import fingerprint


def raw_base_profile(bam, target):
    counts, exclusions = Counter(), Counter()
    with pysam.AlignmentFile(str(bam), 'rb') as reader:
        for read in reader.fetch(target['chrom'], target['pos'] - 1, target['pos']):
            if read.flag & 0xF04:
                exclusions['flags_0xF04'] += 1
                continue
            query = next((q for q, r in read.get_aligned_pairs() if r == target['pos'] - 1), None)
            if query is None:
                exclusions['no_aligned_base'] += 1
                continue
            quality = None if read.query_qualities is None else read.query_qualities[query]
            counts[(read.query_sequence[query].upper(),
                    '20+' if read.mapping_quality >= 20 else 'below20',
                    'unknown' if quality is None else '20+' if quality >= 20 else 'below20',
                    bool(read.is_paired and not read.is_proper_pair))] += 1
    return {'rows': [{'base': base, 'MQ_bin': mq, 'raw_BQ_bin': bq, 'improper_pair': improper, 'reads': n}
                     for (base, mq, bq, improper), n in sorted(counts.items())],
            'excluded_read_counts': dict(exclusions),
            'scope': 'Raw aligned bases with flags 0xF04 excluded, without BAQ or collapsing mates; not independent fragment counts or a genotype.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--vendor-vcf', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    state = json.loads((args.run / 'run-status.json').read_text())
    if state['status'] != 'complete':
        raise ValueError('audit not complete')
    for name in ('evidence.json', 'markdup.bam', 'markdup.bam.bai'):
        if fingerprint(args.run / name)['sha256'] != state['output_hashes'][name]:
            raise ValueError('audit output hash changed: ' + name)
    evidence = json.loads((args.run / 'evidence.json').read_text())
    vendor = fingerprint(args.vendor_vcf)
    if vendor['sha256'] != 'bec64c68d3ad10d05b8071d07cdf87d217013d66d3ba02e24b505035c9d42175':
        raise ValueError('unexpected source vendor VCF')
    targets = evidence['config']['targets']
    result = {'created_at': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'run': state, 'evidence': evidence,
              'vendor_VCF': vendor, 'vendor_calls': target_calls(args.vendor_vcf, targets, evidence['config']['sample']),
              'raw_alignment_profiles': {t['key']: raw_base_profile(args.run / 'markdup.bam', t) for t in targets},
              'exporter': fingerprint(Path(__file__)),
              'scope': 'Only four selected SNPs and aggregate quality counts; no read names, raw reads or genome-wide patient calls.'}
    after = args.vendor_vcf.stat()
    if (after.st_size, after.st_mtime_ns) != (vendor['bytes'], vendor['mtime_ns']):
        raise ValueError('vendor VCF changed during export')
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'output': args.output.name, 'bytes': args.output.stat().st_size,
                      'sha256': fingerprint(args.output)['sha256'], 'targets': len(targets)}))


if __name__ == '__main__':
    main()
