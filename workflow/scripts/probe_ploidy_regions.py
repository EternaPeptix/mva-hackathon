#!/usr/bin/env python3
"""Exercise the actual two-caller regional-ploidy runner on a fictional genome."""
import argparse
import json
from pathlib import Path
import random
import subprocess
import sys

import pysam

from prepare_gatk_compatible_cram import digest, save, stamp


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root, out = args.root.resolve(), args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    rng = random.Random(20260920)
    sequence = ''.join(rng.choice('ACGT') for _ in range(8000))
    reference = out / 'fictional.fa'
    reference.write_text('>control\n' + sequence + '\n')
    pysam.faidx(str(reference))
    with (out / 'fictional.dict').open('w') as stream:
        subprocess.run([str(root / 'tools/bio/bin/samtools'), 'dict', str(reference)], stdout=stream, check=True)
    targets = []
    for position, ploidy, genotype, insertion in [(1001, 2, [0, 1], False), (4001, 1, [1], False), (5001, 1, [1], True)]:
        ref = sequence[position - 1]
        alt = ref + 'AG' if insertion else 'ACGT'[('ACGT'.index(ref) + 1) % 4]
        targets.append({'pos': position, 'ref': ref, 'alt': alt, 'expected_GT': genotype, 'expected_ploidy': ploidy})
    bam = out / 'fictional.bam'
    header = {'HD': {'VN': '1.6', 'SO': 'coordinate'}, 'SQ': [{'SN': 'control', 'LN': 8000}],
              'RG': [{'ID': 'synthetic', 'SM': 'fictional', 'PL': 'ILLUMINA'}]}
    with pysam.AlignmentFile(str(bam), 'wb', header=header) as stream:
        for target in targets:
            for i in range(80):
                start = target['pos'] - 1 - 80 + i // 2
                offset = target['pos'] - 1 - start
                is_alt = target['expected_ploidy'] == 1 or i % 2 == 0
                insertion = len(target['alt']) > 1
                if insertion:
                    query = sequence[start:start + 147]
                    query = query[:offset + 1] + 'AG' + query[offset + 1:]
                    cigar = f'{offset + 1}M2I{146 - offset}M'
                else:
                    query = sequence[start:start + 149]
                    if is_alt:
                        query = query[:offset] + target['alt'] + query[offset + 1:]
                    cigar = '149M'
                read = pysam.AlignedSegment(stream.header)
                read.query_name, read.query_sequence = f'fictional_{target["pos"]}_{i}', query
                read.query_qualities = [40] * len(query)
                read.reference_id, read.reference_start, read.mapping_quality = 0, start, 60
                read.cigarstring, read.flag = cigar, 16 if i % 4 in (0, 1) else 0
                read.set_tag('RG', 'synthetic')
                read.set_tag('NM', 2 if insertion else int(is_alt))
                stream.write(read)
    pysam.index(str(bam))
    config = {'sample': 'fictional', 'expected_read_groups': 1, 'reference_sha256': digest(reference),
              'alignment_sha256': digest(bam), 'regions': [{'chrom': 'control', 'length': 8000,
                'intervals_0based_half_open': [[0, 3000, 2], [3000, 8000, 1]]}],
              'scope': 'Synthetic diploid-to-haploid model transition; fictional reads, no patient evidence.'}
    save(out / 'config.json', config)
    runner = Path(__file__).with_name('call_ploidy_regions.py')
    command = [sys.executable, str(runner), '--root', str(root), '--reference', str(reference),
               '--alignment', str(bam), '--config', str(out / 'config.json'), '--output', str(out / 'run')]
    with (out / 'runner.log').open('w') as stream:
        subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, check=True)
    status = json.loads((out / 'run/status.json').read_text())
    if status['status'] != 'complete':
        raise ValueError('native runner incomplete')
    for name, expected in status['output_sha256'].items():
        if digest(out / 'run' / name) != expected:
            raise ValueError('native output differs')
    observed = []
    for caller in ('gatk', 'bcftools'):
        with pysam.VariantFile(str(out / 'run' / (caller + '.original.vcf.gz'))) as stream:
            calls = list(stream)
        for target in targets:
            matches = [r for r in calls if r.pos == target['pos'] and r.ref == target['ref'] and r.alts == (target['alt'],)]
            if len(matches) != 1:
                raise ValueError('missing or ambiguous native control allele')
            record = matches[0]
            gt = list(record.samples['fictional']['GT'])
            if sorted(gt) != target['expected_GT']:
                raise ValueError('wrong native allele or genotype ploidy')
            observed.append({'caller': caller, **target, 'observed_GT': gt,
                             'DP': record.samples['fictional'].get('DP'), 'GQ': record.samples['fictional'].get('GQ')})
    result = {'completed_at': stamp(), 'status': 'complete', 'command': command, 'control_genotypes': observed,
              'runner_status_sha256': digest(out / 'run/status.json'),
              'files': {str(p.relative_to(root / 'repo')): {'sha256': digest(p), 'bytes': p.stat().st_size}
                        for p in (Path(__file__), runner)},
              'native_tool_sha256': status['source_and_tool_sha256'],
              'scope': 'Both native callers recover the three fictional alleles with the specified diploid/haploid GT shapes. This is a software-contract control, not sensitivity or clinical validation.'}
    save(out / 'summary.json', result)
    print(json.dumps({'status': 'complete', 'matched_allele_genotypes': len(observed), 'summary_sha256': digest(out / 'summary.json')}))


if __name__ == '__main__':
    main()
