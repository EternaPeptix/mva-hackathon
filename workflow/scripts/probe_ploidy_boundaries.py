#!/usr/bin/env python3
"""Native fictional boundary control for GATK regional and fixed ploidy."""
import argparse
import json
from pathlib import Path
import random
import subprocess
import sys

import pysam

from call_ploidy_regions import ploidy_at
from prepare_gatk_compatible_cram import digest, save, stamp


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root, out = args.root.resolve(), args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    rng = random.Random(202609201)
    sequence = ''.join(rng.choice('ACGT') for _ in range(12000))
    reference = out / 'fictional.fa'
    reference.write_text('>control\n' + sequence + '\n')
    pysam.faidx(str(reference))
    with (out / 'fictional.dict').open('w') as stream:
        subprocess.run([str(root / 'tools/bio/bin/samtools'), 'dict', str(reference)], stdout=stream, check=True)
    regions = [{'chrom': 'control', 'length': 12000,
        'intervals_0based_half_open': [[0, 3000, 2], [3000, 6000, 1], [6000, 9000, 2], [9000, 12000, 1]]}]
    targets = []
    for boundary in (3000, 6000, 9000):
        for pos in (boundary - 35, boundary, boundary + 1, boundary + 35):
            p, _ = ploidy_at(regions, 'control', pos - 1)
            ref = sequence[pos - 1]
            targets.append({'pos': pos, 'ref': ref, 'alt': 'ACGT'[('ACGT'.index(ref) + 1) % 4],
                'expected_ploidy': p, 'expected_GT': [1] if p == 1 else [0, 1], 'boundary': boundary})
    header = {'HD': {'VN': '1.6', 'SO': 'coordinate'}, 'SQ': [{'SN': 'control', 'LN': 12000}],
        'RG': [{'ID': 'synthetic', 'SM': 'fictional', 'PL': 'ILLUMINA'}]}
    bam = out / 'fictional.bam'
    with pysam.AlignmentFile(str(bam), 'wb', header=header) as stream:
        for boundary in (3000, 6000, 9000):
            for i in range(200):
                start = boundary - 110 + i // 2
                query = list(sequence[start:start + 149])
                edits = 0
                for target in targets:
                    if start <= target['pos'] - 1 < start + 149 and (target['expected_ploidy'] == 1 or i % 2 == 0):
                        query[target['pos'] - 1 - start] = target['alt']
                        edits += 1
                read = pysam.AlignedSegment(stream.header)
                read.query_name = f'fictional_{boundary}_{i}'
                read.query_sequence = ''.join(query)
                read.query_qualities = [40] * 149
                read.reference_id, read.reference_start, read.mapping_quality = 0, start, 60
                read.cigarstring, read.flag = '149M', 16 if i % 4 in (0, 1) else 0
                read.set_tag('RG', 'synthetic')
                read.set_tag('NM', edits)
                stream.write(read)
    pysam.index(str(bam))
    config = {'sample': 'fictional', 'expected_read_groups': 1,
        'reference_sha256': digest(reference), 'alignment_sha256': digest(bam), 'regions': regions,
        'scope': 'Error-free fictional reads at both sides, last base and first base of three ploidy transitions. Software contract only.'}
    save(out / 'config.json', config)
    save(out / 'targets.json', targets)
    shared = ['--root', str(root), '--reference', str(reference), '--alignment', str(bam), '--config', str(out / 'config.json')]
    regional = Path(__file__).with_name('call_ploidy_regions.py')
    fixed = Path(__file__).with_name('call_fixed_ploidy_regions.py')
    command = [sys.executable, str(regional), *shared, '--output', str(out / 'regional')]
    with (out / 'regional-launch.log').open('w') as stream:
        original = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT)
    original_status = json.loads((out / 'regional/status.json').read_text())
    stages = json.loads((out / 'regional/commands.json').read_text())
    if any(not any(s['stage'] == label and s.get('exit_code') == 0 for s in stages) for label in ('gatk', 'pileup', 'bcftools')):
        raise ValueError('regional native callers did not finish')
    if original.returncode == 0 or original_status['status'] != 'failed' or original_status['error_type'] != 'ValueError':
        raise ValueError('expected strict regional postprocessing failure was not reproduced')
    if 'genotype/model ploidy mismatch' not in (out / 'regional-launch.log').read_text():
        raise ValueError('regional failure was not the expected boundary contract')
    reused = out / 'regional/bcftools.original.vcf.gz'
    fixed_command = [sys.executable, str(fixed), *shared, '--output', str(out / 'fixed'),
        '--bcftools-vcf', str(reused), '--bcftools-sha256', digest(reused)]
    with (out / 'fixed-launch.log').open('w') as stream:
        subprocess.run(fixed_command, stdout=stream, stderr=subprocess.STDOUT, check=True)
    fixed_status = json.loads((out / 'fixed/status.json').read_text())
    if fixed_status['status'] != 'complete':
        raise ValueError('fixed recall incomplete')
    for name, expected in fixed_status['output_sha256'].items():
        if digest(out / 'fixed' / name) != expected:
            raise ValueError('fixed output changed')
    observations, mismatches = [], []
    for mode, caller in [('regional', 'gatk'), ('fixed', 'gatk'), ('fixed', 'bcftools')]:
        with pysam.VariantFile(str(out / mode / (caller + '.original.vcf.gz'))) as reader:
            records = list(reader)
        if len(records) != len(targets):
            raise ValueError('unexpected number of synthetic variant records')
        for target in targets:
            matches = [r for r in records if r.pos == target['pos'] and r.ref == target['ref'] and r.alts == (target['alt'],)]
            if len(matches) != 1:
                raise ValueError('synthetic allele missing or ambiguous')
            values = matches[0].samples['fictional']
            gt = list(values['GT'])
            row = {'mode': mode, 'caller': caller, **target, 'observed_GT': gt,
                'DP': values.get('DP'), 'GQ': values.get('GQ'), 'AD': list(values.get('AD') or []),
                'PL': list(values.get('PL') or [])}
            observations.append(row)
            if mode == 'fixed' and sorted(gt) != target['expected_GT']:
                raise ValueError('fixed caller did not recover expected original genotype')
            if len(gt) != target['expected_ploidy']:
                mismatches.append(row)
    if not mismatches or any(row['mode'] != 'regional' for row in mismatches):
        raise ValueError('boundary failure/recovery contract not demonstrated')
    result = {'completed_at': stamp(), 'status': 'complete', 'regional_expected_failure_exit_code': original.returncode,
        'commands': [command, fixed_command], 'control_genotypes': observations,
        'regional_ploidy_exceptions': mismatches, 'targets': len(targets), 'verified_fixed_caller_genotypes': len(targets) * 2,
        'fixed_status_sha256': digest(out / 'fixed/status.json'),
        'files': {str(p.relative_to(root / 'repo')): {'sha256': digest(p), 'bytes': p.stat().st_size}
            for p in (Path(__file__), regional, fixed, Path(__file__).with_name('prepare_gatk_compatible_cram.py'))},
        'native_tool_sha256': fixed_status['source_and_tool_sha256'],
        'scope': 'Native GATK boundary behavior reproduced and fixed-mode genotype recovery checked against fictional ground truth, including exact interval endpoints. Not clinical validation or a sensitivity benchmark.'}
    save(out / 'summary.json', result)
    print(json.dumps({'status': 'complete', 'targets': len(targets), 'regional_ploidy_exceptions': len(mismatches),
        'verified_fixed_caller_genotypes': len(targets) * 2, 'summary_sha256': digest(out / 'summary.json')}))


if __name__ == '__main__':
    main()
