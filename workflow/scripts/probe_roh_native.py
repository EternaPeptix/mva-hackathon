#!/usr/bin/env python3
"""Run native ROH controls with known genotypes, AF and an unobserved gap."""
import argparse
from collections import Counter
import datetime
import hashlib
import json
import os
from pathlib import Path
import random
import subprocess
import time

import pysam
from segmented_roh import run_segmented


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fixture(path, gap=False):
    rng = random.Random(20260919)
    header = pysam.VariantHeader()
    for chrom in ('1', '2'):
        header.contigs.add(chrom, length=20_000_000)
    header.info.add('POP_AF', 1, 'Float', 'Known population AF in a synthetic experiment')
    header.info.add('AC', 'A', 'Integer', 'Alternate allele count in this single synthetic sample')
    header.info.add('AN', 1, 'Integer', 'Called allele count in this single synthetic sample')
    header.formats.add('GT', 1, 'String', 'Synthetic genotype')
    header.formats.add('PL', 'G', 'Integer', 'Synthetic genotype likelihoods')
    header.add_sample('SYNTH')
    counts = Counter()
    with pysam.VariantFile(str(path), 'w', header=header) as writer:
        for chrom in ('1', '2'):
            for pos in range(5001, 20_000_000, 10000):
                af = rng.choice([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
                if chrom == '1' and 5_000_000 <= pos < 15_000_000:
                    allele = int(rng.random() < af)
                    gt = (allele, allele)
                else:
                    gt = tuple(sorted((int(rng.random() < af), int(rng.random() < af))))
                if gap and chrom == '1' and 8_000_000 <= pos < 12_000_000:
                    counts['unobserved_markers_omitted'] += 1
                    continue
                record = writer.new_record(contig=chrom, start=pos-1, alleles=('A', 'G'))
                record.filter.add('PASS')
                record.info['POP_AF'] = af
                record.info['AC'] = (sum(gt),)
                record.info['AN'] = 2
                record.samples['SYNTH']['GT'] = gt
                record.samples['SYNTH']['PL'] = [(0,60,120),(60,0,60),(120,60,0)][sum(gt)]
                writer.write(record)
                counts[chrom + '_markers'] += 1
                counts['homref' if sum(gt) == 0 else 'heterozygous' if sum(gt) == 1 else 'homalt'] += 1
    return dict(counts)


def run(root, executable, label, source, frequency_args):
    output = root / (label + '.txt')
    command = [str(executable), 'roh', '--skip-indels', '--samples', 'SYNTH',
               '--hw-to-az', '6.7e-8', '--az-to-hw', '5e-9',
               '--output-type', 'sr', '--output', str(output)] + frequency_args + [str(source)]
    started = time.monotonic()
    result = subprocess.run(command, capture_output=True, text=True)
    elapsed = time.monotonic() - started
    (root / (label + '.stderr')).write_text(result.stderr)
    if result.returncode:
        raise RuntimeError('native ROH control failed: ' + label)
    regions = []
    site_counts = Counter()
    for line in output.read_text().splitlines():
        row = line.split('\t')
        if row[0] == 'RG':
            regions.append({'chrom': row[2], 'start_1based': int(row[3]),
                            'end_inclusive': int(row[4]), 'length_bp': int(row[5]),
                            'markers': int(row[6]), 'mean_model_phred': float(row[7])})
        elif row[0] == 'ST':
            site_counts[row[2] + '_' + ('AZ' if row[4] == '1' else 'HW')] += 1
    return {'command': command, 'elapsed_seconds': elapsed, 'exit_code': result.returncode,
            'output_sha256': digest(output), 'stderr_sha256': digest(root / (label + '.stderr')),
            'regions': regions, 'site_states': dict(site_counts)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bcftools', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    os.nice(10)
    started = time.monotonic()
    inputs = {}
    for name,gap in [('observed',False),('gap',True)]:
        p = args.output / (name + '.vcf')
        inputs[name] = {'counts': fixture(p,gap), 'sha256': digest(p)}
    runs = {}
    runs['external_AF'] = run(args.output,args.bcftools,'external_AF',args.output/'observed.vcf',['--AF-tag','POP_AF'])
    runs['self_AC_AN'] = run(args.output,args.bcftools,'self_AC_AN',args.output/'observed.vcf',[])
    runs['external_AF_gap'] = run(args.output,args.bcftools,'external_AF_gap',args.output/'gap.vcf',['--AF-tag','POP_AF'])
    true_regions = runs['external_AF']['regions']
    matching = [r for r in true_regions if r['chrom']=='1'
                and min(r['end_inclusive'],15_000_000)-max(r['start_1based']-1,5_000_000) >= 8_000_000]
    if not matching or any(r['chrom']=='2' and r['length_bp']>=1_000_000 for r in true_regions):
        raise AssertionError('external-frequency positive/control expectation not met')
    split = run_segmented(args.output/'gap.vcf','SYNTH',args.bcftools,args.output/'gap-segmented')
    split_regions = [r for segment in split['segments'] for r in segment['regions']]
    if len(split['gap_boundaries']) != 1 or len(split['segments']) != 3:
        raise AssertionError('expected independent models around the gap and on control chromosome')
    if any(r['chrom']=='1' and r['start_1based']<=8_000_000 and r['end_inclusive']>=12_000_000
           for r in split_regions):
        raise AssertionError('segmented result bridges the unobserved gap')
    if len([r for r in split_regions if r['chrom']=='1' and r['length_bp']>=2_500_000]) != 2:
        raise AssertionError('flanking homozygous regions were not retained')
    result = {'completed_at': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'elapsed_seconds': time.monotonic()-started, 'seed':20260919,
        'bcftools_version': subprocess.check_output([str(args.bcftools),'--version'],text=True),
        'bcftools_sha256': digest(args.bcftools), 'script_sha256': digest(Path(__file__)),
        'truth': {'chrom':'1','start_0based':5_000_000,'end_exclusive':15_000_000,
                  'gap_0based_half_open':[8_000_000,12_000_000]},
        'inputs': inputs, 'runs': runs, 'positive_and_control_expectations_passed': True,
        'segmented_gap_summary_sha256':digest(args.output/'gap-segmented/summary.json'),
        'segmented_gap_regions':split_regions,
        'independent_gap_models_and_flanking_recovery_passed':True,
        'gap_bridged_by_raw_HMM': any(r['chrom']=='1' and r['start_1based']<=8_000_000
            and r['end_inclusive']>=12_000_000 for r in runs['external_AF_gap']['regions']),
        'scope': 'Synthetic likelihood/AF model experiment. No read generation, mapping, clinical sensitivity, patient interval or parental-origin inference.'}
    (args.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
