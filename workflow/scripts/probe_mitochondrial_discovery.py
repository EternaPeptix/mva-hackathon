#!/usr/bin/env python3
"""End-to-end native control: synthetic circular paired reads on public rCRS sequence."""
import argparse
from collections import Counter
import json
from pathlib import Path
import random
import subprocess
import sys

import pysam

from prepare_gatk_compatible_cram import digest, save, stamp
from run_mitochondrial_discovery import sequence


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('root', 'resources', 'output'):
        p.add_argument('--' + name, type=Path, required=True)
    args = p.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    bases = sequence(args.resources / 'Homo_sapiens_assembly38.chrM.fasta')[1]
    samtools, bwa, bcf = [args.root / 'tools/bio/bin' / x for x in ('samtools', 'bwa', 'bcftools')]
    ref = out / 'synthetic-reference.fa'
    ref.write_text('>MT\n' + bases + '\n')
    variants = []
    for pos, kind, fraction in [(50, 'SNV', 0.10), (16540, 'SNV', 0.25), (6200, 'INS', 0.20), (9000, 'SNV', 1.0)]:
        reference = bases[pos - 1]
        alternate = reference + 'GAT' if kind == 'INS' else {'A':'C', 'C':'G', 'G':'T', 'T':'A'}[reference]
        variants.append({'pos': pos, 'ref': reference, 'alt': alternate, 'mixture_fraction': fraction})
    truth = out / 'truth.vcf'
    truth.write_text('##fileformat=VCFv4.2\n##contig=<ID=MT,length=16569>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' +
                     ''.join(f"MT\t{x['pos']}\t.\t{x['ref']}\t{x['alt']}\t.\tPASS\t.\n" for x in sorted(variants,key=lambda x:x['pos'])))
    haplotypes = []
    for optional in [None, 0, 1, 2]:
        selected = [variants[3]] + ([] if optional is None else [variants[optional]])
        hap = bases
        for variant in sorted(selected, key=lambda x:x['pos'], reverse=True):
            j = variant['pos'] - 1
            assert hap[j:j+len(variant['ref'])] == variant['ref']
            hap = hap[:j] + variant['alt'] + hap[j+len(variant['ref']):]
        haplotypes.append(hap)
    rng = random.Random(20260920)
    allocation = Counter()
    r1, r2 = out / 'reads.1.fastq', out / 'reads.2.fastq'
    complement = str.maketrans('ACGTN', 'TGCAN')
    with r1.open('w') as left, r2.open('w') as right:
        for i in range(30000):
            choice = rng.random()
            h = 0 if choice < .45 else (1 if choice < .55 else (2 if choice < .80 else 3))
            allocation[str(h)] += 1
            hap = haplotypes[h]
            start, length = rng.randrange(len(hap)), rng.randrange(300, 551)
            fragment = (hap + hap)[start:start+length]
            first, second = fragment[:150], fragment[-150:].translate(complement)[::-1]
            name = f'SIM:1:FLOW:1:1:{i//500+1}:{i%500+1}'
            left.write('@' + name + '\n' + first + '\n+\n' + 'I'*150 + '\n')
            right.write('@' + name + '\n' + second + '\n+\n' + 'I'*150 + '\n')
    commands = []
    with (out / 'preparation.log').open('w') as log:
        for command in [[samtools, 'faidx', ref], [samtools, 'dict', '-o', ref.with_suffix('.dict'), ref], [bwa, 'index', ref],
                        [bcf, 'norm', '-c', 'e', '-f', ref, '-Ov', '-o', out/'truth.normalized.vcf', truth]]:
            commands.append(list(map(str, command)))
            subprocess.run(commands[-1], stdout=log, stderr=subprocess.STDOUT, check=True)
        command = [bwa, 'mem', '-t', '2', '-Y', '-R', r'@RG\tID:synthetic\tSM:fictional\tPL:ILLUMINA\tLB:synthetic', ref, r1, r2]
        commands.append(list(map(str, command)))
        with (out/'original.sam').open('w') as dest:
            subprocess.run(commands[-1], stdout=dest, stderr=log, check=True)
        bam = out / 'original.bam'
        for command in [[samtools, 'sort', '-@', '2', '-o', bam, out/'original.sam'], [samtools, 'index', bam]]:
            commands.append(list(map(str, command)))
            subprocess.run(commands[-1], stdout=log, stderr=subprocess.STDOUT, check=True)
    config = {'sample': 'fictional', 'expected_read_groups': 1, 'input_contig': 'MT',
              'reference_sha256': digest(ref), 'alignment_sha256': digest(bam),
              'mitochondrial_sequence_sha256': json.loads((args.resources/'manifest.json').read_text())['regular_sequence_sha256'],
              'scope': 'Synthetic circular paired-read control on public reference; not biological replication or clinical sensitivity.'}
    save(out/'config.json', config)
    runner = Path(__file__).with_name('run_mitochondrial_discovery.py')
    command = [sys.executable, runner, '--root', args.root, '--reference', ref, '--alignment', bam,
               '--resources', args.resources, '--config', out/'config.json', '--output', out/'run']
    commands.append(list(map(str, command)))
    save(out/'commands.json', commands)
    with (out/'runner.log').open('w') as log:
        subprocess.run(commands[-1], stdout=log, stderr=subprocess.STDOUT, check=True)
    normalized = out/'observed.normalized.vcf.gz'
    with (out/'normalization.log').open('w') as log:
        subprocess.run(list(map(str, [bcf, 'norm', '-c', 'e', '-f', args.resources/'Homo_sapiens_assembly38.chrM.fasta',
                      '-m', '-any', '-Oz', '-o', normalized, out/'run/masked.assuming-no-contamination.vcf.gz'])), stdout=log, stderr=subprocess.STDOUT, check=True)
    with pysam.VariantFile(str(normalized)) as reader:
        observed = {(x.pos, x.ref, x.alts[0]): {'AF': x.samples['fictional']['AF'][0], 'AD': x.samples['fictional']['AD'],
                    'filters': list(x.filter)} for x in reader}
    with pysam.VariantFile(str(out/'truth.normalized.vcf')) as reader:
        identities = [(x.pos, x.ref, x.alts[0]) for x in reader]
    checks = []
    expected_by_ref_alt = {(v['ref'],v['alt']): v for v in variants}
    for key in identities:
        if key not in observed:
            raise ValueError('expected synthetic allele not recovered: ' + str(key))
        result = observed[key]
        # Match the known single insertion after normalization; SNV positions remain fixed.
        expected = next(v for v in variants if (len(v['alt'])>1 and len(key[2])>1) or
                        (v['pos']==key[0] and v['ref']==key[1] and v['alt']==key[2]))
        if abs(result['AF']-expected['mixture_fraction']) > .10:
            raise ValueError('synthetic allele fraction inconsistent with planted mixture')
        if result['filters'] != ['PASS']:
            raise ValueError('synthetic allele fails implemented filters: ' + str(key))
        checks.append({'pos':key[0], 'ref':key[1], 'alt':key[2], 'planted_mixture_fraction':expected['mixture_fraction'], **result})
    state = json.loads((out/'run/status.json').read_text())
    summary = json.loads((out/'run/summary.json').read_text())
    if state['status'] != 'complete' or summary['rejected_liftover_records'] or summary['coverage']['positions'] != 16569:
        raise ValueError('incomplete synthetic circular run')
    receipt = {'completed_at':stamp(), 'status':'complete', 'read_pairs':30000, 'haplotype_pair_allocation':dict(allocation),
        'controls':checks, 'all_reference_positions_accounted':16569, 'liftover_rejects':0,
        'runner_status_sha256':digest(out/'run/status.json'),
        'source_sha256':{str(x):digest(x) for x in [Path(__file__),runner]},
        'scope':'Four planted alleles recovered after real alignment, reversion, remapping, duplicate marking, calling, liftover and filtering. Error-free synthetic reads are not a clinical sensitivity or NUMT/contamination validation.'}
    save(out/'summary.json', receipt)
    print(json.dumps(receipt))


if __name__ == '__main__':
    main()
