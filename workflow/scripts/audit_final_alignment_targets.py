#!/usr/bin/env python3
"""Bounded four-lane SNP review preserving whole-alignment duplicate marks.

Uses the existing tested SNP evidence/reference-confidence helpers. Extraction
does not realign, re-mark duplicates, fix mate fields, or filter source records.
Raw regional BAMs and VCFs remain in the owned analysis output directory.
"""
import argparse
from collections import Counter
import json
import os
from pathlib import Path
import subprocess
import time

import pysam

from audit_targeted_snvs import (fragment_evidence, target_calls, reference_support,
                                  reference_PL_gap, validate_targets)
from prepare_gatk_compatible_cram import digest, save, stamp
from summarize_targeted_snvs import raw_base_profile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('root', 'alignment-state', 'config', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    os.nice(10)
    root, output = args.root.resolve(), args.output.resolve()
    if output.exists():
        raise ValueError('existing target output preserved')
    provenance = json.loads(args.alignment_state.read_text())
    if provenance.get('status') != 'complete' or provenance.get('expected_lanes') != 4:
        raise ValueError('completed four-lane alignment required')
    source = Path(provenance['cram_path'])
    if source.stat().st_size != provenance['cram_bytes'] or digest(source) != provenance['cram_sha256']:
        raise ValueError('completed alignment differs from its receipt')
    config = json.loads(args.config.read_text())
    reference = root / 'ref/realignment/GRCh38_nochr.fa'
    if digest(reference) != config['reference_sha256']:
        raise ValueError('reference differs from target selection')
    validate_targets(config['targets'], reference)
    output.mkdir(parents=True, exist_ok=False)
    state = {'status': 'running', 'pid': os.getpid(), 'started_at': stamp(),
             'scope': 'Selected SNPs in completed four-lane alignment; same biological sample.',
             'source_cram_sha256': provenance['cram_sha256'],
             'source_state_sha256': digest(args.alignment_state),
             'config_sha256': digest(args.config), 'reference_sha256': config['reference_sha256'],
             'script_sha256': digest(Path(__file__))}
    commands = []
    runtime = root / 'tools/bio'
    samtools, bcftools = runtime / 'bin/samtools', runtime / 'bin/bcftools'
    java = [runtime / 'lib/jvm/bin/java', '-XX:ActiveProcessorCount=2', '-Xmx6g', '-jar',
            root / 'gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar']
    state['tool_and_helper_sha256'] = {str(p): digest(p) for p in [samtools, bcftools, java[0], java[4]] +
        [Path(__file__).with_name(name) for name in ('audit_targeted_snvs.py', 'audit_baseline_reads.py',
                                                   'summarize_targeted_snvs.py', 'prepare_gatk_compatible_cram.py')]}

    def run(command, stage):
        state.update(current_stage=stage)
        save(output / 'status.json', state)
        row = {'stage': stage, 'command': list(map(str, command)), 'started_at': stamp()}
        commands.append(row)
        save(output / 'commands.json', commands)
        with (output / (stage + '.log')).open('w') as log:
            result = subprocess.run(row['command'], stdout=log, stderr=subprocess.STDOUT)
        row.update(exit_code=result.returncode, completed_at=stamp())
        save(output / 'commands.json', commands)
        if result.returncode:
            raise RuntimeError('target stage failed: ' + stage)

    try:
        intervals = []
        for target in sorted(config['targets'], key=lambda r: (r['chrom'], r['pos'])):
            chrom, start, end = target['chrom'], max(0, target['pos'] - 1 - config['slice_padding_bp']), target['pos'] + config['slice_padding_bp']
            if intervals and intervals[-1][0] == chrom and intervals[-1][2] >= start:
                intervals[-1][2] = max(end, intervals[-1][2])
            else:
                intervals.append([chrom, start, end])
        bed = output / 'regions.bed'
        bed.write_text(''.join(f'{chrom}\t{start}\t{end}\n' for chrom, start, end in intervals))
        bam = output / 'targets.bam'
        run([samtools, 'view', '-@', '1', '--no-PG', '-T', reference, '-b', '-M', '-L', bed,
             '--save-counts', output / 'extraction-counts.json', '-o', bam, source], 'extract-preserving-records')
        run([samtools, 'quickcheck', '-v', bam], 'quickcheck')
        run([samtools, 'index', bam], 'index')
        with pysam.AlignmentFile(str(bam)) as reads:
            groups = reads.header.to_dict().get('RG', [])
            if len(groups) != 4 or {r.get('SM') for r in groups} != {config['sample']}:
                raise ValueError('four expected sample read groups absent')
            group_counts = Counter(read.get_tag('RG') for read in reads)
            if set(group_counts) != {r['ID'] for r in groups}:
                raise ValueError('all four lanes must contribute to selected regions')
        evidence = {'selection_config_sha256': digest(args.config), 'targets': config['targets'],
                    'scope': state['scope'], 'source_cram_sha256': provenance['cram_sha256'],
                    'regions_BED_0based_half_open': intervals, 'region_record_counts_by_RG': dict(group_counts),
                    'whole_alignment_duplicate_marks_preserved': True,
                    'reads': {t['key']: fragment_evidence(bam, reference, t, config['read_filters']) for t in config['targets']},
                    'raw_base_profiles': {t['key']: raw_base_profile(bam, t) for t in config['targets']}}
        regions = [f'{chrom}:{start + 1}-{end}' for chrom, start, end in intervals]
        limits = [part for region in regions for part in ('-L', region)]
        run(java + ['HaplotypeCaller', '-R', reference, '-I', bam, '-O', output / 'gatk.g.vcf.gz',
                    '--pair-hmm-implementation', 'LOGLESS_CACHING', '-ERC', 'BP_RESOLUTION'] + limits, 'gatk-reference-confidence')
        run(java + ['GenotypeGVCFs', '-R', reference, '-V', output / 'gatk.g.vcf.gz',
                    '-O', output / 'gatk.allsites.vcf.gz', '--include-non-variant-sites', 'true'] + limits, 'gatk-genotype-allsites')
        run([bcftools, 'mpileup', '-Ob', '-o', output / 'pileup.bcf', '-f', reference, '-r', ','.join(regions),
             '-a', 'AD,DP', '-q', '20', '-Q', '20', '-d', '100000', '--ff', '0xF04', bam], 'bcftools-pileup')
        run([bcftools, 'call', '-m', '-A', '-a', 'GQ', '-Oz', '-o', output / 'bcftools.allsites.vcf.gz', output / 'pileup.bcf'], 'bcftools-allsites')
        targets = output / 'alleles.tsv'
        targets.write_text(''.join(f"{t['chrom']}\t{t['pos']}\t{t['ref']},{t['alt']}\n" for t in config['targets']))
        pysam.tabix_compress(str(targets), str(targets) + '.gz', force=False)
        pysam.tabix_index(str(targets) + '.gz', seq_col=0, start_col=1, end_col=1, zerobased=False)
        run([bcftools, 'call', '-m', '-A', '-C', 'alleles', '-T', str(targets) + '.gz', '-a', 'GQ',
             '-Oz', '-o', output / 'bcftools.targeted.vcf.gz', output / 'pileup.bcf'], 'bcftools-targeted')
        evidence['calls'] = {}
        for name, filename in [('GATK_BP_gVCF', 'gatk.g.vcf.gz'), ('GATK_genotyped', 'gatk.allsites.vcf.gz'),
                               ('bcftools_allsites', 'bcftools.allsites.vcf.gz'), ('bcftools_targeted', 'bcftools.targeted.vcf.gz')]:
            path = output / filename
            run([bcftools, 'index', '-t', '-f', path], name + '-index')
            evidence['calls'][name] = target_calls(path, config['targets'], config['sample'])
        evidence['reference_support'] = {}
        for target in config['targets']:
            key = target['key']
            row = reference_support(evidence['reads'][key], {name: data[key] for name, data in evidence['calls'].items()},
                                    config['reference_support_policy'])
            if row['status'] == 'reference_supported_in_retained_lanes':
                row['status'] = 'reference_supported_in_completed_four_lane_alignment'
            row['meaning'] = 'Selected SNP reference evidence only; no clinical negative or mosaicism/structural exclusion.'
            evidence['reference_support'][key] = row
        evidence['target_reference_PL_gaps'] = {t['key']: [reference_PL_gap(c) for c in
            evidence['calls']['bcftools_targeted'][t['key']]['exact_site_records']] for t in config['targets']}
        evidence['limitations'] = ['Same child and sequencing data, not independent biological replication.',
            'Selected SNP windows; does not establish genome-wide callability, phase or clinical causality.',
            'Native diploid calls without BQSR. Original multiallelic alleles and genotype fields are retained.',
            'Caller depth and filtered fragment counts are different quantities. PL gap is not GQ or causal probability.',
            'The original selection config describes historical lanes; this run explicitly replaces that scope with all four lanes.']
        save(output / 'evidence.json', evidence)
        state.update(status='complete', completed_at=stamp(), commands_sha256=digest(output / 'commands.json'),
                     output_sha256={p.name: digest(p) for p in output.iterdir() if p.is_file() and p.name != 'status.json'})
        save(output / 'status.json', state)
        print(json.dumps({'status': 'complete', 'targets': len(config['targets']), 'output': str(output)}))
    except BaseException as error:
        state.update(status='failed', failed_at=stamp(), error_type=type(error).__name__)
        save(output / 'status.json', state)
        raise


if __name__ == '__main__':
    main()
