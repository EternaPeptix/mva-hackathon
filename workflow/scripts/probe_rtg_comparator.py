#!/usr/bin/env python3
"""Native RTG controls with literal alleles and independently specified outcomes.

These fixtures validate comparison semantics, not biological caller accuracy.
Phased split records test equivalent sequence representation, not phase scoring.
"""
import argparse
import json
import os
from pathlib import Path
import random
import shutil
import time

import pysam

from rtg_small_variant_helpers import comparison_command, count_outputs, execute, fingerprint, runtime, save, stamp


def fixture_reference():
    rng = random.Random(19092026)
    sequence = list(''.join(rng.choices('ACGT', k=2048)))
    for pos, bases in [(101, 'A'), (201, 'GATT'), (301, 'C'), (401, 'GAAAAAAC'),
                       (501, 'AC'), (801, 'GAAAAAAC')]:
        sequence[pos-1:pos-1+len(bases)] = bases
    return ''.join(sequence)


def write_vcf(path, rows, sequence):
    text = ['##fileformat=VCFv4.2', '##contig=<ID=chrControl,length=2048>',
            '##FILTER=<ID=LowQual,Description="Fictional filtered control">',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002']
    for pos, ref, alt, genotype, quality, filt in rows:
        if sequence[pos-1:pos-1+len(ref)] != ref: raise ValueError('literal fixture REF is invalid')
        text.append('\t'.join(map(str, ['chrControl', pos, '.', ref, alt, quality, filt, '.', 'GT', genotype])))
    raw = path.with_suffix(''); raw.write_text('\n'.join(text) + '\n')
    pysam.tabix_compress(str(raw), str(path), force=False)
    pysam.tabix_index(str(path), preset='vcf', force=False)


def run(root, output, config_path):
    output.mkdir(parents=True, exist_ok=False); started = time.monotonic()
    config = json.loads(config_path.read_text()); rtg, rtg_files = runtime(root, config)
    save(output / 'status.json', {'status': 'running', 'pid': os.getpid(), 'started_at': stamp()})
    shutil.copyfile(config_path, output / 'frozen-comparison-config.json')
    shutil.copyfile(Path(__file__), output / 'executed-probe_rtg_comparator.py')
    sequence = fixture_reference(); reference = output / 'control.fasta'
    reference.write_text('>chrControl\n' + '\n'.join(sequence[i:i+60] for i in range(0,len(sequence),60)) + '\n')
    sdf = output / 'control.sdf'; format_command = execute([rtg, 'format', '-o', sdf, reference], output, 'format-control')
    snv = (101, 'A', 'G', '0/1', 30, 'PASS')
    deletion = (201, 'GATT', 'G', '0/1', 30, 'PASS')
    insertion = (301, 'C', 'CTGA', '0/1', 30, 'PASS')
    cases = [
        ('exact_snv_deletion_insertion', [snv,deletion,insertion], [snv,deletion,insertion], [3,3,0,0], [0,2048]),
        ('wrong_alternate_allele', [snv], [(101,'A','T','0/1',30,'PASS')], [0,0,1,1], [0,2048]),
        ('correct_allele_wrong_genotype', [snv], [(101,'A','G','1/1',30,'PASS')], [0,0,1,1], [0,2048]),
        ('correct_homozygous_genotype', [(101,'A','G','1/1',30,'PASS')], [(101,'A','G','1/1',30,'PASS')], [1,1,0,0], [0,2048]),
        ('shifted_repeat_deletion', [(401,'GA','G','0/1',30,'PASS')], [(404,'AA','A','0/1',30,'PASS')], [1,1,0,0], [0,2048]),
        ('mnp_as_two_snvs', [(501,'AC','GT','0/1',30,'PASS')], [(501,'A','G','0|1',30,'PASS'),(502,'C','T','0|1',30,'PASS')], [1,2,0,0], [0,2048]),
        ('boundary_shift_with_full_context', [(805,'AA','A','0/1',30,'PASS')], [(801,'GA','G','0/1',30,'PASS')], [1,1,0,0], [803,820]),
        ('filtered_call_is_excluded', [snv], [(101,'A','G','0/1',30,'LowQual')], [0,0,0,1], [0,2048]),
        ('outside_bed_is_not_scored', [snv], [(101,'A','T','0/1',30,'PASS')], [0,0,0,0], [200,350]),
        ('zero_qual_is_included', [snv], [(101,'A','G','0/1',0,'.')], [1,1,0,0], [0,2048]),
        ('missing_genotype_is_not_a_call', [snv], [(101,'A','G','./.',30,'PASS')], [0,0,0,1], [0,2048]),
        ('spanning_deletion_star_is_not_a_snv', [snv], [(101,'A','*','0/1',30,'PASS')], [0,0,0,1], [0,2048]),
    ]
    # Freeze literal expected outcomes before executing any comparisons.
    save(output / 'literal-cases.json', {'count_order': ['tp_baseline','tp_calls','fp','fn'], 'cases': cases})
    results = []; failures = []
    for name, baseline, calls, expected, interval in cases:
        folder = output / name; folder.mkdir(); b = folder / 'baseline.vcf.gz'; c = folder / 'calls.vcf.gz'
        write_vcf(b, baseline, sequence); write_vcf(c, calls, sequence)
        bed = folder / 'evaluation.bed'; bed.write_text('chrControl\t%d\t%d\n'%tuple(interval))
        command = execute(comparison_command(rtg, config, b, c, sdf, bed, folder/'native','chrControl'), folder, 'compare')
        result = count_outputs(folder / 'native'); observed = [result['counts'][k] for k in ('tp_baseline','tp_calls','fp','fn')]
        passed = observed == expected
        row = {'name': name, 'passed': passed, 'expected': expected, 'observed': observed,
               'command': command, 'comparison': result}
        results.append(row)
        if not passed: failures.append(name)
        save(output / 'progress.json', {'completed_cases': len(results), 'failed_cases': failures})
    summary = {'status': 'complete' if not failures else 'failed', 'completed_at': stamp(),
               'elapsed_seconds': time.monotonic()-started, 'config': fingerprint(config_path),
               'runtime': rtg_files, 'reference': fingerprint(reference), 'format_command': format_command,
               'literal_cases': fingerprint(output/'literal-cases.json'), 'cases': results,
               'cases_passed': len(results)-len(failures), 'cases_total': len(cases), 'failed_cases': failures,
               'code': fingerprint(Path(__file__)),
               'limits': ['Fictional comparator controls; not sequencing or clinical validation',
                          'Genotype-aware comparison does not establish parental or long-range phase',
                          'Equivalent representations can produce different baseline and call record counts']}
    save(output / 'summary.json', summary); save(output / 'status.json', {k: summary[k] for k in ('status','completed_at','cases_passed','cases_total','failed_cases')})
    if failures: raise ValueError('native comparator controls failed: ' + ','.join(failures))
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('root','output','config'): parser.add_argument('--'+key,required=True,type=Path)
    args = parser.parse_args(); result = run(args.root,args.output,args.config)
    print(json.dumps({k: result[k] for k in ('status','cases_passed','cases_total','elapsed_seconds')}))
