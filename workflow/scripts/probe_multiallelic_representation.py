#!/usr/bin/env python3
"""Test the interaction between multiallelic splitting and strict diploid comparison.

Fictional comparator inputs only. The frozen public comparison is not changed.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import time

import pysam

from probe_rtg_comparator import fixture_reference, write_vcf
from rtg_small_variant_helpers import comparison_command, count_outputs, execute, fingerprint, runtime, save, stamp


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('root', 'output', 'config'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Existing experiment preserved')
    root, out = args.root.resolve(), args.output.resolve()
    config = json.loads(args.config.read_text())
    rtg, runtime_files = runtime(root, config)
    out.mkdir(parents=True)
    start = time.monotonic()
    sequence = fixture_reference()
    reference = out / 'reference.fa'
    reference.write_text('>chrControl\n' + sequence + '\n')
    pysam.faidx(str(reference))
    sdf = out / 'reference.sdf'
    execute([rtg, 'format', '-o', sdf, reference], out, 'format')
    bed = out / 'evaluation.bed'
    bed.write_text('chrControl\t0\t2048\n')
    baseline = out / 'baseline.vcf.gz'
    # The two haplotypes delete two or four A bases from the same homopolymer.
    raw = [(401, 'GAAAA', 'GAA,G', '1/2', 30, 'PASS')]
    write_vcf(baseline, raw, sequence)
    split = out / 'split.vcf.gz'
    bcftools = root / 'tools/bio/bin/bcftools'
    normalize = execute([bcftools, 'norm', '--check-ref', 'e', '-f', reference,
                         '-m', '-any', '-Oz', '-o', split, baseline], out, 'split')
    execute([bcftools, 'index', '-t', split], out, 'split-index')
    phased = out / 'split-phased.vcf.gz'
    with pysam.VariantFile(str(split)) as source, pysam.VariantFile(str(phased), 'wz', header=source.header) as dest:
        split_rows = []
        for record in source:
            split_rows.append(str(record).strip())
            record.samples['HG002'].phased = True
            dest.write(record)
    pysam.tabix_index(str(phased), preset='vcf')
    wrong_gt = out / 'wrong-genotype.vcf.gz'
    write_vcf(wrong_gt, [(401, 'GAAAA', 'GAA,G', '1/1', 30, 'PASS')], sequence)
    # Keep the haplotype proof independent of native comparison labels. Apply only
    # each haplotype's nonreference edit; reference assertions are tested by RTG.
    expected_haplotypes = sorted(sequence[:400] + alt + sequence[405:] for alt in ('GAA', 'G'))
    with pysam.VariantFile(str(phased)) as vcf:
        records = list(vcf)
    replayed = []
    for hap in (0, 1):
        edits = [r for r in records if r.samples['HG002']['GT'][hap] != 0]
        if len(edits) != 1:
            raise ValueError('Unexpected split nonreference edit count')
        record = edits[0]
        alt = record.alleles[record.samples['HG002']['GT'][hap]]
        replayed.append(sequence[:record.start] + alt + sequence[record.start + len(record.ref):])
    if sorted(replayed) != expected_haplotypes:
        raise ValueError('Literal intended haplotypes differ after splitting')
    cases = [
        ('unsplit_strict', baseline, False, [1, 1, 0, 0]),
        ('split_unphased_strict', split, False, [0, 0, 2, 1]),
        ('split_phased_strict', phased, False, [0, 0, 2, 1]),
        ('split_unphased_ref_overlap', split, True, [1, 2, 0, 0]),
        ('split_phased_ref_overlap', phased, True, [1, 2, 0, 0]),
        ('wrong_GT_ref_overlap', wrong_gt, True, [0, 0, 1, 1]),
    ]
    save(out / 'prespecified-expectations.json', {
        'interpretation': 'Hypothesized strict/reference-overlap behavior from native documentation, recorded before comparisons.',
        'count_order': ['tp_baseline', 'tp_calls', 'fp', 'fn'],
        'cases': [{'name': n, 'calls': str(p), 'ref_overlap': overlap, 'expected': expected}
                  for n, p, overlap, expected in cases]})
    results = []
    for name, calls, overlap, expected in cases:
        target = out / name
        command = comparison_command(rtg, config, baseline, calls, sdf, bed, target, 'chrControl')
        if overlap:
            command.append('--ref-overlap')
        command_record = execute(command, out, name)
        compared = count_outputs(target)
        observed = [compared['counts'][k] for k in ('tp_baseline', 'tp_calls', 'fp', 'fn')]
        results.append({'name': name, 'ref_overlap': overlap, 'expected': expected,
                        'observed': observed, 'expectation_met': observed == expected,
                        'command': command_record, 'comparison': compared})
    result = {'status': 'complete', 'completed_at': stamp(), 'elapsed_seconds': time.monotonic() - start,
              'all_prespecified_expectations_met': all(r['expectation_met'] for r in results),
              'cases': results, 'config': fingerprint(args.config), 'runtime': runtime_files,
              'code': fingerprint(Path(__file__)), 'normalization': normalize,
              'literal_unsplit_row': raw, 'native_split_rows': split_rows,
              'intended_nonreference_haplotypes_identical': True,
              'haplotype_sha256': sorted(hashlib.sha256(s.encode()).hexdigest() for s in replayed),
              'limits': ['Fictional comparator experiment, not real-read variant calling or patient phase evidence.',
                         'The ref-overlap arms are explicit diagnostic interventions; the frozen public configuration remains unchanged.',
                         'A correct nonreference-edit interpretation does not make conflicting reference assertions in split records disappear.']}
    save(out / 'summary.json', result)
    print(json.dumps({'expectations_met': result['all_prespecified_expectations_met'],
                      'cases': [{'name': r['name'], 'observed': r['observed']} for r in results]}))


if __name__ == '__main__':
    main()
