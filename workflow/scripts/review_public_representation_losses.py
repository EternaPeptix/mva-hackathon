#!/usr/bin/env python3
"""Inspect every lost truth match when retaining multiallelic genotypes.

Bounded raw-call context is descriptive, not independent read/haplotype validation.
"""
import argparse
import json
from pathlib import Path

import pysam

from audit_public_small_variants import compact, fingerprint, key
from rtg_small_variant_helpers import save, stamp
from summarize_public_representation import checked, counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('reconciliation', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Existing loss review preserved')
    receipt = args.reconciliation / 'summary.json'
    reconciliation = json.loads(receipt.read_text())
    if not reconciliation['identical_truth_universe_all_four_comparisons']:
        parser.error('Identical truth universes required')
    summaries = {}
    verified = {}
    for arm, pin in reconciliation['source_comparisons'].items():
        verified[arm + '/summary'] = checked(pin)
        summaries[arm] = json.loads(Path(pin['path']).read_text())
    examples = []
    for method in ('gatk', 'bcftools'):
        found = {}
        for arm, summary in summaries.items():
            pin = summary['callers'][method]['files']['tp_baseline']
            verified[arm + '/' + method + '/truth_matches'] = checked(pin)
            found[arm] = counts(pin['path'])
            verified[arm + '/' + method + '/calls'] = checked(summary['inputs'][method])
        losses = found['frozen_split'] - found['retained_multiallelic']
        expected = reconciliation['truth_match_transitions'][method]['lost_truth_matches']['records']
        if sum(losses.values()) != expected or expected > 20:
            raise ValueError('Loss census differs or exceeds bounded complete review')
        for identity in sorted(losses):
            chrom, pos, ref, alts, gt = identity
            example = {'method': method, 'truth': {'chrom': chrom, 'pos': pos, 'ref': ref, 'alt': alts, 'GT': gt},
                       'copies': losses[identity], 'nearby_calls': {}}
            for arm, summary in summaries.items():
                with pysam.VariantFile(summary['inputs'][method]['path']) as calls:
                    rows = [compact(record) for record in calls.fetch(chrom, max(0, pos - 51), pos + max(50, len(ref)))]
                    if len(rows) > 30:
                        raise ValueError('Context exceeds bounded review; no silent truncation')
                    example['nearby_calls'][arm] = rows
            examples.append(example)
    args.output.mkdir(parents=True)
    summary = {'status': 'complete', 'completed_at': stamp(), 'code': fingerprint(Path(__file__)),
               'reconciliation': fingerprint(receipt), 'verified_inputs': verified,
               'all_lost_truth_matches_reviewed': True, 'examples': examples,
               'limits': ['All losses are included; nearby call context is descriptive and not an independent reassembly or read review.',
                          'Native matching can depend on neighboring records; a same-position genotype difference is not an exhaustive causal attribution.',
                          'These are public HG002 examples, not patient findings.']}
    save(args.output / 'summary.json', summary)
    print(json.dumps({'status': 'complete', 'lost_truth_records_reviewed': sum(e['copies'] for e in examples)}))


if __name__ == '__main__':
    main()
