#!/usr/bin/env python3
"""Reassess retained synthetic SV output without changing its calls or original receipt."""
import argparse
import json
from pathlib import Path

from probe_structural_variants import compare_truth, read_calls, save, stamp
from simulate_sv_reads import sha


def reassess(experiment, output):
    if output.exists():
        raise ValueError('reassessment output already exists')
    original_path = experiment / 'summary.json'
    truth_path = experiment / 'fixture/truth.json'
    calls_path = experiment / 'calls.bcf'
    original = json.loads(original_path.read_text())
    truth = json.loads(truth_path.read_text())
    if original['truth_sha256'] != sha(truth_path):
        raise ValueError('retained truth differs from original experiment')
    if original['output_sha256']['calls.bcf'] != sha(calls_path):
        raise ValueError('retained native calls differ from original experiment')
    raw_calls = json.loads(json.dumps(read_calls(calls_path, 'SYNTH_SV')))
    if raw_calls != original['native_calls']:
        raise ValueError('parsed native calls differ from original receipt')
    comparison = compare_truth(truth, raw_calls)
    result = {
        'reassessed_at': stamp(),
        'original_summary_sha256': sha(original_path),
        'original_probe_source_sha256': original['probe_source_sha256'],
        'current_probe_source_sha256': sha(Path(__file__).with_name('probe_structural_variants.py')),
        'reassessment_source_sha256': sha(Path(__file__)),
        'truth_sha256': sha(truth_path),
        'native_calls_sha256': sha(calls_path),
        'retained_calls_reparsed_identically': True,
        'comparison': comparison,
        'fully_detected_event_ids': [row['event_id'] for row in comparison['events']
                                    if row['site_detected_within_tolerance']
                                    and row['required_connection_types_detected']],
        'scope': 'Corrected endpoint-order/orientation comparison of the original synthetic calls. '
                 'No realignment, recall, truth change, or general performance estimate. '
                 'The original receipt and original comparison remain unchanged.',
    }
    save(output, result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--experiment', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = reassess(args.experiment, args.output)
    print(json.dumps({'reassessed_at': result['reassessed_at'],
                      'fully_detected_event_ids': result['fully_detected_event_ids'],
                      'unmatched_calls': len(result['comparison']['unmatched_calls']),
                      'control_calls': len(result['comparison']['control_chromosome_calls'])}, indent=2))


if __name__ == '__main__':
    main()
