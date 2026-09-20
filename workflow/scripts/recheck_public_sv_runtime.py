#!/usr/bin/env python3
"""Repeat the frozen public SV comparisons using a validated rebuilt runtime.

Reuses original caller outputs and comparison inputs. It does not rerun alignment,
DELLY or SvABA, or establish patient/whole-genome accuracy.
"""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import sys
import time

from build_structural_exclusions import fingerprint, save, stamp
from compare_svaba_public import allele_records, key
from probe_sv_comparator import compare


def require(condition, message):
    if not condition:
        raise ValueError(message)


def check(path, expected):
    got = fingerprint(path)
    require(got['sha256'] == expected['sha256'], 'retained input differs: ' + str(path))
    return got


def signatures(path):
    rows = allele_records(path)
    return Counter(key(row) + (row['decision'], tuple(row['GT'] or ()), row['refined']) for row in rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('original-root', 'runtime-root', 'validation', 'config', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--expected-validation-sha256', required=True)
    args = parser.parse_args(); os.nice(10)
    old, runtime, out = args.original_root.resolve(), args.runtime_root.resolve(), args.output.resolve()
    require(Path(sys.prefix).resolve() == (runtime / 'tools/truvari-5.5.0').resolve(), 'run through the rebuilt comparator environment')
    require(importlib.metadata.version('truvari') == '5.5.0', 'comparator version differs')
    require(fingerprint(args.validation)['sha256'] == args.expected_validation_sha256, 'runtime validation identity differs')
    validation = json.loads(args.validation.read_text())
    require(validation['status'] == 'complete' and validation['wheel_provenance_mode'] == 'published_sources', 'source-built runtime not validated')
    require(sum(r['passed'] for r in validation['tests']) == 20, 'selected runtime controls did not all pass')
    config_identity = fingerprint(args.config)
    require(config_identity['sha256'] == '86b2e380e5dd390b05ee7deb3393171976ebe84f828bb556e3e4f5ef8f667e80', 'frozen comparison policy differs')
    config = json.loads(args.config.read_text())
    fresh_controls = args.validation.parent / 'structural-comparator/summary.json'
    controls = json.loads(fresh_controls.read_text())
    require(controls['all_passed'] and len(controls['cases']) == 6
            and controls['configuration']['sha256'] == config_identity['sha256'], 'current comparator controls differ')
    prior_dirs = {'DELLY': old / 'results/structural-public-comparison-20260919-v2',
                  'SvABA': old / 'results/svaba-public-comparison-20260919-v4'}
    prior = {name: json.loads((path / 'summary.json').read_text()) for name, path in prior_dirs.items()}
    first = prior['DELLY']
    baseline = Path(first['inputs']['baseline']['path'])
    reference = Path(first['inputs']['reference']['path'])
    bed = Path(first['inputs']['benchmark_bed']['path'])
    inputs = {name: check(path, first['inputs'][name]) for name, path in (
        ('baseline', baseline), ('reference', reference), ('benchmark_bed', bed))}
    comp_delly = Path(first['inputs']['comp']['path'])
    inputs['DELLY'] = check(comp_delly, first['inputs']['comp'])
    conversion_path = Path(prior['SvABA']['conversion']['path'])
    check(conversion_path, prior['SvABA']['conversion'])
    conversion = json.loads(conversion_path.read_text())
    comp_svaba = Path(conversion['output_vcf']['path'])
    inputs['SvABA'] = check(comp_svaba, conversion['output_vcf'])
    out.mkdir(parents=True, exist_ok=False); started = time.monotonic()
    state = {'status': 'running', 'pid': os.getpid(), 'started_at': stamp()}
    save(out / 'status.json', state)
    try:
        def run(name, comp):
            folder = out / name; folder.mkdir()
            result = compare(runtime, baseline, comp, reference, bed, folder / 'benchmark', config)
            save(folder / 'comparison.json', result)
            metric_checks = {}
            for stage in ('bench', 'refined'):
                names = ('TP-base', 'TP-comp', 'FP', 'FN', 'base cnt', 'comp cnt')
                before = {key: prior[name]['comparison'][stage][key] for key in names}
                after = {key: result[stage][key] for key in names}
                require(before == after, name + ' public counts changed at ' + stage)
                metric_checks[stage] = after
            records_checked = {}
            for side in ('base', 'comp'):
                filename = 'refine.' + side + '.vcf.gz'
                original = prior_dirs[name] / 'benchmark' / filename
                repeated = folder / 'benchmark' / filename
                before, after = signatures(original), signatures(repeated)
                require(before == after, name + ' allele, genotype, decision or refinement state changed: ' + side)
                records_checked[side] = sum(after.values())
            return {'caller': name, 'count_checks': metric_checks, 'allele_GT_decision_refinement_records_equal': records_checked,
                    'comparison': result, 'prior_summary': fingerprint(prior_dirs[name] / 'summary.json')}
        with ThreadPoolExecutor(max_workers=2) as pool:
            futures = [pool.submit(run, 'DELLY', comp_delly), pool.submit(run, 'SvABA', comp_svaba)]
            results = [future.result() for future in futures]
        result = {'completed_at': stamp(), 'elapsed_seconds': time.monotonic() - started,
                  'script': fingerprint(Path(__file__)), 'native_comparator_code': fingerprint(Path(__file__).with_name('probe_sv_comparator.py')),
                  'runtime_validation': fingerprint(args.validation), 'fresh_controls': fingerprint(fresh_controls),
                  'configuration': config_identity, 'inputs': inputs, 'comparisons': results,
                  'scope': 'Public HG002 chr20 comparator reproduction on retained caller outputs; 802 truth records per comparison.',
                  'limits': ['Original caller outputs reused; no new read alignment, variant discovery, sample or biological replication.',
                             'Record comparison covers allele identity, GT, TP/FP/FN decision and refinement-state multiset; not compressed-file byte equality or every annotation.',
                             'Earlier public sensitivity limitations remain. Matching comparator results do not establish full workflow or clinical accuracy.']}
        save(out / 'summary.json', result)
        state.update(status='complete', completed_at=result['completed_at'], summary_sha256=fingerprint(out / 'summary.json')['sha256'])
        save(out / 'status.json', state)
        print(json.dumps({'status': 'complete', 'comparisons': [{k: row[k] for k in ('caller', 'count_checks', 'allele_GT_decision_refinement_records_equal')} for row in results]}))
    except BaseException as error:
        state.update(status='failed', failed_at=stamp(), error_type=type(error).__name__, error=str(error))
        save(out / 'status.json', state)
        raise


if __name__ == '__main__':
    main()
