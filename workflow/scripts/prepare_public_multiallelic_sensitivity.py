#!/usr/bin/env python3
"""Preserve original multiallelic calls for an exploratory representation sensitivity.

Replays the frozen split normalization and requires identical complete record
payloads. A separate normalization retains multiallelic genotypes. No caller,
reference, matching option, QUAL threshold or original result is changed.
"""
import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import time

import pysam

from rtg_small_variant_helpers import audit_vcf, execute, fingerprint, save, stamp


def record_fingerprint(path):
    payload = hashlib.sha256()
    formats = Counter()
    parsed_qualities = Counter()
    count = 0
    with pysam.VariantFile(str(path)) as vcf:
        if list(vcf.header.samples) != ['HG002']:
            raise ValueError('Unexpected sample')
        for record in vcf:
            if record.contig != 'chr20':
                raise ValueError('Unexpected non-chr20 record')
            rendered = str(record)
            payload.update(rendered.encode())
            sample = record.samples['HG002']
            rendered_qual = rendered.split('\t')[5]
            formats[json.dumps({'FORMAT': dict(sample), 'phased': sample.phased,
                                'serialized_QUAL': rendered_qual, 'FILTER': list(record.filter)}, sort_keys=True)] += 1
            parsed_qualities[(rendered_qual, record.qual)] += 1
            count += 1
    return {'records': count, 'record_payload_sha256': payload.hexdigest()}, formats, parsed_qualities


def verify_metadata(source_records, source_formats, source_quals,
                    retained_records, retained_formats, retained_quals):
    if retained_records['records'] != source_records['records'] or retained_formats != source_formats:
        raise ValueError('Normalization changed record count or sample FORMAT/serialized QUAL/FILTER multiset')
    # htslib renders VCF QUAL with limited significant digits. For example,
    # 10289.64 is rendered as 10289.6; parsing the emitted VCF then produces a
    # different float32. Require exact equality of the rendered tokens, not an
    # arbitrary numeric tolerance that could hide a real score change.
    return {
        'FORMAT_serialized_QUAL_FILTER_multiset_preserved': True,
        'parsed_QUAL_numeric_multiset_records_changed': sum((source_quals - retained_quals).values()),
        'QUAL_policy': 'Exact htslib-rendered QUAL token equality; binary parsed numeric changes from VCF serialization are counted separately. No numeric tolerance or score tuning.',
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('root', 'original-comparison', 'controls', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    root, output = args.root.resolve(), args.output.resolve()
    if output.exists():
        parser.error('Existing representation study preserved')
    original = json.loads((args.original_comparison / 'summary.json').read_text())
    controls = json.loads((args.controls / 'summary.json').read_text())
    if original['status'] != 'complete' or not controls['all_prespecified_expectations_met']:
        parser.error('Completed baseline and successful literal representation controls required')
    pilot = root / 'results/germline-public-pilot-20260919'
    pilot_summary = json.loads((pilot / 'summary.json').read_text())
    native_status = json.loads((pilot / 'calls/status.json').read_text())
    if pilot_summary['status'] != 'complete' or native_status['status'] != 'complete':
        parser.error('Completed original native calling required')
    reference = Path(original['inputs']['reference']['path'])
    if fingerprint(reference)['sha256'] != original['inputs']['reference']['sha256']:
        parser.error('Public reference changed')
    shard = pilot / 'calls/shard-0019'
    completed = json.loads((shard / 'complete.json').read_text())
    output.mkdir(parents=True)
    (output / ('executed-' + Path(__file__).name)).write_bytes(Path(__file__).read_bytes())
    os.nice(10)
    started = time.monotonic()
    state = {'status': 'running', 'pid': os.getpid(), 'started_at': stamp(),
             'scope': 'Exploratory post-comparison sensitivity: left-normalized multiallelic records retained; no new variant calling.'}
    save(output / 'status.json', state)
    try:
        results = {}
        verification = {}
        for method, raw_name in (('gatk', 'gatk.raw.vcf.gz'), ('bcftools', 'bcftools.raw.bcf')):
            frozen_shard = shard / completed[method]['path']
            if fingerprint(frozen_shard)['sha256'] != completed[method]['sha256']:
                raise ValueError('Completed normalized shard changed')
            attempt = frozen_shard.parent
            if fingerprint(attempt / 'commands.json')['sha256'] != completed['commands_sha256']:
                raise ValueError('Native command history changed')
            source = attempt / raw_name
            source_pin = fingerprint(source)
            source_records, source_formats, source_quals = record_fingerprint(source)
            replay = output / (method + '.replayed-split.vcf.gz')
            bcftools = root / 'tools/bio/bin/bcftools'
            split_command = execute([bcftools, 'norm', '--check-ref', 'e', '-f', reference,
                                     '-m', '-any', '-Oz', '-o', replay, source], output, method + '-split-replay')
            frozen_whole = Path(original['inputs'][method]['path'])
            if fingerprint(frozen_whole)['sha256'] != original['inputs'][method]['sha256']:
                raise ValueError('Completed gathered calls changed')
            replay_records, _, _ = record_fingerprint(replay)
            frozen_records, _, _ = record_fingerprint(frozen_whole)
            if replay_records != frozen_records:
                raise ValueError('Raw sources no longer reproduce all frozen split record payloads')
            retained = output / (method + '.normalized.multiallelic.vcf.gz')
            keep_command = execute([bcftools, 'norm', '--check-ref', 'e', '-f', reference,
                                    '-Oz', '-o', retained, source], output, method + '-keep-multiallelic')
            execute([bcftools, 'index', '-t', retained], output, method + '-index')
            retained_records, retained_formats, retained_quals = record_fingerprint(retained)
            metadata_check = verify_metadata(source_records, source_formats, source_quals,
                                             retained_records, retained_formats, retained_quals)
            audit = audit_vcf(retained, reference, 'HG002', 'chr20', require_raw_qual=True)
            results[method] = {**fingerprint(retained), 'input_audit': audit}
            verification[method] = {'original_raw': source_pin, 'original_raw_records': source_records,
                                    'frozen_split_records': frozen_records,
                                    'replayed_split_record_payloads_identical': True,
                                    'retained_multiallelic_records': retained_records,
                                    'metadata_check': metadata_check,
                                    'replay_command': split_command, 'retain_command': keep_command}
        summary = {'status': 'complete', 'completed_at': stamp(), 'elapsed_seconds': time.monotonic() - started,
                   'scope': state['scope'], 'outputs': results, 'verification': verification,
                   'original_comparison': fingerprint(args.original_comparison / 'summary.json'),
                   'representation_controls': fingerprint(args.controls / 'summary.json'),
                   'original_pilot': fingerprint(pilot / 'summary.json'),
                   'config_sha256': original['config']['sha256'], 'code': fingerprint(Path(__file__)),
                   'limits': ['This arm was motivated by observed baseline discordances and is explicitly exploratory.',
                              'Original split-call results remain primary recorded evidence and are not overwritten.',
                              'Reference left-normalization is retained; only multiallelic splitting is omitted.',
                              'FORMAT/rendered QUAL/FILTER checks are multiset checks, not independent proof of every variant normalization.',
                              'VCF serialization can round the underlying numeric QUAL value; exact rendered tokens must agree and affected numeric values are counted. The frozen threshold is zero.',
                              'The native genotype-aware comparator and fixed score/region policy remain unchanged.']}
        save(output / 'summary.json', summary)
        state.update(status='complete', completed_at=summary['completed_at'])
        save(output / 'status.json', state)
        print(json.dumps({'status': 'complete', 'input_records': {m: v['original_raw_records']['records'] for m, v in verification.items()},
                          'split_replay_identical': True, 'elapsed_seconds': summary['elapsed_seconds']}))
    except BaseException as error:
        state.update(status='failed', failed_at=stamp(), error_type=type(error).__name__, error=str(error))
        save(output / 'status.json', state)
        raise


if __name__ == '__main__':
    main()
