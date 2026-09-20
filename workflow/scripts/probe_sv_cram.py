#!/usr/bin/env python3
"""Check synthetic SV evidence and native calls through the intended CRAM format."""
import argparse
from collections import Counter
import itertools
import json
import os
from pathlib import Path
import subprocess
import time

import pysam

from probe_structural_variants import read_calls, save, stamp
from simulate_sv_reads import sha


def read_identity(read):
    # CRAM may regenerate MD/NM. Preserve every other tag, alignment field and base.
    return (read.query_name, read.flag, read.reference_id, read.reference_start,
            read.mapping_quality, read.cigartuples, read.next_reference_id,
            read.next_reference_start, read.template_length, read.query_sequence,
            read.query_qualities, sorted((key, value) for key, value in read.get_tags()
                                         if key not in ('MD', 'NM')))


def run(experiment, output, samtools, delly):
    if output.exists():
        raise ValueError('CRAM comparison output already exists')
    original = json.loads((experiment / 'summary.json').read_text())
    truth_path = experiment / 'fixture/truth.json'
    truth = json.loads(truth_path.read_text())
    reference = experiment / 'fixture/reference.fa'; bam = experiment / 'aligned.bam'
    if truth['scope'] != 'Synthetic diploid simple-SV reference; no patient sequence or observed biological effect.':
        raise ValueError('synthetic fixture required')
    if original['truth_sha256'] != sha(truth_path) or sha(bam) != original['output_sha256']['aligned.bam']:
        raise ValueError('retained experiment input identity mismatch')
    if sha(reference) != truth['files']['reference.fa']['sha256']:
        raise ValueError('synthetic reference identity mismatch')
    if sha(delly) != original['tools']['delly']['sha256'] or sha(samtools) != original['tools']['samtools']['sha256']:
        raise ValueError('native tool differs from original experiment')
    output.mkdir(parents=True, exist_ok=False); start = time.monotonic(); commands = []
    state = {'status': 'running', 'started_at': stamp(), 'pid': os.getpid()}
    save(output / 'status.json', state)

    def command(args):
        log = output / ('step-%02d.log' % len(commands)); tick = time.monotonic()
        with log.open('w') as stream:
            result = subprocess.run(list(map(str, args)), stdout=stream, stderr=subprocess.STDOUT)
        commands.append({'command': list(map(str, args)), 'exit_code': result.returncode,
                         'elapsed_seconds': time.monotonic() - tick, 'log_sha256': sha(log)})
        save(output / 'commands.json', commands)
        if result.returncode:
            raise RuntimeError('native synthetic CRAM stage failed')

    try:
        cram = output / 'aligned.cram'
        command([samtools, 'view', '-C', '-@', '1', '-T', reference, '-o', cram, bam])
        command([samtools, 'index', cram])
        command([samtools, 'quickcheck', '-v', cram])
        counts = Counter()
        with pysam.AlignmentFile(str(bam)) as before, \
                pysam.AlignmentFile(str(cram), reference_filename=str(reference)) as after:
            if before.references != after.references or before.lengths != after.lengths:
                raise ValueError('reference dictionary changed in CRAM')
            if before.header.to_dict()['RG'] != after.header.to_dict()['RG']:
                raise ValueError('read-group metadata changed in CRAM')
            for a, b in itertools.zip_longest(before, after):
                if a is None or b is None or read_identity(a) != read_identity(b):
                    raise ValueError('read evidence changed in CRAM; identities omitted from error')
                counts['records_compared'] += 1
                if a.is_supplementary: counts['supplementary_records_compared'] += 1
                if a.is_duplicate: counts['duplicate_records_compared'] += 1
        native = [row for row in original['commands'] if row['stage'] == 'native structural discovery']
        if len(native) != 1:
            raise ValueError('original discovery command not unique')
        args = list(native[0]['command'])
        if args[0] != str(delly) or args[-1] != str(bam):
            raise ValueError('unexpected original discovery paths')
        args[-1] = str(cram); args[args.index('-o') + 1] = str(output / 'calls.bcf')
        command(args)
        calls = json.loads(json.dumps(read_calls(output / 'calls.bcf', 'SYNTH_SV')))
        if calls != original['native_calls']:
            save(output / 'discordant-calls.json', {'BAM': original['native_calls'], 'CRAM': calls})
            raise ValueError('native SV records differ between BAM and CRAM')
        result = {'completed_at': stamp(), 'elapsed_seconds': time.monotonic() - start,
                  'probe_source_sha256': sha(Path(__file__)),
                  'parser_source_sha256': sha(Path(__file__).with_name('probe_structural_variants.py')),
                  'original_summary_sha256': sha(experiment / 'summary.json'),
                  'truth_sha256': sha(truth_path), 'reference_sha256': sha(reference),
                  'original_bam_sha256': sha(bam), 'commands': commands,
                  'read_identity_counts': dict(counts), 'CRAM_same_read_evidence': True,
                  'compared_fields': 'Every read name, sequence, base quality, flag, coordinate, CIGAR, MQ, '
                                     'mate field and tag except regenerable MD/NM; reference dictionary and RG metadata.',
                  'native_call_count': len(calls), 'native_call_fields_identical': True,
                  'native_calls': calls,
                  'outputs': {name: sha(output / name) for name in
                              ('aligned.cram', 'aligned.cram.crai', 'calls.bcf', 'calls.bcf.csi')},
                  'scope': 'Same synthetic reads and caller options through BAM/CRAM. This validates the '
                           'tested format path, not new biological evidence or patient SV sensitivity.'}
        save(output / 'summary.json', result)
        state.update(status='complete', completed_at=result['completed_at']); save(output / 'status.json', state)
        return result
    except BaseException as error:
        state.update(status='failed', failed_at=stamp(), error_type=type(error).__name__)
        save(output / 'status.json', state); raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('experiment', 'output', 'samtools', 'delly'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args(); os.nice(10)
    result = run(args.experiment, args.output, args.samtools, args.delly)
    print(json.dumps({key: result[key] for key in ('completed_at', 'elapsed_seconds',
                      'read_identity_counts', 'native_call_count', 'native_call_fields_identical')}, indent=2))


if __name__ == '__main__':
    main()
