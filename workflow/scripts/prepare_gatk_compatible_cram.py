#!/usr/bin/env python3
"""Preserve a CRAM 3.1 source and verify a CRAM 3.0 copy for pinned GATK.

This changes file encoding, not alignment or caller parameters. A native verifier
compares every decoded record including full flags, tags, sequence, qualities,
coordinates, CIGAR, mate fields and absolute order. Original inputs and failed
runs stay. Samtools checksum helpers below are only for diagnostic controls:
the pinned version did not detect a synthetic duplicate-flag change.
"""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import subprocess
import shutil
import time

import pysam


def stamp():
    return time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())


def digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def save(path, data):
    temporary = Path(str(path) + '.tmp')
    temporary.write_text(json.dumps(data, indent=2) + '\n')
    os.replace(temporary, path)


def cram_version(path):
    with Path(path).open('rb') as stream:
        prefix = stream.read(6)
    if len(prefix) != 6 or prefix[:4] != b'CRAM':
        raise ValueError('expected CRAM input')
    return list(prefix[4:])


def checksum_rows(path):
    rows = [line.split() for line in Path(path).read_text().splitlines()
            if line.strip() and not line.startswith('#')]
    if not rows or any(len(row) != 11 for row in rows):
        raise ValueError('unexpected full checksum report layout')
    keys = [tuple(row[:2]) for row in rows]
    if len(set(keys)) != len(keys):
        raise ValueError('duplicate checksum groups')
    return sorted(rows)


def compare_checksums(before, after, expected_records):
    first, second = checksum_rows(before), checksum_rows(after)
    if first != second:
        raise ValueError('decoded alignment checksum mismatch; read data omitted')
    all_rows = [row for row in first if row[:2] == ['all', 'all']]
    if len(all_rows) != 1 or int(all_rows[0][2]) != expected_records:
        raise ValueError('record count differs from expected source QC')
    return {'records': expected_records, 'group_rows_compared': len(first),
            'all_record_checksums': dict(zip(
                ['flag_sequence', 'name', 'quality', 'all_aux_tags', 'position',
                 'mapq_cigar', 'mate', 'combined'], all_rows[0][3:])),
            'scope': 'Diagnostic samtools CRC comparison; insufficient for full flag equality.'}


def prepare(args):
    source, reference, output, root = [Path(p).resolve() for p in
                                       (args.source, args.reference, args.output, args.root)]
    if output.exists():
        raise ValueError('existing compatibility output preserved; use a new output')
    if args.expected_records < 1 or not 1 <= args.threads <= 8:
        raise ValueError('invalid expected record count or thread limit')
    if cram_version(source) != [3, 1]:
        raise ValueError('this recovery requires a CRAM 3.1 source')
    if digest(source) != args.expected_source_sha256 or digest(reference) != args.expected_reference_sha256:
        raise ValueError('source or reference checksum mismatch')
    if shutil.disk_usage(output.parent).free < source.stat().st_size * 2 + 1024**3:
        raise ValueError('insufficient space for a preserved compatibility copy')
    samtools = root / 'tools/bio/bin/samtools'
    java = root / 'tools/bio/lib/jvm/bin/java'
    jar = root / 'gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar'
    output.mkdir(parents=True, exist_ok=False)
    with (output / '.compatibility.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        state = {'status': 'running', 'pid': os.getpid(), 'started_at': stamp(),
                 'source': str(source), 'output_directory': str(output),
                 'script_sha256': digest(Path(__file__))}
        commands = []
        save(output / 'status.json', state)

        def run(command, name):
            state.update(current_stage=name, stage_started_at=stamp())
            save(output / 'status.json', state)
            entry = {'stage': name, 'command': list(map(str, command)), 'started_at': stamp()}
            commands.append(entry)
            save(output / 'commands.json', commands)
            with (output / (name + '.log')).open('w') as log:
                result = subprocess.run(entry['command'], stdout=log, stderr=subprocess.STDOUT,
                                        pass_fds=(lock.fileno(),))
            entry.update(exit_code=result.returncode, completed_at=stamp())
            save(output / 'commands.json', commands)
            if result.returncode:
                raise RuntimeError('compatibility stage failed: ' + name)

        try:
            verifier_source = Path(__file__).with_name('verify_alignment_records.cpp')
            verifier = output / 'verify-alignment-records'
            runtime = root / 'tools/bio'
            run(['clang++', '-std=c++17', '-O2', '-Wall', '-Wextra', '-I', runtime / 'include',
                 verifier_source, '-L', runtime / 'lib', '-lhts',
                 '-Wl,-rpath,' + str(runtime / 'lib'), '-o', verifier], 'compile-record-verifier')
            partial = output / (args.sample + '.partial.cram')
            final = output / (args.sample + '.cram')
            run([samtools, 'view', '-@', str(args.threads), '--no-PG', '-T', reference,
                 '-O', 'CRAM,version=3.0,store_md=1,store_nm=1', '--write-index',
                 '--save-counts', output / 'conversion-counts.json', '-o', partial, source],
                'convert-cram30')
            if cram_version(partial) != [3, 0]:
                raise ValueError('converted file is not CRAM 3.0')
            run([samtools, 'quickcheck', '-v', partial], 'quickcheck')
            with pysam.AlignmentFile(str(source), reference_filename=str(reference)) as before, \
                    pysam.AlignmentFile(str(partial), reference_filename=str(reference)) as after:
                if before.header.to_dict() != after.header.to_dict():
                    raise ValueError('alignment header changed during conversion')
                if {row.get('SM') for row in after.header.to_dict().get('RG', [])} != {args.sample}:
                    raise ValueError('unexpected sample read groups')
                after.check_index()
                with pysam.FastaFile(str(reference)) as fasta:
                    if list(zip(after.references, after.lengths)) != list(zip(fasta.references, fasta.lengths)):
                        raise ValueError('reference and alignment dictionaries differ')
            run([verifier, reference, source, partial, str(args.expected_records),
                 str(min(args.threads, 3))], 'compare-every-decoded-record')
            lines = (output / 'compare-every-decoded-record.log').read_text().splitlines()
            verified = json.loads(lines[-1])
            if verified.get('records_compared') != args.expected_records or not verified.get('all_compared_fields_equal'):
                raise ValueError('direct decoded-record verification did not complete')
            run([java, '-XX:ActiveProcessorCount=2', '-Xmx2g', '-jar', jar, 'CountReads',
                 '-R', reference, '-I', partial, '-L', args.probe_region], 'gatk-indexed-read-probe')
            # No caller is allowed to depend on an input that changed during recovery.
            if digest(source) != args.expected_source_sha256:
                raise ValueError('original alignment changed during recovery')
            os.replace(partial, final)
            os.replace(str(partial) + '.crai', str(final) + '.crai')
            summary = {'completed_at': stamp(), 'source': str(source),
                       'source_sha256': args.expected_source_sha256,
                       'source_cram_version': [3, 1], 'cram_path': str(final),
                       'cram_sha256': digest(final), 'cram_bytes': final.stat().st_size,
                       'cram_version': [3, 0], 'index_sha256': digest(str(final) + '.crai'),
                       'reference_sha256': args.expected_reference_sha256,
                       'samtools_sha256': digest(samtools), 'gatk_jar_sha256': digest(jar),
                       'script_sha256': digest(Path(__file__)), 'sample': args.sample,
                       'verifier_source_sha256': digest(verifier_source),
                       'verifier_binary_sha256': digest(verifier),
                       'header_equal': True, 'index_readable': True,
                       'gatk_indexed_probe_region': args.probe_region,
                       'decoded_record_verification': verified,
                       'commands_sha256': digest(output / 'commands.json'),
                       'scope': 'Verified file-format compatibility copy; original alignment retained. '
                                'No realignment, filtering, variant classification or clinical conclusion.'}
            save(output / 'summary.json', summary)
            state.update(status='complete', current_stage='verified compatible alignment',
                         completed_at=summary['completed_at'], summary_sha256=digest(output / 'summary.json'))
            save(output / 'status.json', state)
            return summary
        except BaseException as error:
            state.update(status='failed', failed_at=stamp(), error_type=type(error).__name__)
            save(output / 'status.json', state)
            raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('root', 'source', 'reference', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    for name in ('sample', 'expected-source-sha256', 'expected-reference-sha256', 'probe-region'):
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--expected-records', type=int, required=True)
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--call-output', type=Path,
                        help='After successful verification, run the unchanged diploid caller in a new directory.')
    parser.add_argument('--expected-caller-sha256', help='Required with --call-output.')
    args = parser.parse_args()
    os.nice(10)
    if args.call_output and args.call_output.exists():
        parser.error('existing caller output preserved; choose a new directory')
    caller_path = Path(__file__).with_name('call_genome.py')
    if args.call_output and digest(caller_path) != args.expected_caller_sha256:
        parser.error('caller differs from the required source identity')
    summary = prepare(args)
    if args.call_output:
        if digest(caller_path) != args.expected_caller_sha256 or args.call_output.exists():
            raise ValueError('caller identity or output changed during conversion; preserved without launch')
        # call_genome supplies its own definition, lock and complete/failed receipt.
        command = [args.root / 'tools/bio/bin/python', Path(__file__).with_name('call_genome.py'),
                   '--root', args.root, '--reference', args.reference,
                   '--alignment', summary['cram_path'], '--sample', args.sample,
                   '--output', args.call_output, '--workers', '4']
        launch = {'status': 'running', 'started_at': stamp(), 'command': list(map(str, command)),
                  'caller_sha256': digest(Path(__file__).with_name('call_genome.py')),
                  'compatible_alignment_summary_sha256': digest(args.output / 'summary.json')}
        save(args.output / 'calling.json', launch)
        with (args.output / 'calling.log').open('w') as log:
            result = subprocess.run(launch['command'], stdout=log, stderr=subprocess.STDOUT)
        launch.update(status='complete' if result.returncode == 0 else 'failed',
                      finished_at=stamp(), exit_code=result.returncode)
        save(args.output / 'calling.json', launch)
        if result.returncode:
            raise RuntimeError('native caller failed; compatible alignment and logs retained')
    print(json.dumps({'status': 'complete', 'cram_path': summary['cram_path']}))


if __name__ == '__main__':
    main()
