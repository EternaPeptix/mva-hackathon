#!/usr/bin/env python3
"""Rebuild the four captured Conda environments from cached exact archives and test native callers.

Creates one new output tree. It never changes an existing environment or runs on
patient data. The existing synthetic native-caller test is used without edits.
"""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time


def digest(path):
    result = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def write(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def package_identity(prefix):
    return {data['name']: {k: data[k] for k in ('version', 'build', 'subdir', 'url', 'sha256')}
            for p in (prefix / 'conda-meta').glob('*.json')
            for data in [json.loads(p.read_text())]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', required=True, type=Path)
    parser.add_argument('--export', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    root, export, output = (p.resolve() for p in (args.root, args.export, args.output))
    if output.exists():
        parser.error('Existing rebuild preserved; choose a new output')
    captured = json.loads((export / 'runtime.json').read_text())
    for data in captured['conda'].values():
        if digest(export / data['lock']) != data['lock_sha256']:
            parser.error('A captured lock file changed')
    # This run overlaps a large retained-alignment workflow; do not consume its reserve.
    if shutil.disk_usage(root).free < 350_000_000_000:
        parser.error('Insufficient disk headroom for this additional validation run')
    os.nice(10)
    output.mkdir(parents=True)
    state = {'status': 'running', 'pid': os.getpid(),
             'started_at': datetime.datetime.now(datetime.timezone.utc).isoformat(),
             'export_sha256': digest(export / 'runtime.json'),
             'script_sha256': digest(Path(__file__)), 'environments': {},
             'scope': 'Fresh native macOS arm64 Conda prefixes from captured archives; synthetic reads only.'}
    write(output / 'status.json', state)
    started = time.monotonic()
    env = dict(os.environ, PYTHONNOUSERSITE='1', PYTHONDONTWRITEBYTECODE='1',
               OMP_NUM_THREADS='2', OPENBLAS_NUM_THREADS='2', VECLIB_MAXIMUM_THREADS='2')
    env.pop('PYTHONPATH', None)
    env.pop('CONDA_PREFIX', None)
    try:
        for name, data in captured['conda'].items():
            prefix = output / 'tools' / name
            command = [str(root / 'tools/bin/micromamba'), 'create', '--no-rc', '--no-env',
                       '--offline', '--yes', '--always-copy', '--no-pyc',
                       '--safety-checks', 'enabled', '--extra-safety-checks', '--json',
                       '--root-prefix', str(root / 'tools/mamba'), '--prefix', str(prefix),
                       '--file', str(export / data['lock'])]
            stage_started = time.monotonic()
            write(output / (name + '-command.json'), command)
            with (output / (name + '-create.log')).open('x') as log:
                result = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT)
            if result.returncode:
                raise RuntimeError('Native creation failed: ' + name)
            actual = package_identity(prefix)
            expected = {p['name']: {k: p[k] for k in ('version', 'build', 'subdir', 'url', 'sha256')}
                        for p in data['packages']}
            if actual != expected:
                write(output / (name + '-identity-mismatch.json'), {'actual': actual, 'expected': expected})
                raise RuntimeError('Rebuilt package identities differ: ' + name)
            state['environments'][name] = {'package_count': len(actual), 'all_package_identities_equal': True,
                                            'elapsed_seconds': round(time.monotonic() - stage_started, 3)}
            write(output / 'status.json', state)
            print('[rebuilt] ' + name, flush=True)

        jar_rel = 'gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar'
        jar_pin = next(p for p in captured['runtime_files'] if p['path'] == jar_rel)
        if digest(root / jar_rel) != jar_pin['sha256']:
            raise RuntimeError('GATK changed since export')
        (output / jar_rel).parent.mkdir()
        shutil.copy2(root / jar_rel, output / jar_rel)
        if digest(output / jar_rel) != jar_pin['sha256']:
            raise RuntimeError('Rebuilt-root GATK copy differs')
        native = [str(output / 'tools/bio/bin/python'), '-m', 'unittest', 'discover',
                  '-s', str(root / 'repo/workflow/tests'), '-p', 'test_call_native.py', '-v']
        state['test_source_sha256'] = digest(root / 'repo/workflow/tests/test_call_native.py')
        state['caller_source_sha256'] = digest(root / 'repo/workflow/scripts/call_genome.py')
        state['native_test_command'] = native
        env['MVA_PROJECT_ROOT'] = str(output)
        env['PATH'] = str(output / 'tools/bio/bin') + os.pathsep + '/usr/bin:/bin:/usr/sbin:/sbin'
        write(output / 'status.json', state)
        test_started = time.monotonic()
        with (output / 'native-caller-tests.log').open('x') as log:
            test = subprocess.run(native, env=env, stdout=log, stderr=subprocess.STDOUT)
        test_text = (output / 'native-caller-tests.log').read_text()
        if test.returncode or 'Ran 2 tests' not in test_text or 'skipped' in test_text or not test_text.rstrip().endswith('OK'):
            raise RuntimeError('Native caller tests failed, skipped or did not execute both checks')
        state['native_caller_tests'] = {'tests': 2, 'skipped': 0, 'exit_code': test.returncode,
                                         'elapsed_seconds': round(time.monotonic() - test_started, 3)}
        state['runtime_files'] = []
        for rec in captured['runtime_files']:
            p = output / rec['path']
            if p.is_file():
                new_hash = digest(p)
                state['runtime_files'].append({'path': rec['path'], 'original_sha256': rec['sha256'],
                                               'rebuilt_sha256': new_hash, 'identical_bytes': new_hash == rec['sha256']})
        state.update(status='complete', completed_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                     elapsed_seconds=round(time.monotonic() - started, 3),
                     limits=['This verifies four Conda environments, a pinned GATK copy, and a two-contig synthetic caller experiment.',
                             'Pip environments, public pilot accuracy, patient calls and full scientific reproduction are not validated by these checks.',
                             'Relocated installed executable bytes can differ; package archives remain pinned and compared by identity.'])
        write(output / 'status.json', state)
        print(json.dumps({'status': 'complete', 'output': str(output), 'elapsed_seconds': state['elapsed_seconds']}))
    except BaseException as error:
        state.update(status='failed', error_type=type(error).__name__, error=str(error),
                     failed_at=datetime.datetime.now(datetime.timezone.utc).isoformat())
        write(output / 'status.json', state)
        raise


if __name__ == '__main__':
    main()
