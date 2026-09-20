#!/usr/bin/env python3
"""Restore captured package versions from a verified wheel set and run existing controls."""
import argparse
import datetime
import json
import os
from pathlib import Path
import subprocess
import time

from export_runtime import PROBE, sha256


def save(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def identities(data):
    return sorted((r['name'].lower().replace('_', '-'), r['version'], r['origin'])
                  for r in data['distributions'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('root', 'export', 'rebuild', 'wheelhouse', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--published-source-build', action='store_true',
                        help='Explicitly validate newly compiled wheels, without original-payload equality claims')
    parser.add_argument('--expected-wheelhouse-sha256')
    args = parser.parse_args()
    root, capture_dir, rebuild, wheelhouse, out = (getattr(args, n).resolve()
        for n in ('root', 'export', 'rebuild', 'wheelhouse', 'output'))
    if out.exists():
        parser.error('Existing Python rebuild verification preserved')
    capture = json.loads((capture_dir / 'runtime.json').read_text())
    conda = json.loads((rebuild / 'status.json').read_text())
    wheels = json.loads((wheelhouse / 'status.json').read_text())
    if conda['status'] != 'complete' or wheels['status'] != 'complete':
        parser.error('Completed Conda rebuild and verified wheel preparation required')
    mode = wheels.get('build_mode', 'preserved_wheels')
    if args.published_source_build != (mode == 'published_sources'):
        parser.error('Newly source-built wheels require the explicit published-source-build mode')
    if args.published_source_build:
        if (not args.expected_wheelhouse_sha256
                or sha256(wheelhouse / 'status.json') != args.expected_wheelhouse_sha256):
            parser.error('Pin the exact completed source-build receipt')
        if (len(wheels.get('builds', [])) != 5
                or {r['name'] for r in wheels['builds']} != {'cnvpytor', 'pyabpoa', 'pytabix', 'pywfa', 'whatshap'}
                or len(wheels['wheels']) != 45):
            parser.error('Unexpected published-source build scope')
    elif mode != 'preserved_wheels':
        parser.error('Unknown wheel preparation mode')
    if conda['export_sha256'] != sha256(capture_dir / 'runtime.json') or wheels['export_sha256'] != conda['export_sha256']:
        parser.error('Environment captures disagree')
    py_pin = next(p for p in conda['runtime_files'] if p['path'] == 'tools/bio/bin/python')
    if sha256(rebuild / 'tools/bio/bin/python') != py_pin['rebuilt_sha256']:
        parser.error('Rebuilt base interpreter changed')
    for row in wheels['wheels']:
        if sha256(wheelhouse / 'wheels' / row['filename']) != row['sha256']:
            parser.error('Wheel changed after preparation')
    for name in capture['venvs']:
        if name != 'bio' and (rebuild / 'tools' / name).exists():
            parser.error('A target Python environment already exists')
    os.nice(10)
    out.mkdir(parents=True)
    start = time.monotonic()
    state = {'status': 'running', 'pid': os.getpid(),
             'started_at': datetime.datetime.now(datetime.timezone.utc).isoformat(),
             'script_sha256': sha256(Path(__file__)), 'export_sha256': conda['export_sha256'],
             'wheelhouse_receipt_sha256': sha256(wheelhouse / 'status.json'),
             'wheel_provenance_mode': mode,
             'environments': {}, 'tests': [],
             'scope': 'Fresh Python environments on the rebuilt native base; existing synthetic controls.'}
    env = dict(os.environ, PYTHONNOUSERSITE='1', PYTHONDONTWRITEBYTECODE='1',
               OMP_NUM_THREADS='2', OPENBLAS_NUM_THREADS='2', VECLIB_MAXIMUM_THREADS='2',
               MVA_PROJECT_ROOT=str(rebuild), MVA_BIO_BIN=str(rebuild / 'tools/bio/bin'))
    for key in ('PYTHONPATH', 'PYTHONHOME', 'VIRTUAL_ENV', 'CONDA_PREFIX'):
        env.pop(key, None)
    env['PATH'] = str(rebuild / 'tools/bio/bin') + ':/usr/bin:/bin:/usr/sbin:/sbin'
    env['MPLCONFIGDIR'] = str(out / 'matplotlib')
    save(out / 'status.json', state)

    def run(command, stem):
        save(out / (stem + '-command.json'), command)
        with (out / (stem + '.log')).open('x') as log:
            result = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT)
        if result.returncode:
            raise RuntimeError(stem + ' failed; inspect retained log')
        return (out / (stem + '.log')).read_text()

    try:
        for name, data in capture['venvs'].items():
            if name == 'bio':
                continue
            prefix = rebuild / 'tools' / name
            command = [str(rebuild / 'tools/bio/bin/python'), '-m', 'venv']
            if data['include_system_site_packages']:
                command += ['--copies', '--system-site-packages']
            command.append(str(prefix))
            run(command, name + '-create')
            python = str(prefix / 'bin/python')
            run([python, '-m', 'pip', '--isolated', 'install', '--disable-pip-version-check',
                 '--no-index', '--find-links', str(wheelhouse / 'wheels'), '--no-deps',
                 '--require-hashes', '--ignore-installed', '--no-compile', '-r',
                 str(wheelhouse / (name + '-local.wheels.txt'))], name + '-install')
            probe = subprocess.run([python, '-c', PROBE], capture_output=True, text=True, env=env, check=True)
            installed = json.loads(probe.stdout)
            if identities(installed) != identities(data):
                save(out / (name + '-identity-mismatch.json'), {'observed': identities(data), 'rebuilt': identities(installed)})
                raise RuntimeError('Python package identity or inheritance differs: ' + name)
            check = run([python, '-m', 'pip', '--isolated', 'check'], name + '-dependency-check')
            state['environments'][name] = {'local_distributions': sum(r['origin'] == 'local' for r in installed['distributions']),
                                            'base_distributions': sum(r['origin'] == 'base' for r in installed['distributions']),
                                            'all_version_and_origin_identities_equal': True,
                                            'dependency_check': check.strip()}
            save(out / 'status.json', state)
            print('[rebuilt Python] ' + name, flush=True)

        cases = [('alignment', 'bio', 'test_align_native.py', 1),
                 ('phasing', 'whatshap-2.8', 'test_read_phasing.py', 6),
                 ('copy-number', 'cnvpytor-1.3.2', 'test_copy_number.py', 7)]
        for label, tool, pattern, count in cases:
            command = [str(rebuild / 'tools' / tool / 'bin/python'), '-m', 'unittest', 'discover',
                       '-s', str(root / 'repo/workflow/tests'), '-p', pattern, '-v']
            test_start = time.monotonic()
            result = run(command, label + '-tests')
            if ('Ran ' + str(count) + ' test') not in result or 'skipped' in result or not result.rstrip().endswith('OK'):
                raise RuntimeError('Expected controls did not all run: ' + label)
            state['tests'].append({'name': label, 'passed': count, 'skipped': 0,
                                   'elapsed_seconds': round(time.monotonic() - test_start, 3),
                                   'test_sha256': sha256(root / 'repo/workflow/tests' / pattern)})
            save(out / 'status.json', state)
        comparator = root / 'repo/workflow/scripts/probe_sv_comparator.py'
        config = root / 'repo/workflow/configs/public-sv-pilot-comparison.json'
        controls = out / 'structural-comparator'
        test_start = time.monotonic()
        run([str(rebuild / 'tools/truvari-5.5.0/bin/python'), str(comparator),
             '--root', str(rebuild), '--output', str(controls), '--config', str(config)], 'structural-comparator')
        result = json.loads((controls / 'summary.json').read_text())
        if not result['all_passed'] or len(result['cases']) != 6:
            raise RuntimeError('Native structural comparator controls failed or changed in count')
        state['tests'].append({'name': 'structural comparator with full-context MAFFT refinement', 'passed': 6,
                               'elapsed_seconds': round(time.monotonic() - test_start, 3),
                               'test_sha256': sha256(comparator), 'config_sha256': sha256(config)})
        state.update(status='complete', completed_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                     elapsed_seconds=round(time.monotonic() - start, 3),
                     limits=['These are existing synthetic execution/interpretation controls, not patient or public-benchmark accuracy.',
                             'Copy-number fixtures model aligned read density, not realignment, biological mosaics or clinical sensitivity.',
                             ('New source-built payload differences are recorded separately; no equality to the old cached wheel bytes is claimed.'
                              if args.published_source_build else
                              'Wheel payload equality omits generated files and metadata, as described in wheel preparation.')])
        save(out / 'status.json', state)
        print(json.dumps({'status': 'complete', 'python_environments': len(state['environments']),
                          'controls_passed': sum(r['passed'] for r in state['tests']), 'elapsed_seconds': state['elapsed_seconds']}))
    except BaseException as error:
        state.update(status='failed', failed_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                     error_type=type(error).__name__, error=str(error))
        save(out / 'status.json', state)
        raise


if __name__ == '__main__':
    main()
