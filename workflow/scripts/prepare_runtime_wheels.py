#!/usr/bin/env python3
"""Pin wheels for the captured Python environments and compare their installed code payloads.

Uses already-built matching wheels from pip's cache, then downloads remaining
binary distributions from public PyPI. No environment is installed or changed.
"""
import argparse
import datetime
import email
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import zipfile

from packaging.utils import canonicalize_name, parse_wheel_filename


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write(path, data):
    path.write_text(json.dumps(data, indent=2) + '\n')


def wheel_identity(path):
    name, version, _, _ = parse_wheel_filename(path.name)
    with zipfile.ZipFile(path) as archive:
        members = [n for n in archive.namelist() if n.endswith('.dist-info/METADATA')]
        if len(members) != 1:
            raise ValueError('Unexpected wheel metadata count')
        metadata = email.message_from_bytes(archive.read(members[0]))
        if canonicalize_name(metadata['Name']) != name or metadata['Version'] != str(version):
            raise ValueError('Wheel filename and metadata disagree')
    return name, str(version)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--export', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root, export, out = (p.resolve() for p in (args.root, args.export, args.output))
    if out.exists():
        parser.error('Existing wheel preparation preserved; choose a new output')
    capture = json.loads((export / 'runtime.json').read_text())
    wanted = {}
    for name, env in capture['venvs'].items():
        if name == 'bio':
            continue
        for row in env['distributions']:
            if row['origin'] == 'local':
                key = canonicalize_name(row['name']), row['version']
                wanted.setdefault(key, []).append((name, row))
    os.nice(10)
    out.mkdir(parents=True)
    wheels = out / 'wheels'
    wheels.mkdir()
    state = {'status': 'running', 'pid': os.getpid(),
             'started_at': datetime.datetime.now(datetime.timezone.utc).isoformat(),
             'script_sha256': digest(Path(__file__)),
             'export_sha256': digest(export / 'runtime.json'), 'wheels': [], 'payload_checks': []}
    write(out / 'status.json', state)
    env = dict(os.environ, PYTHONNOUSERSITE='1', PYTHONDONTWRITEBYTECODE='1')
    for key in ('PYTHONPATH', 'PYTHONHOME'):
        env.pop(key, None)
    python = root / 'tools/bio/bin/python'
    try:
        cache = subprocess.run([str(python), '-m', 'pip', '--isolated', 'cache', 'list', '--format=abspath'],
                               capture_output=True, text=True, env=env, check=True)
        cache_records = {}
        for line in cache.stdout.splitlines():
            p = Path(line)
            if p.suffix != '.whl' or not p.is_file():
                continue
            key = wheel_identity(p)
            if key not in wanted:
                continue
            target = wheels / p.name
            if target.exists() and digest(target) != digest(p):
                raise ValueError('Conflicting cached wheels')
            if not target.exists():
                shutil.copy2(p, target)
            origin_path = p.parent / 'origin.json'
            origin = json.loads(origin_path.read_text()) if origin_path.exists() else None
            if origin and origin.get('url', '').startswith('file://'):
                origin['url'] = 'local-archive:' + Path(origin['url']).name
            cache_records[p.name] = {'source': 'pre-existing locally built wheel', 'source_origin': origin}
        pins = out / 'all-local.versions.txt'
        pins.write_text(''.join(name + '==' + version + '\n' for name, version in sorted(wanted)))
        command = [str(python), '-m', 'pip', '--isolated', 'download', '--disable-pip-version-check',
                   '--index-url', 'https://pypi.org/simple', '--only-binary=:all:', '--no-deps',
                   '--find-links', str(wheels), '--dest', str(wheels), '-r', str(pins)]
        write(out / 'download-command.json', command)
        with (out / 'download.log').open('x') as log:
            subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
        available = {}
        for p in sorted(wheels.glob('*.whl')):
            key = wheel_identity(p)
            if key not in wanted or key in available:
                raise ValueError('Unexpected or duplicate wheel identity')
            available[key] = p
            state['wheels'].append({'name': key[0], 'version': key[1], 'filename': p.name,
                                    'bytes': p.stat().st_size, 'sha256': digest(p),
                                    **cache_records.get(p.name, {'source': 'public PyPI binary selected for this native runtime'})})
        if set(available) != set(wanted):
            raise ValueError('Wheel set does not cover every captured local distribution')
        mismatches = []
        for key, owners in sorted(wanted.items()):
            with zipfile.ZipFile(available[key]) as archive:
                for owner, row in owners:
                    site = root / 'tools' / owner / Path(row['metadata_path']).parent
                    checked, ignored = 0, []
                    for member in archive.infolist():
                        if member.is_dir():
                            continue
                        path = PurePosixPath(member.filename)
                        if path.is_absolute() or '..' in path.parts:
                            raise ValueError('Unsafe wheel member path')
                        if path.parts[0].endswith('.dist-info'):
                            continue
                        if path.parts[0].endswith('.data'):
                            if len(path.parts) < 3 or path.parts[1] not in ('purelib', 'platlib'):
                                ignored.append(member.filename)
                                continue
                            path = PurePosixPath(*path.parts[2:])
                        installed = site / str(path)
                        wheel_hash = hashlib.sha256(archive.read(member)).hexdigest()
                        if not installed.is_file() or digest(installed) != wheel_hash:
                            mismatches.append({'environment': owner, 'package': key[0],
                                               'member': member.filename,
                                               'reason': 'missing' if not installed.is_file() else 'different_bytes'})
                        checked += 1
                    state['payload_checks'].append({'environment': owner, 'name': key[0],
                                                     'version': key[1], 'package_payload_files': checked,
                                                     'uncompared_data_or_scripts': ignored})
        write(out / 'payload-mismatches.json', mismatches)
        if mismatches:
            raise ValueError('Wheel payload differs from installed runtime; inspect preserved mismatch details')
        for name, env_data in capture['venvs'].items():
            if name == 'bio':
                continue
            lines = ['# Exact locally verified wheel artifacts; use --no-index --find-links and --require-hashes.']
            for row in env_data['distributions']:
                if row['origin'] != 'local':
                    continue
                key = canonicalize_name(row['name']), row['version']
                lines.append(key[0] + '==' + key[1] + ' --hash=sha256:' + digest(available[key]))
            (out / (name + '-local.wheels.txt')).write_text('\n'.join(lines) + '\n')
        state.update(status='complete', completed_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                     wheel_count=len(available), payload_file_comparisons=sum(p['package_payload_files'] for p in state['payload_checks']),
                     mismatches=0,
                     limits=['Package code and data under site-packages are compared byte for byte to the observed runtime.',
                             'Distribution metadata, generated launchers, bytecode and wheel .data scripts/other destinations are not covered by that equality.',
                             'Cached source-built wheels retain origin records; their hashes are not claims of publisher-built wheel identity.',
                             'This preparation does not yet prove a rebuilt environment runs correctly.'])
        write(out / 'status.json', state)
        print(json.dumps({k: state[k] for k in ('status', 'wheel_count', 'payload_file_comparisons', 'mismatches')}))
    except BaseException as error:
        state.update(status='failed', failed_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                     error_type=type(error).__name__, error=str(error))
        write(out / 'status.json', state)
        raise


if __name__ == '__main__':
    main()
