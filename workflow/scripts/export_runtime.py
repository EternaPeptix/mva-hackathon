#!/usr/bin/env python3
"""Export the observed Mac runtime without modifying its environments or reading case data.

Conda exports pin archive URLs and SHA-256 values. Pip exports pin installed local
versions and describe inherited packages; they are not hash-locked wheel sets.
"""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import subprocess
from urllib.parse import urlsplit


CONDA_ENVS = ('bio', 'delly-2.6.0', 'mafft-7.526', 'svaba-1.2.0')
VENV_ENVS = ('cnvpytor-1.3.2', 'truvari-5.5.0', 'whatshap-2.8')
PROBE = r'''
import hashlib, importlib.metadata as metadata, json, pathlib, sys
prefix = pathlib.Path(sys.prefix).resolve()
base = pathlib.Path(sys.base_prefix).resolve()
rows = []
seen = set()
for dist in metadata.distributions():
    path = pathlib.Path(dist._path).resolve()
    if path in seen:
        continue
    seen.add(path)
    if path.is_relative_to(prefix):
        origin = 'local'
        relative = str(path.relative_to(prefix))
    elif path.is_relative_to(base):
        origin = 'base'
        relative = str(path.relative_to(base))
    else:
        raise RuntimeError('Unexpected distribution outside the selected runtime')
    direct = dist.read_text('direct_url.json')
    wheel = dist.read_text('WHEEL')
    meta = dist.read_text('METADATA') or dist.read_text('PKG-INFO')
    rows.append({'name': dist.metadata['Name'], 'version': dist.version,
                 'origin': origin, 'metadata_path': relative,
                 'installer': (dist.read_text('INSTALLER') or '').strip() or None,
                 'direct_url': json.loads(direct) if direct else None,
                 'wheel_metadata': wheel,
                 'metadata_sha256': hashlib.sha256(meta.encode()).hexdigest() if meta else None})
print(json.dumps({'python': sys.version, 'base_python_version': sys.version.split()[0],
                  'distributions': sorted(rows, key=lambda x: (x['name'].lower(), x['origin']))}))
'''


def sha256(path):
    before = path.stat()
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            value.update(block)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise RuntimeError('File changed during capture: ' + path.name)
    return value.hexdigest()


def write_json(path, data):
    with path.open('x') as stream:
        json.dump(data, stream, indent=2)
        stream.write('\n')


def package_url(value):
    url = urlsplit(value)
    if (url.scheme != 'https' or url.hostname != 'conda.anaconda.org'
            or url.username or url.password or url.query or url.fragment
            or url.path.split('/')[1] not in ('bioconda', 'conda-forge')):
        raise ValueError('Unexpected package origin; inspect locally before exporting')
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    out = args.output.resolve()
    if out.exists():
        parser.error('Existing export preserved; choose a new output')
    if platform.system() != 'Darwin' or platform.machine() != 'arm64':
        parser.error('This export describes the observed native Apple Silicon runtime')
    os.nice(10)
    out.mkdir(parents=True)
    env = dict(os.environ, PYTHONNOUSERSITE='1', PYTHONDONTWRITEBYTECODE='1')
    env.pop('PYTHONPATH', None)
    checked_archives = {}
    manifest = {'captured_at': datetime.datetime.now(datetime.timezone.utc).isoformat(),
                'platform': {'system': platform.system(), 'architecture': platform.machine(),
                             'macos': platform.mac_ver()[0]},
                'script_sha256': sha256(Path(__file__)), 'conda': {}, 'venvs': {},
                'limits': ['Installed Conda metadata supplies archive identity; cached archive hashes are checked separately.',
                           'Pip version pins do not identify original wheels or source builds unless direct provenance exists.',
                           'This export alone is not a successful rebuild or end-to-end scientific reproduction.',
                           'No raw sequencing, genome-wide tables, patient records or account settings are inspected.']}
    for name in CONDA_ENVS:
        prefix = root / 'tools' / name
        records = []
        for source in sorted((prefix / 'conda-meta').glob('*.json')):
            data = json.loads(source.read_text())
            url = package_url(data['url'])
            expected = data['sha256']
            if not re.fullmatch('[0-9a-f]{64}', expected):
                raise ValueError('Missing archive SHA-256 for ' + data['name'])
            archive = Path(data['package_tarball_full_path']).parent / urlsplit(url).path.rsplit('/', 1)[1]
            actual = None
            if archive.is_file():
                if str(archive) not in checked_archives:
                    checked_archives[str(archive)] = sha256(archive)
                actual = checked_archives[str(archive)]
                if actual != expected:
                    raise ValueError('Cached archive differs from installed record: ' + data['name'])
            records.append({key: data[key] for key in ('name', 'version', 'build', 'subdir', 'url', 'sha256', 'md5')})
            records[-1].update(cached_archive_verified=actual == expected,
                               installed_metadata_sha256=sha256(source))
        if not records or len({r['name'] for r in records}) != len(records):
            raise ValueError('Missing or duplicate Conda package records: ' + name)
        records.sort(key=lambda r: r['name'])
        lock = out / (name + '-osx-arm64.explicit.txt')
        with lock.open('x') as stream:
            stream.write('# Native macOS arm64; exact observed Conda archives.\n@EXPLICIT\n')
            for record in records:
                stream.write(record['url'] + '#' + record['sha256'] + '\n')
        manifest['conda'][name] = {'packages': records, 'package_count': len(records),
                                    'lock': lock.name, 'lock_sha256': sha256(lock),
                                    'cached_archives_verified': sum(r['cached_archive_verified'] for r in records)}

    for name in ('bio',) + VENV_ENVS:
        prefix = root / 'tools' / name
        probe = subprocess.run([str(prefix / 'bin/python'), '-c', PROBE], env=env,
                               capture_output=True, text=True, timeout=30, check=True)
        data = json.loads(probe.stdout)
        for row in data['distributions']:
            direct = row['direct_url']
            if direct:
                url = urlsplit(direct['url'])
                if url.scheme == 'file':
                    direct['url'] = 'local-archive:' + Path(url.path).name
                elif url.scheme != 'https' or url.username or url.password or url.query:
                    raise ValueError('Unreviewed Python distribution origin')
        config = prefix / 'pyvenv.cfg'
        include_base = False
        if config.exists():
            settings = dict(line.split(' = ', 1) for line in config.read_text().splitlines() if ' = ' in line)
            include_base = settings['include-system-site-packages'].lower() == 'true'
        data['include_system_site_packages'] = include_base
        data['base_environment'] = 'bio' if name != 'bio' else None
        manifest['venvs'][name] = data
        if name != 'bio':
            local = [row for row in data['distributions'] if row['origin'] == 'local']
            pin = out / (name + '-local.versions.txt')
            with pin.open('x') as stream:
                stream.write('# Observed LOCAL distributions, including pip; version pins, not archive hashes.\n')
                stream.write('# Recreate the documented base-environment inheritance first.\n')
                for row in sorted(local, key=lambda r: r['name'].lower()):
                    stream.write(row['name'] + '==' + row['version'] + '\n')
            data.update(local_pin_file=pin.name, local_pin_sha256=sha256(pin))

    binaries = [root / 'tools' / name / binary for name, binary in (
        ('bin', 'micromamba'), ('bio', 'bin/python'), ('bio', 'bin/bwa'),
        ('bio', 'bin/samtools'), ('bio', 'bin/bcftools'), ('bio', 'lib/jvm/bin/java'),
        ('delly-2.6.0', 'bin/delly'), ('svaba-1.2.0', 'bin/svaba'),
        ('mafft-7.526', 'bin/mafft'), ('rtg-tools-3.13', 'rtg'),
        ('rtg-tools-3.13', 'RTG.jar'), ('rtg-tools-3.13', 'rtg.cfg'))]
    jars = []
    for base in (root / 'gatk-4.6.1.0', root / 'exomiser/exomiser-cli-15.1.0', root / 'lirical'):
        jars.extend(sorted(base.rglob('*.jar')))
    manifest['runtime_files'] = [{'path': str(p.relative_to(root)), 'bytes': p.stat().st_size,
                                  'sha256': sha256(p)} for p in binaries + jars]
    manifest['verified_unique_cached_archives'] = len(checked_archives)
    write_json(out / 'runtime.json', manifest)
    print(json.dumps({'output': str(out), 'conda_environments': len(CONDA_ENVS),
                      'conda_package_records': sum(x['package_count'] for x in manifest['conda'].values()),
                      'venv_environments': len(VENV_ENVS), 'unique_cached_archives': len(checked_archives),
                      'runtime_files': len(manifest['runtime_files'])}))


if __name__ == '__main__':
    main()
