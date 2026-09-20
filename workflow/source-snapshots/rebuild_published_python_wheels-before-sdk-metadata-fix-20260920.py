#!/usr/bin/env python3
"""Obtain pinned publisher wheels and rebuild the five source-only dependencies.

Creates a new build environment and wheelhouse without changing analysis runtimes.
Native files shipped inside a source archive are omitted from its separate build
copy; original archives remain intact. New wheels require separate runtime controls.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import datetime
import email
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import platform
import re
import subprocess
import tarfile
import time
from urllib.parse import urlsplit
from urllib.request import Request, urlopen
import zipfile


def digest(path):
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def canonical(name):
    return re.sub(r'[-_.]+', '-', name).lower()


def native(data):
    return data[:4] in (b'\x7fELF', b'\xfe\xed\xfa\xce', b'\xce\xfa\xed\xfe',
        b'\xfe\xed\xfa\xcf', b'\xcf\xfa\xed\xfe', b'\xca\xfe\xba\xbe', b'\xbe\xba\xfe\xca') or data[:8] == b'!<arch>\n'


def fetch(url, filename, size, expected, folder):
    origin = urlsplit(url)
    require(origin.scheme == 'https' and origin.hostname == 'files.pythonhosted.org'
            and not origin.username and not origin.query and not origin.fragment, 'unreviewed artifact host')
    require(Path(filename).name == filename and filename not in ('', '.', '..'), 'unsafe artifact filename')
    path = folder / filename
    with urlopen(Request(url, headers={'User-Agent': 'MVA-source-build-reproduction/1.0'}), timeout=45) as stream, path.open('xb') as target:
        observed = 0
        while True:
            block = stream.read(1024 * 1024)
            if not block:
                break
            observed += len(block)
            require(observed <= size, 'artifact exceeds expected size')
            target.write(block)
    require(path.stat().st_size == size and digest(path) == expected, 'publisher artifact identity differs')
    return path


def wheel_info(path):
    with zipfile.ZipFile(path) as archive:
        names = archive.namelist()
        require(len(names) == len(set(names)), 'duplicate wheel member')
        for name in names:
            p = PurePosixPath(name)
            require(not p.is_absolute() and '..' not in p.parts, 'unsafe wheel member')
        # Bundled dependencies (for example inside pip/_vendor) can retain their
        # own metadata. Only the top-level dist-info describes this wheel.
        metadata = [name for name in names if name.endswith('.dist-info/METADATA')
                    and len(PurePosixPath(name).parts) == 2]
        wheel = [name for name in names if name.endswith('.dist-info/WHEEL')
                 and len(PurePosixPath(name).parts) == 2]
        require(len(metadata) == len(wheel) == 1, 'wheel metadata count differs')
        record = email.message_from_bytes(archive.read(metadata[0]))
        return {'name': canonical(record['Name']), 'version': record['Version'], 'filename': path.name,
                'bytes': path.stat().st_size, 'sha256': digest(path),
                'wheel_metadata': archive.read(wheel[0]).decode()}


def clean_source(archive_path, folder):
    folder.mkdir()
    omitted, kept, tops = [], [], set()
    with tarfile.open(archive_path, 'r:gz') as archive:
        for member in archive.getmembers():
            name = PurePosixPath(member.name)
            require(not name.is_absolute() and '..' not in name.parts and len(name.parts) >= 1, 'unsafe source member')
            tops.add(name.parts[0])
            if member.isdir():
                (folder / str(name)).mkdir(parents=True, exist_ok=True)
                continue
            require(member.isfile(), 'source link or special member requires review')
            raw = archive.extractfile(member).read()
            identity = {'path': str(name), 'bytes': len(raw), 'sha256': hashlib.sha256(raw).hexdigest()}
            if native(raw) or name.suffix in ('.o', '.so', '.dylib', '.a', '.dll', '.pyd'):
                omitted.append(identity)
                continue
            target = folder / str(name)
            target.parent.mkdir(parents=True, exist_ok=True)
            with target.open('xb') as stream:
                stream.write(raw)
            target.chmod(member.mode & 0o777)
            os.utime(target, (member.mtime, member.mtime))
            kept.append(identity)
    require(len(tops) == 1 and kept, 'source archive root differs')
    return folder / next(iter(tops)), {'retained_source_files': len(kept),
        'retained_source_manifest_sha256': hashlib.sha256(json.dumps(kept, sort_keys=True).encode()).hexdigest(),
        'omitted_prebuilt_members': omitted}


def payloads(path):
    with zipfile.ZipFile(path) as archive:
        return {m.filename: hashlib.sha256(archive.read(m)).hexdigest() for m in archive.infolist()
                if not m.is_dir() and '.dist-info/' not in m.filename}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('source-receipt', 'build-tools', 'runtime-export', 'python', 'original-wheelhouse', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    for name in ('source-receipt', 'build-tools', 'runtime-export', 'python'):
        parser.add_argument('--expected-' + name + '-sha256', required=True)
    args = parser.parse_args()
    require(platform.system() == 'Darwin' and platform.machine() == 'arm64', 'macOS arm64 build required')
    for name in ('source_receipt', 'build_tools', 'runtime_export', 'python'):
        require(digest(getattr(args, name)) == getattr(args, 'expected_' + name + '_sha256'), name + ' identity differs')
    sources = json.loads(args.source_receipt.read_text())
    tools = json.loads(args.build_tools.read_text())
    capture = json.loads(args.runtime_export.read_text())
    require(sources['publisher_wheels'] == 40 and sources['locally_built_wheels_with_publisher_source_identity'] == 5,
            'unexpected runtime artifact scope')
    require(len(sources['sources']) == 45, 'expected 45 distinct package artifacts')
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False)
    for name in ('sources', 'wheels', 'build-dependencies', 'source-copies', 'logs'):
        (out / name).mkdir()
    os.nice(10)
    started = time.monotonic()
    state = {'status': 'running', 'build_mode': 'published_sources', 'pid': os.getpid(),
        'started_at': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'script_sha256': digest(Path(__file__)),
        'export_sha256': digest(args.runtime_export), 'source_receipt_sha256': digest(args.source_receipt),
        'build_tools_sha256': digest(args.build_tools), 'python_sha256': digest(args.python),
        'commands': [], 'builds': [], 'wheels': [], 'build_environment': {}}
    save(out / 'status.json', state)
    env = dict(os.environ)
    for key in ('PYTHONPATH', 'PYTHONHOME', 'VIRTUAL_ENV', 'CONDA_PREFIX', 'CFLAGS', 'CXXFLAGS',
                'LDFLAGS', 'CPPFLAGS', 'ARCHFLAGS', 'SDKROOT', 'MACOSX_DEPLOYMENT_TARGET', '_PYTHON_HOST_PLATFORM'):
        env.pop(key, None)
    env.update(PYTHONNOUSERSITE='1', PYTHONDONTWRITEBYTECODE='1', OMP_NUM_THREADS='1',
        OPENBLAS_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1', MAKEFLAGS='-j2',
        CC='/usr/bin/clang', CXX='/usr/bin/clang++',
        SOURCE_DATE_EPOCH=str(tools['source_date_epoch']), PIP_NO_INDEX='1', PIP_DISABLE_PIP_VERSION_CHECK='1')
    env['PATH'] = str(args.python.resolve().parent) + ':/usr/bin:/bin:/usr/sbin:/sbin'

    def run(command, stem):
        tick = time.monotonic()
        log = out / 'logs' / (stem + '.log')
        save(out / 'logs' / (stem + '-command.json'), list(map(str, command)))
        with log.open('x') as stream:
            result = subprocess.run(list(map(str, command)), env=env, stdout=stream, stderr=subprocess.STDOUT)
        state['commands'].append({'name': stem, 'exit_code': result.returncode,
            'elapsed_seconds': time.monotonic() - tick, 'log_sha256': digest(log)})
        save(out / 'status.json', state)
        require(result.returncode == 0, stem + ' failed; original inputs and build log retained')
        return log.read_text()

    try:
        def download_record(record):
            return fetch(record['publisher_url'], record['publisher_filename'], record['publisher_size'],
                record['publisher_sha256'], out / ('sources' if record['locally_built_wheel'] else 'wheels'))
        with ThreadPoolExecutor(max_workers=3) as pool:
            list(pool.map(download_record, sources['sources']))
        for row in tools['build_wheels']:
            p = fetch(row['url'], row['filename'], row['bytes'], row['sha256'], out / 'build-dependencies')
            observed = wheel_info(p)
            require((observed['name'], observed['version']) == (canonical(row['name']), row['version']), 'build wheel identity differs')
        requirements = out / 'build-requirements.txt'
        requirements.write_text(''.join(canonical(row['name']) + '==' + row['version'] + ' --hash=sha256:' + row['sha256'] + '\n'
                                       for row in tools['build_wheels']))
        build_env = out / 'build-environment'
        run([args.python, '-m', 'venv', '--copies', build_env], 'create-build-environment')
        python = build_env / 'bin/python'
        env['PATH'] = str(build_env / 'bin') + ':/usr/bin:/bin:/usr/sbin:/sbin'
        run([python, '-m', 'pip', '--isolated', 'install', '--no-index', '--no-deps', '--no-compile',
             '--find-links', out / 'build-dependencies', '--require-hashes', '-r', requirements], 'install-build-tools')
        run([python, '-m', 'pip', '--isolated', 'check'], 'build-dependency-check')
        state['build_environment'] = json.loads(run([python, '-c',
            'import sys,platform,importlib.metadata as m,json; print(json.dumps({"python":sys.version,"machine":platform.machine(),"packages":sorted((d.metadata["Name"],d.version) for d in m.distributions())}))'], 'build-environment-identity'))
        require(state['build_environment']['python'].startswith('3.12.'), 'CPython 3.12 build required')
        state['compiler'] = run(['/usr/bin/clang', '--version'], 'compiler-version').strip()
        state['sdk_path'] = run(['/usr/bin/xcrun', '--show-sdk-path'], 'sdk-path').strip()
        state['sdk_version'] = run(['/usr/bin/xcrun', '--show-sdk-version'], 'sdk-version').strip()
        for row in sources['sources']:
            if not row['locally_built_wheel']:
                continue
            state['stage'] = 'building ' + row['name']; save(out / 'status.json', state)
            source_dir, cleaning = clean_source(out / 'sources' / row['publisher_filename'], out / 'source-copies' / row['name'])
            save(out / 'source-copies' / (row['name'] + '-preparation.json'), cleaning)
            before = set((out / 'wheels').glob('*.whl'))
            run([python, '-m', 'pip', '--isolated', 'wheel', '--no-index', '--no-deps', '--no-build-isolation',
                 '--no-cache-dir', '--wheel-dir', out / 'wheels', source_dir], 'build-' + row['name'])
            added = set((out / 'wheels').glob('*.whl')) - before
            require(len(added) == 1, 'source build did not produce exactly one new wheel')
            wheel = next(iter(added)); identity = wheel_info(wheel)
            require((identity['name'], identity['version']) == (canonical(row['name']), row['version']), 'built version differs')
            original = args.original_wheelhouse / row['local_wheel']
            require(digest(original) == row['local_wheel_sha256'], 'original comparison wheel differs')
            old, new = payloads(original), payloads(wheel)
            with zipfile.ZipFile(wheel) as z:
                binaries = [m.filename for m in z.infolist() if not m.is_dir() and native(z.read(m)[:8])]
            state['builds'].append({'name': row['name'], 'version': row['version'],
                'source_sha256': row['publisher_sha256'], 'preparation': cleaning, 'wheel': identity,
                'native_wheel_members': binaries, 'payload_comparison': {
                    'same_bytes': sorted(n for n in old.keys() & new.keys() if old[n] == new[n]),
                    'different_bytes': sorted(n for n in old.keys() & new.keys() if old[n] != new[n]),
                    'removed': sorted(old.keys() - new.keys()), 'added': sorted(new.keys() - old.keys()),
                    'scope': 'Archive package payload excluding dist-info; changes are recorded, not assumed equivalent.'}})
            save(out / 'status.json', state)
            print('[built] ' + row['name'], flush=True)
        available = {}
        for path in sorted((out / 'wheels').glob('*.whl')):
            item = wheel_info(path); key = (item['name'], item['version'])
            require(key not in available, 'duplicate final wheel identity')
            available[key] = item; state['wheels'].append(item)
        expected = {(canonical(r['name']), r['version']) for r in sources['sources']}
        require(set(available) == expected and len(available) == 45, 'final wheel set differs')
        for name, data in capture['venvs'].items():
            if name == 'bio':
                continue
            lines = ['# Newly built source wheels plus exact publisher wheels; macOS arm64 CPython 3.12 only.']
            for row in data['distributions']:
                if row['origin'] == 'local':
                    key = canonical(row['name']), row['version']
                    lines.append(key[0] + '==' + key[1] + ' --hash=sha256:' + available[key]['sha256'])
            (out / (name + '-local.wheels.txt')).write_text('\n'.join(lines) + '\n')
        state.update(status='complete', stage='wheelhouse built; runtime validation required',
            completed_at=datetime.datetime.now(datetime.timezone.utc).isoformat(), elapsed_seconds=time.monotonic() - started,
            limits=['Published sources and build-tool artifacts are pinned; the unknown historical build environment is not reconstructed.',
                    'Prebuilt members omitted only from fresh extracted build copies; original source and original wheels are preserved.',
                    'Five newly generated wheels are not assumed byte-identical to old cached builds; runtime controls are a separate stage.',
                    'Some upstream wheel tags label native content as pure Python; this workflow is restricted to macOS arm64 CPython 3.12.',
                    'A successful source build does not establish patient accuracy, complete scientific reproduction or clinical adequacy.'])
        save(out / 'status.json', state)
        print(json.dumps({'status': state['status'], 'wheels': len(state['wheels']), 'source_builds': len(state['builds'])}))
    except BaseException as error:
        state.update(status='failed', error_type=type(error).__name__, error=str(error),
                     failed_at=datetime.datetime.now(datetime.timezone.utc).isoformat())
        save(out / 'status.json', state)
        raise


if __name__ == '__main__':
    main()
