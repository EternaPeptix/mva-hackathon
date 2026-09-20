#!/usr/bin/env python3
"""One bounded Studio run: verify inputs, prepare databases, run two Exomiser analyses."""
import argparse
import fcntl
import json
import os
from pathlib import Path
import shutil
import stat
import subprocess
import sys
import time
import zipfile

from fetch_resource import fetch, digest


VCF = 'WGS_EX2312012_HGWCNDSX7.vcf.gz'
VCF_SHA256 = 'bec64c68d3ad10d05b8071d07cdf87d217013d66d3ba02e24b505035c9d42175'
DATABASES = [('2512_phenotype', 13132847504), ('2512_hg38', 23557715707)]


def stamp():
    return time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())


def save(path, data):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(data, indent=2) + '\n')
    os.replace(temporary, path)


def extract_database(archive, directory, name):
    target = directory / name
    receipt = directory / (name + '.extracted.json')
    if target.exists():
        if not receipt.is_file():
            raise RuntimeError('database directory has no extraction receipt; preserving it for inspection')
        record = json.loads(receipt.read_text())
        if record['archive_size'] != archive.stat().st_size:
            raise RuntimeError('archive changed since extraction')
        if any(not (directory / item['path']).is_file() or (directory / item['path']).stat().st_size != item['size']
               for item in record['files']):
            raise RuntimeError('extracted database is incomplete')
        return
    staging = directory / ('.extract-' + name + '-' + str(os.getpid()))
    staging.mkdir()
    with zipfile.ZipFile(archive) as package:
        entries = package.infolist()
        expanded = sum(item.file_size for item in entries)
        if shutil.disk_usage(directory).free < expanded + 50_000_000_000:
            raise RuntimeError('insufficient free space with 50 GB reserve')
        for item in entries:
            relative = Path(item.filename)
            if relative.is_absolute() or '..' in relative.parts or not relative.parts or relative.parts[0] != name:
                raise RuntimeError('unexpected archive member path')
            if stat.S_ISLNK(item.external_attr >> 16):
                raise RuntimeError('archive symlinks are not permitted')
        # ZipFile verifies each member CRC while extracting. Publish only the completed tree.
        package.extractall(staging)
        files = [{'path': item.filename, 'size': item.file_size, 'crc32': '%08x' % item.CRC}
                 for item in entries if not item.is_dir()]
    os.replace(staging / name, target)
    staging.rmdir()
    save(receipt, {'archive': archive.name, 'archive_size': archive.stat().st_size,
                   'extracted_at': stamp(), 'expanded_bytes': expanded, 'files': files})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--analysis-only', action='store_true', help='Reuse databases with completed verification/extraction receipts')
    parser.add_argument('--run-id', default=time.strftime('studio-%Y%m%dT%H%M%SZ', time.gmtime()))
    args = parser.parse_args()
    if not args.run_id or any(c not in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_' for c in args.run_id):
        parser.error('run-id must contain only letters, numbers, hyphens, and underscores')
    root = args.root.resolve()
    state_dir = root / 'state'
    state_dir.mkdir(exist_ok=True)
    with (state_dir / 'exomiser-resume.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        os.nice(10)
        state_path = state_dir / 'exomiser-resume.json'
        if state_path.exists():
            (state_dir / ('exomiser-resume-history-%d.json' % time.time_ns())).write_bytes(state_path.read_bytes())
        state = {'status': 'running', 'started_at': stamp(), 'pid': os.getpid(), 'run_id': args.run_id,
                 'root': str(root), 'threads': args.threads, 'analysis_only': args.analysis_only, 'stages': []}

        def stage(name):
            state['current_stage'] = name
            state['updated_at'] = stamp()
            state['stages'].append({'stage': name, 'started_at': stamp()})
            save(state_path, state)
            print('[stage] ' + name, flush=True)

        try:
            stage('verify supplied VCF against publisher SHA-256')
            digest(root / 'data' / VCF, VCF_SHA256)
            save(state_dir / 'vendor-vcf-verification.json', {'filename': VCF, 'publisher_sha256': VCF_SHA256,
                 'verified_at': stamp(), 'source': 'Hugging Face dataset LFS metadata retrieved 2026-09-18'})
            java = root / 'tools/bio/lib/jvm/bin/java'
            python = root / 'tools/bio/bin/python'
            if not java.is_file() or not python.is_file():
                raise RuntimeError('isolated Java/Python runtime is not ready')
            for name, size in DATABASES:
                archive = root / 'exomiser' / (name + '.zip')
                if args.analysis_only:
                    stage('reuse verified ' + name)
                    receipt = json.loads(Path(str(archive) + '.verified.json').read_text())
                    if archive.stat().st_size != size or receipt.get('size') != size or not receipt.get('archive_integrity_checked'):
                        raise RuntimeError('verified archive receipt is missing or inconsistent')
                    if not (root / 'exomiser' / name).is_dir():
                        raise RuntimeError('verified database has not been extracted')
                else:
                    stage('download and verify ' + name)
                    archive = fetch('https://data.monarchinitiative.org/exomiser/data/' + name + '.zip',
                                    root / 'exomiser', expected_size=size)
                stage('extract and verify ' + name)
                extract_database(archive, root / 'exomiser', name)
            script_dir = Path(__file__).resolve().parent
            for label, extra in [('subject-phenotypes', []), ('parental-history-sensitivity', ['--add-hpo', 'HP:0200067'])]:
                stage('Exomiser ' + label)
                output = root / 'results' / args.run_id / label
                command = [str(python), str(script_dir / 'prioritize.py'), '-v', str(root / 'data' / VCF),
                    '-o', str(output), '-e', str(root / 'exomiser/exomiser-cli-15.1.0'), '-d', str(root / 'exomiser'),
                    '-y', str(script_dir.parent / 'configs/proband01-exomiser.yml'),
                    '-t', str(args.threads), '-m', '16g', '--java', str(java)] + extra
                subprocess.run(command, check=True)
            state.update(status='complete', completed_at=stamp(), current_stage='complete')
            save(state_path, state)
            print('[done] both genome-wide ranking runs completed', flush=True)
        except BaseException as error:
            state.update(status='failed', failed_at=stamp(), error_type=type(error).__name__)
            save(state_path, state)
            print('[FAIL] %s; inputs and incomplete work preserved' % type(error).__name__, flush=True)
            raise


if __name__ == '__main__':
    main()
