#!/usr/bin/env python3
"""Verify all challenge FASTQs and produce a four-lane CRAM on the Studio."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time


def stamp():
    return time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as source:
        for block in iter(lambda: source.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def save(path, value):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    os.replace(temporary, path)


def prepare_reference(root, env):
    source = root / 'ref/vendor/GRCh38_nochr.fa'
    folder = root / 'ref/realignment'
    folder.mkdir(exist_ok=True)
    with (folder / '.index.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        reference_hash = digest(source)
        reference = folder / source.name
        if reference.exists():
            if not os.path.samefile(source, reference) and digest(reference) != reference_hash:
                raise RuntimeError('existing realignment reference differs from source')
        else:
            os.link(source, reference)
        source_index = Path(str(source) + '.fai')
        if not Path(str(reference) + '.fai').exists():
            shutil.copy2(source_index, str(reference) + '.fai')
        receipt = folder / 'bwa-index.json'
        suffixes = ['.amb', '.ann', '.bwt', '.pac', '.sa']
        if receipt.is_file():
            prior = json.loads(receipt.read_text())
            if prior['reference_sha256'] != reference_hash:
                raise RuntimeError('reference changed since index construction')
            for suffix in suffixes:
                path = Path(str(reference) + suffix)
                if not path.is_file() or digest(path) != prior['index_sha256'][suffix]:
                    raise RuntimeError('existing index is incomplete or changed')
            return reference
        if any(Path(str(reference) + suffix).exists() for suffix in suffixes):
            raise RuntimeError('existing unverified index files preserved for inspection')
        stage = folder / ('index-building-' + str(os.getpid()))
        stage.mkdir()
        prefix = stage / source.name
        with (root / 'logs/bwa-reference-index.log').open('w') as log:
            subprocess.run(['bwa', 'index', '-p', str(prefix), str(reference)],
                           env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
        hashes = {}
        for suffix in suffixes:
            path = Path(str(prefix) + suffix)
            if not path.is_file() or path.stat().st_size == 0:
                raise RuntimeError('bwa did not produce all expected index files')
            hashes[suffix] = digest(path)
        for suffix in suffixes:
            os.replace(str(prefix) + suffix, str(reference) + suffix)
        stage.rmdir()
        save(receipt, {'completed_at': stamp(), 'reference_sha256': reference_hash,
                       'index_sha256': hashes, 'construction': 'bwa index; native Apple Silicon'})
        return reference


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--threads', type=int, default=8)
    args = parser.parse_args()
    if args.threads < 1 or args.threads > 8:
        parser.error('use between one and eight bwa worker threads')
    root = args.root.resolve()
    scripts = Path(__file__).resolve().parent
    python = root / 'tools/bio/bin/python'
    env = dict(os.environ, PATH=str(python.parent) + os.pathsep + os.environ.get('PATH', ''),
               MVA_PYTHON=str(python), SORT_MEM='1G')
    state_path = root / 'state/raw-reads.json'
    with (root / 'state/raw-reads.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        os.nice(10)
        if state_path.exists():
            save(root / ('state/raw-reads-history-%d.json' % time.time_ns()), json.loads(state_path.read_text()))
        state = {'status': 'running', 'pid': os.getpid(), 'started_at': stamp(),
                 'root': str(root), 'bwa_worker_threads': args.threads, 'expected_lanes': 4,
                 'scope': 'Verify all FASTQs, align all four lanes, mark duplicates, and validate final CRAM. Calling is separate.',
                 'controller_sha256': digest(Path(__file__)), 'stages': []}

        def stage(name):
            state['current_stage'] = name
            state['stages'].append({'stage': name, 'started_at': stamp()})
            save(state_path, state)
            print('[stage] ' + name, flush=True)

        try:
            if shutil.disk_usage(root).free < 300_000_000_000:
                raise RuntimeError('less than 300 GB free; preserve headroom before beginning full alignment')
            stage('verify FASTQs and build reference index')
            # One reference-index process may run alongside network I/O. Alignment waits for both.
            with ThreadPoolExecutor(max_workers=1) as executor:
                index = executor.submit(prepare_reference, root, env)
                with (root / 'logs/fastq-verification.log').open('a') as log:
                    subprocess.run(['bash', str(scripts / 'fetch_fastqs.sh'), str(root / 'data')],
                                   env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
                stage('wait for verified reference index')
                reference = index.result()
            stage('align all four lanes and validate final CRAM')
            output = root / 'align/from-fastq-20260919'
            with (root / 'logs/four-lane-alignment.log').open('a') as log:
                subprocess.run(['bash', str(scripts / '01_align.sh'), '-r', str(reference),
                    '-o', str(output), '-s', 'WGS_EX2312012', '-f', str(root / 'data'),
                    '-t', str(args.threads), '-n', '4'], env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
            cram = output / 'WGS_EX2312012.cram'
            if not cram.is_file() or not (output / 'logs/finalization.complete').is_file():
                raise RuntimeError('alignment returned without a validated final CRAM')
            subprocess.run(['samtools', 'quickcheck', '-v', str(cram)], env=env, check=True)
            stage('record final CRAM checksum')
            state.update(status='complete', completed_at=stamp(), current_stage='complete',
                         cram_path=str(cram), cram_bytes=cram.stat().st_size, cram_sha256=digest(cram),
                         next_stage='Independent variant calling and read-level evidence review remain to be performed.')
            save(state_path, state)
            print('[done] all four lanes represented in a verified CRAM', flush=True)
        except BaseException as error:
            state.update(status='failed', failed_at=stamp(), error_type=type(error).__name__)
            save(state_path, state)
            print('[FAIL] %s; files retained' % type(error).__name__, flush=True)
            raise


if __name__ == '__main__':
    main()
