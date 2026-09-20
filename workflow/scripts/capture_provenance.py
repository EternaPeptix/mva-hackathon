#!/usr/bin/env python3
"""Record content hashes and tool versions on the machine performing the analysis."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import time


def fingerprint(path, root):
    digest = hashlib.sha256()
    before = path.stat()
    with path.open('rb') as source:
        for block in iter(lambda: source.read(8 * 1024 * 1024), b''):
            digest.update(block)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise RuntimeError('input changed while computing its fingerprint')
    return {'path': str(path.relative_to(root)), 'bytes': after.st_size, 'sha256': digest.hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    if args.output.exists():
        parser.error('existing provenance preserved; choose a new output path')
    os.nice(10)
    files = [root / name for name in (
        'data/WGS_EX2312012_HGWCNDSX7.vcf.gz', 'data/WGS_EX2312012_HGWCNDSX7.vcf.gz.tbi',
        'ref/vendor/GRCh38_nochr.fa', 'ref/vendor/GRCh38_nochr.fa.fai',
        'exomiser/2512_phenotype.zip', 'exomiser/2512_hg38.zip',
        'exomiser/exomiser-cli-15.1.0/exomiser-cli-15.1.0.jar',
        'state/runtime-explicit-lock.txt', 'repo/workflow/configs/proband01-exomiser.yml')]
    files += sorted((root / 'exomiser/exomiser-cli-15.1.0/lib').glob('*.jar'))
    files += sorted((root / 'repo/workflow/scripts').glob('*.py'))
    files += sorted((root / 'repo/workflow/scripts').glob('*.sh'))
    if any(not path.is_file() for path in files):
        parser.error('a required input is not complete yet')
    records = [fingerprint(path, root) for path in files]
    java = root / 'tools/bio/lib/jvm/bin/java'
    commands = {'java': [str(java), '-version'], 'python': [str(root/'tools/bio/bin/python'), '--version'],
                'samtools': [str(root/'tools/bio/bin/samtools'), '--version'],
                'bcftools': [str(root/'tools/bio/bin/bcftools'), '--version'],
                'bwa': [str(root/'tools/bio/bin/bwa')]}
    versions = {}
    for name, command in commands.items():
        result = subprocess.run(command, capture_output=True, text=True, timeout=15)
        lines = (result.stdout + result.stderr).strip().splitlines()
        versions[name] = [line for line in lines if 'Version:' in line] if name == 'bwa' else lines[:2]
    record = {'captured_at': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
              'platform': {'architecture': platform.machine(), 'macos': platform.mac_ver()[0]},
              'tools': versions, 'files': records,
              'hash_note': 'These are hashes computed on the Studio. Publisher verification is recorded separately in download receipts.'}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + '.tmp')
    temporary.write_text(json.dumps(record, indent=2) + '\n')
    os.replace(temporary, args.output)
    print('Provenance recorded for %d files' % len(records))


if __name__ == '__main__':
    main()
