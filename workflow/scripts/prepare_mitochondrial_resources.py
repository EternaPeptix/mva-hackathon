#!/usr/bin/env python3
"""Acquire pinned public mitochondrial references; verify rotation and build native indexes."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import urllib.request

from prepare_gatk_compatible_cram import digest, save, stamp
from run_mitochondrial_discovery import sequence, verify_rotation


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('root', 'sources', 'output'):
        p.add_argument('--' + name, type=Path, required=True)
    args = p.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    sources = json.loads(args.sources.read_text())
    for item in sources['sources']:
        if Path(item['name']).name != item['name']:
            raise ValueError('unsafe public resource name')
        with urllib.request.urlopen(item['url'], timeout=40) as response:
            data = response.read(200001)
        if len(data) != item['bytes'] or hashlib.sha256(data).hexdigest() != item['sha256']:
            raise ValueError('public resource identity differs: ' + item['name'])
        (out / item['name']).write_bytes(data)
    sam = args.root / 'tools/bio/bin/samtools'
    bwa = args.root / 'tools/bio/bin/bwa'
    commands = []
    with (out / 'preparation.log').open('w') as log:
        for name in ['Homo_sapiens_assembly38.chrM.fasta', 'Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta', 'NC_012920.1.fasta']:
            ref = out / name
            for command in [[sam, 'faidx', ref], [sam, 'dict', '-o', ref.with_suffix('.dict'), ref]]:
                commands.append(list(map(str, command)))
                subprocess.run(commands[-1], stdout=log, stderr=subprocess.STDOUT, check=True)
            if name != 'NC_012920.1.fasta':
                commands.append(list(map(str, [bwa, 'index', ref])))
                subprocess.run(commands[-1], stdout=log, stderr=subprocess.STDOUT, check=True)
        java = args.root / 'tools/bio/lib/jvm/bin/java'
        jar = args.root / 'gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar'
        commands.append(list(map(str, [java, '-Xmx1g', '-jar', jar, 'IndexFeatureFile', '-I', out / 'blacklist_sites.hg38.chrM.bed'])))
        subprocess.run(commands[-1], stdout=log, stderr=subprocess.STDOUT, check=True)
    regular = sequence(out / 'Homo_sapiens_assembly38.chrM.fasta')[1]
    shifted = sequence(out / 'Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta')[1]
    ncbi = sequence(out / 'NC_012920.1.fasta')[1]
    verify_rotation(regular, shifted)
    if regular != ncbi:
        raise ValueError('Broad mitochondrial reference differs from NCBI accession')
    save(out / 'commands.json', commands)
    manifest = {'completed_at': stamp(), 'sources_config_sha256': digest(args.sources),
                'source_sha256': digest(Path(__file__)), 'regular_sequence_sha256': hashlib.sha256(regular.encode()).hexdigest(),
                'NCBI_accession': 'NC_012920.1', 'NCBI_sequence_identical': True, 'shift_offset': 8000,
                'files': {p.name: {'bytes': p.stat().st_size, 'sha256': digest(p)} for p in out.iterdir() if p.is_file()}}
    save(out / 'manifest.json', manifest)
    print(json.dumps({'status': 'complete', 'files': len(manifest['files']), 'regular_sequence_sha256': manifest['regular_sequence_sha256']}))


if __name__ == '__main__':
    main()
