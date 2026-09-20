#!/usr/bin/env python3
"""Extract bounded candidate windows from retained lane BAMs for compute-host review.

Only indexed regional extraction runs where the source BAMs are stored. Whole-BAM
hashes are not claimed; source metadata/header/index and extracted-byte hashes are
recorded, and sources must remain unchanged during extraction.
"""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def metadata(path):
    stat = path.stat()
    return {'path': str(path.resolve()), 'bytes': stat.st_size, 'mtime_ns': stat.st_mtime_ns}


def intervals_for_targets(targets, padding):
    intervals = []
    for target in sorted(targets, key=lambda x: (int(x['chrom']), x['pos'])):
        chrom, start, end = target['chrom'], max(0, target['pos'] - 1 - padding), target['pos'] - 1 + len(target['ref']) + padding
        if intervals and intervals[-1][0] == chrom and start <= intervals[-1][2]:
            intervals[-1][2] = max(intervals[-1][2], end)
        else:
            intervals.append([chrom, start, end])
    return intervals


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bam-directory', type=Path, required=True)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--samtools', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    config = json.loads(args.config.read_text())
    intervals = intervals_for_targets(config['targets'], config['slice_padding_bp'])
    bed = args.output / 'regions.bed'
    bed.write_text(''.join('\t'.join(map(str, row)) + '\n' for row in intervals))
    commands, files = [], []
    for lane in ('L001', 'L002'):
        source = args.bam_directory / (config['sample'] + '.' + lane + '.bam')
        index = Path(str(source) + '.bai')
        before = metadata(source)
        index_before = metadata(index)
        header = subprocess.run([str(args.samtools), 'view', '-H', str(source)], capture_output=True, check=True).stdout
        destination = args.output / (lane + '.candidate-windows.bam')
        command = [str(args.samtools), 'view', '-@', '1', '-b', '-M', '-L', str(bed), '-o', str(destination), str(source)]
        commands.append(command)
        with (args.output / (lane + '.extraction.log')).open('wb') as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
            subprocess.run([str(args.samtools), 'quickcheck', '-v', str(destination)], stdout=log, stderr=subprocess.STDOUT, check=True)
        if before != metadata(source) or index_before != metadata(index):
            raise ValueError('source changed during extraction')
        files.append({'lane': lane, 'source': before, 'source_index': dict(index_before, sha256=sha(index)),
            'source_header_sha256': hashlib.sha256(header).hexdigest(), 'source_full_BAM_sha256': None,
            'slice_name': destination.name, 'slice_bytes': destination.stat().st_size, 'slice_sha256': sha(destination)})
    result = {'scope': config['scope'], 'status': 'complete',
        'created_at': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'config_sha256': sha(args.config), 'script_sha256': sha(Path(__file__)),
        'samtools': dict(metadata(args.samtools), sha256=sha(args.samtools)),
        'regions_BED_0based_half_open': intervals, 'regions_bed_sha256': sha(bed),
        'files': files, 'commands': commands, 'invocation': [sys.executable] + sys.argv,
        'elapsed_seconds': time.monotonic() - started,
        'limitation': 'Historical source BAMs are identified by metadata/header/index; their complete bytes were not rehashed. Extracted bytes are SHA-256 pinned.'}
    (args.output / 'slice-provenance.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'status': 'complete', 'regions': len(intervals), 'slice_bytes': sum(row['slice_bytes'] for row in files),
                      'elapsed_seconds': result['elapsed_seconds']}))


if __name__ == '__main__':
    main()
