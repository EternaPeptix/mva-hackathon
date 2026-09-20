"""Read-only completed-callset preservation audit; no calling or interpretation."""
import argparse
from collections import Counter
import datetime
import gzip
import hashlib
import itertools
import json
import os
from pathlib import Path
import platform
import shutil
import signal
import subprocess
import time

ROOT = Path('/Users/jeweled/Documents/mva-hackathon')
RUN = ROOT / 'calls/from-fastq-cram30-20260920'
DEFINITION = '6bb6aca4c6e281ccbe345b6f7a4e7a38f340da8a9138b3cd8faec87ccc532267'
TERMINAL = 'd58b058470ab989ac27257699705401a58a6a81c587b09c4f3018ee64a6788e0'


def require(value, message):
    if not value:
        raise ValueError(message)


def stamp():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def save(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write('\n')


def records(path, sample, dictionary):
    """Decompressed original VCF payload, with sample/dictionary header checks."""
    header_contigs, sample_seen = [], False
    with gzip.open(path, 'rt', encoding='utf-8', newline='') as source:
        for line in source:
            require(line.endswith('\n'), 'unterminated VCF line')
            if line.startswith('##contig=<'):
                fields = dict(item.split('=', 1) for item in line.strip()[10:-1].split(',') if '=' in item)
                header_contigs.append((fields['ID'], int(fields['length'])))
            elif line.startswith('#CHROM\t'):
                require(line.rstrip('\n').split('\t')[9:] == [sample], 'sample header differs')
                require(header_contigs == dictionary, 'contig dictionary differs')
                require(not sample_seen, 'duplicate sample header')
                sample_seen = True
            elif not line.startswith('#'):
                require(sample_seen, 'record before sample header')
                yield line.encode('utf-8')
        require(sample_seen, 'missing sample header')


def compare_streams(parts, final, allowed_contigs):
    """Exact concatenation equality including duplicate multiplicity and all fields."""
    count, payload, per_contig = 0, hashlib.sha256(), Counter()
    positions, last_order = {}, (-1, -1)
    order = {name: i for i, (name, _) in enumerate(allowed_contigs)}
    lengths = dict(allowed_contigs)
    for expected, observed in itertools.zip_longest(itertools.chain.from_iterable(parts), final):
        require(expected is not None and observed is not None, 'record count differs')
        require(expected == observed, 'gather changed decoded record payload at ordinal ' + str(count + 1))
        fields = observed.split(b'\t')
        require(len(fields) == 10, 'record column count differs')
        chrom, pos = fields[0].decode(), int(fields[1])
        require(chrom in lengths and 1 <= pos <= lengths[chrom], 'record outside reference')
        coordinate = order[chrom], pos
        require(coordinate >= last_order, 'records not in dictionary/position order')
        last_order = coordinate
        payload.update(observed)
        count += 1
        per_contig[chrom] += 1
        positions.setdefault(chrom, [pos, pos])[1] = pos
    return dict(records=count, payload_sha256=payload.hexdigest(), records_by_contig=dict(per_contig),
                observed_position_bounds=positions)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require(platform.node().split('.')[0] == '512S2', 'Mac Studio512S2 only')
    out = args.output.resolve()
    require(out.is_relative_to(ROOT / 'results/whole-reference-completion-20260920'), 'isolated output scope required')
    require(shutil.disk_usage(ROOT).free >= 40 * 1024**3, 'less than40GiB free')
    out.mkdir(parents=True, exist_ok=False)
    os.nice(10)
    started, pins = time.monotonic(), {}
    def alarm(*unused):
        raise TimeoutError('30minute audit limit')
    signal.signal(signal.SIGALRM, alarm)
    signal.alarm(1800)
    def pin(path, expected=None):
        path = Path(path)
        require(path.is_file() and not path.is_symlink(), 'missing or symlink input')
        st = path.stat()
        digest = hashlib.sha256()
        with path.open('rb') as stream:
            for block in iter(lambda: stream.read(8*1024*1024), b''):
                digest.update(block)
        after = path.stat()
        require((st.st_size, st.st_mtime_ns) == (after.st_size, after.st_mtime_ns), 'input changed while hashing')
        row = dict(bytes=st.st_size, mtime_ns=st.st_mtime_ns, sha256=digest.hexdigest())
        if expected:
            require(row['sha256'] == expected, 'input SHA mismatch: ' + str(path))
        pins[str(path)] = row
        return row
    try:
        import pysam
        pin(Path(__file__))
        pin(RUN / 'status.json', TERMINAL)
        state = json.loads((RUN / 'status.json').read_text())
        require(state['status'] == 'complete' and state['pid'] == 79898 and state['completed_shards'] == list(range(64)), 'not original completed64-section run')
        require(state['expected_shards'] == 64 and state['definition_sha256'] == DEFINITION, 'definition linkage differs')
        pin(RUN / 'run-definition.json', DEFINITION)
        definition = json.loads((RUN / 'run-definition.json').read_text())
        dictionary = [tuple(row) for group in definition['shards'] for row in group]
        require(len(definition['shards']) == 64 and len({c for c, _ in dictionary}) == len(dictionary) == definition['contigs'], 'contig membership differs')
        require(sum(n for _, n in dictionary) == definition['reference_bases'], 'reference size differs')
        fai = ROOT / 'ref/realignment/GRCh38_nochr.fa.fai'
        pin(fai, 'd81f672f40cfabf2510f37e3f790f2c79665023f4ac950bb9cc6496a40cebf32')
        require([(line.split('\t')[0], int(line.split('\t')[1])) for line in fai.read_text().splitlines()] == dictionary, 'reference index dictionary differs')
        completed = []
        for i, group in enumerate(definition['shards']):
            folder = RUN / ('shard-%04d' % i)
            pin(folder / 'complete.json')
            completed.append(json.loads((folder / 'complete.json').read_text()))
        results = {}
        for method in ('gatk', 'bcftools'):
            pin(RUN / (method + '-complete.json'))
            published = json.loads((RUN / (method + '-complete.json')).read_text())
            require(published == state['outputs'][method], 'terminal output receipt differs')
            final = RUN / ('WGS_EX2312012.' + method + '.normalized.raw.vcf.gz')
            require(published['path'] == str(final), 'unexpected published path')
            pin(final, published['sha256'])
            pin(Path(str(final) + '.tbi'))
            paths, shard_counts = [], []
            for i, row in enumerate(completed):
                folder = RUN / ('shard-%04d' % i)
                p = folder / row[method]['path']
                require(p.resolve().is_relative_to(folder.resolve()) and not Path(row[method]['path']).is_absolute(), 'shard path escape')
                pin(p, row[method]['sha256'])
                paths.append(p)
                shard_counts.append(row[method]['counts']['records'])
            list_path = RUN / (method + '-inputs.list')
            pin(list_path)
            require(list_path.read_text().splitlines() == list(map(str, paths)), 'gather membership/order differs')
            def part(index, path):
                names = {c for c, _ in definition['shards'][index]}
                n = 0
                for record in records(path, definition['sample'], dictionary):
                    require(record.split(b'\t', 1)[0].decode() in names, 'record belongs to another section')
                    n += 1
                    yield record
                require(n == shard_counts[index], 'section receipt count differs')
            result = compare_streams([part(i, p) for i, p in enumerate(paths)], records(final, definition['sample'], dictionary), dictionary)
            require(result['records'] == sum(shard_counts) == published['counts']['records'], 'published count differs')
            indexed_payload, indexed_count = hashlib.sha256(), 0
            with pysam.TabixFile(str(final)) as index:
                for chrom, length in dictionary:
                    if chrom not in index.contigs:
                        require(result['records_by_contig'].get(chrom, 0) == 0, 'nonempty contig absent from index')
                        continue
                    n = 0
                    for line in index.fetch(chrom, 0, length):
                        indexed_payload.update((line + '\n').encode())
                        indexed_count += 1
                        n += 1
                    require(n == result['records_by_contig'].get(chrom, 0), 'indexed contig count differs')
            require(indexed_count == result['records'] and indexed_payload.hexdigest() == result['payload_sha256'], 'indexed payload differs from sequential payload')
            result.update(published_sha256=published['sha256'], index_payload_verified=True, exact_shard_concatenation_verified=True)
            results[method] = result
            print(json.dumps(dict(method=method, records=result['records'], status='verified')), flush=True)
        for path, row in pins.items():
            st = Path(path).stat()
            require((st.st_size, st.st_mtime_ns) == (row['bytes'], row['mtime_ns']), 'input changed during audit')
        require(shutil.disk_usage(ROOT).free >= 40 * 1024**3, 'free space limit')
        receipt = dict(status='complete', completed_at=stamp(), elapsed_seconds=time.monotonic()-started,
            host='512S2', original_caller_pid=79898, caller_completed_at=state['completed_at'], sections=64,
            reference_contigs=definition['contigs'], reference_bases=definition['reference_bases'], methods=results,
            freshly_hashed_files=pins, raw_calls_recomputed=False, alignment_opened=False,
            scope='Exact64-section normalized VCF concatenation and indexed/sequential payload preservation; no genotype/quality recoding, clinical validation, phase or causal inference. Reference/alignment content identity is inherited from the completed run and prior audits, not rehashed here.',
            resources=dict(wall_limit_seconds=1800, minimum_free_bytes=40*1024**3, policy='Streaming bounded state; nice+10; alarm terminates only this verifier; no child computation.'))
        save(out / 'verification.json', receipt)
        save(out / 'status.json', dict(status='complete', completed_at=stamp(), verification_sha256=hashlib.sha256((out / 'verification.json').read_bytes()).hexdigest()))
    except BaseException as error:
        save(out / 'status.json', dict(status='failed', failed_at=stamp(), error_type=type(error).__name__, error=str(error), retained_without_retry=True))
        raise
    finally:
        signal.alarm(0)


if __name__ == '__main__':
    main()
