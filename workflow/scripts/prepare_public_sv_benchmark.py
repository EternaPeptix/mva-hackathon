#!/usr/bin/env python3
"""Fetch checksum-pinned public HG002 chr20 pilot inputs outside the repository.

This prepares real public reads, not synthetic reads or patient data. The public
prealigned chromosome subset tests a caller component, not whole-genome alignment
or clinical sensitivity. Keep the complete truth VCF for complex-allele context.
"""
import argparse
import concurrent.futures
import hashlib
import json
from pathlib import Path
import shutil
import time

from fetch_resource import fetch

GIAB = 'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/'
TRUTH = GIAB + 'release/AshkenazimTrio/HG002_NA24385_son/v5.0q/'
REF = GIAB + 'release/references/GRCh38/'
READS = 'https://storage.googleapis.com/deepvariant/case-study-testdata/'
REFERENCE = 'GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz'
ALIGNMENT = 'HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam'
FILES = [
    (REF + REFERENCE, '3a3347eae0893f96ecf495d1c39e2284', 886344618),
    (REF + REFERENCE + '.fai', '5fddbc109c82980f9436aa5c21a57c61', None),
    (REF + REFERENCE + '.gzi', 'cfa1ee11b1ecb29f936578b27014fbe0', None),
    (TRUTH + 'HG002_GRCh38_v5.0q_stvar.vcf.gz', '34dbeeac5eb3fc09e94c5cb8cfe9e24f', 48040012),
    (TRUTH + 'HG002_GRCh38_v5.0q_stvar.vcf.gz.tbi', 'df639c204bf457f22248220602ea183e', None),
    (TRUTH + 'HG002_GRCh38_v5.0q_stvar.benchmark.bed', '55fc11c80ff376c0a8da99d0c69004f3', None),
    (TRUTH + 'HG002_GRCh38_v5.0q_stvar.exclusion_stats.txt', '92bdc427ad4fb4c1e7954d715e276e43', None),
    (READS + ALIGNMENT, '8f1a541a348803e304275359ac52da57', 1093938924),
    (READS + ALIGNMENT + '.bai', '4d4892f4e13adab90745877b3e01c5b6', 196136),
]


def checksums(path):
    before = path.stat(); md5 = hashlib.md5(); sha = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            md5.update(block); sha.update(block)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError('input changed while hashing')
    return {'path': str(path.resolve()), 'bytes': after.st_size,
            'md5': md5.hexdigest(), 'sha256': sha.hexdigest()}


def prepare(output):
    if output.resolve().is_relative_to(Path(__file__).resolve().parents[2]):
        raise ValueError('public benchmark assets must remain outside the repository')
    output.mkdir(parents=True, exist_ok=True)
    if shutil.disk_usage(output).free < 350 * 2**30:
        raise ValueError('preserve at least 350 GiB free before starting these downloads')
    start = time.monotonic()
    def one(item):
        url, expected, size = item
        path = fetch(url, output, expected_size=size)
        identity = checksums(path)
        if identity['md5'] != expected: raise ValueError('public publisher MD5 mismatch: ' + path.name)
        return {**identity, 'url': url, 'publisher_md5': expected, 'publisher_md5_matches': True}
    results = []
    # Two independent HTTP transfers; do not change the active bio environment.
    with concurrent.futures.ThreadPoolExecutor(2) as pool:
        futures = [pool.submit(one, item) for item in FILES]
        for future in futures: results.append(future.result())
    receipt = {'scope': 'public_HG002_chr20_prealigned_Illumina_35x_pilot',
               'status': 'inputs_verified', 'elapsed_seconds': time.monotonic() - start,
               'completed_at': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
               'files': results, 'preparation_code': checksums(Path(__file__)),
               'checksum_sources': [GIAB + 'release/references/checksum.md5', TRUTH + 'checksum.md5',
                                    'Google Cloud object x-goog-hash md5 headers, retrieved 2026-09-19'],
               'limitations': ['This prealigned chr20 subset cannot validate raw-read alignment or genome-wide performance.',
                               'Interchromosomal rearrangements cannot be assessed from a chr20-only BAM.',
                               'The v5.0q draft benchmark requires paired BED/VCF use and complex-allele-aware comparisons.',
                               'Reference sequence and alignment metadata still require examination before calling.']}
    (output / 'input-receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    return receipt


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args(); result = prepare(args.output)
    print(json.dumps({'status': result['status'], 'files': len(result['files']),
                      'elapsed_seconds': result['elapsed_seconds']}))
