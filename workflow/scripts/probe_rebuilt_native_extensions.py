#!/usr/bin/env python3
"""Exercise rebuilt native alignment and indexed-query dependencies on artificial data."""
import argparse
import datetime
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import subprocess
import sys

import pyabpoa
import pysam
import pywfa.align
from pywfa import WavefrontAligner
import tabix


def require(condition, message):
    if not condition:
        raise ValueError(message)


def file_info(path):
    path = Path(path).resolve()
    return {'path': str(path), 'bytes': path.stat().st_size,
            'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
            'file_description': subprocess.check_output(['/usr/bin/file', str(path)], text=True).strip()}


def validate_alignment(pattern, text, cigar):
    p, t = 0, 0
    for op, n in cigar:
        require(n > 0, 'nonpositive CIGAR length')
        if op in (0, 7):
            require(pattern[p:p+n] == text[t:t+n], 'match operation changes bases')
            p += n; t += n
        elif op == 8:
            require(all(a != b for a, b in zip(pattern[p:p+n], text[t:t+n])), 'mismatch operation contains match')
            p += n; t += n
        elif op == 1:
            t += n
        elif op == 2:
            p += n
        else:
            raise ValueError('unreviewed end-to-end CIGAR operation')
        require(p <= len(pattern) and t <= len(text), 'alignment overruns input')
    require(p == len(pattern) and t == len(text), 'alignment omits input bases')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require(platform.system() == 'Darwin' and platform.machine() == 'arm64' and sys.version_info[:2] == (3, 12),
            'this probe covers macOS arm64 CPython 3.12')
    args.output.mkdir(parents=True, exist_ok=False)
    modules = {'pyabpoa': pyabpoa, 'pywfa': pywfa.align, 'pytabix': tabix}
    native_files = {}
    for name, module in modules.items():
        path = Path(module.__file__).resolve()
        require(path.is_relative_to(Path(sys.prefix).resolve()), 'extension imported outside selected environment')
        native_files[name] = file_info(path)
        require('arm64' in native_files[name]['file_description'], 'native extension is not reported as arm64')
    cases = []
    pattern = 'ACGTCAGTACCTGATCGTAC'
    sequences = [pattern, pattern[:10] + 'G' + pattern[10:], pattern]
    msa = pyabpoa.msa_aligner().msa(sequences, out_cons=True, out_msa=True)
    require(msa.n_seq == 3 and msa.n_cons == 1 and msa.cons_seq == [pattern], 'known majority consensus differs')
    require([row.replace('-', '') for row in msa.msa_seq[:3]] == sequences, 'MSA changes input sequences')
    require(len(set(map(len, msa.msa_seq))) == 1, 'MSA rows differ in width')
    cases.append({'name': 'pyabpoa majority consensus with one insertion', 'inputs': sequences,
                  'consensus': msa.cons_seq, 'msa': msa.msa_seq, 'passed': True})
    examples = [('exact', pattern, 0), ('substitution', pattern[:10] + 'A' + pattern[11:], -4),
                ('insertion', pattern[:10] + 'GGG' + pattern[10:], -12),
                ('deletion', pattern[:10] + pattern[13:], -12)]
    for name, text, expected_score in examples:
        aligner = WavefrontAligner(pattern, span='end-to-end', heuristic=None)
        result = aligner(text, clip_cigar=False)
        require(result.status == 0 and result.score == expected_score, 'known affine alignment score differs: ' + name)
        validate_alignment(pattern, text, result.cigartuples)
        cases.append({'name': 'pywfa ' + name, 'pattern': pattern, 'text': text,
                      'score': result.score, 'cigar': result.cigartuples, 'passed': True})
    bed = args.output / 'artificial.bed'
    rows = [['1', '10', '20', 'a'], ['1', '20', '30', 'b'], ['2', '0', '5', 'c']]
    bed.write_text(''.join('\t'.join(r) + '\n' for r in rows))
    compressed = args.output / 'artificial.bed.gz'
    pysam.tabix_compress(str(bed), str(compressed), force=False)
    pysam.tabix_index(str(compressed), preset='bed', force=False)
    reader = tabix.open(str(compressed))
    checks = [('1', 19, 20, [rows[0]]), ('1', 20, 21, [rows[1]]),
              ('2', 0, 1, [rows[2]]), ('1', 30, 31, [])]
    observations = []
    for chrom, start, end, expected in checks:
        actual = list(reader.query(chrom, start, end))
        require(actual == expected, 'indexed overlap or coordinate-zero result differs')
        observations.append({'chrom': chrom, 'start_0based': start, 'end_exclusive': end, 'rows': actual})
    cases.append({'name': 'pytabix indexed half-open interval boundaries', 'queries': observations, 'passed': True})
    result = {'completed_at': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'python': sys.version, 'prefix': sys.prefix, 'native_files': native_files,
              'versions': {name: importlib.metadata.version(name) for name in ('pyabpoa', 'pywfa', 'pytabix', 'pysam')},
              'cases': cases, 'all_passed': True,
              'scope': 'Six bounded artificial-data checks exercise three imported native extensions; no genome-wide accuracy or clinical sensitivity claim.'}
    (args.output / 'summary.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'all_passed': True, 'cases': len(cases), 'versions': result['versions']}))


if __name__ == '__main__':
    main()
