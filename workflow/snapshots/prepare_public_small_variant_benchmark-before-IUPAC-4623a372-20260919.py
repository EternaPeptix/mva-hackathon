#!/usr/bin/env python3
"""Verify and prepare the pinned public chr20 truth, BED and RTG reference.

No BED clipping of truth alleles, renaming, normalization or genotype editing is
performed. A complete reference round trip checks all 195 sequences in the SDF.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import time

import pysam

from rtg_small_variant_helpers import audit_vcf, execute, fingerprint, runtime, save, stamp

REFERENCE_SHA = '9cce8b926416dd96b152deea85188495b75f7ac8d634cc723a017067be8702b7'
CHR20_SHA = '2442701218d3b2f7e2d7cacd78774e1628fd48a7a13b58bee53ab8b7cb3d1af4'


def sequence_manifest(path):
    rows = []
    with pysam.FastxFile(str(path)) as fasta:
        for record in fasta:
            sequence = record.sequence.upper()
            rows.append({'name': record.name, 'bases': len(sequence),
                         'uppercase_sha256': hashlib.sha256(sequence.encode('ascii')).hexdigest()})
    if len({r['name'] for r in rows}) != len(rows): raise ValueError('duplicate sequence name')
    return rows


def run(root, output, config_path):
    if output.resolve().is_relative_to(Path(__file__).resolve().parents[2]):
        raise ValueError('large public artifacts must remain outside the repository')
    output.mkdir(parents=True, exist_ok=False); start = time.monotonic()
    state = {'status': 'running', 'pid': os.getpid(), 'started_at': stamp()}; save(output / 'status.json', state)
    try:
        config = json.loads(config_path.read_text()); rtg, rtg_files = runtime(root, config)
        inputs = root / 'ref/giab-small-variants-20260919/inputs'
        sources = {name: fingerprint(inputs / name) for name in config['truth_sha256']}
        for name, expected in config['truth_sha256'].items():
            if sources[name]['sha256'] != expected: raise ValueError('public truth input changed: ' + name)
        reference = root / 'ref/giab-sv-20260919/assembly-reference/GRCh38_no_alt.fa'
        reference_file = fingerprint(reference)
        if reference_file['sha256'] != REFERENCE_SHA: raise ValueError('public reference changed')
        source_sequences = sequence_manifest(reference)
        if len(source_sequences) != 195: raise ValueError('unexpected full reference dictionary')
        chr20 = [r for r in source_sequences if r['name'] == 'chr20']
        if len(chr20) != 1 or chr20[0]['uppercase_sha256'] != CHR20_SHA or chr20[0]['bases'] != 64444167:
            raise ValueError('public chr20 sequence changed')
        config_copy = output / 'frozen-comparison-config.json'; shutil.copyfile(config_path, config_copy)
        shutil.copyfile(Path(__file__), output / 'executed-prepare_public_small_variant_benchmark.py')
        truth = output / 'HG002.chr20.v4.2.1.full-context.vcf.gz'
        source_vcf = inputs / 'HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz'
        n = 0
        with pysam.VariantFile(str(source_vcf)) as original, pysam.VariantFile(str(truth), 'wz', header=original.header) as target:
            for record in original.fetch('chr20'):
                target.write(record); n += 1
        pysam.tabix_index(str(truth), preset='vcf', force=False)
        truth_audit = audit_vcf(truth, reference, 'HG002', 'chr20')
        if truth_audit['counts']['records'] != n: raise ValueError('truth extraction changed record count')
        # Independent serialization checks verify the extraction changed no fields.
        with pysam.VariantFile(str(source_vcf)) as original, pysam.VariantFile(str(truth)) as target:
            for a, b in zip(original.fetch('chr20'), target, strict=True):
                if str(a) != str(b): raise ValueError('truth record changed during extraction')
        source_bed = inputs / 'HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed'
        rows = []; lines = []; prior_end = -1
        for line in source_bed.read_text().splitlines():
            fields = line.split('\t')
            if fields[0] != 'chr20': continue
            if len(fields) != 3: raise ValueError('unexpected public BED columns')
            left, right = map(int, fields[1:])
            if not 0 <= left < right <= 64444167 or left < prior_end:
                raise ValueError('invalid or overlapping public BED intervals')
            rows.append([left, right]); lines.append(line); prior_end = right
        if not rows: raise ValueError('chr20 benchmark BED is empty')
        bed = output / 'official-chr20-v4.2.1.bed'; bed.write_text('\n'.join(lines) + '\n')
        state['stage'] = 'format and round-trip complete public reference'; save(output / 'status.json', state)
        sdf = output / 'reference.sdf'; exported = output / 'exported-reference.fasta'
        commands = [execute([rtg, 'format', '-o', sdf, reference], output, 'format-reference'),
                    execute([rtg, 'sdfstats', '--lengths', sdf], output, 'reference-sdfstats'),
                    execute([rtg, 'sdf2fasta', '-i', sdf, '-o', exported, '-Z', '-l', '60'], output, 'roundtrip-reference')]
        exported_sequences = sequence_manifest(exported)
        if exported_sequences != source_sequences: raise ValueError('SDF round-trip altered dictionary or sequence')
        sdf_files = {str(p.relative_to(sdf)): fingerprint(p) for p in sorted(sdf.rglob('*')) if p.is_file()}
        result = {'status': 'complete', 'completed_at': stamp(), 'elapsed_seconds': time.monotonic() - start,
                  'scope': config['scope'], 'config': fingerprint(config_copy), 'sources': sources,
                  'download_receipt': fingerprint(inputs / 'input-receipt.json'),
                  'publisher_integrity_note': config['publisher_integrity_note'],
                  'runtime': rtg_files, 'reference': reference_file, 'reference_sequences': source_sequences,
                  'SDF_roundtrip_all_195_sequences_identical': True, 'sdf_files': sdf_files,
                  'exported_reference': fingerprint(exported), 'commands': commands,
                  'truth': fingerprint(truth), 'truth_index': fingerprint(Path(str(truth) + '.tbi')),
                  'truth_audit': truth_audit, 'truth_all_fields_preserved': True,
                  'bed': fingerprint(bed), 'bed_intervals': len(rows), 'bed_bases': sum(b-a for a,b in rows),
                  'bed_fraction_of_chr20': sum(b-a for a,b in rows)/64444167,
                  'code': fingerprint(Path(__file__))}
        save(output / 'summary.json', result)
        state.update(status='complete', completed_at=result['completed_at']); save(output / 'status.json', state)
        return result
    except BaseException as error:
        state.update(status='failed', failed_at=stamp(), error_type=type(error).__name__, error=str(error))
        save(output / 'status.json', state); raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('root', 'output', 'config'): parser.add_argument('--' + key, type=Path, required=True)
    args = parser.parse_args()
    result = run(args.root, args.output, args.config)
    print(json.dumps({k: result[k] for k in ('status', 'completed_at', 'bed_bases', 'bed_intervals', 'truth_audit')}))
