#!/usr/bin/env python3
"""Additive regional-ploidy discovery; preserve raw multi-allele calls and all prior runs."""
import argparse
from collections import Counter
import json
import os
from pathlib import Path
import subprocess
import time

import pysam

from prepare_gatk_compatible_cram import digest, save, stamp


def validate_regions(config, lengths):
    regions = config['regions']
    if not regions or len({r['chrom'] for r in regions}) != len(regions):
        raise ValueError('unique nonempty contig selection required')
    for region in regions:
        chrom, end = region['chrom'], 0
        if lengths.get(chrom) != region['length']:
            raise ValueError('contig length differs')
        for interval in region['intervals_0based_half_open']:
            start, stop, ploidy = interval
            if start != end or stop <= start or ploidy not in (1, 2):
                raise ValueError('ploidy intervals must cover selected contig without gaps or overlaps')
            end = stop
        if end != region['length']:
            raise ValueError('incomplete contig ploidy definition')
    return regions


def ploidy_at(regions, chrom, position_0based):
    matches = [(p, end) for r in regions if r['chrom'] == chrom
               for start, end, p in r['intervals_0based_half_open'] if start <= position_0based < end]
    if len(matches) != 1:
        raise ValueError('position outside unique ploidy definition')
    return matches[0]


def inspect_calls(path, regions, sample):
    counts, boundaries, missing = Counter(), [], 0
    with pysam.VariantFile(str(path)) as reader:
        if list(reader.header.samples) != [sample]:
            raise ValueError('unexpected sample')
        for record in reader:
            expected, end = ploidy_at(regions, record.contig, record.start)
            gt = record.samples[sample].get('GT')
            if not gt or None in gt:
                missing += 1
            if gt and len(gt) != expected:
                raise ValueError(f'genotype/model ploidy mismatch at {record.contig}:{record.pos}')
            counts[str(expected)] += 1
            if record.stop > end:
                boundaries.append({'chrom': record.contig, 'pos': record.pos, 'ref': record.ref,
                                   'alts': list(record.alts or []), 'boundary_end_0based': end})
    return {'records_by_expected_ploidy': dict(counts), 'missing_or_partial_GT_records': missing,
            'reference_span_crosses_ploidy_boundary': boundaries,
            'scope': 'Genotype shape and interval accounting only, not truth or clinical validation.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('root', 'reference', 'alignment', 'config', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    os.nice(10)
    root, out = args.root.resolve(), args.output.resolve()
    if out.exists():
        raise ValueError('existing ploidy run preserved')
    started = time.monotonic()
    config = json.loads(args.config.read_text())
    if digest(args.reference) != config['reference_sha256'] or digest(args.alignment) != config['alignment_sha256']:
        raise ValueError('input identity differs from frozen ploidy model')
    with pysam.FastaFile(str(args.reference)) as fasta, pysam.AlignmentFile(str(args.alignment), reference_filename=str(args.reference)) as reads:
        lengths = dict(zip(fasta.references, fasta.lengths))
        if list(zip(reads.references, reads.lengths)) != list(lengths.items()):
            raise ValueError('reference/alignment dictionary differs')
        reads.check_index()
        groups = reads.header.to_dict().get('RG', [])
        if {g.get('SM') for g in groups} != {config['sample']} or len(groups) != config['expected_read_groups']:
            raise ValueError('sample/read groups differ')
        regions = validate_regions(config, lengths)
        for masked in config.get('required_N_masks', []):
            sequence = fasta.fetch(masked['chrom'], masked['start'], masked['end']).upper()
            if len(sequence) != masked['end'] - masked['start'] or set(sequence) != {'N'}:
                raise ValueError('required reference mask differs')
    out.mkdir(parents=True, exist_ok=False)
    runtime = root / 'tools/bio'
    bcf, samtools = runtime / 'bin/bcftools', runtime / 'bin/samtools'
    java = [runtime / 'lib/jvm/bin/java', '-XX:ActiveProcessorCount=2', '-Xmx6g', '-jar',
            root / 'gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar']
    state = {'status': 'running', 'pid': os.getpid(), 'started_at': stamp(), 'scope': config['scope'],
             'config': config, 'source_and_tool_sha256': {str(p): digest(p) for p in
                 (Path(__file__), args.config, samtools, bcf, java[0], java[4], Path(__file__).with_name('prepare_gatk_compatible_cram.py'))}}
    commands = []

    def run(command, label):
        state.update(current_stage=label)
        save(out / 'status.json', state)
        row = {'stage': label, 'command': list(map(str, command)), 'started_at': stamp()}
        commands.append(row)
        save(out / 'commands.json', commands)
        with (out / (label + '.log')).open('w') as log:
            p = subprocess.run(row['command'], stdout=log, stderr=subprocess.STDOUT)
        row.update(exit_code=p.returncode, completed_at=stamp())
        save(out / 'commands.json', commands)
        if p.returncode:
            raise RuntimeError('ploidy stage failed: ' + label)

    try:
        interval_bed, gatk_bed, bcf_ploidy, samples = [out / name for name in ('regions.bed', 'gatk-ploidy.bed', 'bcftools-ploidy.txt', 'sample-model.txt')]
        interval_bed.write_text(''.join(f"{r['chrom']}\t0\t{r['length']}\n" for r in regions))
        gatk_bed.write_text(''.join(f"{r['chrom']}\t{start}\t{end}\t{p}\n" for r in regions for start, end, p in r['intervals_0based_half_open']))
        bcf_ploidy.write_text(''.join(f"{r['chrom']}\t{start + 1}\t{end}\tDEPTH_MODEL\t{p}\n" for r in regions for start, end, p in r['intervals_0based_half_open']))
        samples.write_text(config['sample'] + '\tDEPTH_MODEL\n')
        run([samtools, 'quickcheck', '-v', args.alignment], 'quickcheck')
        gatk = out / 'gatk.original.vcf.gz'
        run(java + ['HaplotypeCaller', '-R', args.reference, '-I', args.alignment, '-O', gatk,
            '-L', interval_bed, '--sample-ploidy', '2', '--ploidy-regions', gatk_bed,
            '--pair-hmm-implementation', 'LOGLESS_CACHING'], 'gatk')
        run([bcf, 'mpileup', '-Ob', '-o', out / 'pileup.bcf', '-f', args.reference, '-R', interval_bed,
             '-a', 'AD,DP,SP', '-q', '20', '-Q', '20', '-d', '10000', '--ff', '0xF04', args.alignment], 'pileup')
        original_bcf = out / 'bcftools.original.vcf.gz'
        run([bcf, 'call', '-mv', '-A', '-a', 'GQ', '--ploidy-file', bcf_ploidy, '--samples-file', samples,
             '-Oz', '-o', original_bcf, out / 'pileup.bcf'], 'bcftools')
        evidence = {}
        for name, raw in [('gatk', gatk), ('bcftools', original_bcf)]:
            run([bcf, 'index', '-t', '-f', raw], name + '-original-index')
            evidence[name] = inspect_calls(raw, regions, config['sample'])
            normalized = out / (name + '.normalized.vcf.gz')
            run([bcf, 'norm', '-c', 'e', '-f', args.reference, '-m', '-any', '--old-rec-tag', 'SF_ORIG',
                 '-Oz', '-o', normalized, raw], name + '-normalize')
            run([bcf, 'index', '-t', normalized], name + '-normalized-index')
        save(out / 'evidence.json', {'config': config, 'original_call_accounting': evidence,
            'limitations': ['Ploidy assumptions are a computational model, not clinical karyotype or a mosaicism test.',
                'Raw original multi-allele calls are authoritative for GT and allele indices; normalized records are identities.',
                'Original-call REF spans crossing a ploidy boundary require individual review before integration.',
                'Variant-only absence does not establish callability or reference genotype. No BQSR or clinical filtering.',
                'Prior all-contig diploid discovery and original source alignments remain unchanged.']})
        state.update(status='complete', completed_at=stamp(), elapsed_seconds=time.monotonic() - started,
                     output_sha256={p.name: digest(p) for p in out.iterdir() if p.is_file() and p.name != 'status.json'})
        save(out / 'status.json', state)
        print(json.dumps({'status': 'complete', 'output': str(out), 'original_calls': evidence}))
    except BaseException as error:
        state.update(status='failed', failed_at=stamp(), error_type=type(error).__name__)
        save(out / 'status.json', state)
        raise


if __name__ == '__main__':
    main()
