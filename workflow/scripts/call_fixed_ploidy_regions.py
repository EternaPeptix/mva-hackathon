#!/usr/bin/env python3
"""Recall each ploidy class independently and reuse an identified bcftools VCF.

GATK 4.6.1.0 resolves --ploidy-regions by the maximum ploidy overlapping an
active region, which can differ from the ploidy at a variant's start. These
independent constant-ploidy calls avoid that contract. Every caller output is
retained; ownership is assigned by the original record's start, before any
normalization. Spanning variants remain explicitly unresolved.
"""
import argparse
from collections import Counter
import heapq
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import time

import pysam

from call_ploidy_regions import inspect_calls, ploidy_at, validate_regions
from prepare_gatk_compatible_cram import digest, save, stamp


def record_contract(record, sample, ploidy):
    values = record.samples[sample]
    gt = values.get('GT')
    if not gt or None in gt or len(gt) != ploidy:
        raise ValueError(f'constant-ploidy genotype contract failed at {record.contig}:{record.pos}')
    if any(a < 0 or a >= len(record.alleles) for a in gt):
        raise ValueError('genotype allele index outside original allele list')
    likelihoods = values.get('PL')
    expected = math.comb(len(record.alleles) + ploidy - 1, ploidy)
    if likelihoods is None or len(likelihoods) != expected:
        raise ValueError('original genotype likelihood cardinality differs from ploidy/alleles')


def merge_fixed_calls(paths, output, regions, sample, excluded_path):
    """Keep the unique model owner; never recode GT or collapse original ALTs."""
    streams = {p: pysam.VariantFile(str(path)) for p, path in paths.items()}
    try:
        first = streams[min(streams)]
        header = first.header.copy()
        header.add_meta('INFO', items=[('ID', 'SF_FIXED_PLOIDY'), ('Number', 1),
            ('Type', 'Integer'), ('Description', 'Constant ploidy of retained original caller output; original GT/alleles preserved')])
        for p, path in sorted(paths.items()):
            header.add_meta('SF_FIXED_PLOIDY_SOURCE', items=[('ID', f'P{p}'),
                ('File', Path(path).name), ('SHA256', digest(path))])
        order = {chrom: i for i, chrom in enumerate(header.contigs)}
        for reader in streams.values():
            if list(reader.header.samples) != [sample] or list(reader.header.contigs) != list(header.contigs):
                raise ValueError('fixed-call sample or sequence dictionary differs')
        counts = Counter()
        spans = []
        previous = {}
        heap = []
        serial = 0

        def advance(p):
            nonlocal serial
            record = next(streams[p], None)
            if record is not None:
                position = (order[record.contig], record.start)
                if p in previous and position < previous[p]:
                    raise ValueError('constant-ploidy source is not coordinate sorted')
                previous[p] = position
                serial += 1
                heapq.heappush(heap, (*position, serial, p, record))

        for p in streams:
            advance(p)
        seen_at_position = set()
        last_position = None
        with pysam.VariantFile(str(output), 'wz', header=header) as writer, Path(excluded_path).open('w') as excluded:
            while heap:
                contig_order, start, _, p, record = heapq.heappop(heap)
                record_contract(record, sample, p)
                expected, boundary = ploidy_at(regions, record.contig, record.start)
                counts[f'P{p}_raw_records'] += 1
                if expected != p:
                    counts[f'P{p}_outside_start_ownership'] += 1
                    excluded.write(json.dumps({'source_ploidy': p, 'expected_ploidy': expected,
                        'chrom': record.contig, 'pos': record.pos, 'ref': record.ref,
                        'alts': list(record.alts or []), 'GT': list(record.samples[sample]['GT']),
                        'reason': 'record start belongs to the other independent ploidy call; original output retained'}) + '\n')
                else:
                    position = (contig_order, start)
                    if position != last_position:
                        seen_at_position.clear()
                        last_position = position
                    key = (record.ref, tuple(record.alts or []))
                    if key in seen_at_position:
                        raise ValueError('duplicate original allele record in owned calls')
                    seen_at_position.add(key)
                    if record.stop > boundary:
                        spans.append({'chrom': record.contig, 'pos': record.pos, 'ref': record.ref,
                            'alts': list(record.alts or []), 'boundary_end_0based': boundary,
                            'disposition': 'retained for individual review; start ownership does not resolve spanning biology'})
                    record.translate(header)
                    record.info['SF_FIXED_PLOIDY'] = p
                    writer.write(record)
                    counts[f'P{p}_owned_records'] += 1
                advance(p)
        return {'record_counts': dict(counts), 'reference_span_crosses_ploidy_boundary': spans,
            'scope': 'Original multi-allele genotype and likelihood preservation with explicit start ownership; not truth validation.'}
    finally:
        for reader in streams.values():
            reader.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('root', 'reference', 'alignment', 'config', 'output', 'bcftools-vcf'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--bcftools-sha256', required=True)
    args = parser.parse_args()
    root, out = args.root.resolve(), args.output.resolve()
    if out.exists():
        raise ValueError('existing fixed-ploidy output preserved')
    out.mkdir(parents=True, exist_ok=False)
    os.nice(10)
    started = time.monotonic()
    state = {'status': 'running', 'pid': os.getpid(), 'started_at': stamp(), 'current_stage': 'validate-inputs'}
    save(out / 'status.json', state)
    commands = []

    def run(command, label):
        state.update(current_stage=label)
        save(out / 'status.json', state)
        row = {'stage': label, 'command': list(map(str, command)), 'started_at': stamp()}
        commands.append(row)
        save(out / 'commands.json', commands)
        with (out / (label + '.log')).open('w') as stream:
            result = subprocess.run(row['command'], stdout=stream, stderr=subprocess.STDOUT)
        row.update(exit_code=result.returncode, completed_at=stamp())
        save(out / 'commands.json', commands)
        if result.returncode:
            raise RuntimeError('fixed-ploidy stage failed: ' + label)

    try:
        config = json.loads(args.config.read_text())
        if digest(args.reference) != config['reference_sha256'] or digest(args.alignment) != config['alignment_sha256']:
            raise ValueError('frozen reference or alignment identity differs')
        if digest(args.bcftools_vcf) != args.bcftools_sha256:
            raise ValueError('reused bcftools identity differs')
        with pysam.FastaFile(str(args.reference)) as fasta, pysam.AlignmentFile(str(args.alignment), reference_filename=str(args.reference)) as reads:
            lengths = dict(zip(fasta.references, fasta.lengths))
            if list(zip(reads.references, reads.lengths)) != list(lengths.items()):
                raise ValueError('alignment/reference dictionary differs')
            reads.check_index()
            groups = reads.header.to_dict().get('RG', [])
            if len(groups) != config['expected_read_groups'] or {g.get('SM') for g in groups} != {config['sample']}:
                raise ValueError('sample/read groups differ')
            regions = validate_regions(config, lengths)
            for mask in config.get('required_N_masks', []):
                seq = fasta.fetch(mask['chrom'], mask['start'], mask['end']).upper()
                if len(seq) != mask['end'] - mask['start'] or set(seq) != {'N'}:
                    raise ValueError('mandatory N mask differs')
        reused_accounting = inspect_calls(args.bcftools_vcf, regions, config['sample'])
        runtime = root / 'tools/bio'
        bcf, samtools = runtime / 'bin/bcftools', runtime / 'bin/samtools'
        java = [runtime / 'lib/jvm/bin/java', '-XX:ActiveProcessorCount=2', '-Xmx6g', '-jar',
                root / 'gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar']
        sources = (Path(__file__), args.config, Path(__file__).with_name('call_ploidy_regions.py'),
            Path(__file__).with_name('prepare_gatk_compatible_cram.py'), bcf, samtools, java[0], java[4])
        identities = {str(p): digest(p) for p in sources}
        state.update(config=config, source_and_tool_sha256=identities,
            reused_bcftools={'path': str(args.bcftools_vcf.resolve()), 'sha256': args.bcftools_sha256,
                'accounting': reused_accounting},
            scope='Complete selected-region GATK recall with one constant-ploidy process per class. Existing bcftools calls reused byte-for-byte. No genotype coercion or previous-run mutation.')
        save(out / 'status.json', state)
        run([samtools, 'quickcheck', '-v', args.alignment], 'quickcheck')
        fixed_paths = {}
        coverage = []
        for p in sorted({p for r in regions for _, _, p in r['intervals_0based_half_open']}):
            rows = [(r['chrom'], start, end) for r in regions for start, end, q in r['intervals_0based_half_open'] if p == q]
            intervals = out / f'ploidy-{p}.bed'
            intervals.write_text(''.join(f'{chrom}\t{start}\t{end}\n' for chrom, start, end in rows))
            coverage.append({'ploidy': p, 'intervals': rows, 'selected_bases': sum(end - start for _, start, end in rows)})
            fixed_paths[p] = out / f'gatk.ploidy-{p}.original.vcf.gz'
            run(java + ['HaplotypeCaller', '-R', args.reference, '-I', args.alignment,
                '-O', fixed_paths[p], '-L', intervals, '--sample-ploidy', str(p),
                '--interval-padding', '0', '--interval-merging-rule', 'OVERLAPPING_ONLY',
                '--pair-hmm-implementation', 'LOGLESS_CACHING',
                '--assembly-region-out', out / f'gatk.ploidy-{p}.assembly-regions.igv'], f'gatk-ploidy-{p}')
        raw = out / 'gatk.original.vcf.gz'
        ownership = merge_fixed_calls(fixed_paths, raw, regions, config['sample'], out / 'gatk.outside-ownership.jsonl')
        reused = out / 'bcftools.original.vcf.gz'
        shutil.copyfile(args.bcftools_vcf, reused)
        if digest(reused) != args.bcftools_sha256:
            raise ValueError('copied bcftools output differs')
        evidence = {}
        for name, path in [('gatk', raw), ('bcftools', reused)]:
            evidence[name] = inspect_calls(path, regions, config['sample'])
            run([bcf, 'index', '-t', path], name + '-original-index')
            normalized = out / (name + '.normalized.vcf.gz')
            run([bcf, 'norm', '-c', 'e', '-f', args.reference, '-m', '-any', '--old-rec-tag', 'SF_ORIG',
                '-Oz', '-o', normalized, path], name + '-normalize')
            run([bcf, 'index', '-t', normalized], name + '-normalized-index')
        for path, expected in identities.items():
            if digest(path) != expected:
                raise ValueError('source or native tool changed during recall')
        save(out / 'evidence.json', {'completed_at': stamp(), 'coverage': coverage,
            'gatk_start_ownership': ownership, 'original_call_accounting': evidence,
            'limitations': ['Computational ploidy model only; metadata remains unchanged.',
                'Original per-ploidy VCFs preserve all emitted calls, including any outside start ownership.',
                'REF spans crossing a ploidy boundary require individual review before integration.',
                'Normalized VCFs provide allele identities; original records provide GT and allele indices.',
                'Variant-only absence is not reference evidence or callability. No clinical filtering or sensitivity claim.']})
        state.update(status='complete', completed_at=stamp(), elapsed_seconds=time.monotonic() - started,
            output_sha256={p.name: digest(p) for p in out.iterdir() if p.is_file() and p.name != 'status.json'})
        save(out / 'status.json', state)
        print(json.dumps({'status': 'complete', 'output': str(out), 'accounting': evidence}))
    except BaseException as error:
        state.update(status='failed', failed_at=stamp(), error_type=type(error).__name__, error=str(error))
        save(out / 'status.json', state)
        raise


if __name__ == '__main__':
    main()
