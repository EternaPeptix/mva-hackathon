#!/usr/bin/env python3
"""Run bounded, read-only-input WhatsHap comparisons and preserve phase uncertainty."""
import argparse
from collections import Counter, defaultdict
import datetime
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

import pysam

PHASE_FIELDS = ('PS', 'PGT', 'PID', 'HP', 'PQ')


def fingerprint(path):
    path = Path(path).resolve()
    before = path.stat()
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(block)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError('input changed while hashing')
    return {'path': str(path), 'bytes': after.st_size, 'mtime_ns': after.st_mtime_ns,
            'sha256': digest.hexdigest()}


def genotype_signature(record, sample_name):
    gt = record.samples[sample_name].get('GT')
    unordered = None if gt is None else tuple(sorted(gt, key=lambda x: -1 if x is None else x))
    return (record.chrom, record.pos, record.ref, record.alts, unordered)


def evidence_signature(record, sample_name):
    """Preserve non-phase sample evidence through input preparation and phasing."""
    return tuple(sorted((key, value) for key, value in record.samples[sample_name].items()
                        if key not in PHASE_FIELDS and key != 'GT'))


def ps_phase(record, sample_name, alt_index):
    """Use this run's PS/GT only; legacy PGT/PID must not override fresh phasing."""
    sample = record.samples[sample_name]
    gt, ps = sample.get('GT'), sample.get('PS')
    if (not sample.phased or not gt or len(gt) != 2 or None in gt
            or gt.count(alt_index) != 1 or ps in (None, '.', '')):
        return None
    return (record.chrom, str(ps), gt.index(alt_index))


def relation(first, second):
    if first is None or second is None or first[:2] != second[:2]:
        return 'unknown'
    return 'cis' if first[2] == second[2] else 'trans'


def prepare_vcf(source, output, reference, sample_name, region):
    counts = Counter()
    signatures = []
    with pysam.VariantFile(str(source)) as reader, pysam.FastaFile(str(reference)) as fasta:
        if list(reader.header.samples) != [sample_name]:
            raise ValueError('expected exactly the configured sample in VCF')
        with pysam.VariantFile(str(output), 'w', header=reader.header) as writer:
            for record in reader.fetch(region['chrom'], region['start'] - 1, region['end']):
                if record.pos < region['start'] or record.pos > region['end']:
                    continue
                counts['source_records_in_region'] += 1
                if set(record.filter.keys()) - {'PASS'}:
                    counts['excluded_filter_records'] += 1
                    continue
                if fasta.fetch(record.chrom, record.start, record.start + len(record.ref)).upper() != record.ref:
                    raise ValueError('regional variant reference mismatch')
                sample = record.samples[sample_name]
                signatures.append(genotype_signature(record, sample_name))
                sample.phased = False
                # Remove old fields, rather than assigning None to String values:
                # pysam 0.24 can round-trip those empty strings as NUL bytes.
                for tag in PHASE_FIELDS:
                    if tag in record.format:
                        del record.format[tag]
                writer.write(record)
                counts['retained_records'] += 1
    if not signatures:
        raise ValueError('no eligible records in region')
    return dict(counts), signatures


def summarize(input_vcf, phased_vcf, sample_name, targets):
    for path in (input_vcf, phased_vcf):
        if b'\x00' in Path(path).read_bytes():
            raise ValueError('NUL byte in regional VCF')
    with pysam.VariantFile(str(input_vcf)) as reader:
        initial = list(reader)
        before = [genotype_signature(r, sample_name) for r in initial]
        evidence_before = [evidence_signature(r, sample_name) for r in initial]
    phases, target_rows = [], []
    blocks = defaultdict(list)
    het_positions = []
    with pysam.VariantFile(str(phased_vcf)) as reader:
        if list(reader.header.samples) != [sample_name]:
            raise ValueError('phaser changed sample identity')
        records = list(reader)
        after = [genotype_signature(r, sample_name) for r in records]
        if before != after:
            raise ValueError('phaser changed variant records or unphased genotypes')
        if evidence_before != [evidence_signature(r, sample_name) for r in records]:
            raise ValueError('phaser changed non-phase sample evidence')
        for r in records:
            sample = r.samples[sample_name]
            gt = sample.get('GT')
            if not gt or len(gt) != 2 or None in gt or gt[0] == gt[1]:
                continue
            het_positions.append(r.pos)
            ps = sample.get('PS')
            if sample.phased and ps not in (None, '.', ''):
                blocks[(r.chrom, str(ps))].append(r.pos)
        for target in targets:
            found = [r for r in records if r.chrom == target['chrom'] and r.pos == target['pos']
                     and r.ref == target['ref'] and target['alt'] in (r.alts or ())]
            phase = None
            row = {'target': target, 'matching_records': len(found)}
            if len(found) == 1:
                record = found[0]
                sample = record.samples[sample_name]
                phase = ps_phase(record, sample_name, record.alts.index(target['alt']) + 1)
                row.update(GT=sample.get('GT'), phased=sample.phased, PS=sample.get('PS'),
                           alternate_haplotype_index=phase[2] if phase else None)
                positions = blocks.get(phase[:2], []) if phase else []
                row.update(block_heterozygous_records=len(positions),
                           block_span_bp=max(positions) - min(positions) if positions else None)
            target_rows.append(row)
            phases.append(phase)
    endpoints = sorted(t['pos'] for t in targets)
    between = sorted({p for p in het_positions if endpoints[0] <= p <= endpoints[-1]})
    spans = [max(v) - min(v) for v in blocks.values()]
    return {'records_and_unordered_genotypes_preserved': True,
            'non_phase_sample_evidence_preserved': True, 'record_count': len(after),
            'heterozygous_records': len(het_positions), 'phased_heterozygous_records': sum(map(len, blocks.values())),
            'phase_blocks': len(blocks), 'maximum_block_span_bp': max(spans, default=0),
            'heterozygous_records_strictly_between_candidates': sum(endpoints[0] < p < endpoints[-1] for p in het_positions),
            'maximum_gap_between_intervening_heterozygous_positions_bp':
                max((b - a for a, b in zip(between, between[1:])), default=None),
            'targets': target_rows, 'candidate_pair_relation': relation(*phases)}


def phase_command(executable, vcf, bam, reference, sample, contig, output, read_list, only_snvs=False):
    cmd = [str(executable), 'phase', '--reference', str(reference), '--sample', sample,
           '--chromosome', contig, '--tag', 'PS', '--mapping-quality', '20',
           '--no-genetic-haplotyping', '--output-read-list', str(read_list), '-o', str(output)]
    if only_snvs:
        cmd.append('--only-snvs')
    return cmd + [str(vcf), str(bam)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root, output = args.root.resolve(), args.output.resolve()
    config = json.loads(args.config.read_text())
    if len(config['targets']) != 2:
        parser.error('exactly two target alleles are required')
    output.mkdir(parents=True, exist_ok=False)
    os.nice(10)
    started = time.monotonic()
    state = {'status': 'running', 'started_at': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'commands': []}

    def save_state():
        (output / 'status.json').write_text(json.dumps(state, indent=2) + '\n')

    def run(cmd, label):
        state['commands'].append(cmd)
        state['current_stage'] = label
        save_state()
        with (output / (label + '.log')).open('w') as log:
            subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True,
                           env=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1'))

    try:
        reference = root / config['reference']
        bam = root / config['bam']
        executable = root / config['whatshap']
        samtools = root / 'tools/bio/bin/samtools'
        inputs = {'reference': fingerprint(reference), 'reference_index': fingerprint(Path(str(reference) + '.fai')),
                  'bam': fingerprint(bam),
                  'bam_index': fingerprint(Path(str(bam) + '.bai')),
                  'config': fingerprint(args.config), 'script': fingerprint(Path(__file__))}
        if inputs['reference']['sha256'] != config['reference_sha256']:
            raise ValueError('reference identity differs from configuration')
        version = subprocess.run([str(executable), '--version'], capture_output=True, text=True, check=True).stdout.strip()
        if version != config['whatshap_version']:
            raise ValueError('unreviewed phaser version')
        with pysam.AlignmentFile(str(bam), 'rb') as reader, pysam.FastaFile(str(reference)) as fasta:
            if dict(zip(reader.references, reader.lengths)) != dict(zip(fasta.references, fasta.lengths)):
                raise ValueError('alignment/reference dictionary mismatch')
            groups = reader.header.to_dict().get('RG', [])
            if not groups or any(g.get('SM') != config['sample'] for g in groups):
                raise ValueError('alignment sample metadata differs from configuration')
        filtered = output / 'filtered.bam'
        run([str(samtools), 'view', '-b', '-F', '0xF04', '-q', '20', '-o', str(filtered), str(bam)], 'filter-bam')
        run([str(samtools), 'index', str(filtered)], 'index-bam')
        runs = []
        for name, source_relative in config['vcfs'].items():
            source = root / source_relative
            inputs[name + '_vcf'] = fingerprint(source)
            for extension in ('.tbi', '.csi'):
                index = Path(str(source) + extension)
                if index.exists():
                    inputs[name + '_index' + extension] = fingerprint(index)
            vcf = output / (name + '.input.vcf')
            selection, _ = prepare_vcf(source, vcf, reference, config['sample'], config['region'])
            for only_snvs in (False, True):
                label = name + ('-snvs' if only_snvs else '-all')
                phased = output / (label + '.vcf')
                cmd = phase_command(executable, vcf, filtered, reference, config['sample'],
                                    config['region']['chrom'], phased, output / (label + '.reads.tsv'), only_snvs)
                t0 = time.monotonic()
                run(cmd, label)
                elapsed = time.monotonic() - t0
                result = summarize(vcf, phased, config['sample'], config['targets'])
                result.update(name=label, only_snvs=only_snvs, selection=selection,
                              phase_elapsed_seconds=elapsed, input_vcf=fingerprint(vcf), output_vcf=fingerprint(phased))
                runs.append(result)
                run([str(executable), 'stats', '--sample', config['sample'], '--chromosome', config['region']['chrom'],
                     '--tsv', str(output / (label + '.stats.tsv')), '--block-list', str(output / (label + '.blocks.tsv')),
                     str(phased)], label + '-stats')
        for value in inputs.values():
            stat = Path(value['path']).stat()
            if (stat.st_size, stat.st_mtime_ns) != (value['bytes'], value['mtime_ns']):
                raise ValueError('input changed during phasing')
        evidence = {'scope': config['scope'], 'inputs': inputs, 'whatshap_version': version,
                    'region': config['region'], 'filtered_bam': fingerprint(filtered), 'runs': runs,
                    'read_policy': {'excluded_flags': '0xF04', 'minimum_mapping_quality': 20,
                                    'hard_minimum_base_quality': None,
                                    'base_quality_handling': 'WhatsHap realignment/quality model; no BQ20 hard filter',
                                    'supplementary_used': False, 'duplicates_used': False,
                                    'statistical_population_or_pedigree_phasing': False,
                                    'prior_VCF_phase_fields_removed': list(PHASE_FIELDS)},
                    'limitations': ['Same two retained lanes and DNA sample; not new biological replication.',
                                    'Bounded region and caller-provided genotypes; uncalled markers/classes remain outside scope.',
                                    'A phase block is software-derived linkage evidence, not parental origin or clinical classification.',
                                    'Separate or missing PS blocks leave the pair unknown; no phase is inferred from phenotype.'],
                    'completed_at': datetime.datetime.now(datetime.timezone.utc).isoformat(),
                    'elapsed_seconds': time.monotonic() - started, 'commands': state['commands']}
        (output / 'evidence.json').write_text(json.dumps(evidence, indent=2) + '\n')
        state.update(status='complete', current_stage='complete', completed_at=evidence['completed_at'],
                     elapsed_seconds=evidence['elapsed_seconds'], evidence_sha256=fingerprint(output / 'evidence.json')['sha256'])
        save_state()
        print(json.dumps({'runs': [{'name': x['name'], 'phase': x['candidate_pair_relation'],
                                   'blocks': x['phase_blocks']} for x in runs], 'elapsed_seconds': evidence['elapsed_seconds']}))
    except BaseException as error:
        state.update(status='failed', error_type=type(error).__name__)
        save_state()
        raise


if __name__ == '__main__':
    main()
