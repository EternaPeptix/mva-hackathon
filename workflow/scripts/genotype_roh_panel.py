#!/usr/bin/env python3
"""Explicit marker genotyping and exploratory, gap-separated ROH sensitivity.

All alignment subsets, site-level genotypes and model outputs stay outside Git.
This research workflow does not assign ancestry, UPD, clinical status or EPCR.
"""
import argparse
from collections import Counter
import datetime
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import time

import pysam

from panel_marker_evidence import assess_marker
from segmented_roh import POPULATION_TAGS, run_segmented


PATIENT_SCOPE = 'verified_four_lane_autosomal_marker_genotyping'
PANEL_FILES = ('panel.vcf.gz', 'panel.vcf.gz.tbi', 'sites.bed',
               'alleles.tsv.gz', 'alleles.tsv.gz.tbi')


def stamp():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def fingerprint(path):
    path = Path(path).resolve()
    before = path.stat()
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(block)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError('input changed while hashing')
    return {'path': str(path), 'bytes': after.st_size, 'mtime_ns': after.st_mtime_ns,
            'sha256': h.hexdigest()}


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def validate_config(config):
    if config['scope'] not in ('synthetic', PATIENT_SCOPE):
        raise ValueError('unreviewed genotyping scope')
    chroms = config['chromosomes']
    if not chroms or len(chroms) != len(set(chroms)) or any(c not in map(str, range(1, 23)) for c in chroms):
        raise ValueError('distinct explicit autosomes required')
    if config['scope'] == PATIENT_SCOPE and chroms != list(map(str, range(1, 23))):
        raise ValueError('patient analysis requires all 22 autosomes')
    if config['excluded_flags'] != 0xF04 or config['exclude_unknown_MQ_255'] is not True:
        raise ValueError('unreviewed read flag or unknown mapping-quality policy')
    if (not 1 <= config['minimum_mapping_quality'] <= 254
            or not 1 <= config['minimum_base_quality'] <= 93
            or config['pileup_maximum_depth'] <= config['evidence']['maximum_depth']):
        raise ValueError('invalid read-quality or pileup depth policy')
    evidence = config['evidence']
    if (not 1 <= evidence['minimum_depth'] <= evidence['maximum_depth']
            or evidence['nonpanel_minimum_reads'] < 1
            or not 0 < evidence['nonpanel_minimum_fraction'] <= 1):
        raise ValueError('invalid marker-evidence thresholds')
    if config['population_profiles'] != list(POPULATION_TAGS):
        raise ValueError('pooled and five explicit population sensitivities required')
    gaps = config['pooled_gap_sensitivity_bp']
    if (not gaps or gaps != sorted(set(gaps)) or any(type(n) is not int or n < 1 for n in gaps)
            or config['primary_maximum_gap_bp'] not in gaps or config['minimum_segment_markers'] < 2):
        raise ValueError('invalid gap or segment policy')


def verify_completion(config, alignment, identity, completion_path):
    if config['scope'] == 'synthetic':
        return {'scope': 'Synthetic input, no patient completion claim'}
    if completion_path is None:
        raise ValueError('verified final four-lane completion receipt required')
    receipt = fingerprint(completion_path)
    state = json.loads(completion_path.read_text())
    if (state.get('status') != 'complete' or state.get('expected_lanes') != 4
            or alignment.suffix != '.cram'
            or Path(state.get('cram_path', '')).resolve() != alignment.resolve()
            or state.get('cram_bytes') != identity['bytes']
            or state.get('cram_sha256') != identity['sha256']):
        raise ValueError('final four-lane completion receipt does not match input')
    return receipt


def load_panel(panel, reference, config, summary):
    """Verify every selected public allele and every representation, independently of GT."""
    chroms = config['chromosomes']
    if summary['config']['reference_sha256'] != config['reference_sha256']:
        raise ValueError('population panel uses a different reference')
    if summary['config']['chromosomes'] != chroms:
        raise ValueError('population panel chromosome scope differs')
    order = {c: i for i, c in enumerate(chroms)}
    counts = Counter()
    rows = []
    previous = None
    with pysam.FastaFile(str(reference)) as fasta, pysam.VariantFile(str(panel / 'panel.vcf.gz')) as source:
        if list(source.header.samples) or list(source.header.contigs) != chroms:
            raise ValueError('independent sites-only panel with exact chromosome order required')
        for c in chroms:
            if source.header.contigs[c].length != fasta.get_reference_length(c):
                raise ValueError('panel/reference chromosome length mismatch')
        for row in source:
            key = (order[row.contig], row.pos)
            if previous is not None and key <= previous:
                raise ValueError('panel positions must be strictly ordered and unique')
            if (len(row.alleles) != 2 or row.ref == row.alts[0]
                    or any(len(a) != 1 or a not in 'ACGT' for a in row.alleles)
                    or tuple(row.filter) != ('PASS',)):
                raise ValueError('selected PASS biallelic SNP required')
            if fasta.fetch(row.contig, row.pos - 1, row.pos).upper() != row.ref:
                raise ValueError('panel REF differs from alignment reference')
            for tag in POPULATION_TAGS:
                af = row.info.get(tag)
                if type(af) not in (int, float) or not math.isfinite(af) or not 0 < af < 1:
                    raise ValueError('missing or unusable external frequency profile')
            rows.append(row.copy()); counts[row.contig] += 1; previous = key
        panel_header = source.header.copy()
    if (dict(counts) != summary['selected_markers_by_chromosome']
            or len(rows) != summary['selected_markers'] or set(counts) != set(chroms)):
        raise ValueError('panel marker counts differ from preparation receipt')
    expected_bed = ''.join(f'{r.contig}\t{r.pos-1}\t{r.pos}\n' for r in rows)
    expected_alleles = ''.join(f'{r.contig}\t{r.pos}\t{r.ref},{r.alts[0]}\n' for r in rows)
    with gzip.open(panel / 'alleles.tsv.gz', 'rt') as stream:
        actual_alleles = stream.read()
    if (panel / 'sites.bed').read_text() != expected_bed or actual_alleles != expected_alleles:
        raise ValueError('population VCF, BED and allele-table representations disagree')
    for name in ('panel.vcf.gz', 'alleles.tsv.gz'):
        with pysam.TabixFile(str(panel / name)) as indexed:
            if set(indexed.contigs) != set(chroms):
                raise ValueError('panel index contigs differ')
            for chrom in chroms:
                if sum(1 for _ in indexed.fetch(chrom)) != counts[chrom]:
                    raise ValueError('panel index marker counts differ')
    return rows, panel_header


def genotype_label(gt):
    if gt is None or not gt or any(x is None for x in gt):
        return 'missing'
    if len(gt) != 2 or any(x not in (0, 1) for x in gt):
        raise ValueError('unexpected called genotype in diploid biallelic panel')
    return ('homref', 'heterozygous', 'homalt')[sum(gt)]


def record_evidence(record, sample):
    return {'alleles': record.alleles,
            'sample': {key: record.samples[sample].get(key) for key in ('GT', 'AD', 'DP', 'GQ', 'PL')}}


def qualify_calls(rows, panel_header, pileup_path, calls_path, sample, policy, output):
    """Join by the exact ordered marker key; never substitute another population allele."""
    keys = [(r.contig, r.pos) for r in rows]
    key_order = {key: i for i, key in enumerate(keys)}
    pileups = {}
    previous = -1
    with pysam.VariantFile(str(pileup_path)) as source:
        if list(source.header.samples) != [sample]:
            raise ValueError('original pileup sample identity mismatch')
        for record in source:
            key = (record.contig, record.pos)
            if key not in key_order or key_order[key] <= previous:
                raise ValueError('original pileup has unexpected, duplicate or reordered markers')
            i = key_order[key]
            if record.ref != rows[i].ref:
                raise ValueError('original pileup reference allele mismatch')
            data = record_evidence(record, sample)
            ad, dp = data['sample']['AD'], data['sample']['DP']
            if ad is not None and dp is not None and all(type(x) is int and x >= 0 for x in ad) and sum(ad) > dp:
                raise ValueError('original allele depths exceed reported depth')
            pileups[key] = data; previous = i
    totals = Counter(); per_chrom = {}; exclusions = Counter(); depths = Counter()
    with pysam.VariantFile(str(calls_path)) as source:
        if list(source.header.samples) != [sample]:
            raise ValueError('constrained call sample identity mismatch')
        header = source.header.copy()
        for tag in panel_header.info:
            if tag in header.info:
                raise ValueError('population annotation unexpectedly already exists in caller output')
            field = panel_header.info[tag]
            header.info.add(tag, field.number, field.type, field.description)
        header.info.add('MVA_USABLE', 0, 'Flag', 'Meets configured exploratory marker evidence rules; not clinical status')
        header.info.add('MVA_EXCLUDE', '.', 'String', 'Explicit exclusion reasons; excluded is not reference')
        raw = output / 'annotated-all-markers.vcf.gz'
        qualified = output / 'qualified-markers.vcf.gz'
        with pysam.VariantFile(str(raw), 'wz', header=header) as all_writer, \
                pysam.VariantFile(str(qualified), 'wz', header=header) as usable_writer, \
                gzip.open(output / 'marker-evidence.jsonl.gz', 'wt') as evidence_log:
            for i, record in enumerate(source):
                if i >= len(rows) or (record.contig, record.pos) != keys[i]:
                    raise ValueError('constrained caller omitted, reordered or added a requested marker')
                panel_record = rows[i]
                original = pileups.get(keys[i]); call = record_evidence(record, sample)
                decision = assess_marker(panel_record.alleles, original, call, **policy)
                category = genotype_label(call['sample']['GT'])
                chrom_counts = per_chrom.setdefault(record.contig, Counter())
                for counter in (totals, chrom_counts):
                    counter['requested_and_returned_markers'] += 1
                    counter['raw_genotype_' + category] += 1
                    if call['sample']['GQ'] is None: counter['raw_GQ_missing'] += 1
                    counter['qualified_markers' if decision['usable_for_declared_ROH_analysis'] else 'excluded_markers'] += 1
                exclusions.update(decision['exclusion_reasons'])
                depths[str(call['sample']['DP'])] += 1
                annotated = record.copy(); annotated.translate(header)
                for tag in panel_header.info:
                    annotated.info[tag] = panel_record.info[tag]
                if decision['usable_for_declared_ROH_analysis']:
                    annotated.info['MVA_USABLE'] = True
                    usable_writer.write(annotated)
                else:
                    annotated.info['MVA_EXCLUDE'] = tuple(decision['exclusion_reasons'])
                all_writer.write(annotated)
                evidence_log.write(json.dumps({'chrom': record.contig, 'pos_1based': record.pos,
                    'selected_alleles': panel_record.alleles,
                    'population_frequencies': {tag: panel_record.info[tag] for tag in POPULATION_TAGS},
                    'original_pileup': original, 'raw_constrained_call': call, 'qualification': decision},
                    allow_nan=False) + '\n')
    if totals['requested_and_returned_markers'] != len(rows):
        raise ValueError('constrained caller returned fewer than all requested markers')
    for path in (raw, qualified):
        pysam.tabix_index(str(path), preset='vcf', force=False)
    return {'counts': dict(totals), 'per_chromosome': {c: dict(n) for c, n in per_chrom.items()},
            'exclusion_reason_counts_nonexclusive': dict(exclusions), 'raw_depth_histogram': dict(depths),
            'original_pileup_records': len(pileups)}


def run_analysis(alignment, reference, panel, config, bcftools, samtools, output, completion_state=None):
    validate_config(config)
    if config['scope'] == PATIENT_SCOPE and output.resolve().is_relative_to(Path(__file__).resolve().parents[2]):
        raise ValueError('patient genotype and read artifacts must remain outside the repository')
    output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    state = {'status': 'running', 'started_at': stamp(), 'pid': os.getpid(), 'stage': 'input provenance'}
    save(output / 'status.json', state)
    commands = []
    try:
        index_suffix = '.crai' if alignment.suffix == '.cram' else '.bai'
        inputs = {name: fingerprint(path) for name, path in {
            'alignment': alignment, 'alignment_index': Path(str(alignment) + index_suffix),
            'reference': reference, 'reference_index': Path(str(reference) + '.fai'),
            'panel_summary': panel / 'summary.json', 'bcftools': bcftools, 'samtools': samtools,
            'runner': Path(__file__), 'marker_filter': Path(__file__).with_name('panel_marker_evidence.py'),
            'segmented_roh': Path(__file__).with_name('segmented_roh.py')}.items()}
        for name in ('reference', 'panel_summary', 'bcftools', 'samtools'):
            if inputs[name]['sha256'] != config[name + '_sha256']:
                raise ValueError(name + ' differs from pinned configuration')
        completion = verify_completion(config, alignment, inputs['alignment'], completion_state)
        if config['scope'] != 'synthetic': inputs['alignment_completion'] = completion
        summary = json.loads((panel / 'summary.json').read_text())
        for name in PANEL_FILES:
            identity = fingerprint(panel / name)
            if (identity['bytes'] != summary['outputs'][name]['bytes']
                    or identity['sha256'] != summary['outputs'][name]['sha256']):
                raise ValueError('panel output differs from preparation receipt: ' + name)
            inputs['panel_' + name] = identity
        versions = {}
        for name, executable in [('bcftools', bcftools), ('samtools', samtools)]:
            versions[name] = subprocess.check_output([str(executable), '--version'], text=True)
            if versions[name].splitlines()[0] != name + ' 1.24':
                raise ValueError('unreviewed native executable version')
        rows, panel_header = load_panel(panel, reference, config, summary)
        with pysam.FastaFile(str(reference)) as fasta, pysam.AlignmentFile(str(alignment),
                reference_filename=str(reference)) as reader:
            if list(zip(fasta.references, fasta.lengths)) != list(zip(reader.references, reader.lengths)):
                raise ValueError('alignment/reference dictionaries differ')
            if not reader.has_index() or reader.header.to_dict().get('HD', {}).get('SO') != 'coordinate':
                raise ValueError('indexed coordinate-sorted alignment required')
            groups = reader.header.to_dict().get('RG', [])
            if (not groups or any(g.get('SM') != config['sample'] or not g.get('ID') for g in groups)
                    or len({g['ID'] for g in groups}) != len(groups)):
                raise ValueError('alignment sample or read-group metadata mismatch')
            if config['scope'] == PATIENT_SCOPE and len(groups) != 4:
                raise ValueError('four sequencing-lane read groups required')
            read_groups = {g['ID'] for g in groups}
        def command(stage, args):
            state.update(stage=stage); save(output / 'status.json', state)
            log = output / ('step-%02d.log' % len(commands))
            tick = time.monotonic()
            with log.open('w') as stream:
                result = subprocess.run(list(map(str, args)), stdout=stream, stderr=subprocess.STDOUT)
            commands.append({'stage': stage, 'command': list(map(str, args)), 'exit_code': result.returncode,
                             'elapsed_seconds': time.monotonic() - tick, 'log': fingerprint(log)})
            save(output / 'commands.json', commands)
            if result.returncode:
                raise RuntimeError('native stage failed: ' + stage)
        subset = output / 'marker-reads.bam'
        command('extract eligible marker-overlapping reads', [samtools, 'view', '-b', '-M', '-L', panel / 'sites.bed',
            '-T', reference, '-q', config['minimum_mapping_quality'], '-F', '0xF04', '-e', 'mapq != 255',
            '--save-counts', output / 'extraction-counts.json', '-o', subset, alignment])
        command('index marker read subset', [samtools, 'index', subset])
        read_counts = Counter(); rg_counts = Counter()
        with pysam.AlignmentFile(str(subset)) as reader:
            for read in reader:
                if (read.flag & config['excluded_flags'] or read.mapping_quality == 255
                        or read.mapping_quality < config['minimum_mapping_quality']):
                    raise ValueError('extracted read violates quality/flag policy')
                if not read.has_tag('RG') or read.get_tag('RG') not in read_groups:
                    raise ValueError('eligible marker read has missing or unknown read group')
                read_counts[read.reference_name] += 1; rg_counts[read.get_tag('RG')] += 1
        extraction = json.loads((output / 'extraction-counts.json').read_text())
        required = ('records_processed', 'records_filter_accepted', 'records_filter_rejected')
        if (any(type(extraction.get(k)) is not int or extraction[k] < 0 for k in required)
                or extraction['records_processed'] != extraction['records_filter_accepted'] + extraction['records_filter_rejected']
                or extraction['records_filter_accepted'] != sum(read_counts.values())):
            raise ValueError('native extraction/read count reconciliation failed')
        command('original marker pileup', [bcftools, 'mpileup', '-Ob', '-o', output / 'pileup.bcf',
            '-f', reference, '-R', panel / 'sites.bed', '-a', 'AD,DP,SP',
            '-q', config['minimum_mapping_quality'], '-Q', config['minimum_base_quality'],
            '-d', config['pileup_maximum_depth'], '--ff', '0xF04', '--skip-indels', subset])
        command('explicit diploid panel genotypes', [bcftools, 'call', '-m', '--keep-alts', '--insert-missed',
            '--constrain', 'alleles', '--targets-file', panel / 'alleles.tsv.gz', '--ploidy', '2',
            '--annotate', 'GQ', '-Ob', '-o', output / 'calls.bcf', output / 'pileup.bcf'])
        state.update(stage='population annotation and original-evidence qualification'); save(output / 'status.json', state)
        qualification = qualify_calls(rows, panel_header, output / 'pileup.bcf', output / 'calls.bcf',
                                      config['sample'], config['evidence'], output)
        models = []
        scenarios = [(tag, config['primary_maximum_gap_bp']) for tag in config['population_profiles']]
        scenarios += [('POP_AF', gap) for gap in config['pooled_gap_sensitivity_bp'] if gap != config['primary_maximum_gap_bp']]
        if qualification['counts'].get('qualified_markers', 0):
            for tag, gap in scenarios:
                name = tag + '-gap-' + str(gap)
                state.update(stage='independent homozygosity models', current_model=name); save(output / 'status.json', state)
                folder = output / name
                result = run_segmented(output / 'qualified-markers.vcf.gz', config['sample'], bcftools,
                                       folder, gap, config['minimum_segment_markers'], tag)
                models.append({'population_af_tag': tag, 'maximum_gap_bp': gap, 'status': result['status'],
                    'counts': result['counts'], 'segments': len(result['segments']),
                    'gap_boundaries': len(result['gap_boundaries']),
                    'regions': [r for segment in result['segments'] for r in segment['regions']],
                    'summary': fingerprint(folder / 'summary.json')})
        for identity in inputs.values():
            st = Path(identity['path']).stat()
            if (st.st_size, st.st_mtime_ns) != (identity['bytes'], identity['mtime_ns']):
                raise ValueError('input changed during genotyping or model analysis')
        output_names = ['marker-reads.bam', 'marker-reads.bam.bai', 'extraction-counts.json', 'pileup.bcf', 'calls.bcf',
                        'annotated-all-markers.vcf.gz', 'annotated-all-markers.vcf.gz.tbi', 'qualified-markers.vcf.gz',
                        'qualified-markers.vcf.gz.tbi', 'marker-evidence.jsonl.gz', 'commands.json']
        result = {'started_at': state['started_at'], 'completed_at': stamp(), 'elapsed_seconds': time.monotonic() - started,
            'config': config, 'inputs': inputs, 'alignment_completion': completion, 'runtime': versions,
            'commands': commands, 'marker_preparation_receipt_validated': True,
            'read_extraction': {'native_counts': extraction, 'accepted_by_chromosome': dict(read_counts),
                'accepted_by_read_group': dict(rg_counts), 'scope': 'Indexed marker-overlap extraction only; not all genome reads or independent molecules.'},
            'marker_qualification': qualification, 'homozygosity_models': models,
            'analysis_status': 'completed_requested_models' if models else 'no_qualified_markers_no_models',
            'outputs': {name: fingerprint(output / name) for name in output_names},
            'limitations': ['Mapped-read-to-marker workflow; clinical sensitivity and alignment truth are not established by synthetic controls.',
                'Public accessibility and selected-marker retention are not base-level patient callability.',
                'Read quality, depth, third-allele and model-gap thresholds are research choices requiring sensitivity review.',
                'Population-frequency profiles do not assign ancestry; spacing does not establish linkage independence.',
                'No parental origin, UPD exclusion, copy neutrality, clinical classification, EPCR or drug response is inferred.']}
        save(output / 'summary.json', result)
        state.update(status='complete', stage='complete', completed_at=result['completed_at'],
                     summary_sha256=fingerprint(output / 'summary.json')['sha256'])
        save(output / 'status.json', state)
        return result
    except BaseException as error:
        state.update(status='failed', failed_at=stamp(), error_type=type(error).__name__)
        save(output / 'status.json', state)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('alignment', 'reference', 'panel', 'config', 'bcftools', 'samtools', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--completion-state', type=Path)
    args = parser.parse_args(); os.nice(10)
    config = json.loads(args.config.read_text())
    result = run_analysis(args.alignment, args.reference, args.panel, config, args.bcftools, args.samtools,
                          args.output, args.completion_state)
    print(json.dumps({'completed_at': result['completed_at'], 'elapsed_seconds': result['elapsed_seconds'],
                      'marker_counts': result['marker_qualification']['counts'],
                      'model_comparisons': len(result['homozygosity_models']),
                      'analysis_status': result['analysis_status']}, indent=2))


if __name__ == '__main__':
    main()
