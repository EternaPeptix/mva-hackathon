#!/usr/bin/env python3
"""Retained synthetic mapped-read experiments for the integrated marker workflow."""
import argparse
from collections import Counter
import json
import os
from pathlib import Path
import random
import time

import pysam

from genotype_roh_panel import PANEL_FILES, fingerprint, run_analysis, save, stamp
from segmented_roh import POPULATION_TAGS


def create_fixture(root, mode, bcftools, samtools):
    root.mkdir(parents=True, exist_ok=False)
    rows = []
    if mode in ('controls', 'empty'):
        lengths = {'1': 15000}
        cases = [('homref', {'A': 30}, 40, 60, 0), ('heterozygous', {'A': 15, 'G': 15}, 40, 60, 0),
                 ('homalt', {'G': 30}, 40, 60, 0), ('no_coverage', {}, 40, 60, 0),
                 ('low_depth', {'A': 3}, 40, 60, 0), ('third_allele_mixture', {'A': 20, 'C': 10}, 40, 60, 0),
                 ('third_allele_only', {'C': 30}, 40, 60, 0), ('low_base_quality', {'A': 30}, 10, 60, 0),
                 ('duplicate_only', {'A': 30}, 40, 60, 0x400), ('unknown_mapping_quality', {'A': 30}, 40, 255, 0),
                 ('low_mapping_quality', {'A': 30}, 40, 10, 0), ('secondary_only', {'A': 30}, 40, 60, 0x100),
                 ('supplementary_only', {'A': 30}, 40, 60, 0x800), ('QC_failed_only', {'A': 30}, 40, 60, 0x200)]
        for i, (label, bases, bq, mq, flag) in enumerate(cases):
            rows.append({'chrom': '1', 'pos': 1001 + i * 1000, 'label': label,
                         'bases': {} if mode == 'empty' else bases, 'BQ': bq, 'MQ': mq, 'flag': flag,
                         'AF': 0.3, 'GT': {'homref': (0, 0), 'heterozygous': (0, 1), 'homalt': (1, 1)}.get(label)})
    elif mode == 'planted':
        lengths = {'1': 20000000, '2': 20000000}
        rng = random.Random(20260919)
        for chrom in lengths:
            for pos in range(5001, lengths[chrom], 10000):
                af = rng.choice([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
                if chrom == '1' and 5000000 <= pos < 15000000:
                    a = int(rng.random() < af); gt = (a, a)
                else:
                    gt = tuple(sorted((int(rng.random() < af), int(rng.random() < af))))
                gap = chrom == '1' and 8000000 <= pos < 12000000
                counts = [{'A': 30}, {'A': 15, 'G': 15}, {'G': 30}][sum(gt)]
                rows.append({'chrom': chrom, 'pos': pos, 'label': 'missing_gap' if gap else 'known_genotype',
                             'bases': {} if gap else counts, 'BQ': 40, 'MQ': 60, 'flag': 0, 'AF': af, 'GT': gt})
    else:
        raise ValueError('unknown synthetic experiment')
    reference = root / 'synthetic.fa'
    with reference.open('w') as stream:
        for chrom, length in lengths.items():
            stream.write('>' + chrom + '\n' + 'A' * length + '\n')
    pysam.faidx(str(reference))
    header = {'HD': {'VN': '1.6', 'SO': 'coordinate'}, 'SQ': [{'SN': c, 'LN': n} for c, n in lengths.items()],
              'RG': [{'ID': 'SYNTH.1', 'SM': 'SYNTH', 'PL': 'ILLUMINA'}]}
    alignment = root / 'synthetic.bam'; count = 0
    with pysam.AlignmentFile(str(alignment), 'wb', header=header) as writer:
        for row in rows:
            for base, depth in row['bases'].items():
                for i in range(depth):
                    read = pysam.AlignedSegment(writer.header)
                    read.query_name = f'synthetic_{count}'
                    sequence = list('A' * 150); sequence[75] = base
                    read.query_sequence = ''.join(sequence); read.query_qualities = [row['BQ']] * 150
                    read.reference_id = list(lengths).index(row['chrom']); read.reference_start = row['pos'] - 76
                    read.mapping_quality = row['MQ']; read.cigarstring = '150M'; read.flag = row['flag'] | (16 if i % 2 else 0)
                    read.set_tag('RG', 'SYNTH.1'); writer.write(read); count += 1
    pysam.index(str(alignment))
    panel = root / 'panel'; panel.mkdir()
    header = pysam.VariantHeader()
    for chrom, length in lengths.items(): header.contigs.add(chrom, length=length)
    for tag in POPULATION_TAGS: header.info.add(tag, 1, 'Float', 'Known synthetic independent frequency ' + tag)
    header.info.add('POP_AN', 1, 'Integer', 'Synthetic reference-population allele count')
    header.info.add('POP_NS', 1, 'Integer', 'Synthetic reference-population sample count')
    with pysam.VariantFile(str(panel / 'panel.vcf.gz'), 'wz', header=header) as writer:
        for data in rows:
            r = writer.new_record(contig=data['chrom'], start=data['pos'] - 1, alleles=('A', 'G'))
            r.filter.add('PASS')
            for i, tag in enumerate(POPULATION_TAGS): r.info[tag] = data['AF'] + (i * 0.01 if mode != 'planted' else 0)
            r.info['POP_AN'] = 5000; r.info['POP_NS'] = 2500; writer.write(r)
    pysam.tabix_index(str(panel / 'panel.vcf.gz'), preset='vcf')
    (panel / 'sites.bed').write_text(''.join(f"{r['chrom']}\t{r['pos']-1}\t{r['pos']}\n" for r in rows))
    with pysam.BGZFile(str(panel / 'alleles.tsv.gz'), 'wb') as stream:
        stream.write(''.join(f"{r['chrom']}\t{r['pos']}\tA,G\n" for r in rows).encode())
    pysam.tabix_index(str(panel / 'alleles.tsv.gz'), seq_col=0, start_col=1, end_col=1)
    panel_summary = {'config': {'reference_sha256': fingerprint(reference)['sha256'], 'chromosomes': list(lengths)},
        'selected_markers': len(rows), 'selected_markers_by_chromosome': dict(Counter(r['chrom'] for r in rows)),
        'outputs': {n: fingerprint(panel / n) for n in PANEL_FILES}, 'scope': 'Synthetic known-frequency fixture only'}
    save(panel / 'summary.json', panel_summary)
    config = json.loads(Path(__file__).parents[1].joinpath('configs/roh-genotyping.json').read_text())
    config.update(scope='synthetic', sample='SYNTH', chromosomes=list(lengths),
                  reference_sha256=fingerprint(reference)['sha256'],
                  panel_summary_sha256=fingerprint(panel / 'summary.json')['sha256'],
                  bcftools_sha256=fingerprint(bcftools)['sha256'], samtools_sha256=fingerprint(samtools)['sha256'],
                  minimum_segment_markers=20 if mode == 'planted' else 2)
    save(root / 'config.json', config)
    save(root / 'truth.json', {'mode': mode, 'read_records_written': count, 'markers': rows,
        'scope': 'Constructed mapped reads, independent markers and known genotype truth; no alignment/duplicate-marking accuracy claim.'})
    return rows, alignment, reference, panel, config


def check_result(rows, mode, result, output):
    with pysam.VariantFile(str(output / 'annotated-all-markers.vcf.gz')) as source:
        observed = list(source)
    if len(observed) != len(rows): raise AssertionError('synthetic marker accounting failed')
    usable = 0
    for truth, record in zip(rows, observed):
        if (record.contig, record.pos) != (truth['chrom'], truth['pos']):
            raise AssertionError('synthetic marker identity changed')
        expected_usable = bool(truth['bases']) and (mode == 'planted' or truth['label'] in ('homref', 'heterozygous', 'homalt'))
        if bool(record.info.get('MVA_USABLE', False)) != expected_usable:
            raise AssertionError('synthetic marker qualification error: ' + truth['label'])
        if expected_usable and record.samples['SYNTH']['GT'] != truth['GT']:
            raise AssertionError('known synthetic genotype not recovered')
        if not truth['bases'] and any(a is not None for a in record.samples['SYNTH']['GT']):
            raise AssertionError('unobserved site converted to genotype')
        usable += expected_usable
    if result['marker_qualification']['counts'].get('qualified_markers', 0) != usable:
        raise AssertionError('qualified marker reconciliation failed')
    if mode == 'empty':
        if result['homozygosity_models'] or result['analysis_status'] != 'no_qualified_markers_no_models':
            raise AssertionError('empty evidence treated as an analyzed negative')
    else:
        if len(result['homozygosity_models']) != 8: raise AssertionError('required sensitivity scenario missing')
        for model in result['homozygosity_models']:
            if model['counts']['input_markers'] != usable: raise AssertionError('model marker count changed')
            if mode == 'planted':
                regions = model['regions']
                if any(r['chrom'] == '2' for r in regions): raise AssertionError('false control chromosome region')
                if any(r['start_1based'] <= 8000000 and r['end_inclusive'] >= 12000000 for r in regions):
                    raise AssertionError('model bridged unobserved gap')
                if len([r for r in regions if r['chrom'] == '1' and r['length_bp'] >= 2500000]) != 2:
                    raise AssertionError('flanking planted regions not recovered')
    return {'requested_markers': len(rows), 'qualified_markers': usable,
            'all_known_qualified_genotypes_correct': True, 'unobserved_sites_remain_missing': True,
            'expected_exclusions_verified': True, 'model_comparisons': len(result['homozygosity_models']),
            'planted_regions_and_control_checked': mode == 'planted',
            'empty_evidence_explicitly_unanalyzed': mode == 'empty'}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mode', choices=('controls', 'empty', 'planted'), required=True)
    for name in ('bcftools', 'samtools', 'output'): p.add_argument('--' + name, type=Path, required=True)
    a = p.parse_args(); a.output.mkdir(parents=True, exist_ok=False); os.nice(10)
    started = time.monotonic()
    rows, alignment, reference, panel, config = create_fixture(a.output / 'fixture', a.mode, a.bcftools, a.samtools)
    result = run_analysis(alignment, reference, panel, config, a.bcftools, a.samtools, a.output / 'analysis')
    checks = check_result(rows, a.mode, result, a.output / 'analysis')
    validation = {'completed_at': stamp(), 'elapsed_seconds_including_fixture': time.monotonic() - started,
        'mode': a.mode, 'checks': checks, 'probe_source': fingerprint(Path(__file__)),
        'pipeline_summary': fingerprint(a.output / 'analysis/summary.json'),
        'scope': 'Synthetic mapped-read-to-marker-to-model validation. No patient genotype, clinical sensitivity, ancestry, copy neutrality or parental origin is established.'}
    save(a.output / 'validation.json', validation)
    print(json.dumps(validation, indent=2))


if __name__ == '__main__': main()
