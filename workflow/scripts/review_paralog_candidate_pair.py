#!/usr/bin/env python3
"""Review a nearby pair without confusing local phase with genomic-locus assignment.

Retain raw regional calls and reads, count fragment-level phase at several mapping
thresholds, and remap the observed primary reads/mates individually to the complete
reference. BWA's reported alternatives are a sensitivity analysis, not an exhaustive
homology search or independent evidence that a variant belongs to a particular locus.
Patient sequences and detailed output must stay outside the source repository.
"""
import argparse
from collections import Counter, defaultdict
import datetime
import json
import os
from pathlib import Path
import subprocess
import time

import pysam

from audit_candidate_phase import haplotypes, spanning_sequence
from review_cnv_context import fingerprint, require, save, sha, verify


LABELS = ('0|0', '1|0', '0|1', '1|1', 'other', 'ambiguous', 'mate_conflict')


def merge_intervals(intervals):
    """Merge overlapping/adjacent zero-based, half-open intervals on one contig."""
    out = []
    for start, end in sorted(intervals):
        require(0 <= start < end, 'Invalid interval')
        if out and start <= out[-1][1]:
            out[-1][1] = max(end, out[-1][1])
        else:
            out.append([start, end])
    return out


def phase_observation(read, start, end, expected, minimum_base_quality):
    sequence = spanning_sequence(read, start, end, minimum_base_quality)
    labels = expected.get(sequence, []) if sequence is not None else []
    return None if sequence is None else labels[0] if len(labels) == 1 else 'ambiguous' if labels else 'other'


def summarize_phase(rows, minimum_mapping_quality):
    groups = defaultdict(list)
    excluded = Counter()
    for row in rows:
        reasons = []
        if row['flag'] & 0xF04:
            reasons.append('excluded_flag')
        if row['mapping_quality'] == 255 or row['mapping_quality'] < minimum_mapping_quality:
            reasons.append('mapping_quality')
        if row['phase_label'] is None:
            reasons.append('not_spanning_or_low_base_quality_or_reference_skip')
        if reasons:
            excluded.update(reasons)
        else:
            groups[(row['RG'], row['query_name'])].append(row)
    fragments = []
    for (rg, name), observations in sorted(groups.items()):
        labels = {r['phase_label'] for r in observations}
        label = next(iter(labels)) if len(labels) == 1 else 'mate_conflict'
        strands = {r['reverse'] for r in observations}
        fragments.append({'RG': rg, 'query_name': name, 'label': label,
                          'strand': 'both' if len(strands) == 2 else 'reverse' if True in strands else 'forward',
                          'observation_ids': [r['id'] for r in observations]})
    counts = Counter(f['label'] for f in fragments)
    details = {}
    for label in LABELS:
        fs = [f for f in fragments if f['label'] == label]
        ids = {i for f in fs for i in f['observation_ids']}
        observed = [r for r in rows if r['id'] in ids]
        details[label] = {'fragments': len(fs), 'by_read_group': dict(Counter(f['RG'] for f in fs)),
                          'strand': dict(Counter(f['strand'] for f in fs)),
                          'spanning_observation_MQ': dict(Counter(r['mapping_quality'] for r in observed)),
                          'spanning_observation_NM': dict(Counter(r['tags'].get('NM') for r in observed)),
                          'spanning_observation_CIGAR': dict(Counter(r['cigar'] for r in observed))}
    return {'minimum_mapping_quality': minimum_mapping_quality,
            'counts': {label: counts[label] for label in LABELS}, 'details': details,
            'eligible_fragments': len(fragments), 'excluded_record_reasons_nonexclusive': dict(excluded),
            'fragments': fragments}


def record_details(read, identifier):
    require(read.has_tag('RG'), 'Missing read group')
    return {'id': identifier, 'RG': read.get_tag('RG'), 'query_name': read.query_name,
            'read_number': 1 if read.is_read1 else 2 if read.is_read2 else 0,
            'flag': read.flag, 'chrom': read.reference_name, 'start_0based': read.reference_start,
            'end_exclusive': read.reference_end, 'mapping_quality': read.mapping_quality,
            'cigar': read.cigarstring, 'reverse': read.is_reverse,
            'mate_chrom': read.next_reference_name, 'mate_start_0based': read.next_reference_start,
            'template_length': read.template_length, 'tags': dict(read.tags)}


def remap_details(read):
    return {'chrom': read.reference_name, 'start_0based': read.reference_start,
            'end_exclusive': read.reference_end, 'flag': read.flag, 'mapping_quality': read.mapping_quality,
            'cigar': read.cigarstring, 'query_alignment_length': read.query_alignment_length,
            'tags': dict(read.tags)}


def within_locus(row, locus):
    return (row['chrom'] == locus['chrom'] and row['end_exclusive'] is not None
            and row['start_0based'] < locus['end_exclusive']
            and row['end_exclusive'] > locus['start_0based'])


def compare_reported_mappings(rows, locus):
    """Keep ties and missing alternatives; no implied infinite gap or unique mapping."""
    scored = [r for r in rows if not r['flag'] & (4 | 2048) and 'AS' in r['tags']]
    inside = [r for r in scored if within_locus(r, locus)]
    outside = [r for r in scored if not within_locus(r, locus)]
    gene_score = max((r['tags']['AS'] for r in inside), default=None)
    other_score = max((r['tags']['AS'] for r in outside), default=None)
    primary = [r for r in rows if not r['flag'] & (256 | 2048)]
    require(len(primary) == 1, 'Expected exactly one primary remap per query')
    gap = gene_score - other_score if gene_score is not None and other_score is not None else None
    return {'best_reported_gene_AS': gene_score, 'best_reported_other_AS': other_score,
            'gene_minus_other_AS': gap, 'primary_MQ': primary[0]['mapping_quality'],
            'primary_overlaps_gene': within_locus(primary[0], locus),
            'best_reported_other_locations': [r for r in outside if r['tags']['AS'] == other_score]}


def variant_details(record):
    return {'chrom': record.contig, 'pos': record.pos, 'end': record.stop, 'id': record.id,
            'ref': record.ref, 'alts': record.alts, 'qual': record.qual,
            'filter': list(record.filter), 'info': dict(record.info),
            'samples': {s: dict(record.samples[s]) | {'GT_phased': record.samples[s].phased} for s in record.samples}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('root', 'config', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    output = args.output.resolve()
    require(not output.exists() and not output.is_relative_to(Path(__file__).resolve().parents[2]),
            'Use a new output directory outside the source tree')
    config = json.loads(args.config.read_text())
    output.mkdir(parents=True)
    tick = time.monotonic()
    os.nice(10)
    state = {'status': 'running', 'pid': os.getpid(), 'started_at': datetime.datetime.now(datetime.timezone.utc).isoformat()}
    save(output / 'status.json', state)
    try:
        native_dir = root / config['structural_directory']
        require(sha(native_dir / 'summary.json') == config['structural_summary_sha256'], 'Structural summary differs')
        native = json.loads((native_dir / 'summary.json').read_text())
        require(native['analysis_status'] == 'complete_with_calls', 'Structural input incomplete')
        verify(native_dir / 'provenance.json', native['provenance'])
        proof = json.loads((native_dir / 'provenance.json').read_text())
        monitored = {}
        inputs = {}
        for name in ('alignment', 'alignment_index', 'reference', 'reference_index', 'samtools'):
            info = proof['inputs'][name]
            path = Path(info['path'])
            monitored[str(path)] = (path.stat().st_size, path.stat().st_mtime_ns)
            verify(path, info)
            inputs[name] = info
        source = Path(inputs['alignment']['path'])
        reference = Path(inputs['reference']['path'])
        samtools = Path(inputs['samtools']['path'])
        prefix = root / config['bwa_reference']
        index_info = root / config['bwa_index_record']
        index_proof = json.loads(index_info.read_text())
        require(sha(prefix) == index_proof['reference_sha256'] == inputs['reference']['sha256'], 'BWA reference differs')
        inputs['bwa_reference'] = fingerprint(prefix)
        monitored[str(prefix)] = (prefix.stat().st_size, prefix.stat().st_mtime_ns)
        inputs['bwa_index_record'] = fingerprint(index_info)
        inputs['bwa_indices'] = {}
        for suffix, expected in index_proof['index_sha256'].items():
            path = Path(str(prefix) + suffix)
            require(sha(path) == expected, 'BWA index differs: ' + suffix)
            inputs['bwa_indices'][suffix] = fingerprint(path)
            monitored[str(path)] = (path.stat().st_size, path.stat().st_mtime_ns)
        bwa = root / config['bwa_binary']
        inputs['bwa'] = fingerprint(bwa)
        for name in ('calls.bcf', 'calls.bcf.csi'):
            verify(native_dir / name, native['outputs'][name])
        with pysam.VariantFile(str(native_dir / 'calls.bcf')) as reader:
            structural = [r.copy() for r in reader if r.id in config['structural_record_ids']]
        require(sorted(r.id for r in structural) == sorted(config['structural_record_ids']), 'Structural IDs missing or duplicated')
        chrom = config['chrom']
        require(all(r.contig == chrom for r in structural), 'Other-contig structural record')
        save(output / 'structural-records.json', [variant_details(r) for r in structural])

        shard = root / config['small_variant_shard']
        require(sha(shard / 'complete.json') == config['small_variant_completion_sha256'], 'Small-variant completion differs')
        completion = json.loads((shard / 'complete.json').read_text())
        inputs['small_variant_completion'] = fingerprint(shard / 'complete.json')
        region = config['variant_window_0based']
        variants = {}
        inputs['small_variant_calls'] = {}
        for caller in ('gatk', 'bcftools'):
            path = shard / completion[caller]['path']
            require(sha(path) == completion[caller]['sha256'], 'Small-variant call file differs')
            inputs['small_variant_calls'][caller] = fingerprint(path)
            inputs['small_variant_calls'][caller + '_index'] = fingerprint(Path(str(path) + '.tbi'))
            with pysam.VariantFile(str(path)) as reader:
                variants[caller] = [variant_details(r) for r in reader.fetch(chrom, *region)]
        save(output / 'regional-variants.json', variants)
        sites = [tuple(s) for s in config['sites']]
        require(len(sites) == 2 and sites[0][0] + len(sites[0][1]) <= sites[1][0], 'Expected ordered nonoverlapping pair')
        start = sites[0][0] - 2
        end = sites[-1][0] - 1 + len(sites[-1][1]) + 2
        centers = sorted({sites[0][0] + 1} | {p for r in structural for p in (r.pos, r.stop)})
        radius = config['endpoint_read_window_radius']
        intervals = merge_intervals([[max(0, p-radius-1), p+radius] for p in centers])
        bed = output / 'regions.bed'
        bed.write_text(''.join(f'{chrom}\t{a}\t{b}\n' for a,b in intervals))
        bam_path = output / 'regional-reads.bam'
        commands = []
        def run(command, stdout=None):
            number = len(commands)
            commands.append(command)
            save(output / 'commands.json', commands)
            with (output / f'command-{number}.log').open('w') as log:
                subprocess.run(command, stdout=stdout if stdout is not None else log, stderr=log, check=True)
        run([str(samtools), 'view', '-@', '1', '--no-PG', '-T', str(reference), '-b', '-M', '-L', str(bed), '-o', str(bam_path), str(source)])
        run([str(samtools), 'quickcheck', '-v', str(bam_path)])
        run([str(samtools), 'index', str(bam_path)])
        queries = {}
        reference_queries = []
        with pysam.FastaFile(str(reference)) as fa:
            expected = haplotypes(fa.fetch(chrom, start, end).upper(), start, sites)
            for p in centers:
                for n in config['reference_query_radii']:
                    require(p-n-1 >= 0 and p+n <= fa.get_reference_length(chrom), 'Reference window out of bounds')
                    identifier = f'reference_{p}_{n}'
                    queries[identifier] = fa.fetch(chrom, p-n-1, p+n)
                    reference_queries.append({'id': identifier, 'chrom': chrom, 'start_0based': p-n-1, 'end_exclusive': p+n})
        observations = []
        observed_mates = []
        with pysam.AlignmentFile(str(bam_path), 'rb') as bam, pysam.AlignmentFile(str(source), 'rc', reference_filename=str(reference)) as mates:
            header = bam.header.to_dict()
            declared = {r['ID'] for r in header.get('RG', [])}
            for i, read in enumerate(bam.fetch(chrom, start, end)):
                row = record_details(read, 'read_' + str(i))
                require(row['RG'] in declared, 'Undeclared read group')
                row['phase_label'] = phase_observation(read, start, end, expected, config['minimum_base_quality'])
                row['single_end_remap_query'] = None
                if not read.is_secondary and not read.is_supplementary:
                    require(read.query_sequence is not None, 'Primary read lacks sequence')
                    queries[row['id']] = read.get_forward_sequence()
                    row['single_end_remap_query'] = row['id']
                    if not read.is_paired or read.mate_is_unmapped:
                        row['mate_status'] = 'unpaired_or_mate_unmapped'
                    else:
                        try:
                            mate = mates.mate(read)
                        except ValueError as error:
                            row['mate_status'] = 'lookup_failed'
                            row['mate_error'] = str(error)
                        else:
                            require(mate.query_name == read.query_name and mate.has_tag('RG')
                                    and mate.get_tag('RG') == read.get_tag('RG')
                                    and mate.is_read1 != read.is_read1 and mate.is_read2 != read.is_read2
                                    and not mate.is_secondary and not mate.is_supplementary,
                                    'Retrieved mate identity is inconsistent')
                            require(mate.reference_id == read.next_reference_id and mate.reference_start == read.next_reference_start,
                                    'Retrieved mate coordinate is inconsistent')
                            identifier = row['id'] + '_mate'
                            queries[identifier] = mate.get_forward_sequence()
                            mate_row = record_details(mate, identifier)
                            mate_row['source_observation_id'] = row['id']
                            observed_mates.append(mate_row)
                            row['mate_status'] = 'retrieved_and_identity_checked'
                            row['mate_id'] = identifier
                observations.append(row)
        phase = {str(mq): summarize_phase(observations, mq) for mq in config['mapping_quality_thresholds']}
        fasta_path = output / 'remap-queries.fa'
        require(all(s and set(s.upper()) <= set('ACGTN') for s in queries.values()), 'Invalid query sequence')
        fasta_path.write_text(''.join(f'>{name}\n{sequence}\n' for name,sequence in queries.items()))
        sam_path = output / 'single-end-remap.sam'
        with sam_path.open('w') as sam:
            run([str(bwa), 'mem', '-a', '-Y', '-t', '2', str(prefix), str(fasta_path)], stdout=sam)
        remaps = defaultdict(list)
        with pysam.AlignmentFile(str(sam_path), 'r') as sam:
            for read in sam:
                require(read.query_name in queries, 'Unknown remap query')
                remaps[read.query_name].append(remap_details(read))
        require(set(remaps) == set(queries), 'Remap query census differs')
        comparison = {name: compare_reported_mappings(rows, config['gene_locus']) for name,rows in remaps.items()}
        base_phase = phase[str(config['primary_mapping_quality'])]
        remap_summary = {}
        for label in LABELS:
            fragments = [f for f in base_phase['fragments'] if f['label'] == label]
            ids = {i for f in fragments for i in f['observation_ids']}
            found = [comparison[i] for i in ids if i in comparison]
            remap_summary[label] = {
                'qualified_original_fragments': len(fragments), 'qualified_spanning_observations': len(ids),
                'remapped_spanning_observations': len(found),
                'single_end_primary_MQ': dict(Counter(r['primary_MQ'] for r in found)),
                'single_end_gene_minus_other_AS': dict(Counter('unreported_comparison' if r['gene_minus_other_AS'] is None else str(r['gene_minus_other_AS']) for r in found)),
                'fragments_with_any_spanning_single_end_MQ_at_least_20': sum(any(comparison[i]['primary_MQ'] != 255 and comparison[i]['primary_MQ'] >= 20 for i in f['observation_ids'] if i in comparison) for f in fragments),
                'fragments_with_any_spanning_primary_outside_gene': sum(any(not comparison[i]['primary_overlaps_gene'] for i in f['observation_ids'] if i in comparison) for f in fragments),
            }
        save(output / 'read-review.json', {'header': header, 'phase_window_0based': [start, end],
             'expected_haplotypes': dict(expected), 'observations': observations,
             'observed_mates': observed_mates, 'phase_thresholds': phase})
        save(output / 'single-end-remap.json', {'reference_queries': reference_queries,
             'reported_alignments': dict(remaps), 'comparisons': comparison, 'summary': remap_summary})
        require(monitored == {p:(Path(p).stat().st_size,Path(p).stat().st_mtime_ns) for p in monitored}, 'Input changed during review')
        summary = {'completed_at': datetime.datetime.now(datetime.timezone.utc).isoformat(),
                   'elapsed_seconds': time.monotonic()-tick, 'configuration': config,
                   'inputs': inputs | {'native_summary': fingerprint(native_dir/'summary.json'),
                                      'native_provenance': fingerprint(native_dir/'provenance.json'),
                                      'script': fingerprint(Path(__file__)), 'config': fingerprint(args.config),
                                      'phase_helper': fingerprint(Path(__file__).with_name('audit_candidate_phase.py')),
                                      'io_helper': fingerprint(Path(__file__).with_name('review_cnv_context.py'))},
                   'pysam_version': pysam.__version__, 'whole_alignment_reference_and_BWA_indices_SHA256_rechecked': True,
                   'regional_variant_records': {k:len(v) for k,v in variants.items()},
                   'structural_records': [variant_details(r) for r in structural],
                   'phase': {mq:{k:v for k,v in value.items() if k != 'fragments'} for mq,value in phase.items()},
                   'candidate_overlapping_alignment_records': len(observations),
                   'mate_statuses': dict(Counter(r.get('mate_status','secondary_or_supplementary_not_queried') for r in observations)),
                   'single_end_remap_summary': remap_summary,
                   'reference_remap_queries': [r | {'reported_alignments': remaps[r['id']]} for r in reference_queries],
                   'outputs': {p.name:fingerprint(p) for p in output.iterdir() if p.is_file() and p.name not in ('summary.json','status.json')},
                   'limitations': [
                       'Post hoc review of one pair and selected structural calls; not a complete FANCD2 gene, other-variant, or diagnostic exclusion.',
                       'Phase requires one read spanning both sites, BQ at least the configured threshold, existing duplicate marks, flags 0xF04 excluded, MQ255 excluded; no BAQ.',
                       'Single-end BWA remapping reuses the same reads and reference and removes paired-end information; score/MQ changes do not prove a mapping artifact or correct origin.',
                       'BWA -a retains aligner-reported alternatives; this is not an exhaustive homology search. Absent alternatives and absent caller records are not proof of uniqueness or a reference genotype.',
                       'Read/mate queries may repeat the same molecule; only RG/QNAME fragment tallies are counted as phase evidence. No independent samples or parental phase.',
                       'Structural records and reference-window homology are descriptive context; no complete SV allele, breakpoint validation, copy-number change or clinical classification inferred.'
                   ]}
        save(output / 'summary.json', summary)
        state.update(status='complete', completed_at=summary['completed_at'], summary_sha256=sha(output/'summary.json'))
        save(output / 'status.json', state)
        print(json.dumps({'elapsed_seconds': summary['elapsed_seconds'], 'phase': {k:v['counts'] for k,v in phase.items()},
                          'single_end_remap_summary': remap_summary, 'mate_statuses': summary['mate_statuses']}))
    except BaseException as error:
        state.update(status='failed', error_type=type(error).__name__, error=str(error))
        save(output / 'status.json', state)
        raise


if __name__ == '__main__':
    main()
