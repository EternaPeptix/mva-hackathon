#!/usr/bin/env python3
"""Review selected indels in completed four-lane reads without changing source records.

Exact local sequence counts are descriptive, not a calibrated genotyper. All
original multi-allele caller genotypes are retained when normalizing identities.
"""
import argparse
from collections import Counter, defaultdict
import json
import os
from pathlib import Path
import re
import subprocess
import time

import pysam

from audit_baseline_reads import collapse_fragments, exact_calls
from audit_candidate_phase import spanning_sequence
from prepare_gatk_compatible_cram import digest, save, stamp


def window_hypotheses(fasta, target, left, right):
    if left < 1 or right < 1:
        raise ValueError('require aligned flanks on both sides')
    chrom, pos, ref = (target[k] for k in ('chrom', 'pos', 'ref'))
    start, end = pos - 1 - left, pos - 1 + len(ref) + right
    if start < 0 or end > fasta.get_reference_length(chrom):
        raise ValueError('window outside reference')
    sequence = fasta.fetch(chrom, start, end).upper()
    if sequence[left:left + len(ref)] != ref:
        raise ValueError('target reference mismatch')
    hypotheses = {'REF': sequence}
    for label, alt in [('ALT', target['alt'])] + list(target.get('other_hypotheses', {}).items()):
        if label in hypotheses or label in ('other', 'mate_conflict', 'ambiguous'):
            raise ValueError('duplicate or reserved hypothesis label')
        if not re.fullmatch('[ACGT]+', alt) or alt == ref:
            raise ValueError('require a distinct concrete alternate sequence')
        hypotheses[label] = sequence[:left] + alt + sequence[left + len(ref):]
    return start, end, hypotheses


def window_evidence(bam, reference, target, left, right, minimum_mapping_quality=20):
    with pysam.FastaFile(str(reference)) as fasta:
        start, end, hypotheses = window_hypotheses(fasta, target, left, right)
    lookup = defaultdict(list)
    for label, sequence in hypotheses.items():
        lookup[sequence].append(label)
    fragments, excluded, features = defaultdict(list), Counter(), defaultdict(Counter)
    by_lane = defaultdict(lambda: defaultdict(list))
    other_sequences = defaultdict(set)
    with pysam.AlignmentFile(str(bam), 'rb') as alignment:
        for read in alignment.fetch(target['chrom'], start, end):
            if read.flag & 0xF04:
                excluded['excluded_flags'] += 1
                continue
            if read.mapping_quality == 255 or read.mapping_quality < minimum_mapping_quality:
                excluded['low_or_unknown_mapping_quality'] += 1
                continue
            sequence = spanning_sequence(read, start, end, 20)
            if sequence is None:
                excluded['not_spanning_or_low_quality_or_reference_skip'] += 1
                continue
            labels = lookup.get(sequence, [])
            label = labels[0] if len(labels) == 1 else 'ambiguous' if labels else 'other'
            if not read.has_tag('RG'):
                raise ValueError('eligible record lacks read group')
            key = (read.get_tag('RG'), read.query_name)
            observation = (label, read.is_reverse)
            fragments[key].append(observation)
            by_lane[key[0]][key].append(observation)
            if label == 'other':
                other_sequences[sequence].add(key)
            features[label]['read_observations'] += 1
            features[label]['soft_clipped_observations'] += any(op == 4 for op, _ in read.cigartuples or [])
            features[label]['NM_at_least_5_observations'] += read.has_tag('NM') and read.get_tag('NM') >= 5
    labels = list(hypotheses) + ['ambiguous']
    result = collapse_fragments(fragments, labels)
    result.update(key=target['key'], window_0based_half_open=[target['chrom'], start, end],
                  flanks_bp=[left, right], hypotheses=hypotheses,
                  per_lane={lane: collapse_fragments(rows, labels) for lane, rows in by_lane.items()},
                  excluded_read_counts=dict(excluded), raw_read_features=dict(features),
                  other_sequence_length_fragment_counts=dict(Counter(
                      len(sequence) for sequence, keys in other_sequences.items() for _ in keys)),
                  filters={'minimum_mapping_quality': minimum_mapping_quality, 'MQ255_excluded': True,
                           'minimum_raw_base_quality_entire_window': 20, 'exclude_flags': '0xF04',
                           'BAQ': False, 'one_fragment_per_RG_and_query_name': True},
                  limitation='Exact single-edit hypotheses; nearby differences and incomplete windows remain unresolved. Window sensitivities overlap and must not be added.')
    return result


def anchor_mapping(bam, target):
    qualities, flags = Counter(), Counter()
    with pysam.AlignmentFile(str(bam), 'rb') as alignment:
        for read in alignment.fetch(target['chrom'], target['pos'] - 1, target['pos']):
            if read.flag & 0xF04:
                flags[hex(read.flag & 0xF04)] += 1
            else:
                qualities[read.mapping_quality] += 1
    return {'key': target['key'], 'nonexcluded_reads_by_MQ': dict(qualities),
            'excluded_read_flags': dict(flags),
            'scope': 'Alignment overlaps anchor before BQ/window tests; mates not collapsed, not a genotype.'}


def call_context(path, targets, sample, padding):
    result = {}
    with pysam.VariantFile(str(path)) as reader:
        if list(reader.header.samples) != [sample]:
            raise ValueError('caller sample differs')
        for target in targets:
            rows = []
            for record in reader.fetch(target['chrom'], max(0, target['pos'] - 1 - padding),
                                       target['pos'] - 1 + len(target['ref']) + padding):
                values = record.samples[sample]
                rows.append({'chrom': record.contig, 'pos': record.pos, 'stop': record.stop,
                             'ref': record.ref, 'alts': list(record.alts or []),
                             'FILTER': list(record.filter), 'QUAL': record.qual,
                             'phased': values.phased,
                             **{field: values.get(field) for field in ('GT', 'DP', 'GQ', 'AD', 'PL', 'PGT', 'PID', 'PS')}})
            result[target['key']] = rows
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('root', 'alignment-state', 'config', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    os.nice(10)
    root, output = args.root.resolve(), args.output.resolve()
    if output.exists():
        raise ValueError('existing results preserved')
    started = time.monotonic()
    provenance = json.loads(args.alignment_state.read_text())
    if provenance.get('status') != 'complete' or provenance.get('expected_lanes') != 4:
        raise ValueError('completed four-lane alignment required')
    source = Path(provenance['cram_path'])
    if source.stat().st_size != provenance['cram_bytes'] or digest(source) != provenance['cram_sha256']:
        raise ValueError('source alignment differs from receipt')
    config = json.loads(args.config.read_text())
    reference = root / 'ref/realignment/GRCh38_nochr.fa'
    if digest(reference) != config['reference_sha256']:
        raise ValueError('reference differs')
    with pysam.FastaFile(str(reference)) as fasta:
        for t in config['targets']:
            if t['key'] != f"{t['chrom']}:{t['pos']}:{t['ref']}:{t['alt']}" or len(t['ref']) == len(t['alt']):
                raise ValueError('require correctly keyed indel targets')
            for left, right in config['window_flanks_bp']:
                window_hypotheses(fasta, t, left, right)
    output.mkdir(parents=True, exist_ok=False)
    runtime = root / 'tools/bio'
    samtools, bcftools = runtime / 'bin/samtools', runtime / 'bin/bcftools'
    java = [runtime / 'lib/jvm/bin/java', '-XX:ActiveProcessorCount=2', '-Xmx6g', '-jar',
            root / 'gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar']
    paths = [Path(__file__), args.config, args.alignment_state, samtools, bcftools, java[0], java[4]]
    paths += [Path(__file__).with_name(name) for name in ('audit_baseline_reads.py', 'audit_candidate_phase.py',
               'screen_secondary_findings.py', 'candidate_read_validation.py', 'prepare_gatk_compatible_cram.py')]
    state = {'status': 'running', 'pid': os.getpid(), 'started_at': stamp(),
             'scope': 'Selected indels in complete four-lane alignment of the same biological sample.',
             'source_cram_sha256': provenance['cram_sha256'], 'reference_sha256': config['reference_sha256'],
             'source_and_tool_sha256': {str(p): digest(p) for p in paths}}
    commands = []

    def run(command, stage):
        state.update(current_stage=stage)
        save(output / 'status.json', state)
        row = {'stage': stage, 'command': list(map(str, command)), 'started_at': stamp()}
        commands.append(row)
        save(output / 'commands.json', commands)
        with (output / (stage + '.log')).open('w') as stream:
            result = subprocess.run(row['command'], stdout=stream, stderr=subprocess.STDOUT)
        row.update(exit_code=result.returncode, completed_at=stamp())
        save(output / 'commands.json', commands)
        if result.returncode:
            raise RuntimeError('stage failed: ' + stage)

    try:
        intervals = []
        for t in sorted(config['targets'], key=lambda t: (t['chrom'], t['pos'])):
            start, end = max(0, t['pos'] - 1 - config['slice_padding_bp']), t['pos'] - 1 + len(t['ref']) + config['slice_padding_bp']
            if intervals and intervals[-1][0] == t['chrom'] and intervals[-1][2] >= start:
                intervals[-1][2] = max(end, intervals[-1][2])
            else:
                intervals.append([t['chrom'], start, end])
        bed = output / 'regions.bed'
        bed.write_text(''.join(f'{chrom}\t{start}\t{end}\n' for chrom, start, end in intervals))
        bam = output / 'targets.bam'
        run([samtools, 'view', '-@', '1', '--no-PG', '-T', reference, '-b', '-M', '-L', bed,
             '--save-counts', output / 'extraction-counts.json', '-o', bam, source], 'extract-preserving-records')
        run([samtools, 'quickcheck', '-v', bam], 'quickcheck')
        run([samtools, 'index', bam], 'index')
        with pysam.AlignmentFile(str(bam)) as reads:
            groups = reads.header.to_dict().get('RG', [])
            counts = Counter(read.get_tag('RG') for read in reads)
            if len(groups) != 4 or {r.get('SM') for r in groups} != {config['sample']} or set(counts) != {r['ID'] for r in groups}:
                raise ValueError('four expected sample lanes must contribute')
        evidence = {'config': config, 'scope': state['scope'], 'region_records_by_RG': dict(counts),
                    'whole_alignment_duplicate_marks_preserved': True,
                    'source_cram_sha256': provenance['cram_sha256'],
                    'regions_BED_0based_half_open': intervals,
                    'indel_windows': [window_evidence(bam, reference, t, left, right, mq)
                        for t in config['targets'] for left, right in config['window_flanks_bp']
                        for mq in config['mapping_quality_thresholds']],
                    'anchor_mapping': [anchor_mapping(bam, t) for t in config['targets']]}
        regions = [f'{chrom}:{start + 1}-{end}' for chrom, start, end in intervals]
        limits = [item for region in regions for item in ('-L', region)]
        run(java + ['HaplotypeCaller', '-R', reference, '-I', bam, '-O', output / 'gatk.g.vcf.gz',
                    '--sample-ploidy', '2', '--pair-hmm-implementation', 'LOGLESS_CACHING', '-ERC', 'BP_RESOLUTION'] + limits,
            'gatk-reference-confidence')
        run(java + ['GenotypeGVCFs', '-R', reference, '-V', output / 'gatk.g.vcf.gz',
                    '-O', output / 'gatk.allsites.vcf.gz', '--include-non-variant-sites', 'true'] + limits, 'gatk-genotype-allsites')
        run([bcftools, 'mpileup', '-Ob', '-o', output / 'pileup.bcf', '-f', reference, '-r', ','.join(regions),
             '-a', 'AD,DP', '-q', '20', '-Q', '20', '-d', '100000', '--ff', '0xF04', bam], 'bcftools-pileup')
        run([bcftools, 'call', '-m', '-A', '-a', 'GQ', '-Oz', '-o', output / 'bcftools.allsites.vcf.gz',
             output / 'pileup.bcf'], 'bcftools-allsites')
        evidence['exact_normalized_calls'], evidence['caller_context'] = {}, {}
        for name in ('gatk', 'bcftools'):
            source_vcf = output / (name + '.allsites.vcf.gz')
            normalized = output / (name + '.normalized.vcf.gz')
            run([bcftools, 'index', '-t', '-f', source_vcf], name + '-raw-index')
            run([bcftools, 'norm', '-f', reference, '-c', 'e', '-m', '-any', '--old-rec-tag', 'SF_ORIG',
                 '-Oz', '-o', normalized, source_vcf], name + '-normalize')
            run([bcftools, 'index', '-t', normalized], name + '-normalized-index')
            evidence['exact_normalized_calls'][name] = exact_calls(source_vcf, normalized, config['targets'], config['sample'])
            evidence['caller_context'][name] = call_context(source_vcf, config['targets'], config['sample'], 50)
        run([bcftools, 'index', '-t', '-f', output / 'gatk.g.vcf.gz'], 'gatk-gvcf-index')
        evidence['gVCF_anchor_records'] = call_context(output / 'gatk.g.vcf.gz', config['targets'], config['sample'], 0)
        evidence['limitations'] = ['Selected indels in the same biological sample; not independent replication.',
            'Exact windows retain unmatched, incomplete and conflicting observations, not inferred reference genotypes.',
            'Window size and MQ sensitivity counts overlap and cannot be added; no BAQ for direct windows.',
            'Caller AD, direct fragments and raw mapping profiles measure different quantities.',
            'Missing exact variant records and reference calls at the anchor do not exclude the indel.',
            'No locus validation, clinical classification, transcript validation or causal probability is implied.']
        save(output / 'evidence.json', evidence)
        state.update(status='complete', completed_at=stamp(), elapsed_seconds=time.monotonic() - started,
                     output_sha256={p.name: digest(p) for p in output.iterdir() if p.is_file() and p.name != 'status.json'})
        save(output / 'status.json', state)
        print(json.dumps({'status': 'complete', 'targets': len(config['targets']), 'output': str(output)}))
    except BaseException as error:
        state.update(status='failed', failed_at=stamp(), error_type=type(error).__name__)
        save(output / 'status.json', state)
        raise


if __name__ == '__main__':
    main()
