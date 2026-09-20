#!/usr/bin/env python3
"""Native adaptation of pinned Broad mitochondrial calling, with explicit unfinished gates.

Uses only reads admitted by the upstream chrM subset filters. This is NOT a
whole-genome NUMT recruitment method or a contamination/clinical validation.
"""
import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

import pysam

from prepare_gatk_compatible_cram import digest, save, stamp


def sequence(path):
    with pysam.FastaFile(str(path)) as handle:
        if len(handle.references) != 1:
            raise ValueError('single mitochondrial contig required')
        return handle.references[0], handle.fetch(handle.references[0]).upper()


def verify_rotation(original, shifted, offset=8000):
    if len(original) != 16569 or shifted != original[offset:] + original[:offset]:
        raise ValueError('unexpected mitochondrial sequence or circular rotation')
    return True


def original_position(shifted_position, length=16569, offset=8000):
    if not 1 <= shifted_position <= length:
        raise ValueError('position outside shifted contig')
    return (shifted_position - 1 + offset) % length + 1


def reference_partition():
    original = set(range(576, 16025))
    lifted = {original_position(p) for p in range(8025, 9145)}
    if original & lifted or original | lifted != set(range(1, 16570)):
        raise ValueError('calling/coverage partitions do not cover each base once')
    return original, lifted


def primary_identity(path):
    rows, flags, groups = {}, Counter(), Counter()
    with pysam.AlignmentFile(str(path), check_sq=False) as handle:
        for read in handle:
            flags[str(read.flag)] += 1
            if read.is_secondary or read.is_supplementary:
                continue
            rg = read.get_tag('RG')
            groups[rg] += 1
            key = (rg, read.query_name, bool(read.is_read1), bool(read.is_read2))
            if key in rows:
                raise ValueError('duplicate primary read identity')
            if not read.is_paired or read.is_read1 == read.is_read2:
                raise ValueError('upstream paired-read method requires paired primary reads')
            bases = read.get_forward_sequence()
            qualities = read.get_forward_qualities()
            if bases is None or qualities is None:
                raise ValueError('missing read sequence/qualities')
            rows[key] = hashlib.sha256(bases.encode() + b'\x00' + bytes(qualities)).hexdigest()
    pairs = Counter((rg, name) for rg, name, _, _ in rows)
    if any(n != 2 for n in pairs.values()):
        raise ValueError('incomplete primary pair after upstream subset filters')
    names = Counter(name for rg, name in pairs)
    if any(n != 1 for n in names.values()):
        raise ValueError('QNAME collision across read groups would confuse native pairing')
    return rows, {'primary_records': len(rows), 'complete_pairs': len(pairs),
                  'flags': dict(flags), 'read_groups': dict(groups)}


def load_depth(path, length=16569):
    values = {}
    with path.open() as handle:
        for line in handle:
            chrom, pos, count = line.rstrip().split('\t')
            pos, count = int(pos), int(count)
            if chrom != 'chrM' or pos in values or count < 0:
                raise ValueError('invalid or duplicate depth position')
            values[pos] = count
    if set(values) != set(range(1, length + 1)):
        raise ValueError('missing or extra depth positions')
    return values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('root', 'reference', 'alignment', 'resources', 'config', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    os.nice(10)
    root, out = args.root.resolve(), args.output.resolve()
    if out.exists():
        raise ValueError('existing mitochondrial run preserved')
    start = time.monotonic()
    config = json.loads(args.config.read_text())
    for p, key in [(args.reference, 'reference_sha256'), (args.alignment, 'alignment_sha256')]:
        if digest(p) != config[key]:
            raise ValueError('frozen input identity differs: ' + key)
    resources = json.loads((args.resources / 'manifest.json').read_text())
    for name, info in resources['files'].items():
        if digest(args.resources / name) != info['sha256']:
            raise ValueError('resource identity differs: ' + name)
    regular = args.resources / 'Homo_sapiens_assembly38.chrM.fasta'
    shifted = args.resources / 'Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta'
    name, bases = sequence(regular)
    shifted_name, shifted_bases = sequence(shifted)
    verify_rotation(bases, shifted_bases)
    if name != 'chrM' or shifted_name != name:
        raise ValueError('unexpected resource contig')
    if hashlib.sha256(bases.encode()).hexdigest() != config['mitochondrial_sequence_sha256']:
        raise ValueError('mitochondrial sequence differs')
    with pysam.FastaFile(str(args.reference)) as genome, pysam.AlignmentFile(str(args.alignment), reference_filename=str(args.reference)) as reads:
        if genome.fetch(config['input_contig']).upper() != bases:
            raise ValueError('input mitochondrion differs from versioned resources')
        if list(zip(genome.references, genome.lengths)) != list(zip(reads.references, reads.lengths)):
            raise ValueError('alignment/reference dictionary differs')
        reads.check_index()
        groups = reads.header.to_dict().get('RG', [])
        if len(groups) != config['expected_read_groups'] or {g['SM'] for g in groups} != {config['sample']}:
            raise ValueError('sample/read groups differ')
    reference_partition()
    out.mkdir(parents=True, exist_ok=False)
    samtools, bwa = [root / 'tools/bio/bin' / n for n in ('samtools', 'bwa')]
    java = [root / 'tools/bio/lib/jvm/bin/java', '-XX:ActiveProcessorCount=2', '-Xmx6g', '-jar',
            root / 'gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar']
    state = {'status': 'running', 'pid': os.getpid(), 'started_at': stamp(), 'config': config,
             'source_and_tool_sha256': {str(p): digest(p) for p in
                 [Path(__file__), args.config, args.resources / 'manifest.json', samtools, bwa, java[0], java[4],
                  Path(__file__).with_name('prepare_gatk_compatible_cram.py')]}}
    commands = []

    def run(command, label, stdout=None):
        state['current_stage'] = label
        save(out / 'status.json', state)
        row = {'stage': label, 'command': list(map(str, command)), 'started_at': stamp()}
        commands.append(row)
        save(out / 'commands.json', commands)
        with (out / (label + '.log')).open('w') as log:
            if stdout is None:
                code = subprocess.run(row['command'], stdout=log, stderr=subprocess.STDOUT).returncode
            else:
                with stdout.open('w') as dest:
                    code = subprocess.run(row['command'], stdout=dest, stderr=log).returncode
        row.update(exit_code=code, completed_at=stamp())
        save(out / 'commands.json', commands)
        if code:
            raise RuntimeError('native mitochondrial stage failed: ' + label)

    try:
        run([samtools, 'quickcheck', '-v', args.alignment], 'quickcheck')
        subset, reverted, fastq = [out / n for n in ('subset.bam', 'reverted.bam', 'interleaved.fastq')]
        run(java + ['PrintReads', '-R', args.reference, '-I', args.alignment, '-L', config['input_contig'],
                    '--read-filter', 'MateOnSameContigOrNoMappedMateReadFilter',
                    '--read-filter', 'MateUnmappedAndUnmappedReadFilter', '-O', subset], 'subset')
        before, subset_counts = primary_identity(subset)
        run(java + ['RevertSam', '--INPUT', subset, '--OUTPUT', reverted, '--OUTPUT_BY_READGROUP', 'false',
                    '--VALIDATION_STRINGENCY', 'LENIENT', '--ATTRIBUTE_TO_CLEAR', 'FT', '--ATTRIBUTE_TO_CLEAR', 'CO',
                    '--SORT_ORDER', 'queryname', '--RESTORE_ORIGINAL_QUALITIES', 'false', '--RESTORE_HARDCLIPS', 'false'], 'revert')
        after, reverted_counts = primary_identity(reverted)
        if before != after:
            raise ValueError('reversion changed a primary read, quality or identity')
        save(out / 'read-reversion.json', {'before': subset_counts, 'after': reverted_counts,
             'all_primary_forward_sequences_and_qualities_identical': True,
             'scope': 'Subset admitted by upstream filters, not every possible mitochondrial/NUMT read.'})
        del before, after
        run(java + ['SamToFastq', '--INPUT', reverted, '--FASTQ', fastq, '--INTERLEAVE', 'true', '--NON_PF', 'true'], 'fastq')
        alignments, calls, stats = {}, {}, {}
        for arm, ref, interval in [('regular', regular, 'chrM:576-16024'), ('shifted', shifted, 'chrM:8025-9144')]:
            sam, merged, marked, bam = [out / (arm + suffix) for suffix in ('.sam', '.merged.bam', '.marked.bam', '.bam')]
            run([bwa, 'mem', '-K', '100000000', '-p', '-v', '3', '-t', '2', '-Y', ref, fastq], arm + '-bwa', sam)
            run(java + ['MergeBamAlignment', '--ALIGNED_BAM', sam, '--UNMAPPED_BAM', reverted, '--OUTPUT', merged,
                        '--REFERENCE_SEQUENCE', ref, '--VALIDATION_STRINGENCY', 'LENIENT', '--EXPECTED_ORIENTATIONS', 'FR',
                        '--ATTRIBUTES_TO_RETAIN', 'X0', '--ATTRIBUTES_TO_REMOVE', 'NM', '--ATTRIBUTES_TO_REMOVE', 'MD',
                        '--PAIRED_RUN', 'true', '--SORT_ORDER', 'queryname', '--ALIGNED_READS_ONLY', 'false',
                        '--CLIP_ADAPTERS', 'false', '--ADD_MATE_CIGAR', 'true', '--MAX_INSERTIONS_OR_DELETIONS', '-1',
                        '--PRIMARY_ALIGNMENT_STRATEGY', 'MostDistant', '--UNMAPPED_READ_STRATEGY', 'COPY_TO_TAG',
                        '--ALIGNER_PROPER_PAIR_FLAGS', 'true', '--UNMAP_CONTAMINANT_READS', 'true', '--ADD_PG_TAG_TO_READS', 'false'], arm + '-merge')
            run(java + ['MarkDuplicates', '--INPUT', merged, '--OUTPUT', marked, '--METRICS_FILE', out / (arm + '.duplicates.txt'),
                        '--VALIDATION_STRINGENCY', 'LENIENT', '--OPTICAL_DUPLICATE_PIXEL_DISTANCE', '2500',
                        '--ASSUME_SORT_ORDER', 'queryname', '--CLEAR_DT', 'false', '--ADD_PG_TAG_TO_READS', 'false'], arm + '-duplicates')
            run([samtools, 'sort', '-@', '2', '-o', bam, marked], arm + '-sort')
            run([samtools, 'index', bam], arm + '-index')
            current, counts = primary_identity(bam)
            expected, _ = primary_identity(reverted)
            # MergeBamAlignment can hard-clip contaminant alignments. Fail rather than silently accepting changed primary data.
            if current != expected:
                raise ValueError('realignment changed primary read identity/sequence/quality')
            save(out / (arm + '.read-counts.json'), counts)
            del current, expected
            raw = out / (arm + '.raw.vcf.gz')
            run(java + ['Mutect2', '-R', ref, '-I', bam, '-L', interval, '-O', raw,
                        '--read-filter', 'MateOnSameContigOrNoMappedMateReadFilter',
                        '--read-filter', 'MateUnmappedAndUnmappedReadFilter', '--annotation', 'StrandBiasBySample',
                        '--mitochondria-mode', 'true', '--max-reads-per-alignment-start', '75', '--max-mnp-distance', '0',
                        '--pair-hmm-implementation', 'LOGLESS_CACHING'], arm + '-mutect2')
            depth = out / (arm + '.depth.tsv')
            run([samtools, 'depth', '-aa', '-q', '20', '-Q', '20', '-G', '0xF04', '-s', bam], arm + '-depth', depth)
            alignments[arm], calls[arm], stats[arm] = bam, raw, Path(str(raw) + '.stats')
        lifted, rejected, merged_vcf = [out / n for n in ('shifted.lifted.vcf.gz', 'shifted.rejected.vcf.gz', 'merged.raw.vcf.gz')]
        run(java + ['LiftoverVcf', '--INPUT', calls['shifted'], '--OUTPUT', lifted, '--REFERENCE_SEQUENCE', regular,
                    '--CHAIN', args.resources / 'ShiftBack.chain', '--REJECT', rejected], 'liftover')
        run(java + ['MergeVcfs', '--INPUT', calls['regular'], '--INPUT', lifted, '--OUTPUT', merged_vcf], 'merge-vcfs')
        merged_stats = out / 'merged.stats'
        run(java + ['MergeMutectStats', '--stats', stats['regular'], '--stats', stats['shifted'], '-O', merged_stats], 'merge-stats')
        filtered = out / 'filtered.assuming-no-contamination.vcf.gz'
        run(java + ['FilterMutectCalls', '-R', regular, '-V', merged_vcf, '-O', filtered, '--stats', merged_stats,
                    '--mitochondria-mode', 'true', '--max-alt-allele-count', '4', '--min-allele-fraction', '0',
                    '--contamination-estimate', '0'], 'initial-filter')
        masked = out / 'masked.assuming-no-contamination.vcf.gz'
        run(java + ['VariantFiltration', '-V', filtered, '-O', masked, '--apply-allele-specific-filters',
                    '--mask', args.resources / 'blacklist_sites.hg38.chrM.bed', '--mask-name', 'blacklisted_site'], 'mask')
        depths = {arm: load_depth(out / (arm + '.depth.tsv')) for arm in ('regular', 'shifted')}
        original_region, control_region = reference_partition()
        combined = {p: depths['regular'][p] for p in original_region}
        combined.update({original_position(p): depths['shifted'][p] for p in range(8025, 9145)})
        if set(combined) != set(range(1, 16570)):
            raise ValueError('combined coverage missingness')
        with (out / 'combined.depth.tsv').open('w') as handle:
            handle.write('chrom\tpos\tref\tMQ20_BQ20_depth\tarm\n')
            for p in sorted(combined):
                handle.write(f'chrM\t{p}\t{bases[p-1]}\t{combined[p]}\t' + ('regular' if p in original_region else 'shifted') + '\n')
        variants = []
        with pysam.VariantFile(str(masked)) as reader:
            if list(reader.header.samples) != [config['sample']]:
                raise ValueError('unexpected final sample')
            for record in reader:
                if record.ref != bases[record.start:record.stop]:
                    raise ValueError('final reference allele differs')
                sample = record.samples[config['sample']]
                variants.append({'chrom': record.contig, 'pos': record.pos, 'ref': record.ref, 'alts': list(record.alts or []),
                    'filters': list(record.filter), 'GT': sample.get('GT'), 'AD': sample.get('AD'), 'DP': sample.get('DP'),
                    'AF': sample.get('AF'), 'F1R2': sample.get('F1R2'), 'F2R1': sample.get('F2R1')})
        with pysam.VariantFile(str(rejected)) as reader:
            rejected_rows = [{'chrom': x.contig, 'pos': x.pos, 'ref': x.ref, 'alts': x.alts, 'filters': list(x.filter)} for x in reader]
        concrete_depth = [combined[p] for p in sorted(combined) if bases[p-1] in 'ACGT']
        summary = {'completed_at': stamp(), 'scope': config['scope'], 'sample': config['sample'],
            'subset_read_counts': subset_counts, 'variants': variants, 'rejected_liftover_records': rejected_rows,
            'coverage': {'positions': len(combined), 'ACGT_positions': len(concrete_depth),
                'mean_concrete_depth': sum(concrete_depth)/len(concrete_depth),
                'concrete_bases_at_least': {str(n): sum(x >= n for x in concrete_depth) for n in [1, 10, 20, 100, 1000]}},
            'limits': ['Contamination was not estimated. Filter output assumes zero and is explicitly provisional.',
                'No dedicated NUMT recruitment or nuclear competition validation; mt-only remapping can increase ambiguous-read MQ.',
                'Subset/reversion and mitochondrial duplicate marking differ from original whole-genome alignment; both are preserved.',
                'Caller AF and filtered per-base depth use different read/quality models; neither is molecular tissue heteroplasmy truth.',
                'Liftover rejects and original multi-allele records remain visible; no absent variant is a clinical negative.',
                'Reference is the verified versioned rCRS sequence, not a guarantee of full variant/structural sensitivity.']}
        save(out / 'summary.json', summary)
        state.update(status='complete', completed_at=stamp(), elapsed_seconds=time.monotonic()-start,
                     output_sha256={p.name: digest(p) for p in out.iterdir() if p.is_file() and p.name != 'status.json'})
        save(out / 'status.json', state)
        print(json.dumps({'status': 'complete', 'variants': len(variants), 'rejected_liftover_records': len(rejected_rows)}))
    except BaseException as error:
        state.update(status='failed', failed_at=stamp(), error_type=type(error).__name__)
        save(out / 'status.json', state)
        raise


if __name__ == '__main__':
    main()
