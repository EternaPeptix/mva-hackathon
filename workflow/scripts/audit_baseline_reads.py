#!/usr/bin/env python3
"""Bounded retained-lane read review; export aggregate evidence, not read identities."""
import argparse
from collections import Counter, defaultdict
import itertools
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import pysam

from audit_candidate_phase import spanning_sequence
from candidate_read_validation import pileup_evidence
from screen_secondary_findings import fingerprint, original_call


def collapse_fragments(fragments, labels):
    counts, strands = Counter(), defaultdict(Counter)
    for observations in fragments.values():
        values = {label for label, reverse in observations}
        label = next(iter(values)) if len(values) == 1 else 'mate_conflict'
        orientations = {reverse for _, reverse in observations}
        direction = 'both' if len(orientations) > 1 else 'reverse' if True in orientations else 'forward'
        counts[label] += 1
        strands[label][direction] += 1
    return {'counts': {label: counts[label] for label in labels + ['other', 'mate_conflict']},
            'strand_counts': dict(strands), 'eligible_fragments': len(fragments)}


def count_snv_configuration(bam, reference, sites, minimum_mapping_quality=20, minimum_base_quality=20):
    if not 2 <= len(sites) <= 8 or len({s[0] for s in sites}) != 1:
        raise ValueError('require two to eight SNPs on one contig')
    if any(len(s[2]) != 1 or len(s[3]) != 1 or s[2] == s[3] for s in sites):
        raise ValueError('SNV configuration requires distinct single-base alleles')
    if len({s[1] for s in sites}) != len(sites):
        raise ValueError('duplicate SNP position')
    contig = sites[0][0]
    with pysam.FastaFile(str(reference)) as fasta:
        if any(fasta.fetch(chrom, pos - 1, pos).upper() != ref for chrom, pos, ref, alt in sites):
            raise ValueError('SNP reference mismatch')
    positions = {site[1] - 1 for site in sites}
    fragments, excluded = defaultdict(list), Counter()
    with pysam.AlignmentFile(str(bam), 'rb') as alignment:
        for read in alignment.fetch(contig, min(positions), max(positions) + 1):
            if read.flag & 0xF04 or read.mapping_quality < minimum_mapping_quality:
                excluded['flags_or_mapping_quality'] += 1
                continue
            if read.query_qualities is None or any(op == 3 for op, length in read.cigartuples or []):
                excluded['missing_quality_or_reference_skip'] += 1
                continue
            mapping = {r: q for q, r in read.get_aligned_pairs() if r in positions}
            if any(mapping.get(pos) is None for pos in positions):
                excluded['not_covering_all_sites_or_deletion'] += 1
                continue
            if any(read.query_qualities[mapping[pos]] < minimum_base_quality for pos in positions):
                excluded['low_target_base_quality'] += 1
                continue
            bases = [read.query_sequence[mapping[pos - 1]].upper() for _, pos, _, _ in sites]
            bits = ['0' if base == ref else '1' if base == alt else '?' for base, (_, _, ref, alt) in zip(bases, sites)]
            label = '|'.join(bits) if '?' not in bits else 'other'
            group = read.get_tag('RG') if read.has_tag('RG') else ''
            fragments[(group, read.query_name)].append((label, read.is_reverse))
    labels = ['|'.join(bits) for bits in itertools.product('01', repeat=len(sites))]
    return dict(collapse_fragments(fragments, labels), excluded_read_counts=dict(excluded), sites=sites,
        method='Aligned-base co-occurrence at specified SNPs; intervening sequence need not match reference.',
        filters={'minimum_mapping_quality': minimum_mapping_quality, 'minimum_target_base_quality': minimum_base_quality,
                 'exclude_flags': '0xF04', 'one_fragment_per_RG_and_query_name': True,
                 'one_read_must_cover_all_sites': True, 'BAQ_applied': False},
        limitation='Co-occurrence does not validate locus assignment in homologous regions or establish clinical significance.')


def indel_window_evidence(bam, reference, site):
    contig, pos, ref, alt = site
    start, end = pos - 2, pos - 1 + len(ref) + 2
    with pysam.FastaFile(str(reference)) as fasta:
        expected_ref = fasta.fetch(contig, start, end).upper()
    if expected_ref[1:1 + len(ref)] != ref:
        raise ValueError('indel reference mismatch')
    expected_alt = expected_ref[:1] + alt + expected_ref[1 + len(ref):]
    if expected_ref == expected_alt:
        raise ValueError('identical indel hypotheses')
    fragments, excluded = defaultdict(list), Counter()
    with pysam.AlignmentFile(str(bam), 'rb') as alignment:
        for read in alignment.fetch(contig, start, end):
            if read.flag & 0xF04 or read.mapping_quality < 20:
                excluded['flags_or_mapping_quality'] += 1
                continue
            sequence = spanning_sequence(read, start, end, 20)
            if sequence is None:
                excluded['not_spanning_or_low_quality_or_refskip'] += 1
                continue
            label = 'REF' if sequence == expected_ref else 'ALT' if sequence == expected_alt else 'other'
            group = read.get_tag('RG') if read.has_tag('RG') else ''
            fragments[(group, read.query_name)].append((label, read.is_reverse))
    return dict(collapse_fragments(fragments, ['REF', 'ALT']), site=site, excluded_read_counts=dict(excluded),
        method='Exact CIGAR-derived local allele sequence including flanks, MQ/BQ >=20, flags 0xF04 excluded, no BAQ.',
        limitation='Nearby changes or equivalent indel alignments can be unmatched; this does not establish clinical effect.')


def raw_snv_mapping_profile(bam, site):
    contig, pos, ref, alt = site
    counts, excluded, alt_features = defaultdict(Counter), Counter(), Counter()
    with pysam.AlignmentFile(str(bam), 'rb') as alignment:
        for read in alignment.fetch(contig, pos - 1, pos):
            if read.flag & 0xF04:
                excluded['flags'] += 1
                continue
            query = next((q for q, r in read.get_aligned_pairs() if r == pos - 1), None)
            if query is None or read.query_qualities is None or read.query_qualities[query] < 20:
                excluded['no_aligned_base_or_low_base_quality'] += 1
                continue
            base = read.query_sequence[query].upper()
            label = 'REF' if base == ref else 'ALT' if base == alt else 'other'
            mq = read.mapping_quality
            bucket = '0' if mq == 0 else '1-19' if mq < 20 else '20-39' if mq < 40 else '40-59' if mq < 60 else '60+'
            counts[label][bucket] += 1
            if label == 'ALT':
                alt_features['reads'] += 1
                alt_features['soft_clipped_reads'] += any(op == 4 for op, length in read.cigartuples or [])
                if read.has_tag('NM') and read.get_tag('NM') >= 5:
                    alt_features['NM_at_least_5_reads'] += 1
    return {'site': site, 'read_counts_by_mapping_quality': dict(counts), 'excluded_read_counts': dict(excluded),
            'ALT_read_features': dict(alt_features),
            'scope': 'Raw per-read QC profile, BQ>=20 with excluded flags 0xF04; mates are not collapsed and BAQ is not applied. Not a genotype.'}


def exact_calls(raw, normalized, targets, sample):
    result = []
    with pysam.VariantFile(str(raw)) as source, pysam.VariantFile(str(normalized)) as reader:
        if list(source.header.samples) != [sample]:
            raise ValueError('unexpected caller sample')
        for row in targets:
            matches = [original_call(record, source, sample) for record in reader.fetch(row['chrom'], row['pos'] - 1, row['pos'])
                       if record.pos == row['pos'] and record.ref == row['ref'] and record.alts == (row['alt'],)]
            result.append({'key': row['key'], 'exact_normalized_calls': matches,
                           'absence_is_not_reference_or_callability_proof': True})
    return result


def run(args):
    os.nice(10)
    started = time.monotonic()
    root, output = args.root.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    config = json.loads(args.config.read_text())
    manifest = json.loads((args.slices / 'slice-provenance.json').read_text())
    state = {'status': 'running', 'pid': os.getpid(), 'started_at': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
             'scope': config['scope'], 'commands': [], 'invocation': [sys.executable] + sys.argv}
    def save():
        (output / 'run-status.json').write_text(json.dumps(state, indent=2) + '\n')
    save()
    try:
        paths = {'config': args.config, 'slice_manifest': args.slices / 'slice-provenance.json',
                 'reference': root / 'ref/vendor/GRCh38_nochr.fa',
                 'reference_index': root / 'ref/vendor/GRCh38_nochr.fa.fai',
                 'reference_dictionary': root / 'ref/vendor/GRCh38_nochr.dict',
                 'script': Path(__file__), 'allele_helpers': Path(__file__).with_name('screen_secondary_findings.py'),
                 'pileup_helpers': Path(__file__).with_name('candidate_read_validation.py'),
                 'sequence_helpers': Path(__file__).with_name('audit_candidate_phase.py'),
                 'samtools': root / 'tools/bio/bin/samtools', 'bcftools': root / 'tools/bio/bin/bcftools',
                 'java': root / 'tools/bio/lib/jvm/bin/java',
                 'gatk_jar': root / 'gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar'}
        inputs = {name: fingerprint(path) for name, path in paths.items()}
        state['inputs'] = inputs
        if inputs['config']['sha256'] != manifest['config_sha256'] or manifest['status'] != 'complete':
            raise ValueError('slice/config identity mismatch')
        if inputs['reference']['sha256'] != config['reference_sha256']:
            raise ValueError('wrong reference')
        reference = paths['reference']
        with pysam.FastaFile(str(reference)) as fasta:
            dictionary = dict(zip(fasta.references, fasta.lengths))
            for row in config['targets']:
                if fasta.fetch(row['chrom'], row['pos'] - 1, row['pos'] - 1 + len(row['ref'])).upper() != row['ref']:
                    raise ValueError('target reference mismatch')
        slices = []
        for row in manifest['files']:
            path = args.slices / row['slice_name']
            inputs[row['slice_name']] = fingerprint(path)
            if inputs[row['slice_name']]['sha256'] != row['slice_sha256']:
                raise ValueError('slice checksum mismatch')
            with pysam.AlignmentFile(str(path), 'rb') as alignment:
                if dict(zip(alignment.references, alignment.lengths)) != dictionary:
                    raise ValueError('slice/reference dictionaries differ')
                if {rg['SM'] for rg in alignment.header.to_dict().get('RG', [])} != {config['sample']}:
                    raise ValueError('slice read-group sample differs')
            slices.append(str(path))
        env = dict(os.environ, PATH=str(root / 'tools/bio/bin') + os.pathsep + os.environ.get('PATH', ''))
        def execute(command, label):
            state['current_stage'] = label
            state['commands'].append(command)
            save()
            with (output / (label + '.log')).open('w') as log:
                subprocess.run(list(map(str, command)), env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
        execute(['samtools', 'merge', '-@', '1', '-o', str(output / 'merged.bam')] + slices, 'merge')
        for command, label in [(['samtools', 'collate', '-@', '1', '-o', str(output / 'collated.bam'), str(output / 'merged.bam')], 'collate'),
            (['samtools', 'fixmate', '-@', '1', '-m', str(output / 'collated.bam'), str(output / 'fixmate.bam')], 'fixmate'),
            (['samtools', 'sort', '-@', '1', '-m', '128M', '-o', str(output / 'sorted.bam'), str(output / 'fixmate.bam')], 'sort'),
            (['samtools', 'markdup', '-@', '1', '-f', str(output / 'duplicate-marking.txt'), str(output / 'sorted.bam'), str(output / 'markdup.bam')], 'markdup'),
            (['samtools', 'quickcheck', '-v', str(output / 'markdup.bam')], 'quickcheck'),
            (['samtools', 'index', str(output / 'markdup.bam')], 'index')]:
            execute(command, label)
        bam = output / 'markdup.bam'
        tuples = [(r['chrom'], r['pos'], r['ref'], r['alt']) for r in config['targets']]
        snps = [s for s in tuples if len(s[2]) == len(s[3]) == 1]
        by_key = dict(zip([r['key'] for r in config['targets']], tuples))
        evidence = {'config': config, 'source_slices': manifest,
            'SNP_BAQ_pileup': pileup_evidence(bam, reference, snps),
            'SNP_raw_mapping_profiles': [raw_snv_mapping_profile(bam, site) for site in snps],
            'indel_windows': [indel_window_evidence(bam, reference, site) for site in tuples if site not in snps],
            'SNP_pair_cooccurrence': [count_snv_configuration(bam, reference, [by_key[k] for k in pair]) for pair in config['snv_phase_pairs']]}
        cluster_keys = sorted({k for pair in config['snv_phase_pairs'] for k in pair}, key=lambda k: int(k.split(':')[1]))
        evidence['SERPINA1_four_SNP_cooccurrence'] = count_snv_configuration(bam, reference, [by_key[k] for k in cluster_keys])
        regions = [f'{chrom}:{start + 1}-{end}' for chrom, start, end in manifest['regions_BED_0based_half_open']]
        gatk = output / 'gatk.raw.vcf.gz'
        execute([str(paths['java']), '-XX:ActiveProcessorCount=4', '-Xmx8g', '-jar', str(paths['gatk_jar']),
                 'HaplotypeCaller', '-R', str(reference), '-I', str(bam), '-O', str(gatk),
                 '--pair-hmm-implementation', 'LOGLESS_CACHING'] + [x for region in regions for x in ('-L', region)], 'gatk')
        execute(['bcftools', 'mpileup', '-Ob', '-o', str(output / 'pileup.bcf'), '-f', str(reference), '-r', ','.join(regions),
                 '-a', 'AD,DP,SP', '-q', '20', '-Q', '20', '-d', '100000', '--ff', '0xF04', str(bam)], 'bcftools-pileup')
        bcf = output / 'bcftools.raw.vcf.gz'
        execute(['bcftools', 'call', '-mv', '-Oz', '-o', str(bcf), str(output / 'pileup.bcf')], 'bcftools-call')
        for name, source in [('gatk', gatk), ('bcftools', bcf)]:
            execute(['bcftools', 'index', '-t', '-f', str(source)], name + '-raw-index')
            normalized = output / (name + '.normalized.vcf.gz')
            execute(['bcftools', 'norm', '-f', str(reference), '-c', 'e', '-m', '-any', '--old-rec-tag', 'SF_ORIG',
                     '-Oz', '-o', str(normalized), str(source)], name + '-normalize')
            execute(['bcftools', 'index', '-t', str(normalized)], name + '-normalized-index')
            evidence[name + '_calls'] = exact_calls(source, normalized, config['targets'], config['sample'])
        evidence['limitations'] = ['Two historical lanes, not new four-lane data; raw regional discovery calls, not clinical accuracy.',
            'Variant-only caller absence is not a reference genotype or evidence of callability.',
            'Direct SNP co-occurrence uses raw aligned target bases; neighboring indels and locus mapping require separate review.',
            'Pileup uses BAQ and overlap handling; direct phase/indel windows do not use BAQ. These counts are different measurements.',
            'Read support does not validate transcript relevance, functional effect, population rarity or disease causality.']
        (output / 'evidence.json').write_text(json.dumps(evidence, indent=2) + '\n')
        for record in inputs.values():
            stat = Path(record['path']).stat()
            if (stat.st_size, stat.st_mtime_ns) != (record['bytes'], record['mtime_ns']):
                raise ValueError('input changed during review')
        state['output_hashes'] = {p.name: fingerprint(p)['sha256'] for p in output.iterdir() if p.is_file() and p.name != 'run-status.json'}
        state.update(status='complete', completed_at=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
                     elapsed_seconds=time.monotonic() - started, current_stage='complete')
        print(json.dumps({'status': 'complete', 'elapsed_seconds': state['elapsed_seconds'], 'targets': len(tuples),
                          'SERPINA1_four_SNP_cooccurrence': evidence['SERPINA1_four_SNP_cooccurrence']}))
    except BaseException as error:
        state.update(status='failed', error=str(error), elapsed_seconds=time.monotonic() - started)
        raise
    finally:
        save()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--slices', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args())


if __name__ == '__main__':
    main()
