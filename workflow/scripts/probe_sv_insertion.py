#!/usr/bin/env python3
"""Target a known synthetic insertion with two native small-variant callers.

This is a post hoc, same-read comparison. It is not blind discovery, independent
biological replication, general insertion sensitivity, or a patient analysis.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import time

import pysam

from probe_structural_variants import save, stamp
from simulate_sv_reads import sha


def local_haplotype_matches(sequence, window_start, cut, inserted, position, ref, alt):
    """Compare sequences after applying each allele, allowing equivalent left shifts."""
    offset = position - 1 - window_start
    truth_offset = cut - window_start
    if not (0 <= offset <= len(sequence) - len(ref) and 0 <= truth_offset <= len(sequence)):
        return False
    if sequence[offset:offset + len(ref)] != ref or not alt or alt[0] == '<' or any(c in alt for c in '[]*'):
        return False
    observed = sequence[:offset] + alt + sequence[offset + len(ref):]
    expected = sequence[:truth_offset] + inserted + sequence[truth_offset:]
    return observed == expected


def run(experiment, output, samtools, bcftools, java, gatk):
    if output.exists():
        raise ValueError('comparison output already exists')
    original = json.loads((experiment / 'summary.json').read_text())
    truth_path = experiment / 'fixture/truth.json'
    truth = json.loads(truth_path.read_text())
    reference = experiment / 'fixture/reference.fa'
    alignment = experiment / 'aligned.bam'
    if truth['scope'] != 'Synthetic diploid simple-SV reference; no patient sequence or observed biological effect.':
        raise ValueError('this probe requires the synthetic fixture')
    if sha(truth_path) != original['truth_sha256'] or sha(reference) != truth['files']['reference.fa']['sha256']:
        raise ValueError('synthetic truth/reference identity changed')
    if sha(alignment) != original['output_sha256']['aligned.bam']:
        raise ValueError('aligned synthetic reads differ from original experiment')
    with pysam.AlignmentFile(str(alignment)) as reader:
        if {row['SM'] for row in reader.header.to_dict()['RG']} != {'SYNTH_SV'}:
            raise ValueError('expected synthetic sample identity missing')
    insertions = [row for row in truth['events_on_one_haplotype'] if row['type'] == 'INS']
    if len(insertions) != 1:
        raise ValueError('exactly one known synthetic insertion required')
    event = insertions[0]
    chrom = event['chrom']; cut = event['left_cut_0based']
    left = cut - 1000; right = cut + 1000
    with pysam.FastaFile(str(reference)) as genome:
        sequence = genome.fetch(chrom, left, right).upper()
        reference_lengths = dict(zip(genome.references, genome.lengths))
    output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic(); commands = []
    state = {'status': 'running', 'started_at': stamp(), 'pid': os.getpid()}
    save(output / 'status.json', state)
    tools = {name: {'path': str(path), 'sha256': sha(path)} for name, path in
             [('samtools', samtools), ('bcftools', bcftools), ('java', java), ('gatk', gatk)]}

    def command(args):
        log = output / ('step-%02d.log' % len(commands)); tick = time.monotonic()
        with log.open('w') as stream:
            result = subprocess.run(list(map(str, args)), stdout=stream, stderr=subprocess.STDOUT)
        commands.append({'command': list(map(str, args)), 'exit_code': result.returncode,
                         'elapsed_seconds': time.monotonic() - tick, 'log_sha256': sha(log)})
        save(output / 'commands.json', commands)
        if result.returncode:
            raise RuntimeError('native synthetic comparison failed; inspect retained log')

    try:
        dictionary = reference.with_suffix('.dict')
        if not dictionary.exists():
            command([samtools, 'dict', '-o', dictionary, reference])
        header = pysam.AlignmentHeader.from_text(dictionary.read_text()).to_dict()
        if {row['SN']: row['LN'] for row in header['SQ']} != reference_lengths:
            raise ValueError('reference dictionary differs from FASTA')
        region = '%s:%d-%d' % (chrom, left, right)
        command([bcftools, 'mpileup', '-Ob', '-o', output / 'pileup.bcf', '-f', reference,
                 '-r', region, '-a', 'AD,DP', '-q', '20', '-Q', '20', alignment])
        command([bcftools, 'call', '-mv', '-Oz', '-o', output / 'bcftools.vcf.gz', output / 'pileup.bcf'])
        command([bcftools, 'index', output / 'bcftools.vcf.gz'])
        command([java, '-Xmx2g', '-jar', gatk, 'HaplotypeCaller', '-R', reference, '-I', alignment,
                 '-O', output / 'gatk.vcf.gz', '-L', region, '--native-pair-hmm-threads', '2'])
        callers = {}
        for caller in ('bcftools', 'gatk'):
            rows = []
            with pysam.VariantFile(str(output / (caller + '.vcf.gz'))) as reader:
                if list(reader.header.samples) != ['SYNTH_SV']:
                    raise ValueError('caller sample identity mismatch')
                for record in reader:
                    genotype = record.samples['SYNTH_SV']
                    for index, alt in enumerate(record.alts or (), 1):
                        rows.append({'chrom': record.contig, 'pos': record.pos, 'ref': record.ref,
                                     'alt': alt, 'allele_index': index,
                                     'length_delta': len(alt) - len(record.ref),
                                     'GT': genotype.get('GT'), 'AD': genotype.get('AD'),
                                     'GQ': genotype.get('GQ'), 'filter': list(record.filter),
                                     'exact_local_haplotype_match': record.contig == chrom and
                                     local_haplotype_matches(sequence, left, cut, event['inserted_sequence'],
                                                             record.pos, record.ref, alt)})
            callers[caller] = rows
        for identity in tools.values():
            if sha(Path(identity['path'])) != identity['sha256']:
                raise ValueError('native executable changed during comparison')
        result = {'completed_at': stamp(), 'elapsed_seconds': time.monotonic() - started,
                  'probe_source_sha256': sha(Path(__file__)),
                  'original_summary_sha256': sha(experiment / 'summary.json'),
                  'truth_sha256': sha(truth_path), 'reference_sha256': sha(reference),
                  'alignment_sha256': sha(alignment), 'tools': tools, 'commands': commands,
                  'known_truth_region_1based_inclusive': region,
                  'comparison_window_0based_half_open': [left, right], 'callers': callers,
                  'outputs': {name: sha(output / name) for name in
                              ('bcftools.vcf.gz', 'bcftools.vcf.gz.csi', 'gatk.vcf.gz', 'gatk.vcf.gz.tbi')},
                  'scope': __doc__}
        save(output / 'summary.json', result)
        state.update(status='complete', completed_at=result['completed_at'])
        save(output / 'status.json', state)
        return result
    except BaseException as error:
        state.update(status='failed', failed_at=stamp(), error_type=type(error).__name__)
        save(output / 'status.json', state)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('experiment', 'output', 'samtools', 'bcftools', 'java', 'gatk'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args(); os.nice(10)
    result = run(args.experiment, args.output, args.samtools, args.bcftools, args.java, args.gatk)
    print(json.dumps({'completed_at': result['completed_at'], 'elapsed_seconds': result['elapsed_seconds'],
                      'callers': result['callers']}, indent=2))


if __name__ == '__main__':
    main()
