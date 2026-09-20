#!/usr/bin/env python3
"""Inspect selected original caller alternatives without rerunning a genotyper.

Reuse qualified BAQ/fragment and exact-window helpers. Same-position alternate
alleles remain separate hypotheses, and caller records retain native ploidy.
Raw regional reads and case-level outputs stay in the owned output directory.
"""
import argparse
from collections import Counter
import json
import os
from pathlib import Path
import re
import subprocess

import pysam

from audit_final_alignment_indels import anchor_mapping, call_context, window_evidence
from audit_targeted_snvs import fragment_evidence
from prepare_gatk_compatible_cram import digest, save, stamp
from summarize_targeted_snvs import raw_base_profile


def validate_selection(targets, reference):
    if not targets:
        raise ValueError('No selected hypotheses')
    keys, refs = set(), {}
    with pysam.FastaFile(str(reference)) as fasta:
        for row in targets:
            chrom, pos, ref, alt = (row[k] for k in ('chrom', 'pos', 'ref', 'alt'))
            if (not isinstance(pos, int) or pos < 1 or ref == alt
                    or not re.fullmatch('[ACGT]+', ref) or not re.fullmatch('[ACGT]+', alt)):
                raise ValueError('Require distinct concrete alleles and a positive position')
            key = f'{chrom}:{pos}:{ref}:{alt}'
            if row['key'] != key or key in keys:
                raise ValueError('Incorrect or duplicate hypothesis key')
            locus = (chrom, pos)
            if locus in refs and refs[locus] != ref:
                raise ValueError('Same-position reference representations differ')
            if fasta.fetch(chrom, pos-1, pos-1+len(ref)).upper() != ref:
                raise ValueError('Selected REF differs from reference')
            keys.add(key)
            refs[locus] = ref


def require_known_mapping_quality(bam):
    """The reused SNP helper has no MQ255 exclusion; fail before counting it."""
    with pysam.AlignmentFile(str(bam), 'rb') as stream:
        for read in stream:
            if not read.flag & 0xF04 and read.mapping_quality == 255:
                raise ValueError('Unknown MQ255 requires a separately qualified counting path')


def merged_regions(targets, padding, reference):
    regions = []
    with pysam.FastaFile(str(reference)) as fasta:
        for row in sorted(targets, key=lambda x: (x['chrom'], x['pos'])):
            chrom = row['chrom']
            start = max(0, row['pos']-1-padding)
            end = min(fasta.get_reference_length(chrom), row['pos']-1+len(row['ref'])+padding)
            if regions and regions[-1][0] == chrom and regions[-1][2] >= start:
                regions[-1][2] = max(regions[-1][2], end)
            else:
                regions.append([chrom, start, end])
    return regions


def run(args):
    root, output = args.root.resolve(), args.output.resolve()
    if output.exists() or output.is_relative_to(Path(__file__).resolve().parents[2]):
        raise ValueError('Use a new output directory outside the source tree')
    config = json.loads(args.config.read_text())
    output.mkdir(parents=True)
    os.nice(10)
    state = {'status': 'running', 'pid': os.getpid(), 'started_at': stamp(),
             'current_stage': 'verify input identities', 'config_sha256': digest(args.config),
             'script_sha256': digest(Path(__file__))}
    save(output/'status.json', state)
    commands = []
    try:
        inputs, monitored = {}, {}
        for name, row in config['inputs'].items():
            p = root/row['path']
            if not p.is_file() or p.is_symlink() or digest(p) != row['sha256']:
                raise ValueError('Input identity differs: '+name)
            stat = p.stat()
            inputs[name] = {'path': str(p), 'bytes': stat.st_size, 'sha256': row['sha256']}
            monitored[str(p)] = (stat.st_size, stat.st_mtime_ns)
        for name, expected in config['helper_sha256'].items():
            p = Path(__file__).with_name(name)
            if digest(p) != expected:
                raise ValueError('Helper identity differs: '+name)
            monitored[str(p)] = (p.stat().st_size, p.stat().st_mtime_ns)
        alignment, reference = (Path(inputs[k]['path']) for k in ('alignment', 'reference'))
        samtools = Path(inputs['samtools']['path'])
        targets = config['targets']
        validate_selection(targets, reference)
        regions = merged_regions(targets, config['slice_padding_bp'], reference)
        bed = output/'regions.bed'
        bed.write_text(''.join(f'{c}\t{a}\t{b}\n' for c, a, b in regions))
        bam = output/'targets.bam'

        def command(argv, label):
            state['current_stage'] = label
            save(output/'status.json', state)
            row = {'command': list(map(str, argv)), 'stage': label, 'started_at': stamp()}
            commands.append(row)
            save(output/'commands.json', commands)
            with (output/(label+'.log')).open('w') as log:
                result = subprocess.run(row['command'], stdout=log, stderr=subprocess.STDOUT)
            row.update(exit_code=result.returncode, completed_at=stamp())
            save(output/'commands.json', commands)
            if result.returncode:
                raise RuntimeError('Command failed: '+label)

        command([samtools, 'view', '-@', '1', '--no-PG', '-T', reference, '-b', '-M', '-L', bed,
                 '--save-counts', output/'extraction-counts.json', '-o', bam, alignment], 'extract-selected-regions')
        command([samtools, 'quickcheck', '-v', bam], 'quickcheck')
        command([samtools, 'index', bam], 'index')
        with pysam.AlignmentFile(str(bam), 'rb') as stream:
            groups = stream.header.to_dict().get('RG', [])
            if len(groups) != config['expected_read_groups'] or {x.get('SM') for x in groups} != {config['sample']}:
                raise ValueError('Unexpected sample/read-group header')
            counts = Counter(read.get_tag('RG') for read in stream)
            if set(counts) != {x['ID'] for x in groups}:
                raise ValueError('Expected lane contributions absent')
        require_known_mapping_quality(bam)
        result = {'config': config, 'inputs': inputs, 'regions_BED_0based_half_open': regions,
                  'extracted_record_counts_by_RG': dict(counts), 'snvs': [], 'indels': [],
                  'anchor_mapping': [anchor_mapping(bam, t) for t in targets],
                  'original_call_context': {name: call_context(Path(inputs[name]['path']), targets,
                      config['sample'], config['caller_context_padding_bp']) for name in config['caller_inputs']}}
        for t in targets:
            state['current_stage'] = 'count '+t['key']
            save(output/'status.json', state)
            if len(t['ref']) == len(t['alt']) == 1:
                result['snvs'].append({'target': t, 'raw_profile': raw_base_profile(bam, t),
                    'quality_sensitivities': [fragment_evidence(bam, reference, t,
                        {'minimum_mapping_quality': mq, 'minimum_base_quality': 20,
                         'exclude_flags': '0xF04'}) for mq in config['mapping_quality_thresholds']]})
            else:
                result['indels'].append({'target': t, 'window_sensitivities': [
                    window_evidence(bam, reference, t, left, right, mq)
                    for left, right in config['window_flanks_bp']
                    for mq in config['mapping_quality_thresholds']]})
        result['limitations'] = [
            'Selected differences from the same sequencing experiment; not independent biological replication.',
            'No new genotype calls, genotype coercion, clinical classification or causal probability.',
            'Same-position alternate hypotheses and overlapping quality/window sensitivities must not be added.',
            'SNP evidence recomputes BAQ and excludes orphan pairs; exact-window evidence uses raw BQ and retains eligible improper pairs.',
            'A negative exact-window match or absent caller record is not an explicit reference call.',
            'Read support does not independently establish locus assignment in homologous or repetitive regions.',
            'Native unfiltered indexed extraction retains the input alignment fields; no independent full-record extraction comparison is claimed.']
        save(output/'evidence.json', result)
        for name, expected in monitored.items():
            stat = Path(name).stat()
            if (stat.st_size, stat.st_mtime_ns) != expected:
                raise ValueError('Input changed during evidence collection: '+name)
        state.update(status='complete', completed_at=stamp(), inputs=inputs,
                     helper_sha256=config['helper_sha256'],
                     output_sha256={p.name: digest(p) for p in output.iterdir()
                                    if p.is_file() and p.name != 'status.json'})
        save(output/'status.json', state)
        return {'status': 'complete', 'hypotheses': len(targets), 'output': str(output)}
    except BaseException as error:
        state.update(status='failed', failed_at=stamp(), error_type=type(error).__name__, error=str(error))
        save(output/'status.json', state)
        raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('root', 'config', 'output'):
        parser.add_argument('--'+name, required=True, type=Path)
    print(json.dumps(run(parser.parse_args())))
