#!/usr/bin/env python3
"""Compare literal local sequence hypotheses with original mapped fragments.

Constructed sequences and their annotations are hypotheses, never native calls.
Caller lineages sharing one sequence never multiply the observed fragment count.
"""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time

import pysam

from native_coding_records import source_context
from joint_haplotype_helpers import (enumerate_haplotypes, minimal_replacement,
                                    aggregate_fragment_evidence)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def stamp():
    return time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())


def save(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')


def fingerprint(path, expected=None):
    path = Path(path).resolve()
    before = path.stat()
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(block)
    after = path.stat()
    require((before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns),
            'input changed while hashing: ' + str(path))
    result = dict(path=str(path), bytes=after.st_size, mtime_ns=after.st_mtime_ns,
                  sha256=digest.hexdigest())
    if expected:
        require(result['sha256'] == expected['sha256'], 'input SHA differs: ' + str(path))
        if 'bytes' in expected:
            require(result['bytes'] == expected['bytes'], 'input size differs: ' + str(path))
    return result


def collect_native_context(config, pins, targets):
    """Scan entire raw sources: no index mutation or missed upstream REF span."""
    selected = {target['id']: [] for target in targets}
    census = []
    flank = max(config['window_flanks_bp'])
    for source in config['caller_sources']:
        pin = pins[source['input']]
        payload, count = hashlib.sha256(), 0
        with pysam.VariantFile(pin['path']) as reader:
            require(list(reader.header.samples) == [config['sample']], 'native sample mismatch')
            for ordinal, record in enumerate(reader, 1):
                count += 1
                payload.update(str(record).encode())
                begin, end = record.pos - 1, max(record.stop, record.pos - 1 + len(record.ref))
                for target in targets:
                    if (record.contig != target['chrom'] or begin >= target['core_end0'] + flank
                            or end <= target['core_start0'] - flank):
                        continue
                    context = source_context(record, config['sample'], source['caller'], pin['sha256'], ordinal)
                    selected[target['id']].append(dict(context=context, reported_end0=record.stop))
        require(count == source['expected_records'], 'native source census count differs')
        require(payload.hexdigest() == source['expected_rendered_record_payload_sha256'],
                'native source rendered payload differs')
        census.append(dict(input=source['input'], caller=source['caller'], records=count,
                           rendered_record_payload_sha256=payload.hexdigest()))
    for target in targets:
        contexts = selected[target['id']]
        identities = [row['context']['source_record_id'] for row in contexts]
        require(len(identities) == len(set(identities)), 'duplicate source record in selected context')
        require(set(target.get('required_source_record_ids', [])) <= set(identities),
                'selected target source record absent from regional census')
    return selected, census


def window_edits(contexts, start, end):
    edits, omitted = [], []
    for entry in contexts:
        context = entry['context']
        begin, stop = context['pos'] - 1, context['pos'] - 1 + len(context['ref'])
        if begin >= end or stop <= start:
            continue
        for index, alt in enumerate(context['alts'], 1):
            item = dict(caller=context['caller'], source_record_id=context['source_record_id'],
                source_alt_index=index, chrom=context['chrom'], pos=context['pos'],
                ref=context['ref'], alt=alt, called=index in (context['FORMAT'].get('GT') or []))
            reason = ('nonconcrete_ALT_or_REF' if not all(re.fullmatch('[ACGT]+', value or '')
                      for value in (context['ref'], alt)) else 'crosses_window_boundary'
                      if begin < start or stop > end else None)
            if reason:
                omitted.append(dict(item, reason=reason))
            else:
                edits.append(item)
    return edits, omitted


def alignment_payload(reader, regions=None):
    digest, count = hashlib.sha256(), 0
    # Consume each indexed iterator before creating the next one on this handle.
    for region in [None] if regions is None else regions:
        stream = reader if region is None else reader.fetch(*region)
        for record in stream:
            digest.update((record.to_string() + '\n').encode())
            count += 1
    return dict(records=count, rendered_record_payload_sha256=digest.hexdigest())


def nonconsequence_payload(path):
    digest, ids = hashlib.sha256(), []
    with pysam.VariantFile(str(path)) as reader:
        require(not list(reader.header.samples), 'hypothesis VCF must be sites-only')
        for record in reader:
            fields = str(record).rstrip('\n').split('\t')
            fields[7] = ';'.join(x for x in fields[7].split(';') if x.split('=', 1)[0] != 'BCSQ') or '.'
            digest.update(('\t'.join(fields) + '\n').encode())
            ids.append(record.id)
    return dict(ids=ids, payload_sha256=digest.hexdigest())


def annotate_hypotheses(hypotheses, reference, gff, bcftools, output, command):
    """Sites-only constructed alleles; REF hypotheses have no variant record."""
    vcfin = output / 'constructed-hypotheses.vcf'
    with pysam.FastaFile(str(reference)) as fasta:
        header = pysam.VariantHeader()
        header.add_meta('hypothesis_scope', value='Constructed sequence edits, not observed genotype calls or clinical findings')
        for chrom, length in zip(fasta.references, fasta.lengths):
            header.contigs.add(chrom, length=length)
        order = {name: i for i, name in enumerate(fasta.references)}
        ordered = sorted(hypotheses.items(), key=lambda item: (
            order[item[1]['chrom']], (item[1]['constructed_variant'] or {'pos': item[1]['start0'] + 1})['pos'], item[0]))
        with pysam.VariantFile(str(vcfin), 'w', header=header) as writer:
            for identifier, row in ordered:
                variant = row['constructed_variant']
                if variant is None:
                    row['joint_consequences'] = []
                    row['annotation_status'] = 'reference_sequence_no_variant_record'
                    continue
                record = writer.new_record(contig=variant['chrom'], start=variant['pos'] - 1,
                    alleles=(variant['ref'], variant['alt']), id=identifier)
                writer.write(record)
    before_sort = output / 'constructed-normalized-before-sort.bcf'
    normalized, annotated = output / 'constructed-normalized.bcf', output / 'constructed-annotated.bcf'
    command([bcftools, 'norm', '--check-ref', 'e', '-f', reference, '-Ob', '-o', before_sort, vcfin], 'normalize-hypotheses')
    command([bcftools, 'sort', '-Ob', '-o', normalized, before_sort], 'sort-hypotheses')
    command([bcftools, 'csq', '--local-csq', '--samples', '-', '-f', reference, '-g', gff,
             '-Ob', '-o', annotated, normalized], 'annotate-hypotheses')
    before, after = nonconsequence_payload(normalized), nonconsequence_payload(annotated)
    require(before == after, 'hypothesis annotation dropped or changed a nonconsequence field')
    expected = {identifier for identifier, row in hypotheses.items() if row['constructed_variant'] is not None}
    require(len(before['ids']) == len(set(before['ids'])) and set(before['ids']) == expected,
            'hypothesis record ID membership changed during normalization')
    with pysam.FastaFile(str(reference)) as fasta, pysam.VariantFile(str(annotated)) as reader:
        for record in reader:
            row = hypotheses[record.id]
            require(record.contig == row['chrom'] and len(record.alts) == 1, 'hypothesis chromosome/allele shape changed')
            start = min(record.pos - 1, row['start0'])
            end = max(record.pos - 1 + len(record.ref), row['end0'])
            sequence = fasta.fetch(record.contig, start, end).upper()
            offset = record.pos - 1 - start
            require(sequence[offset:offset + len(record.ref)] == record.ref, 'normalized hypothesis REF mismatch')
            actual = sequence[:offset] + record.alts[0] + sequence[offset + len(record.ref):]
            wanted = sequence[:row['start0'] - start] + row['sequence'] + sequence[row['end0'] - start:]
            require(actual == wanted, 'normalized hypothesis changes its full sequence')
            row['normalized_constructed_variant'] = dict(chrom=record.contig, pos=record.pos,
                ref=record.ref, alt=record.alts[0])
            row['joint_consequences'] = list(record.info.get('BCSQ', ()))
            row['annotation_status'] = 'local_constructed_hypothesis_annotation'
    return dict(constructed_variant_records=len(expected), annotation_payload_preserved=True,
                every_normalized_constructed_allele_sequence_verified=True)


def run(args):
    root, output = args.root.resolve(), args.output.resolve()
    require(not output.exists() and not output.is_relative_to(Path(__file__).resolve().parents[2]),
            'use a new output directory outside the source tree')
    output.mkdir(parents=True)
    os.nice(10)
    started = time.monotonic()
    state = dict(status='running', pid=os.getpid(), started_at=stamp(), stage='verify inputs')
    save(output / 'status.json', state)
    pins, commands = {}, []
    try:
        config = json.loads(args.config.read_text())
        pins['config'] = fingerprint(args.config)
        for name, expected in config['inputs'].items():
            pins[name] = fingerprint(root / expected['path'], expected)
        for name in ('review_joint_haplotypes.py', 'joint_haplotype_helpers.py', 'native_coding_records.py'):
            pins['source/' + name] = fingerprint(Path(__file__).with_name(name),
                {'sha256': config['source_sha256'][name]})
        require(config['raw_base_qualities_no_BAQ'] is True, 'unsupported base-quality interpretation')
        require(config['window_flanks_bp'] and all(type(x) is int and x > 0 for x in config['window_flanks_bp']),
                'invalid positive window flanks')
        reference, alignment = (Path(pins[name]['path']) for name in ('reference', 'alignment'))
        samtools, bcftools, gff = (Path(pins[name]['path']) for name in ('samtools', 'bcftools', 'gff'))
        targets = config['targets']
        require(targets and len({x['id'] for x in targets}) == len(targets), 'empty or duplicate target selection')
        callers = sorted({row['caller'] for row in config['caller_sources']})
        require(callers, 'caller sources required')
        regions = []
        with pysam.FastaFile(str(reference)) as fasta:
            order = {name: i for i, name in enumerate(fasta.references)}
            for target in targets:
                require(re.fullmatch('[A-Za-z0-9_.-]+', target['id']) is not None, 'invalid target identifier')
                require(0 <= target['core_start0'] < target['core_end0'] <= fasta.get_reference_length(target['chrom']),
                        'invalid target core interval')
                flank = max(config['window_flanks_bp'])
                require(target['core_start0'] >= flank and target['core_end0'] + flank <= fasta.get_reference_length(target['chrom']),
                        'target lacks reference flanks')
                regions.append([target['chrom'], target['core_start0'] - flank, target['core_end0'] + flank])
        regions.sort(key=lambda x: (order[x[0]], x[1], x[2]))
        merged = []
        for chrom, begin, end in regions:
            if merged and merged[-1][0] == chrom and begin <= merged[-1][2]:
                merged[-1][2] = max(end, merged[-1][2])
            else:
                merged.append([chrom, begin, end])
        state['stage'] = 'recover complete nearby native context'
        save(output / 'status.json', state)
        contexts, census = collect_native_context(config, pins, targets)
        save(output / 'regional-native-context.json', dict(targets=contexts, source_census=census))

        def command(argv, label):
            state['stage'] = label
            save(output / 'status.json', state)
            argv = list(map(str, argv))
            entry = dict(command=argv, label=label, started_at=stamp())
            commands.append(entry)
            save(output / 'commands.json', commands)
            with (output / (label + '.log')).open('wb') as log:
                result = subprocess.run(argv, stdout=log, stderr=subprocess.STDOUT)
            entry.update(exit_code=result.returncode, completed_at=stamp())
            save(output / 'commands.json', commands)
            require(result.returncode == 0, 'native command failed: ' + label)

        bed, bam = output / 'regions.bed', output / 'original-regional-reads.bam'
        bed.write_text(''.join('%s\t%d\t%d\n' % tuple(region) for region in merged))
        command([samtools, 'view', '-@', '1', '--no-PG', '-T', reference, '-b', '-M', '-L', bed,
                 '--save-counts', output / 'extraction-counts.json', '-o', bam, alignment], 'extract-original-reads')
        command([samtools, 'quickcheck', '-v', bam], 'quickcheck')
        command([samtools, 'index', bam], 'index')
        with pysam.AlignmentFile(str(alignment), reference_filename=str(reference)) as original, pysam.AlignmentFile(str(bam), 'rb') as extracted:
            require(original.header.to_dict()['SQ'] == extracted.header.to_dict()['SQ'] and
                    original.header.to_dict().get('RG') == extracted.header.to_dict().get('RG'), 'extracted reference or read-group header changed')
            groups = {row['ID']: row for row in original.header.to_dict().get('RG', [])}
            require(len(groups) == config['expected_read_groups'] and
                    {row.get('SM') for row in groups.values()} == {config['sample']}, 'unexpected original sample/read groups')
            before, after = alignment_payload(original, merged), alignment_payload(extracted)
            require(before == after, 'decoded original and extracted ordered records differ')
        hypotheses, windows = {}, []
        with pysam.FastaFile(str(reference)) as fasta, pysam.AlignmentFile(str(bam), 'rb') as reads:
            for target in targets:
                state['stage'] = 'compare literal haplotypes for ' + target['id']
                save(output / 'status.json', state)
                for flank in config['window_flanks_bp']:
                    begin, end = target['core_start0'] - flank, target['core_end0'] + flank
                    sequence = fasta.fetch(target['chrom'], begin, end).upper()
                    edits, omitted = window_edits(contexts[target['id']], begin, end)
                    family = enumerate_haplotypes(target['chrom'], begin, sequence, edits,
                        callers=callers, max_states=config['maximum_hypothesis_states'],
                        max_haplotype_bases=config['maximum_haplotype_bases'])
                    identifiers, expected_sequences = [], {}
                    for index, hypothesis in enumerate(family['haplotypes']):
                        identifier = '%s_F%d_H%04d' % (target['id'], flank, index)
                        variant = minimal_replacement(target['chrom'], begin, sequence, hypothesis['sequence'])
                        hypotheses[identifier] = dict(hypothesis, id=identifier, target_id=target['id'],
                            chrom=target['chrom'], start0=begin, end0=end, constructed_variant=variant)
                        identifiers.append(identifier)
                        require(hypothesis['sequence'] not in expected_sequences, 'duplicate enumerated sequence class')
                        expected_sequences[hypothesis['sequence']] = [identifier]
                    observed = list(reads.fetch(target['chrom'], begin, end))
                    require(all(read.has_tag('RG') and read.get_tag('RG') in groups for read in observed),
                            'unknown or missing original read group')
                    arms = [dict(minimum_mapping_quality=mq, minimum_base_quality=bq,
                            evidence=aggregate_fragment_evidence(observed, begin, end, bq, mq,
                                expected_sequences=expected_sequences))
                            for mq in config['mapping_quality_thresholds'] for bq in config['base_quality_thresholds']]
                    windows.append(dict(target_id=target['id'], chrom=target['chrom'], start0=begin, end0=end,
                        flank_bp=flank, reference_sequence=sequence, included_native_edits=edits,
                        omitted_overlapping_native_alleles=omitted, hypothesis_ids=identifiers,
                        enumeration_metadata={key: value for key, value in family.items() if key != 'haplotypes'},
                        original_overlapping_read_records=len(observed), arms=arms))
        annotation = annotate_hypotheses(hypotheses, reference, gff, bcftools, output, command)
        evidence = dict(scope='Post hoc local joint-sequence comparison of selected original native variants; no genotyper or parental phasing.',
            config=config, inputs=pins, targets=targets, native_source_census=census,
            original_extraction_payload_check=before, regions_BED_0based_half_open=merged,
            hypotheses=hypotheses, windows=windows, annotation_check=annotation,
            limitations=[
                'Caller combinations include all concrete alternatives, even uncalled ones, as hypothetical sequences only; they are not native genotype or haplotype calls.',
                'Identical sequences have multiple caller/source lineages but one observed sequence/fragment count. Sequence algebra is not biological phase.',
                'Exact spanning sequences use raw base qualities without BAQ, exclude flags0xF04 and unknown MQ255, and retain eligible improper pairs.',
                'Fragments are grouped by original read group and query name. Window/quality arms reuse the same molecules and must not be added.',
                'Unmatched sequences, conflicting mates, missing span and omitted nonconcrete/window-crossing hypotheses remain explicit; absence is not a reference call.',
                'Local mapping/read support does not independently establish locus origin, gene-level or parental phase, loss of function, diagnosis, EPCR or treatment benefit.',
                'Annotations describe constructed combined alleles. Transcript relevance, full coding haplotypes and disease-gene interpretation need separate review.'])
        save(output / 'evidence.json', evidence)
        for name, pin in pins.items():
            stat = Path(pin['path']).stat()
            require((stat.st_size, stat.st_mtime_ns) == (pin['bytes'], pin['mtime_ns']), 'input changed during review: ' + name)
        summary = dict(status='complete', completed_at=stamp(), elapsed_seconds=time.monotonic() - started,
            inputs=pins, selected_targets=[x['id'] for x in targets], windows=len(windows),
            quality_arms=sum(len(x['arms']) for x in windows), sequence_hypotheses=len(hypotheses),
            extraction_payload_check=before, annotation_check=annotation, limitations=evidence['limitations'],
            outputs={p.name: fingerprint(p) for p in sorted(output.iterdir()) if p.is_file() and p.name not in ('status.json', 'summary.json')},
            runtime=dict(python=os.sys.version, pysam=pysam.__version__))
        save(output / 'summary.json', summary)
        state.update(status='complete', completed_at=summary['completed_at'], stage='completed',
                     summary_sha256=fingerprint(output / 'summary.json')['sha256'])
        save(output / 'status.json', state)
        return summary
    except BaseException as error:
        state.update(status='failed', failed_at=stamp(), error_type=type(error).__name__, error=str(error))
        save(output / 'status.json', state)
        raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('root', 'config', 'output'):
        parser.add_argument('--' + name, required=True, type=Path)
    run(parser.parse_args())
