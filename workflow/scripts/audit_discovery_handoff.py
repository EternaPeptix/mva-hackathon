#!/usr/bin/env python3
"""Audit completed autosomal discovery shards before downstream interpretation.

Read-only on calling outputs. No allele/genotype/quality recoding, filtering,
merging, annotation or ranking is performed. Completed shards are snapshotted;
unfinished chromosomes remain explicit and are never treated as negative results.
"""
import argparse
from collections import Counter
from contextlib import contextmanager
import hashlib
import itertools
import json
import os
from pathlib import Path
import re
import subprocess
import time

import pysam


def stamp():
    return time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())


def save(path, value):
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')
    os.replace(temp, path)


def fingerprint(path, expected=None):
    path = Path(path).resolve()
    before = path.stat()
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(chunk)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError('input changed while hashing: ' + str(path))
    pin = dict(path=str(path), bytes=after.st_size, mtime_ns=after.st_mtime_ns,
               sha256=digest.hexdigest())
    if expected is not None and pin['sha256'] != expected:
        raise ValueError('input checksum mismatch: ' + str(path))
    return pin


def contained(parent, relative):
    path = (parent / relative).resolve()
    if Path(relative).is_absolute() or not path.is_relative_to(parent.resolve()):
        raise ValueError('shard output path escapes its directory')
    return path


def select_shards(run, definition, autosomes, require_complete):
    groups = definition['shards']
    flattened = [pair for group in groups for pair in group]
    names = [pair[0] for pair in flattened]
    if (not groups or any(not group for group in groups) or len(names) != len(set(names))
            or any(not isinstance(n, str) or not isinstance(size, int) or size <= 0
                   for n, size in flattened)
            or definition['contigs'] != len(names)
            or definition['reference_bases'] != sum(size for _, size in flattened)
            or definition['ploidy'] != 2):
        raise ValueError('unsupported or inconsistent discovery definition')
    if len(autosomes) != 22 or len(set(autosomes)) != 22 or not set(autosomes) <= set(names):
        raise ValueError('exactly 22 declared autosomes must occur once in the reference')
    selected, pending = [], []
    for index, group in enumerate(groups):
        group_names = [name for name, _ in group]
        relevant = set(group_names) & set(autosomes)
        if not relevant:
            continue
        if set(group_names) != relevant:
            raise ValueError('autosome/nonautosome mixed shard requires a separate reviewed method')
        receipt = run / ('shard-%04d' % index) / 'complete.json'
        if receipt.is_file():
            selected.append((index, group, receipt))
        else:
            pending.extend(group_names)
    if require_complete and pending:
        raise ValueError('autosomal calling is incomplete')
    if not selected:
        raise ValueError('no completed autosomal shards')
    return selected, pending


def census(reader, sample, contigs, indexed=False, quality=None):
    """Hash every rendered record, plus bounded native-field/missingness counts."""
    if list(reader.header.samples) != [sample]:
        raise ValueError('unexpected sample header')
    lengths = dict(contigs)
    if indexed and reader.index is None:
        raise ValueError('normalized VCF index missing or unreadable')
    counts, filters, ploidies, failures = Counter(records=0), Counter(), Counter(), Counter()
    fields = {key: Counter() for key in ('DP', 'GQ', 'AD', 'PL')}
    payload = hashlib.sha256()
    records = (itertools.chain.from_iterable(reader.fetch(name, 0, length) for name, length in contigs)
               if indexed else reader)
    for record in records:
        if record.contig not in lengths or not 1 <= record.pos <= lengths[record.contig]:
            raise ValueError('record outside declared complete contig')
        payload.update(str(record).encode())
        counts['records'] += 1
        call = record.samples[sample]
        alts = record.alts or ()
        gt = tuple(call.get('GT') or ())
        if any(x is not None and not 0 <= x <= len(alts) for x in gt):
            raise ValueError('invalid original genotype allele index')
        ploidies[str(len(gt))] += 1
        counts['multiallelic_records'] += len(alts) > 1
        counts['incomplete_GT_records'] += not gt or None in gt
        counts['two_distinct_nonreference_diploid_GT_records'] += (
            len(gt) == 2 and None not in gt and 0 not in gt and gt[0] != gt[1])
        filters[';'.join(record.filter) or '.'] += 1
        for field in fields:
            value = call.get(field)
            present = value is not None and (not isinstance(value, tuple) or any(x is not None for x in value))
            fields[field]['present' if present else 'missing'] += 1
            if isinstance(value, tuple) and any(x is None for x in value) and present:
                fields[field]['partially_missing'] += 1
        counts['retained_ALT_alleles'] += len(alts)
        called = {x for x in gt if x is not None and x > 0}
        counts['GT_positive_ALT_alleles'] += len(called)
        counts['retained_ALT_alleles_absent_from_GT'] += len(alts) - len(called)
        concrete = {i for i, alt in enumerate(alts, 1)
                    if re.fullmatch('[ACGT]+', record.ref) and re.fullmatch('[ACGT]+', alt)}
        counts['GT_positive_concrete_ALT_alleles'] += len(called & concrete)
        counts['GT_positive_nonconcrete_ALT_alleles'] += len(called - concrete)
        if quality is not None:
            reasons = []
            if not called:
                reasons.append('no_called_alternate')
            if not gt or None in gt:
                reasons.append('partial_or_missing_genotype')
            if len(gt) != 2:
                reasons.append('not_diploid')
            if list(record.filter) != ['PASS']:
                reasons.append('record_not_explicit_PASS')
            for field, threshold, label in (
                    ('DP', quality['minimum_DP'], 'depth'),
                    ('GQ', quality['minimum_GQ'], 'genotype_quality')):
                value = call.get(field)
                if value is None or value < threshold:
                    reasons.append(label + '_missing_or_below_threshold')
            failures.update(reasons)
            counts['records_passing_existing_preannotation_quality_gate'] += not reasons
            counts['GT_positive_concrete_ALT_alleles_passing_existing_preannotation_quality_gate'] += (
                len(called & concrete) if not reasons else 0)
    return dict(counts=dict(counts), filters=dict(filters), called_ploidies=dict(ploidies),
                FORMAT_availability={key: dict(value) for key, value in fields.items()},
                existing_gate_nonexclusive_record_failure_counts=dict(failures),
                rendered_record_payload_sha256=payload.hexdigest())


@contextmanager
def normalized_replay(command, log):
    with log.open('wb') as stderr:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=stderr)
        try:
            with pysam.VariantFile(process.stdout) as reader:
                yield reader
            process.stdout.close()
            if process.wait() != 0:
                raise RuntimeError('native normalization replay failed; see ' + str(log))
        finally:
            if process.poll() is None:
                process.terminate()
                try:
                    process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()


def audit(args):
    run, output = args.run.resolve(), args.output.resolve()
    reference, bcftools = args.reference.resolve(), args.bcftools.resolve()
    output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    state = dict(status='running', started_at=stamp(), pid=os.getpid(), completed_shard_ids=[])
    save(output / 'status.json', state)
    try:
        pins = {}
        def pin(label, path, expected=None):
            value = fingerprint(path, expected)
            pins[label] = value
            return value
        pin('definition', run / 'run-definition.json', args.expected_definition_sha256)
        definition = json.loads((run / 'run-definition.json').read_text())
        config_pin = pin('existing_quality_config', args.quality_config)
        config = json.loads(Path(config_pin['path']).read_text())
        quality = config['quality']
        if quality.get('explicit_PASS') is not True or quality.get('fully_called_diploid') is not True:
            raise ValueError('unsupported existing quality contract')
        if any(not isinstance(quality.get(key), int) or quality[key] < 0
               for key in ('minimum_DP', 'minimum_GQ')):
            raise ValueError('invalid existing quality thresholds')
        if config['sample'] != definition['sample']:
            raise ValueError('quality configuration sample differs from calling sample')
        selected, pending = select_shards(run, definition, args.autosomes.split(','), args.require_complete)
        pin('reference', reference, definition['reference_sha256'])
        pin('reference_index', Path(str(reference) + '.fai'))
        with pysam.FastaFile(str(reference)) as fasta:
            if list(zip(fasta.references, fasta.lengths)) != [tuple(pair) for group in definition['shards'] for pair in group]:
                raise ValueError('reference index differs from discovery definition')
        pin('bcftools', bcftools)
        pin('audit_source', Path(__file__))
        # Freeze receipt membership before work; later completions belong to a new audit.
        frozen = []
        for index, contigs, receipt in selected:
            pin('shard-%04d/complete' % index, receipt)
            frozen.append((index, contigs, receipt.parent, json.loads(receipt.read_text())))
        results = []
        for index, contigs, shard, completed in frozen:
            methods = {}
            for method, raw_name in (('gatk', 'gatk.raw.vcf.gz'), ('bcftools', 'bcftools.raw.bcf')):
                label = 'shard-%04d/%s' % (index, method)
                normalized = contained(shard, completed[method]['path'])
                attempt = normalized.parent
                pin(label + '/normalized', normalized, completed[method]['sha256'])
                pin(label + '/normalized_index', Path(str(normalized) + '.tbi'))
                command_file = attempt / 'commands.json'
                pin(label + '/commands', command_file, completed['commands_sha256'])
                source = attempt / raw_name
                expected_command = ['bcftools', 'norm', '--check-ref', 'e', '-f', str(reference),
                                    '-m', '-any', '-Oz', '-o', str(normalized), str(source)]
                if json.loads(command_file.read_text()).count(expected_command) != 1:
                    raise ValueError('unrecognized frozen normalization command')
                raw_pin = pin(label + '/original_raw_newly_fingerprinted', source)
                with pysam.VariantFile(str(source)) as reader:
                    raw = census(reader, definition['sample'], contigs, quality=quality)
                with pysam.VariantFile(str(normalized)) as reader:
                    recorded = census(reader, definition['sample'], contigs, indexed=True)
                if recorded['counts']['records'] != completed[method]['counts']['records']:
                    raise ValueError('completed shard record count mismatch')
                replay_command = [str(bcftools), 'norm', '--check-ref', 'e', '-f', str(reference),
                                  '-m', '-any', '-Ov', str(source)]
                with normalized_replay(replay_command, output / ('shard-%04d-%s-replay.log' % (index, method))) as reader:
                    replay = census(reader, definition['sample'], contigs)
                if (replay['counts']['records'], replay['rendered_record_payload_sha256']) != (
                        recorded['counts']['records'], recorded['rendered_record_payload_sha256']):
                    raise ValueError('original raw source does not reproduce frozen split record payloads')
                methods[method] = dict(original_raw=raw, original_raw_pin=raw_pin,
                    normalized_split=recorded, replayed_split_payload_identical=True,
                    replay_command=replay_command,
                    raw_provenance='Newly fingerprinted retained raw source; verified by complete-record split replay against the previously hashed normalized shard.')
            results.append(dict(shard_id=index, contigs=contigs, methods=methods))
            state['completed_shard_ids'].append(index)
            save(output / 'status.json', state)
            print(json.dumps(dict(verified_shard=index, completed=len(results), selected=len(frozen))), flush=True)
        for label, value in pins.items():
            current = Path(value['path']).stat()
            if (current.st_size, current.st_mtime_ns) != (value['bytes'], value['mtime_ns']):
                raise ValueError('input changed during handoff audit: ' + label)
        summary = dict(status='complete', completed_at=stamp(), elapsed_seconds=time.monotonic() - started,
            scope='Snapshot of completed autosomal raw discovery shards; no calling, rewriting, ranking or clinical inference.',
            sample=definition['sample'], selected_autosomes=[c for _, cs, _ in selected for c, _ in cs],
            pending_autosomes=pending, all_autosomes_complete=not pending,
            existing_preannotation_quality_gate=quality, inputs=pins, shards=results,
            runtime=dict(python=os.sys.version, pysam=pysam.__version__,
                         bcftools_version=subprocess.check_output([str(bcftools), '--version'], text=True).strip()),
            limitations=[
                'Only completed autosomal shards in the initial snapshot are audited; absence of a shard/record is not reference evidence.',
                'Raw allele lists, GT, AD and PL remain authoritative. Normalized split rows alone can conceal incompatible second alleles.',
                'The existing quality gate is counted before annotation on all raw records. These are not numbers of causal or coding candidates.',
                'FILTER dot means not evaluated; it is not PASS. Missing GQ is not filled with zero or inferred from PL.',
                'No new quality policy, genotype threshold, caller concordance or callability estimate is established.',
                'Replay uses the installed bcftools and decoder. It validates handoff consistency, not independent biological correctness.',
                'Raw-file identities are recorded now; only normalized files and command histories had prior completed-shard hashes.',
                'The audit does not replace the separate X/Y, mitochondrial, CNV, structural or mosaicism analyses.'])
        save(output / 'summary.json', summary)
        state.update(status='complete', completed_at=summary['completed_at'], summary_sha256=fingerprint(output / 'summary.json')['sha256'])
        save(output / 'status.json', state)
        return summary
    except BaseException as error:
        state.update(status='failed', failed_at=stamp(), error_type=type(error).__name__, error=str(error))
        save(output / 'status.json', state)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('run', 'reference', 'bcftools', 'quality-config', 'output'):
        parser.add_argument('--' + key, type=Path, required=True)
    parser.add_argument('--expected-definition-sha256', required=True)
    parser.add_argument('--autosomes', default=','.join(map(str, range(1, 23))))
    parser.add_argument('--require-complete', action='store_true')
    args = parser.parse_args()
    os.nice(10)
    audit(args)


if __name__ == '__main__':
    main()
