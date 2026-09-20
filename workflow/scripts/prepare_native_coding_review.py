#!/usr/bin/env python3
"""Prepare caller-separated coding inventories from a verified handoff snapshot.

Original source ordinal and full native genotype are authoritative. All concrete
protein-altering annotation rows are retained regardless of FILTER/DP/GQ, with
descriptive missingness and structural-model restrictions kept explicit.
"""
import argparse
from collections import Counter, defaultdict
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

import pysam

from audit_discovery_handoff import fingerprint, save, stamp
from rare_biallelic_baseline import effects_by_gene, frequency_state
from native_coding_records import source_context, original_alt_index, allele_call, review_axes, build_models


ORDINAL = 'NC_SOURCE_ORD'
ORIGIN = 'NC_ORIG'


def write_jsonl(path, rows):
    with gzip.open(path, 'wt') as stream:
        for row in rows:
            stream.write(json.dumps(row, sort_keys=True, allow_nan=False) + '\n')


def tag_source(source, destination, sample, contigs, expected_census):
    """Only a derived copy receives an ordinal; full original payload is checked."""
    payload = hashlib.sha256()
    count = 0
    with pysam.VariantFile(str(source)) as reader:
        if list(reader.header.samples) != [sample]:
            raise ValueError('native source sample differs from handoff')
        if ORDINAL in reader.header.info or ORIGIN in reader.header.info or 'BCSQ' in reader.header.info or 'BCSQ' in reader.header.formats:
            raise ValueError('reserved source-origin INFO tag already exists')
        original_header = str(reader.header)
        header = reader.header.copy()
        header.info.add(ORDINAL, 1, 'Integer', 'One-based record ordinal in the pinned native source; not a genotype or quality field')
        with pysam.VariantFile(str(destination), 'wb', header=header) as writer:
            for record in reader:
                if record.contig not in contigs:
                    raise ValueError('native source contains an out-of-snapshot contig')
                count += 1
                payload.update(str(record).encode())
                copied = record.copy()
                copied.translate(header)
                copied.info[ORDINAL] = count
                writer.write(copied)
    if (count, payload.hexdigest()) != (expected_census['counts']['records'], expected_census['rendered_record_payload_sha256']):
        raise ValueError('native source payload differs from the completed handoff')
    replay, copied_count = hashlib.sha256(), 0
    with pysam.VariantFile(str(destination)) as reader:
        for record in reader:
            copied_count += 1
            if record.info[ORDINAL] != copied_count:
                raise ValueError('derived original ordinal changed')
            del record.info[ORDINAL]
            replay.update(str(record).encode())
    if copied_count != count or replay.hexdigest() != payload.hexdigest():
        raise ValueError('derived source copy changed original rendered record payload')
    return dict(records=count, rendered_record_payload_sha256=payload.hexdigest(),
                derived_copy_payload_preserved=True), original_header


def normalize_annotate(bcftools, reference, gff, tagged, output, directory):
    normalized = directory / 'normalized-before-csq.bcf'
    commands = [
        [str(bcftools), 'norm', '--check-ref', 'e', '-f', str(reference), '-m', '-any',
         '--old-rec-tag', ORIGIN, '-Ob', '-o', str(normalized), str(tagged)],
        [str(bcftools), 'csq', '--local-csq', '--samples', '-', '-f', str(reference),
         '-g', str(gff), '-Ob', '-o', str(output), str(normalized)],
    ]
    for index, command in enumerate(commands):
        with (directory / ('annotation-%d.log' % index)).open('wb') as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
    return commands, normalized


def annotation_record_payload(path, sample, remove_consequence=False):
    """Only BCSQ annotations may differ; compare all other rendered fields."""
    count, payload = 0, hashlib.sha256()
    with pysam.VariantFile(str(path)) as reader:
        if list(reader.header.samples) != [sample]:
            raise ValueError('sample changed while annotating records')
        for record in reader:
            parts = str(record).rstrip('\n').split('\t')
            if remove_consequence:
                parts[7] = ';'.join(value for value in parts[7].split(';')
                                    if value.split('=', 1)[0] != 'BCSQ') or '.'
                fields = parts[8].split(':')
                if 'BCSQ' in fields:
                    index = fields.index('BCSQ')
                    parts[8] = ':'.join(value for i, value in enumerate(fields) if i != index)
                    for column in range(9, len(parts)):
                        parts[column] = ':'.join(value for i, value in enumerate(parts[column].split(':')) if i != index)
            payload.update(('\t'.join(parts) + '\n').encode())
            count += 1
    return dict(records=count, rendered_nonconsequence_record_payload_sha256=payload.hexdigest())


def verify_allele_haplotype(fasta, normalized, call):
    """Check an allele's exact sequence edit, independently of the origin tag."""
    original = call['source_allele']
    if normalized.contig != original['chrom']:
        raise ValueError('normalized/source chromosome mismatch')
    start = min(normalized.pos, original['pos']) - 1
    end = max(normalized.pos - 1 + len(normalized.ref), original['pos'] - 1 + len(original['ref']))
    sequence = fasta.fetch(normalized.contig, start, end).upper()
    if len(sequence) != end - start:
        raise ValueError('source/normalized allele extends beyond reference')
    def apply(position, reference, alternate):
        offset = position - 1 - start
        if sequence[offset:offset + len(reference)] != reference:
            raise ValueError('source/normalized REF mismatch')
        return sequence[:offset] + alternate + sequence[offset + len(reference):]
    if apply(original['pos'], original['ref'], original['alt']) != apply(normalized.pos, normalized.ref, normalized.alts[0]):
        raise ValueError('normalized allele is not sequence-equivalent to the indexed original ALT')


def prepare_shard(source, source_pin, shard, method, sample, config, bcftools, reference, gff, directory):
    directory.mkdir()
    contigs = [name for name, _ in shard['contigs']]
    expected = shard['methods'][method]['original_raw']
    tagged = directory / 'ordinal-tagged-original.bcf'
    copied, header = tag_source(source, tagged, sample, contigs, expected)
    (directory / 'original-header.vcf.txt').write_text(header)
    annotated = directory / 'alleles-local-csq.bcf'
    commands, normalized = normalize_annotate(bcftools, reference, gff, tagged, annotated, directory)
    before_annotation = annotation_record_payload(normalized, sample)
    after_annotation = annotation_record_payload(annotated, sample, remove_consequence=True)
    if before_annotation['records'] != shard['methods'][method]['normalized_split']['counts']['records']:
        raise ValueError('normalized record count differs from completed handoff')
    if before_annotation != after_annotation:
        raise ValueError('annotation dropped or changed original nonconsequence record fields')
    counters, selected, nonconcrete = Counter(), defaultdict(list), []
    with pysam.VariantFile(str(annotated)) as reader:
        if list(reader.header.samples) != [sample]:
            raise ValueError('annotation output sample changed')
        if 'BCSQ' not in reader.header.info or 'Local consequence' not in reader.header.info['BCSQ'].description:
            raise ValueError('expected local, not inferred haplotype, consequences')
        for record in reader:
            counters['normalized_allele_records'] += 1
            ordinal = record.info.get(ORDINAL)
            if (not isinstance(ordinal, int) or not 1 <= ordinal <= copied['records']
                    or record.contig not in contigs or len(record.alts or ()) != 1):
                raise ValueError('annotation has invalid original ordinal, contig or allele shape')
            if any(not allele or any(base not in 'ACGT' for base in allele) for allele in record.alleles):
                counters['nonconcrete_normalized_allele_records'] += 1
                nonconcrete.append(dict(source_ordinal=ordinal, normalized_record=str(record),
                                        reason='nonconcrete allele retained for representation review'))
                continue
            if 'BCSQ' not in record.info:
                counters['concrete_records_without_BCSQ'] += 1
                continue
            genes = effects_by_gene(record.info['BCSQ'], config)
            if not genes:
                counters['concrete_annotated_records_without_selected_protein_effect'] += 1
                continue
            counters['protein_altering_allele_records'] += 1
            selected[ordinal].append((record.copy(), genes))
    # Recover originals by immutable ordinal in a second sequential source pass.
    # This avoids an index mutation and disambiguates same-position duplicate calls.
    variants, recovered, seen_ids = [], set(), set()
    with pysam.VariantFile(str(source)) as reader, pysam.FastaFile(str(reference)) as fasta:
        for ordinal, record in enumerate(reader, 1):
            if ordinal not in selected:
                continue
            context = source_context(record, sample, method, source_pin['sha256'], ordinal)
            recovered.add(ordinal)
            for normalized, genes in selected[ordinal]:
                index = original_alt_index(normalized, context, origin_tag=ORIGIN)
                call = allele_call(context, index)
                verify_allele_haplotype(fasta, normalized, call)
                key = '%s:%d:%s:%s' % (normalized.contig, normalized.pos, normalized.ref, normalized.alts[0])
                identity = (call['source_record_id'], index, key)
                if identity in seen_ids:
                    raise ValueError('duplicate source-allele/normalized-row identity from annotation')
                seen_ids.add(identity)
                variants.append(dict(key=key, allele=dict(chrom=normalized.contig, pos=normalized.pos,
                    ref=normalized.ref, alt=normalized.alts[0]), source_call=call, genes=genes,
                    normalized_record=str(normalized), shard_id=shard['shard_id']))
    if recovered != set(selected) or len(variants) != counters['protein_altering_allele_records']:
        raise ValueError('selected coding rows do not reconcile with recovered original records')
    write_jsonl(directory / 'nonconcrete-representation-review.jsonl.gz', nonconcrete)
    return variants, dict(shard_id=shard['shard_id'], contigs=shard['contigs'], source=source_pin,
        original_copy_check=copied, annotation_preservation_check=before_annotation,
        counters=dict(counters), selected_original_records=len(recovered),
        commands=commands, scope='All selected protein-altering concrete allele rows, including uncalled and quality-missing records; original genotypes recovered by source ordinal.')


def attach_frequencies(variants, response_file):
    keys = {variant['key'] for variant in variants}
    rows = {}
    with response_file.open() as stream:
        for line in stream:
            row = json.loads(line)
            if row['key'] not in keys or row['key'] in rows:
                raise ValueError('unexpected or duplicate frequency response')
            rows[row['key']] = frequency_state(row)
    if set(rows) != keys:
        raise ValueError('missing frequency response; absence is not an observed zero')
    for variant in variants:
        variant['frequency'] = rows[variant['key']]
    return len(rows)


def run(args):
    root, output = args.root.resolve(), args.output.resolve()
    reference, gff = args.reference.resolve(), args.gff.resolve()
    output.mkdir(parents=True, exist_ok=False)
    state = dict(status='running', started_at=stamp(), pid=os.getpid(), caller=args.caller, stage='verify snapshot inputs')
    save(output / 'status.json', state)
    started = time.monotonic()
    pins = {}
    def pin(name, path, expected=None):
        value = fingerprint(path, expected)
        pins[name] = value
        return value
    try:
        pin('handoff', args.handoff, args.expected_handoff_sha256)
        handoff = json.loads(args.handoff.read_text())
        if handoff['status'] != 'complete':
            raise ValueError('completed handoff required')
        selected, pending = handoff['selected_autosomes'], handoff['pending_autosomes']
        if (len(selected) != len(set(selected)) or len(pending) != len(set(pending))
                or set(selected) & set(pending) or set(selected + pending) != set(map(str, range(1, 23)))
                or handoff['all_autosomes_complete'] != (not pending)):
            raise ValueError('inconsistent autosomal snapshot scope')
        if pending and not args.allow_partial_snapshot:
            raise ValueError('incomplete autosomal snapshot requires explicit --allow-partial-snapshot')
        if [name for shard in handoff['shards'] for name, _ in shard['contigs']] != selected:
            raise ValueError('handoff shard membership differs from selected autosomes')
        pin('reference', reference, handoff['inputs']['reference']['sha256'])
        pin('reference_index', Path(str(reference) + '.fai'), handoff['inputs']['reference_index']['sha256'])
        pin('gff', gff, args.expected_gff_sha256)
        pin('config', args.config)
        config = json.loads(args.config.read_text())
        if config['sample'] != handoff['sample'] or config['quality_handling']['mode'] != 'descriptive_only_no_FILTER_DP_GQ_exclusion':
            raise ValueError('unsupported sample or quality-handling contract')
        if config['quality_handling']['historical_diagnostic_thresholds'] != {'DP': 10, 'GQ': 20}:
            raise ValueError('unsupported historical quality diagnostic bins')
        scripts = Path(__file__).resolve().parent
        for name in ('prepare_native_coding_review.py', 'native_coding_records.py', 'audit_discovery_handoff.py',
                     'rare_biallelic_baseline.py', 'screen_secondary_findings.py', 'AlleleFrequencyProbe.java'):
            pin('source/' + name, scripts / name)
        bcftools = root / 'tools/bio/bin/bcftools'
        pin('bcftools', bcftools)
        sources = []
        for shard in handoff['shards']:
            method = shard['methods'][args.caller]
            if not method['replayed_split_payload_identical']:
                raise ValueError('source lacks completed split replay')
            source = method['original_raw_pin']
            pinned = pin('native_source/shard-%04d' % shard['shard_id'], source['path'], source['sha256'])
            sources.append((shard, Path(source['path']), pinned))
        state.update(selected_autosomes=selected, pending_autosomes=pending, inputs=pins, config=config)
        save(output / 'status.json', state)
        variants, shard_results = [], []
        for shard, source, source_pin in sources:
            state['stage'] = 'annotate and recover original context for shard-%04d' % shard['shard_id']
            save(output / 'status.json', state)
            rows, checked = prepare_shard(source, source_pin, shard, args.caller, handoff['sample'],
                config, bcftools, reference, gff, output / ('shard-%04d' % shard['shard_id']))
            variants.extend(rows)
            shard_results.append(checked)
            print(json.dumps(dict(completed_shard=shard['shard_id'], coding_allele_rows=len(rows))), flush=True)
        key_counts = Counter(row['key'] for row in variants)
        collisions = sorted(key for key, count in key_counts.items() if count > 1)
        for variant in variants:
            duplicate = key_counts[variant['key']] > 1
            variant['duplicate_normalized_key'] = duplicate
            variant['review_axes'] = review_axes(variant['source_call'], duplicate=duplicate)
        save(output / 'normalized-key-collisions.json', dict(keys=collisions,
             interpretation='All colliding source identities are retained in the allele inventory and kept out of ordinary recessive models.'))
        state['stage'] = 'external allele-frequency counts'
        save(output / 'status.json', state)
        queries = output / 'frequency-queries.txt'
        queries.write_text(''.join(key.replace(':', '\t') + '\n' for key in sorted(key_counts)))
        if args.frequency_responses:
            if not args.expected_frequency_sha256:
                raise ValueError('precomputed frequency responses require their expected SHA-256')
            frequency_file = args.frequency_responses.resolve()
            pin('provided_frequency_responses', frequency_file, args.expected_frequency_sha256)
            frequency_provenance = dict(mode='provided_exact_response_file',
                limitation='A matching file hash establishes identity, not the provenance or accuracy of provided counts.')
        else:
            java = root / 'tools/bio/lib/jvm/bin/java'
            jar = root / 'lirical/lirical-cli-2.4.1/lirical-cli-2.4.1.jar'
            database = root / 'exomiser/2512_hg38/2512_hg38_variants.mv.db'
            for label, path in (('java', java), ('frequency_database', database), ('lirical_jar', jar)):
                pin(label, path)
            frequency_file = output / 'frequencies.jsonl'
            command = [str(java), '-Xmx4g', '-XX:ActiveProcessorCount=2', '-cp', str(jar),
                       str(scripts / 'AlleleFrequencyProbe.java'), str(database)]
            with queries.open('rb') as inp, frequency_file.open('wb') as out, (output / 'frequency.log').open('wb') as log:
                subprocess.run(command, stdin=inp, stdout=out, stderr=log, check=True)
            frequency_provenance = dict(mode='pinned_local_database_query', command=command)
        frequency_responses = attach_frequencies(variants, frequency_file)
        write_jsonl(output / 'coding-alleles.jsonl.gz', variants)
        limits = [config['frequency']['primary_maximum_AF_fraction']] + config['frequency']['sensitivity_maximum_AF_fractions']
        profiles = {}
        for maximum in limits:
            models, review, counts = build_models(variants, maximum)
            write_jsonl(output / ('descriptive-models-AF%s.jsonl.gz' % maximum), models)
            write_jsonl(output / ('pair-representation-phase-review-AF%s.jsonl.gz' % maximum), review)
            profiles[str(maximum)] = dict(models=len(models), retained_review_pairs=len(review), counters=counts,
                genes=sorted({model['gene'] for model in models}),
                interpretation='Descriptive genotype/consequence combinations within this caller and snapshot; not diagnoses, probability ranks or validated loss-of-function findings.')
        for name, value in pins.items():
            current = Path(value['path']).stat()
            if (current.st_size, current.st_mtime_ns) != (value['bytes'], value['mtime_ns']):
                raise ValueError('input changed during coding preparation: ' + name)
        outputs = {str(path.relative_to(output)): fingerprint(path) for path in sorted(output.rglob('*'))
                   if path.is_file() and path.name not in ('status.json', 'summary.json')}
        summary = dict(status='complete', completed_at=stamp(), elapsed_seconds=time.monotonic() - started,
            caller=args.caller, sample=handoff['sample'], selected_autosomes=selected, pending_autosomes=pending,
            all_autosomes_complete=not pending, scope='Caller-separated protein-altering allele inventory with original genotype context and descriptive recessive models; no phenotype scoring or EPCR.',
            coding_allele_rows=len(variants), unique_normalized_keys=len(key_counts), duplicate_normalized_keys=len(collisions),
            source_records_with_coding_rows=len({v['source_call']['source_record_id'] for v in variants}),
            frequency_responses=frequency_responses, frequency_provenance=frequency_provenance,
            review_axis_counts={axis: dict(Counter(json.dumps(v['review_axes'][axis], sort_keys=True) for v in variants))
                               for axis in ('filter_state', 'genotype_state', 'allele_participation', 'structural_model_eligible')},
            profiles=profiles, shards=shard_results, inputs=pins, outputs=outputs, config=config,
            runtime=dict(python=os.sys.version, pysam=pysam.__version__), limitations=config['limitations'])
        save(output / 'summary.json', summary)
        state.update(status='complete', completed_at=summary['completed_at'], stage='completed',
                     summary_sha256=fingerprint(output / 'summary.json')['sha256'])
        save(output / 'status.json', state)
        return summary
    except BaseException as error:
        state.update(status='failed', failed_at=stamp(), error_type=type(error).__name__, error=str(error))
        save(output / 'status.json', state)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('root', 'handoff', 'reference', 'gff', 'config', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--expected-handoff-sha256', required=True)
    parser.add_argument('--expected-gff-sha256', required=True)
    parser.add_argument('--caller', choices=('gatk', 'bcftools'), required=True)
    parser.add_argument('--allow-partial-snapshot', action='store_true')
    parser.add_argument('--frequency-responses', type=Path)
    parser.add_argument('--expected-frequency-sha256')
    args = parser.parse_args()
    os.nice(10)
    run(args)


if __name__ == '__main__':
    main()
