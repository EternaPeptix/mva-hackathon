#!/usr/bin/env python3
"""Exact dated secondary-allele inventory with complete native source provenance.

No caller, clinical classification, quality eligibility gate or coding filter.
"""
import argparse
import datetime
import re
from collections import Counter, defaultdict
import gzip
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import platform
import sys
import time

import pysam

AUTOSOMES = [str(n) for n in range(1, 23)]
CALLERS = ['gatk', 'bcftools']
FROZEN_PUBLIC = {'clinvar': '47ea3626cd90154d82304afe3b8f717aab2557a1330f4ea2c809892386e08d49',
    'clinvar_index': 'b220738d2c086062f136600cdd31a95bec54586890cb7d92868ad78bdc632c03',
    'panel': '8089748ab8645336eedab56e4f548f7c902e2f7d163fa3a559b49c9fbd93151d'}
NATIVE_ORD = 'NS_SOURCE_ORD'
PUBLIC_ORD = 'NS_PUBLIC_ORD'
HELPERS = {'screen_secondary_findings.py': '34b7a4bebfb48d1c460fd579d862f18ff802f9ed56c8a39a5475882bc3e1802d',
           'native_coding_records.py': 'c395dd130a2577e7b25ba0ac62d600d7cb64c9240f3eb218536eff22c5c44cfe'}


def require(test, message):
    if not test:
        raise ValueError(message)


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def same(a, b, label):
    require(canonical(a) == canonical(b), label + ' differs')


def parse(text):
    def pairs(items):
        out = {}
        for key, value in items:
            require(key not in out, 'duplicate JSON key')
            out[key] = value
        return out
    def invalid(value):
        raise ValueError('nonfinite JSON: ' + value)
    return json.loads(text, object_pairs_hook=pairs, parse_constant=invalid)


def stamp():
    return time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')


class Pins:
    def __init__(self):
        self.files = {}

    def get(self, path, expected=None):
        path = Path(path).resolve()
        if expected is not None and 'path' in expected:
            same(str(path), str(Path(expected['path']).resolve()), 'expected pin path')
        key = str(path)
        if key not in self.files:
            before = path.stat(); digest = hashlib.sha256()
            with path.open('rb') as stream:
                for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
                    digest.update(block)
            after = path.stat()
            require((before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns), 'file changed during hash')
            self.files[key] = dict(path=key, bytes=after.st_size, mtime_ns=after.st_mtime_ns, sha256=digest.hexdigest())
        got = self.files[key]
        if expected is not None:
            same(got['sha256'], expected['sha256'], 'expected SHA256')
            if 'bytes' in expected:
                same(got['bytes'], expected['bytes'], 'expected bytes')
        return got

    def read(self, expected):
        self.get(expected['path'], expected)
        return parse(Path(expected['path']).read_text())

    def stable(self):
        for p in self.files.values():
            st = Path(p['path']).stat()
            require((st.st_size, st.st_mtime_ns) == (p['bytes'], p['mtime_ns']), 'input changed during review')


def load_helpers(pins):
    modules = {}
    for name, digest in HELPERS.items():
        path = Path(__file__).parent / 'vendor' / name
        pins.get(path, {'sha256': digest})
        spec = importlib.util.spec_from_file_location('native_secondary_' + path.stem, path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        modules[path.stem] = module
    return modules['screen_secondary_findings'], modules['native_coding_records']


def bind(path, manifest, pins, root=None):
    index = {}
    for value in manifest:
        p = Path(value['path'])
        if not p.is_absolute():
            require(root is not None, 'relative proof without root')
            p = root / p
        key = str(p.resolve())
        require(key not in index, 'duplicate proof path')
        index[key] = dict(value, path=key)
    key = str(Path(path).resolve())
    require(key in index, 'actual path absent from proof')
    return pins.get(path, index[key])


def validate_regions(model, fasta):
    regions = model['regions']
    same([r['chrom'] for r in regions], ['X', 'Y'], 'fixed X/Y region scope')
    for region in regions:
        same(region['length'], fasta.get_reference_length(region['chrom']), 'region/reference length')
        end = 0
        for start, stop, ploidy in region['intervals_0based_half_open']:
            require(all(type(v) is int for v in (start, stop, ploidy)), 'noninteger interval')
            require(start == end and stop > start and ploidy in (1, 2), 'invalid ploidy tiling')
            if region['chrom'] == 'Y':
                same(ploidy, 1, 'Y ploidy')
            end = stop
        same(end, region['length'], 'complete ploidy interval coverage')
    require({p for _, _, p in regions[0]['intervals_0based_half_open']} == {1, 2}, 'X must declare PAR and non-PAR model')
    require(model['required_N_masks'], 'fixed Y masks absent')
    for mask in model['required_N_masks']:
        require(mask['chrom'] == 'Y' and type(mask['start']) is int and type(mask['end']) is int
                and 0 <= mask['start'] < mask['end'] <= fasta.get_reference_length('Y'), 'invalid Y mask')
        seq = fasta.fetch('Y', mask['start'], mask['end']).upper()
        require(len(seq) == mask['end'] - mask['start'] and set(seq) == {'N'}, 'required Y mask differs')
    return regions


def expected_ploidy(record, source, regions):
    require(record.contig in source['contigs'], 'native record outside declared source scope')
    if source['kind'] == 'autosome':
        return 2
    intervals = [(p, end) for r in regions if r['chrom'] == record.contig
                 for start, end, p in r['intervals_0based_half_open'] if start <= record.start < end]
    require(len(intervals) == 1, 'record outside fixed ploidy model')
    p, end = intervals[0]
    require(record.stop <= end, 'native REF span crosses ploidy boundary')
    return p


def genotype_shape(context, expected):
    fields = context['FORMAT']; gt = fields.get('GT') or []
    informative = any(value is not None for value in gt)
    if informative:
        same(len(gt), expected, 'informative native GT/model ploidy')
    state = ('absent' if 'GT' not in fields else 'fully_missing' if not informative
             else 'partial' if None in gt else 'fully_called')
    return dict(state=state, expected_model_ploidy=expected, encoded_GT_slots=len(gt),
                observed_GT_ploidy=len(gt) if informative else None,
                interpretation='Missing GT encoding is retained; encoded slot count does not establish ploidy.')


def prepare_sources(config, inputs, pins, fasta):
    root = Path(config['root']).resolve()
    handoff, verification = pins.read(inputs['handoff']), pins.read(inputs['autosome_verification'])
    same(handoff['status'], 'complete', 'handoff completion')
    same(verification['status'], 'complete', 'native verification completion')
    for scope in (handoff, verification['scope']):
        same(scope['selected_autosomes'], AUTOSOMES, 'all22 source scope')
        same(scope['pending_autosomes'], [], 'pending source scope')
        same(scope['all_autosomes_complete'], True, 'full-autosome completion')
    same(handoff['sample'], config['sample'], 'handoff sample')
    bind(inputs['handoff']['path'], verification['files'], pins)
    for label in ('reference', 'reference_index', 'bcftools'):
        pins.get(inputs[label]['path'], handoff['inputs'][label])
    shards = handoff['shards']
    same([s['shard_id'] for s in shards], list(range(22)), 'all22 source-shard membership')
    same([c for s in shards for c, _ in s['contigs']], AUTOSOMES, 'all22 source chromosome membership')
    sources = []
    for shard in shards:
        require(len(shard['contigs']) == 1, 'one autosome per frozen shard required')
        chrom, length = shard['contigs'][0]
        same(length, fasta.get_reference_length(chrom), 'autosome/reference extent')
        same(sorted(shard['methods']), sorted(CALLERS), 'source caller set')
        for caller in CALLERS:
            item = shard['methods'][caller]
            same(item['replayed_split_payload_identical'], True, 'original split-replay proof')
            pin = pins.get(item['original_raw_pin']['path'], item['original_raw_pin'])
            bind(pin['path'], verification['files'], pins)
            sources.append(dict(caller=caller, kind='autosome', contigs=[chrom], pin=pin,
                                expected_census=item['original_raw'], id=caller + '-' + chrom))
    status = pins.read(inputs['fixed_xy_status'])
    review = pins.read(inputs['fixed_xy_review'])
    proof = pins.read(inputs['fixed_xy_verification'])
    same(status['status'], 'complete', 'fixed X/Y completion')
    same(review['status'], 'complete', 'fixed X/Y independent review')
    same(review['fixed_run_status_sha256'], inputs['fixed_xy_status']['sha256'], 'fixed status/review binding')
    # Legacy fixed-X/Y verification has no status field: successful evidence is
    # its terminal timestamp, matching executed source, complete output count/map,
    # explicit source-preservation result and linked completed review.
    reviewed = datetime.datetime.fromisoformat(review['verified_at'].replace('Z', '+00:00'))
    verified = datetime.datetime.fromisoformat(proof['verified_at'].replace('Z', '+00:00'))
    require(reviewed.tzinfo is not None and verified.tzinfo is not None and verified >= reviewed, 'fixed verification terminal timestamp')
    require(re.fullmatch('[0-9a-f]{64}', review['script_sha256']) is not None, 'fixed review source digest')
    same(proof['script_sha256'], review['script_sha256'], 'fixed verification source identity')
    same(review['fixed_outputs_verified'], len(status['output_sha256']), 'all fixed outputs verified')
    require(type(review['original_files_preserved']) is int and review['original_files_preserved'] > 0, 'fixed original preservation evidence')
    same(review['bcftools_original_reused_byte_for_byte'], True, 'fixed reused source verification')
    same(proof['fixed_outputs'], status['output_sha256'], 'fixed output proof map')
    bind(inputs['fixed_xy_review']['path'], proof['files'], pins, root)
    same(status['config']['sample'], config['sample'], 'fixed X/Y sample')
    same(status['config']['reference_sha256'], inputs['reference']['sha256'], 'fixed reference identity')
    regions = validate_regions(status['config'], fasta)
    xy = Path(inputs['fixed_xy_status']['path']).resolve().parent
    for caller in CALLERS:
        name = caller + '.original.vcf.gz'
        pin = pins.get(xy / name, {'sha256': status['output_sha256'][name]})
        pins.get(xy / (name + '.tbi'), {'sha256': status['output_sha256'][name + '.tbi']})
        counts = review['raw_record_accounting']['fixed_gatk' if caller == 'gatk' else 'bcftools']
        sources.append(dict(caller=caller, kind='fixed_xy', contigs=['X', 'Y'], pin=pin,
                            expected_census={'counts': {'records': counts['records']}}, id=caller + '-fixed-xy'))
    require(len({s['pin']['path'] for s in sources}) == len(sources), 'duplicate source paths')
    for caller in CALLERS:
        chosen = [s for s in sources if s['caller'] == caller]
        same([c for s in chosen for c in s['contigs']], AUTOSOMES + ['X', 'Y'], 'complete distinct caller scope')
    return sources, regions


def quality(call, native, shape):
    axes = native.review_axes(call)
    return {k: axes[k] for k in ('filter_state', 'filter_values', 'DP', 'GQ', 'allele_participation', 'nonconcrete_in_GT')} | dict(
        genotype=shape, called_ALT_present=call['alt_copies'] > 0,
        fully_called_ALT_present=call['fully_called'] and call['alt_copies'] > 0,
        historical_bins_are_eligibility_gate=False)


def literal_equivalence(fasta, normalized, original):
    same(normalized.contig, original['chrom'], 'normalization chromosome')
    start = min(normalized.pos, original['pos']) - 1
    end = max(normalized.pos - 1 + len(normalized.ref), original['pos'] - 1 + len(original['ref']))
    sequence = fasta.fetch(normalized.contig, start, end).upper()
    require(len(sequence) == end - start, 'allele extends beyond reference')
    def edited(pos, ref, alt):
        offset = pos - 1 - start
        require(sequence[offset:offset + len(ref)] == ref, 'source/normalized REF mismatch')
        return sequence[:offset] + alt + sequence[offset + len(ref):]
    same(edited(original['pos'], original['ref'], original['alt']),
         edited(normalized.pos, normalized.ref, normalized.alts[0]), 'literal normalized edit')


def public_lookup(config, inputs, directory, sf, native, fasta):
    genes = sf.load_panel(inputs['panel']['path'], config['expected_genes'])
    originals, counts = {}, Counter()
    counts.update({k: 0 for k in ('original_records','panel_gene_records','selected_source_records',
        'selected_biallelic_source_records','selected_multiallelic_source_records','normalized_lineages',
        'signal_nonprimary_or_nonconcrete_records','selected_nonconcrete_normalized_lineages')})
    selected = directory / 'public-selected.vcf.gz'
    with pysam.VariantFile(inputs['clinvar']['path']) as source:
        require('##reference=GRCh38\n' in str(source.header), 'ClinVar assembly declaration differs')
        require(PUBLIC_ORD not in source.header.info and 'SF_ORIG' not in source.header.info, 'reserved public origin tag')
        dates = [line.split('=',1)[1] for line in str(source.header).splitlines() if line.startswith('##fileDate=')]
        same(dates, [config['public_snapshot']['clinvar_fileDate']], 'dated ClinVar header')
        (directory / 'original-clinvar-header.txt').write_text(str(source.header))
        header = source.header.copy()
        header.info.add(PUBLIC_ORD, 1, 'Integer', 'Original ordinal in pinned dated ClinVar file')
        for chrom in AUTOSOMES + ['X', 'Y']:
            if chrom not in header.contigs:
                header.contigs.add(chrom, length=fasta.get_reference_length(chrom))
        with pysam.VariantFile(str(selected), 'wz', header=header) as writer:
            for ordinal, record in enumerate(source, 1):
                counts['original_records'] += 1
                if not sf.gene_symbols(record.info) & genes.keys():
                    continue
                counts['panel_gene_records'] += 1
                stratum = sf.germline_stratum(record.info)
                if stratum is None:
                    continue
                if record.contig not in AUTOSOMES + ['X', 'Y'] or re.fullmatch('[ACGT]+', record.ref) is None or not any(re.fullmatch('[ACGT]+', alt) for alt in (record.alts or ())):
                    counts['signal_nonprimary_or_nonconcrete_records'] += 1
                    continue
                originals[ordinal] = dict(source_record_id='clinvar:' + inputs['clinvar']['sha256'] + ':' + str(ordinal),
                    source_sha256=inputs['clinvar']['sha256'], source_ordinal=ordinal,
                    chrom=record.contig, pos=record.pos, ref=record.ref, alts=list(record.alts),
                    rendered_record=str(record), INFO=native._json_value(dict(record.info)),
                    id=record.id, QUAL=native._json_value(record.qual), FILTER=list(record.filter))
                copy = record.copy(); copy.translate(header); copy.info[PUBLIC_ORD] = ordinal; writer.write(copy)
                counts['selected_source_records'] += 1
                counts['selected_biallelic_source_records' if len(record.alts) == 1 else 'selected_multiallelic_source_records'] += 1
    lookup = defaultdict(list)
    public_lineages = set()
    with sf.normalized_stream(inputs['bcftools']['path'], inputs['reference']['path'], selected,
                              directory / 'public-normalization.log') as (reader, commands):
        for record in reader:
            ordinal = record.info.get(PUBLIC_ORD)
            require(type(ordinal) is int and ordinal in originals, 'missing public original ordinal')
            context = originals[ordinal]
            if not sf.concrete_alleles(record):
                counts['selected_nonconcrete_normalized_lineages'] += 1
                continue
            index = native.original_alt_index(record, context, 'SF_ORIG')
            original = dict(chrom=context['chrom'], pos=context['pos'], ref=context['ref'], alt=context['alts'][index - 1])
            literal_equivalence(fasta, record, original)
            source_info = {k: tuple(v) if isinstance(v, list) else v for k, v in context['INFO'].items()}
            source_assertion = sf.germline_stratum(source_info)
            require(source_assertion is not None, 'selected original public assertion lost')
            source_stratum, stars = source_assertion
            normalized_assertion = sf.germline_stratum(record.info)
            ambiguous = len(context['alts']) != 1
            if not ambiguous:
                same(normalized_assertion, source_assertion, 'biallelic public assertion normalization')
            stratum = 'multiallelic_assertion_mapping_review' if ambiguous else source_stratum
            identity = (ordinal, index, sf.allele_key(record))
            require(identity not in public_lineages, 'duplicate same public normalization lineage')
            public_lineages.add(identity)
            lookup[sf.allele_key(record)].append(dict(source=context, source_alt_index=index,
                normalized_record=str(record), normalized_allele=dict(chrom=record.contig, pos=record.pos, ref=record.ref, alt=record.alts[0]),
                stratum=stratum, source_record_stratum=source_stratum, stars=stars,
                normalized_record_stratum=normalized_assertion[0] if normalized_assertion else None,
                normalized_record_stars=normalized_assertion[1] if normalized_assertion else None,
                assertion_mapping=dict(status='ambiguous_multiallelic_source' if ambiguous else 'unambiguous_biallelic_source',
                    original_CLNALLE=context['INFO'].get('CLNALLE'), normalized_CLNALLE=native._json_value(record.info['CLNALLE']) if 'CLNALLE' in record.info else None),
                genes=sorted(sf.gene_symbols(record.info) & genes.keys()),
                classification=sf.values(record.info, 'CLNSIG'), review_status=sf.values(record.info, 'CLNREVSTAT'),
                public_original_multiallelic=len(context['alts']) > 1,
                interpretation='Exact dated assertion context; condition/variant/inheritance reportability remains pending.'))
            counts['normalized_lineages'] += 1
    counts['unique_normalized_identities'] = len(lookup)
    return lookup, dict(counts), commands, genes


def tagged_source(source, output, sample, regions, native):
    digest, count = hashlib.sha256(), 0
    with pysam.VariantFile(source['pin']['path']) as reader:
        same(list(reader.header.samples), [sample], 'native source sample')
        require(not {NATIVE_ORD, PUBLIC_ORD, 'SF_ORIG'} & set(reader.header.info), 'reserved native origin tag')
        header = reader.header.copy()
        header.info.add(NATIVE_ORD, 1, 'Integer', 'Original ordinal in byte-pinned native source')
        with pysam.VariantFile(str(output), 'wb', header=header) as writer:
            for count, record in enumerate(reader, 1):
                p = expected_ploidy(record, source, regions)
                context = native.source_context(record, sample, source['caller'], source['pin']['sha256'], count)
                genotype_shape(context, p)
                digest.update(str(record).encode())
                copy = record.copy(); copy.translate(header); copy.info[NATIVE_ORD] = count; writer.write(copy)
    expected = source['expected_census']
    same(count, expected['counts']['records'], 'complete source record census')
    if 'rendered_record_payload_sha256' in expected:
        same(digest.hexdigest(), expected['rendered_record_payload_sha256'], 'complete source payload census')
    replay, replay_count = hashlib.sha256(), 0
    with pysam.VariantFile(str(output)) as reader:
        for replay_count, record in enumerate(reader, 1):
            same(record.info[NATIVE_ORD], replay_count, 'copied source ordinal')
            del record.info[NATIVE_ORD]
            replay.update(str(record).encode())
    same((replay_count, replay.hexdigest()), (count, digest.hexdigest()), 'full copied source payload')
    return dict(records=count, rendered_record_payload_sha256=digest.hexdigest(), copied_original_payload_preserved=True)


def scan_source(source, sample, regions, lookup, inputs, directory, sf, native, fasta):
    directory.mkdir()
    tagged = directory / 'ordinal-tagged-original.bcf'
    census = tagged_source(source, tagged, sample, regions, native)
    selected, counts = defaultdict(list), Counter()
    with sf.normalized_stream(inputs['bcftools']['path'], inputs['reference']['path'], tagged,
                              directory / 'normalization.log') as (reader, commands):
        for record in reader:
            counts['normalized_records'] += 1
            expected_ploidy(record, source, regions)
            ordinal = record.info.get(NATIVE_ORD)
            require(type(ordinal) is int and 1 <= ordinal <= census['records'], 'invalid native normalization ordinal')
            if not sf.concrete_alleles(record):
                counts['nonconcrete_normalized_records'] += 1
                continue
            key = sf.allele_key(record)
            if key in lookup:
                selected[ordinal].append(record.copy())
    rows, recovered, identities = [], set(), set()
    with pysam.VariantFile(source['pin']['path']) as reader:
        for ordinal, record in enumerate(reader, 1):
            if ordinal not in selected:
                continue
            context = native.source_context(record, sample, source['caller'], source['pin']['sha256'], ordinal)
            shape = genotype_shape(context, expected_ploidy(record, source, regions))
            recovered.add(ordinal)
            for normalized in selected[ordinal]:
                index = native.original_alt_index(normalized, context, 'SF_ORIG')
                call = native.allele_call(context, index)
                same(expected_ploidy(normalized, source, regions), shape['expected_model_ploidy'], 'normalization retains ploidy interval class')
                literal_equivalence(fasta, normalized, call['source_allele'])
                key = sf.allele_key(normalized)
                identity = (call['source_record_id'], index, key)
                require(identity not in identities, 'duplicate same source-ALT normalization lineage')
                identities.add(identity)
                rows.append(dict(caller=source['caller'], source_scope=source['id'],
                    normalized_allele=dict(zip(('chrom', 'pos', 'ref', 'alt'), key)),
                    normalized_record=str(normalized), source_call=call, quality=quality(call, native, shape),
                    clinvar=lookup[key], reportability='pending_disease_variant_and_inheritance_review'))
    same(sorted(recovered), sorted(selected), 'all selected native originals recovered')
    counts['exact_matching_occurrences'] = len(rows)
    return rows, dict(source=source, census=census, counts=dict(counts), commands=commands)


def summarize(rows):
    groups = {}
    for row in rows:
        key = canonical(dict(caller=row['caller'], allele=row['normalized_allele']))
        groups.setdefault(key, []).append(row)
    identities, callers = [], {c: Counter() for c in CALLERS}
    for key, occurrences in sorted(groups.items()):
        description = parse(key); caller = description['caller']
        any_alt = any(r['quality']['called_ALT_present'] for r in occurrences)
        complete_alt = any(r['quality']['fully_called_ALT_present'] for r in occurrences)
        lineages = [dict(source_record_id=r['source_call']['source_record_id'], source_alt_index=r['source_call']['source_alt_index']) for r in occurrences]
        identity_metadata = {'source_record_id', 'source_namespace', 'caller', 'sample', 'source_sha256', 'source_ordinal'}
        signatures = {canonical({k:v for k,v in r['source_call']['native_context'].items() if k not in identity_metadata})
                      for r in occurrences}
        identities.append(description | dict(source_lineages=lineages, source_occurrences=len(occurrences),
            duplicate_normalized_identity=len(occurrences) > 1, differing_native_evidence=len(signatures) > 1,
            called_ALT_present=any_alt, fully_called_ALT_present=complete_alt,
            interpretation='No consensus genotype is assigned; every source occurrence is retained.'))
        count = callers[caller]; count['unique_exact_matching_identities'] += 1
        count['exact_matching_occurrences'] += len(occurrences)
        count['unique_called_ALT_identities'] += int(any_alt)
        count['unique_fully_called_ALT_identities'] += int(complete_alt)
        count['unique_identities_without_called_ALT'] += int(not any_alt)
        count['duplicate_normalized_identities'] += int(len(occurrences) > 1)
        if any_alt:
            count.update('called_ALT_unique_identity_stratum/' + value for value in {c['stratum'] for r in occurrences for c in r['clinvar']})
    fields = ('unique_exact_matching_identities','exact_matching_occurrences','unique_called_ALT_identities',
              'unique_fully_called_ALT_identities','unique_identities_without_called_ALT','duplicate_normalized_identities')
    return identities, {c: {field: count[field] for field in fields} | dict(count) for c, count in callers.items()}


def frozen_public_config(config):
    same(config['expected_genes'], 84, 'frozen panel count')
    snapshot = config['public_snapshot']
    same(snapshot['clinvar_fileDate'], '2026-09-13', 'frozen literal ClinVar header date')
    same(snapshot['clinvar_snapshot_label'], '20260913', 'frozen ClinVar snapshot label')
    same(snapshot['clinvar_resource_filename'], 'clinvar_20260913.vcf.gz', 'frozen ClinVar resource filename')
    same(Path(config['inputs']['clinvar']['path']).name, snapshot['clinvar_resource_filename'], 'ClinVar source filename')
    for label, digest in FROZEN_PUBLIC.items():
        same(config['inputs'][label]['sha256'], digest, 'frozen public input ' + label)


def run(config_path, config_sha256, output):
    require(not output.exists(), 'existing output directory preserved')
    pins = Pins(); config = pins.read(dict(path=str(config_path), sha256=config_sha256))
    same(config['schema_version'], 1, 'configuration schema')
    same(config['ready'], True, 'configuration ready flag')
    same(config['quality_policy'], 'descriptive_native_axes_no_eligibility_gate', 'quality policy')
    require(type(config['fixture_only']) is bool, 'fixture scope must be explicit boolean')
    require(type(config['expected_genes']) is int and config['expected_genes'] > 0, 'invalid panel size')
    inputs = config['inputs']
    if not config['fixture_only']:
        frozen_public_config(config)
    required = {'reference','reference_index','bcftools','clinvar','clinvar_index','panel','handoff',
                'autosome_verification','fixed_xy_status','fixed_xy_review','fixed_xy_verification'}
    same(sorted(inputs), sorted(required), 'complete input pin set')
    for item in inputs.values():
        require(Path(item['path']).is_absolute() and type(item['bytes']) is int and item['bytes'] >= 0, 'absolute input path/size required')
        pins.get(item['path'], item)
    same(str(Path(inputs['reference_index']['path']).resolve()), str(Path(inputs['reference']['path'] + '.fai').resolve()), 'reference index path')
    same(str(Path(inputs['clinvar_index']['path']).resolve()), str(Path(inputs['clinvar']['path'] + '.tbi').resolve()), 'ClinVar index path')
    sf, native = load_helpers(pins)
    pins.get(__file__)
    output.mkdir(parents=True)
    state = dict(status='running', started_at=stamp(), fixture_only=config['fixture_only'],
                 scope='Exact dated secondary-allele inventory; no clinical classification, quality eligibility or callability assertion.')
    save(output / 'status.json', state)
    try:
        with pysam.FastaFile(inputs['reference']['path']) as fasta:
            sources, regions = prepare_sources(config, inputs, pins, fasta)
            lookup, public_counts, commands, genes = public_lookup(config, inputs, output, sf, native, fasta)
            rows, censuses = [], []
            for source in sources:
                found, census = scan_source(source, config['sample'], regions, lookup, inputs, output / source['id'], sf, native, fasta)
                rows.extend(found); censuses.append(census)
        with gzip.open(output / 'exact-matching-occurrences.jsonl.gz', 'wt') as stream:
            for row in rows:
                stream.write(canonical(row) + '\n')
        identities, counts = summarize(rows)
        save(output / 'exact-identity-groups.json', identities)
        save(output / 'source-censuses.json', censuses)
        quality_counts = {caller: {axis: dict(Counter((row['quality'][axis]['state'] if axis in ('DP','GQ') else row['quality'][axis])
            for row in rows if row['caller'] == caller)) for axis in ('filter_state','DP','GQ')} for caller in CALLERS}
        result = dict(status='complete', completed_at=stamp(), fixture_only=config['fixture_only'], sample=config['sample'],
            public_snapshot=config['public_snapshot'], panel_genes=sorted(genes), public_counts=public_counts,
            counts_by_caller=counts, exact_occurrence_quality_axes=quality_counts,
            source_scope={c: AUTOSOMES + ['X','Y'] for c in CALLERS}, source_files=len(sources),
            public_normalization_commands=commands, helper_sha256=HELPERS,
            runtime=dict(host=platform.node(), python=sys.version, pysam=pysam.__version__),
            inputs=pins.files, limitations=[
                'Exact matching occurrences and unique called-ALT identities are separate counts; zero of one is not zero of the other.',
                'Every duplicate source occurrence remains visible. No consensus GT, parental phase or cross-caller genotype agreement is inferred.',
                'No exact called-ALT match means only none in the explicit native sources and dated resources; not a clinical negative, callability or absent-record reference genotype.',
                'Noncoding concrete exact matches are retained, but novel, structural, repeat and mosaic detection or complete panel coverage is not established.',
                'Quality fields are native descriptive axes; no PASS/GQ eligibility gate, imputation or clinical sensitivity claim.',
                'ClinVar strata are dated assertion context; condition, variant restrictions and inheritance/reportability remain unreviewed.'])
        pins.stable()
        save(output / 'summary.json', result)
        output_pins = [Pins().get(p) for p in sorted(output.rglob('*')) if p.is_file() and p.name != 'status.json']
        state.update(status='complete', completed_at=result['completed_at'], outputs=output_pins,
                     summary_sha256=Pins().get(output/'summary.json')['sha256'])
        save(output / 'status.json', state)
        return result
    except BaseException as error:
        state.update(status='failed', failed_at=stamp(), error_type=type(error).__name__, error=str(error))
        save(output / 'status.json', state)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--expected-config-sha256', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); os.nice(10)
    result = run(args.config.resolve(), args.expected_config_sha256, args.output.resolve())
    print(json.dumps({'status':result['status'], 'counts_by_caller':result['counts_by_caller']}))


if __name__ == '__main__':
    main()
