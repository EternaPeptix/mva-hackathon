#!/usr/bin/env python3
"""Preserve native calls and describe coding models without a quality gate.

DP 10 and GQ 20 are historical diagnostic bins only. Models describe source
genotypes and consequences; they are neither diagnoses nor probability ranks.
This module deliberately has no filesystem, caller or annotation side effects.
"""
from collections import Counter, defaultdict
from copy import deepcopy
import itertools
import json
import math
import re


def _json_value(value):
    if isinstance(value, (tuple, list)):
        return [_json_value(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if value is None or isinstance(value, (str, int, bool)):
        return value
    if isinstance(value, float) and math.isfinite(value):
        return value
    raise ValueError('unsupported or nonfinite native field value')


def _concrete(value):
    return isinstance(value, str) and re.fullmatch('[ACGT]+', value) is not None


def source_context(record, sample, caller, source_sha256, ordinal):
    """Record ordinal is a nonnegative sequential position in a pinned source."""
    if not isinstance(caller, str) or not re.fullmatch('[A-Za-z0-9_.-]+', caller):
        raise ValueError('invalid caller namespace')
    if not re.fullmatch('[0-9a-f]{64}', source_sha256):
        raise ValueError('source SHA-256 must be lowercase hexadecimal')
    if type(ordinal) is not int or ordinal < 0:
        raise ValueError('invalid source record ordinal')
    if sample not in record.samples:
        raise ValueError('sample absent from native record')
    native = record.samples[sample]
    fields = {key: _json_value(value) for key, value in native.items()}
    gt = fields.get('GT') or []
    alts = list(record.alts or ())
    if any(value is not None and (type(value) is not int or not 0 <= value <= len(alts))
           for value in gt):
        raise ValueError('invalid native GT allele index')
    return {
        'source_record_id': '%s:%s:%d' % (caller, source_sha256, ordinal),
        'source_namespace': json.dumps([caller, sample, source_sha256], separators=(',', ':')),
        'caller': caller, 'sample': sample, 'source_sha256': source_sha256,
        'source_ordinal': ordinal, 'chrom': record.contig, 'pos': record.pos,
        'ref': record.ref, 'alts': alts, 'FORMAT': fields,
        'GT_phased': bool(native.phased), 'FILTER': list(record.filter),
        'QUAL': _json_value(record.qual), 'rendered_record': str(record),
    }


def original_alt_index(normalized_record, context, origin_tag='NC_ORIG'):
    """Verify the old-record tag, without interpreting split GT/AD/PL as native."""
    alts = normalized_record.alts or ()
    if len(alts) != 1:
        raise ValueError('expected one normalized ALT')
    marker = normalized_record.info.get(origin_tag)
    if isinstance(marker, tuple):
        # VCF/BCF decoders can split a String field at the commas in the
        # original ALT list, even when its header declares Number=1. Rejoin
        # only actual strings; the complete marker is validated below.
        if not marker or any(not isinstance(part, str) for part in marker):
            raise ValueError('invalid original-record tag tuple')
        marker = ','.join(marker)
    if marker is None:
        if (normalized_record.contig, normalized_record.pos, normalized_record.ref, list(alts)) != (
                context['chrom'], context['pos'], context['ref'], context['alts']):
            raise ValueError('untagged normalized record differs from biallelic source')
        return 1
    if not isinstance(marker, str) or not marker:
        raise ValueError('unsupported original-record tag type')
    fields = marker.split('|')
    if len(fields) == 4 and ',' not in fields[3]:
        fields.append('1')
    if len(fields) != 5:
        raise ValueError('unsupported original-record tag')
    chrom, position, ref, alternatives, index = fields
    try:
        position, index = int(position), int(index)
    except ValueError as error:
        raise ValueError('invalid original-record tag integer') from error
    if (chrom, position, ref, alternatives.split(',')) != (
            context['chrom'], context['pos'], context['ref'], context['alts']):
        raise ValueError('original-record tag does not match native source')
    if normalized_record.contig != chrom or not 1 <= index <= len(context['alts']):
        raise ValueError('invalid original ALT index or chromosome')
    return index


def allele_call(context, index):
    if type(index) is not int or not 1 <= index <= len(context['alts']):
        raise ValueError('invalid source ALT index')
    fields = context['FORMAT']
    gt = deepcopy(fields.get('GT') or [])
    ad = deepcopy(fields.get('AD'))
    separator = '|' if context['GT_phased'] else '/'
    result = {
        'source_allele': {'chrom': context['chrom'], 'pos': context['pos'],
                          'ref': context['ref'], 'alt': context['alts'][index - 1]},
        'source_all_alts': deepcopy(context['alts']), 'source_alt_index': index,
        'GT': separator.join('.' if value is None else str(value) for value in gt),
        'GT_values': gt, 'GT_phased': context['GT_phased'],
        'alt_copies': gt.count(index), 'called_ploidy': len(gt),
        'fully_called': bool(gt) and None not in gt,
        'alt_depth': ad[index] if isinstance(ad, list) and len(ad) > index else None,
        'FILTER': deepcopy(context['FILTER']), 'QUAL': context['QUAL'],
        'source_record_id': context['source_record_id'], 'caller': context['caller'],
        'source_namespace': context['source_namespace'], 'native_context': deepcopy(context),
    }
    result.update({key: deepcopy(fields.get(key)) for key in ('DP', 'GQ', 'AD', 'PL', 'PGT', 'PID', 'PS')})
    return result


def _numeric_axis(call, name, threshold):
    value = call[name]
    present = name in call['native_context']['FORMAT']
    if value is None:
        state = 'missing'
    elif type(value) not in (int, float) or not math.isfinite(value) or value < 0:
        state = 'invalid_review'
    elif value == 0:
        state = 'explicit_zero'
    elif value < threshold:
        state = 'below_historical_threshold'
    else:
        state = 'at_or_above_historical_threshold'
    return {'state': state, 'value': value, 'historical_threshold': threshold,
            'FORMAT_presence': 'absent' if not present else 'present_missing' if value is None else 'present_value',
            'interpretation': 'Descriptive historical bin only; not an eligibility or validated quality threshold.'}


def review_axes(call, duplicate=False):
    gt, alternatives = call['GT_values'], call['source_all_alts']
    reasons = []
    if not call['fully_called']:
        genotype = 'partial_or_missing'
        reasons.append('partial_or_missing_GT')
    elif len(gt) == 1:
        genotype = 'haploid'
    elif len(gt) == 2:
        genotype = 'fully_called_diploid'
    else:
        genotype = 'other_ploidy'
    if len(gt) != 2:
        reasons.append('not_diploid')
    nonconcrete = sorted({value for value in gt if value is not None and value > 0
                         and not _concrete(alternatives[value - 1])})
    if nonconcrete:
        reasons.append('nonconcrete_allele_in_GT')
    if not _concrete(call['source_allele']['ref']) or not _concrete(call['source_allele']['alt']):
        reasons.append('nonconcrete_reviewed_allele')
    if call['alt_copies'] == 0:
        reasons.append('ALT_not_observed_in_GT')
    if duplicate:
        reasons.append('duplicate_normalized_key')
    filters = call['FILTER']
    return {
        'filter_state': 'explicit_PASS' if filters == ['PASS'] else 'unfiltered' if not filters else 'named_filters',
        'filter_values': deepcopy(filters), 'DP': _numeric_axis(call, 'DP', 10),
        'GQ': _numeric_axis(call, 'GQ', 20), 'genotype_state': genotype,
        'allele_participation': 'observed_in_GT' if call['alt_copies'] else (
            'not_observed_in_partial_or_missing_GT' if not call['fully_called'] else 'ALT_absent_from_GT'),
        'nonconcrete_in_GT': bool(nonconcrete), 'nonconcrete_GT_allele_indices': nonconcrete,
        'duplicate_normalized_key': bool(duplicate),
        'structural_model_eligible': not reasons, 'structural_ineligibility_reasons': reasons,
    }


def _phase_tokens(call):
    """Keep PS and PID independently; agreement must be checked at pair level."""
    gt, index = call['GT_values'], call['source_alt_index']
    tokens, issues = {}, []
    if len(gt) != 2 or None in gt or gt.count(index) != 1:
        return tokens, issues
    namespace = call['source_namespace']
    pgt, pid, ps = call['PGT'], call['PID'], call['PS']
    if pgt not in (None, '.', '') or pid not in (None, '.', ''):
        if isinstance(pgt, str) and re.fullmatch(r'\d+\|\d+', pgt) and pid not in (None, '.', ''):
            phased = [int(value) for value in pgt.split('|')]
            if sorted(phased) == sorted(gt):
                tokens['PID'] = (namespace, str(pid), phased.index(index))
            else:
                issues.append('PGT_does_not_match_GT')
        else:
            issues.append('incomplete_or_invalid_PGT_PID')
    if call['GT_phased'] and ps not in (None, '.', ''):
        tokens['PS'] = (namespace, str(ps), gt.index(index))
    return tokens, issues


def _phase_relation(first, second):
    if first['source_record_id'] == second['source_record_id']:
        if (first['GT_values'] != second['GT_values'] or first['source_all_alts'] != second['source_all_alts']
                or first['source_alt_index'] == second['source_alt_index']
                or first['alt_copies'] != 1 or second['alt_copies'] != 1):
            raise ValueError('inconsistent allele rows sharing one native record')
        return 'opposite_alleles_in_source_GT', {'basis': 'single fully called diploid native GT',
                                              'issues': [], 'parental_phase_established': False}
    x, xi = _phase_tokens(first)
    y, yi = _phase_tokens(second)
    relations = {}
    for kind in ('PS', 'PID'):
        if kind in x and kind in y and x[kind][:2] == y[kind][:2]:
            relations[kind] = 'cis' if x[kind][2] == y[kind][2] else 'opposite'
    issues = sorted(set(xi + yi))
    if len(set(relations.values())) > 1:
        issues.append('conflicting_PS_PID_pair_relations')
    detail = {'compared_source_block_relations': relations, 'issues': issues,
              'parental_phase_established': False}
    if issues or not relations:
        return 'unknown', detail
    return ('cis_in_source_block' if next(iter(relations.values())) == 'cis'
            else 'opposite_alleles_in_source_block'), detail


def _overlapping(first, second):
    x, y = first['allele'], second['allele']
    return x['chrom'] == y['chrom'] and max(x['pos'], y['pos']) < min(
        x['pos'] + len(x['ref']), y['pos'] + len(y['ref']))


def _frequency_eligible(variant, maximum):
    frequency = variant['frequency']
    value = frequency['max_AF_fraction']
    status = frequency['status']
    if status not in ('observed', 'unresolved'):
        raise ValueError('unsupported frequency state')
    if status == 'unresolved':
        if value is not None:
            raise ValueError('unresolved frequency must remain null')
        return True
    if type(value) not in (int, float) or not math.isfinite(value) or not 0 <= value <= 1:
        raise ValueError('observed frequency requires a valid numeric fraction')
    return value <= maximum


def _model_row(gene, first, second, phase, phase_evidence):
    shared = sorted(set(first['genes'][gene]) & set(second['genes'][gene]))
    if not shared:
        raise ValueError('model without shared affected transcript')
    totals = {}
    for transcript in shared:
        groups = [row['genes'][gene][transcript]['group'] for row in (first, second)]
        if any(type(group) is not int or group not in (1, 2) for group in groups):
            raise ValueError('unsupported coding consequence group')
        totals[transcript] = sum(groups)
    best = max(totals.values())
    unique = [first] if first is second else [first, second]
    return {
        'caller': first['source_call']['caller'], 'gene': gene, 'chrom': first['allele']['chrom'],
        'tier': 5 - best, 'tier_interpretation': 'Descriptive consequence tier; not pathogenicity or probability.',
        'frequency_group': 'observed_rare' if all(row['frequency']['status'] == 'observed' for row in unique)
                           else 'frequency_unresolved',
        'phase': phase, 'phase_evidence': phase_evidence,
        'allele_keys': [row['key'] for row in unique],
        'source_record_ids': [row['source_call']['source_record_id'] for row in unique],
        'shared_transcripts_at_best_tier': [name for name in shared if totals[name] == best],
        'shared_transcript_consequence_totals': totals,
        'allele_evidence': [{key: deepcopy(row[key]) for key in ('key', 'allele', 'source_call', 'genes', 'frequency')}
                            | {'review_axes': review_axes(row['source_call'], bool(row.get('duplicate_normalized_key')))}
                           for row in unique],
    }


def build_models(variants, maximum_af):
    """Return descriptive models, retained pair-review rows and accounting.

No FILTER, DP, GQ, QUAL or allele-balance threshold determines eligibility.
Each caller remains separate, including when two callers share an allele key.
"""
    if type(maximum_af) not in (int, float) or not math.isfinite(maximum_af) or not 0 <= maximum_af <= 1:
        raise ValueError('invalid maximum AF fraction')
    groups, counters = defaultdict(list), Counter()
    identities = set()
    for variant in variants:
        counters['input_allele_rows'] += 1
        call = variant['source_call']
        identity = (call['source_record_id'], call['source_alt_index'])
        if identity in identities:
            raise ValueError('duplicate allele-row identity; aggregate annotations before modelling')
        identities.add(identity)
        axes = review_axes(call, bool(variant.get('duplicate_normalized_key')))
        if not axes['structural_model_eligible']:
            counters['structurally_ineligible_allele_rows'] += 1
            counters.update('structural_reason_' + reason for reason in axes['structural_ineligibility_reasons'])
            continue
        if not _frequency_eligible(variant, maximum_af):
            counters['observed_frequency_above_threshold_allele_rows'] += 1
            continue
        counters['retained_' + variant['frequency']['status'] + '_allele_rows'] += 1
        for gene in variant['genes']:
            groups[(call['caller'], variant['allele']['chrom'], gene)].append(variant)
    models, reviews = [], []
    for (_, _, gene), rows in groups.items():
        homo = [row for row in rows if row['source_call']['alt_copies'] == 2]
        hetero = [row for row in rows if row['source_call']['alt_copies'] == 1]
        for row in homo:
            models.append(_model_row(gene, row, row, 'homozygous_ALT_in_source_GT',
                                     {'basis': 'single fully called diploid native GT', 'issues': [],
                                      'parental_phase_established': False}))
        for first, second in itertools.combinations(hetero, 2):
            if not set(first['genes'][gene]) & set(second['genes'][gene]):
                counters['heterozygous_pairs_without_shared_transcript'] += 1
                continue
            phase, evidence = _phase_relation(first['source_call'], second['source_call'])
            row = _model_row(gene, first, second, phase, evidence)
            reasons = []
            if phase == 'cis_in_source_block':
                reasons.append('explicit_source_cis')
            if _overlapping(first, second) and phase != 'opposite_alleles_in_source_GT':
                reasons.append('overlapping_representations')
            if reasons:
                row['review_reasons'] = reasons
                reviews.append(row)
                counters.update('review_pair_' + reason for reason in reasons)
            else:
                models.append(row)
    order = lambda row: (row['caller'], row['chrom'], row['gene'], row['allele_keys'], row['source_record_ids'])
    models.sort(key=order)
    reviews.sort(key=order)
    counters['descriptive_models'] = len(models)
    counters['retained_review_pairs'] = len(reviews)
    return models, reviews, dict(counters)
