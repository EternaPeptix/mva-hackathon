#!/usr/bin/env python3
"""Pure sequence-hypothesis and raw-quality read/fragment evidence helpers.

Enumerated sequences are hypotheses, including uncalled retained ALT alleles.
Neither enumeration nor an absent spanning observation asserts a genotype.
"""
from collections import Counter, defaultdict
from copy import deepcopy
import itertools
import re


EXCLUDED_FLAGS = 0xF04


def _dna(sequence, name, allow_empty=False):
    if not isinstance(sequence, str) or (not sequence and not allow_empty) or (
            sequence and re.fullmatch('[ACGT]+', sequence) is None):
        raise ValueError(name + ' must be concrete uppercase sequence')


def _integer(value, name, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError('invalid ' + name)


def enumerate_haplotypes(chrom, start0, reference, edits, max_states=256,
                         callers=None, max_haplotype_bases=100000):
    """Enumerate caller-specific compatible subsets and collapse exact sequences.

Each edit has caller, source_record_id, source_alt_index, pos (one-based), ref,
alt and called (bool). At most one ALT from each native source record is chosen.
Any intersecting REF spans within a chosen combination are refused. Caller
combinations never mix, and every retained source-label combination is preserved.
The state cap bounds all attempted Cartesian combinations, including rejected
overlaps, before enumeration; the base cap bounds each materialized sequence.
"""
    _integer(start0, 'window start')
    _integer(max_states, 'maximum states', 1)
    _integer(max_haplotype_bases, 'maximum haplotype bases', 1)
    _dna(reference, 'reference window')
    if not isinstance(chrom, str) or not chrom or len(reference) > max_haplotype_bases:
        raise ValueError('invalid chromosome or reference window exceeds base cap')
    edits = list(edits)
    declared = set(callers or ())
    if any(not isinstance(caller, str) or not caller for caller in declared):
        raise ValueError('invalid explicit caller')
    by_caller = defaultdict(lambda: defaultdict(list))
    seen = set()
    for edit in edits:
        caller, record_id = edit['caller'], edit['source_record_id']
        if not isinstance(caller, str) or not caller or not isinstance(record_id, str) or not record_id:
            raise ValueError('invalid edit source identity')
        if callers is not None and caller not in declared:
            raise ValueError('edit caller absent from explicit caller list')
        declared.add(caller)
        _integer(edit['source_alt_index'], 'source ALT index', 1)
        _integer(edit['pos'], 'edit position', 1)
        if type(edit['called']) is not bool:
            raise ValueError('edit called field must explicitly be boolean')
        _dna(edit['ref'], 'edit REF')
        _dna(edit['alt'], 'edit ALT')
        if edit.get('chrom', chrom) != chrom:
            raise ValueError('edit belongs to another chromosome')
        offset = edit['pos'] - 1 - start0
        if offset < 0 or offset + len(edit['ref']) > len(reference):
            raise ValueError('edit REF extends outside the explicit window')
        if reference[offset:offset + len(edit['ref'])] != edit['ref']:
            raise ValueError('edit REF differs from the explicit reference')
        identity = (caller, record_id, edit['source_alt_index'])
        if identity in seen:
            raise ValueError('duplicate source ALT identity')
        seen.add(identity)
        by_caller[caller][record_id].append(deepcopy(edit))
    if not declared:
        raise ValueError('at least one explicit or inferred caller is required')
    estimated = 0
    for caller in declared:
        states = 1
        for options in by_caller[caller].values():
            states *= 1 + len(options)
            if states > max_states:
                raise ValueError('hypothesis state cap exceeded')
        estimated += states
        if estimated > max_states:
            raise ValueError('hypothesis state cap exceeded')
    sequences, rejected_combinations = {}, []
    counters = Counter(attempted_combinations=0, compatible_combinations=0,
                       conflicting_REF_span_combinations=0)
    for caller in sorted(declared):
        choices = [[None] + sorted(options, key=lambda edit: edit['source_alt_index'])
                   for _, options in sorted(by_caller[caller].items())]
        for combination in itertools.product(*choices):
            counters['attempted_combinations'] += 1
            chosen = sorted((edit for edit in combination if edit is not None),
                            key=lambda edit: (edit['pos'], edit['source_record_id'], edit['source_alt_index']))
            if any(first['pos'] - 1 + len(first['ref']) > second['pos'] - 1
                   for first, second in zip(chosen, chosen[1:])):
                counters['conflicting_REF_span_combinations'] += 1
                rejected_combinations.append({'caller': caller, 'edits': deepcopy(chosen),
                                              'reason': 'conflicting_reference_spans'})
                continue
            size = len(reference) + sum(len(edit['alt']) - len(edit['ref']) for edit in chosen)
            if size > max_haplotype_bases:
                raise ValueError('haplotype sequence exceeds base cap')
            pieces, cursor = [], 0
            for edit in chosen:
                offset = edit['pos'] - 1 - start0
                pieces.extend((reference[cursor:offset], edit['alt']))
                cursor = offset + len(edit['ref'])
            pieces.append(reference[cursor:])
            sequence = ''.join(pieces)
            counters['compatible_combinations'] += 1
            row = sequences.setdefault(sequence, {'sequence': sequence, 'source_combinations': []})
            row['source_combinations'].append({'caller': caller, 'edits': deepcopy(chosen),
                                               'reference_only': not chosen,
                                               'contains_uncalled_ALT': any(not edit['called'] for edit in chosen)})
    counters['distinct_sequences'] = len(sequences)
    return {'chrom': chrom, 'start0': start0, 'end0': start0 + len(reference),
            'reference': reference, 'callers': sorted(declared),
            'haplotypes': [sequences[sequence] for sequence in sorted(sequences)],
            'rejected_combinations': rejected_combinations,
            'counters': dict(counters),
            'limits': {'max_states': max_states, 'max_haplotype_bases': max_haplotype_bases},
            'interpretation': 'Compatible sequence hypotheses, including retained uncalled ALTs; not inferred genotypes or phased haplotypes.'}


def minimal_replacement(chrom, start0, reference, alternate, *, left_anchor=None, right_anchor=None):
    """Derive a parsimonious concrete replacement and verify its complete edit.

External anchors, when supplied, are the immediate single reference bases before
or after the explicit window. A whole-window deletion requires an available
external anchor. No unavailable base or contig context is invented. Parsimony
here is not repeat-aware whole-reference left normalization.
"""
    _integer(start0, 'window start')
    _dna(reference, 'reference window')
    _dna(alternate, 'alternate window', allow_empty=True)
    if not isinstance(chrom, str) or not chrom:
        raise ValueError('invalid chromosome')
    for name, anchor in (('left anchor', left_anchor), ('right anchor', right_anchor)):
        if anchor is not None:
            _dna(anchor, name)
            if len(anchor) != 1:
                raise ValueError(name + ' must be one reference base')
    if start0 == 0 and left_anchor is not None:
        raise ValueError('left anchor before contig start')
    if reference == alternate:
        return None
    prefix = 0
    while prefix < min(len(reference), len(alternate)) and reference[prefix] == alternate[prefix]:
        prefix += 1
    suffix = 0
    while (suffix < min(len(reference) - prefix, len(alternate) - prefix)
           and reference[len(reference) - 1 - suffix] == alternate[len(alternate) - 1 - suffix]):
        suffix += 1
    ref_end = len(reference) - suffix
    alt_end = len(alternate) - suffix
    ref, alt = reference[prefix:ref_end], alternate[prefix:alt_end]
    position0 = start0 + prefix
    if not ref or not alt:
        if prefix:
            anchor = reference[prefix - 1]
            ref, alt, position0 = anchor + ref, anchor + alt, position0 - 1
        elif left_anchor is not None:
            ref, alt, position0 = left_anchor + ref, left_anchor + alt, start0 - 1
        elif suffix:
            anchor = reference[ref_end]
            ref, alt = ref + anchor, alt + anchor
        elif right_anchor is not None:
            ref, alt = ref + right_anchor, alt + right_anchor
        else:
            raise ValueError('unanchorable whole-window or entire-contig deletion')
    _dna(ref, 'derived REF')
    _dna(alt, 'derived ALT')
    extended = (left_anchor or '') + reference + (right_anchor or '')
    expected = (left_anchor or '') + alternate + (right_anchor or '')
    extended_start = start0 - int(left_anchor is not None)
    offset = position0 - extended_start
    if offset < 0 or extended[offset:offset + len(ref)] != ref:
        raise ValueError('derived allele REF or coordinate does not match window')
    observed = extended[:offset] + alt + extended[offset + len(ref):]
    if observed != expected:
        raise ValueError('derived allele does not recover alternate window')
    return {'chrom': chrom, 'pos': position0 + 1, 'ref': ref, 'alt': alt}


def extract_spanning_read(read, start0, end0, min_bq):
    """Extract an exact window sequence using raw BQ and explicit CIGAR rules.

Both boundary reference positions must align to query bases. Insertions exactly
at either edge are retained as boundary-insertion rejections, avoiding an implicit
ownership convention. Internal insertions/deletions are represented literally;
reference skips, missing sequence/qualities and malformed CIGAR are rejected.
Flag/MQ/RG/QNAME eligibility is enforced by aggregate_fragment_evidence().
"""
    _integer(start0, 'window start')
    _integer(end0, 'window end', start0 + 1)
    _integer(min_bq, 'minimum raw base quality')
    result = {'status': 'rejected', 'reason': None, 'sequence': None,
              'query_indices': [], 'base_qualities': [], 'minimum_raw_BQ': None,
              'BQ_threshold': min_bq, 'BAQ_applied': False,
              'boundary_policy': 'aligned query base at each reference edge; edge insertions rejected'}
    def reject(reason):
        result['reason'] = reason
        return result
    sequence, qualities, cigar = read.query_sequence, read.query_qualities, read.cigartuples
    if sequence is None:
        return reject('missing_query_sequence')
    if qualities is None:
        return reject('missing_base_qualities')
    if len(sequence) != len(qualities):
        return reject('sequence_quality_length_mismatch')
    if not cigar or read.reference_start is None or read.reference_start < 0:
        return reject('missing_aligned_CIGAR')
    nonclips = [index for index, (operation, _) in enumerate(cigar) if operation not in (4, 5)]
    if not nonclips:
        return reject('missing_aligned_CIGAR')
    if any(operation in (4, 5) for operation, _ in cigar[nonclips[0]:nonclips[-1] + 1]):
        return reject('internal_soft_or_hard_clip')
    if read.reference_start > start0:
        return reject('does_not_span_window_start')
    query, position = 0, read.reference_start
    indices, boundaries = [], set()
    refskip = edge_insertion = False
    for operation, length in cigar:
        if type(length) is not int or length <= 0 or operation not in range(9):
            return reject('invalid_CIGAR_operation')
        if operation in (0, 7, 8):
            left, right = max(position, start0), min(position + length, end0)
            if left < right:
                indices.extend(range(query + left - position, query + right - position))
                if left <= start0 < right:
                    boundaries.add(start0)
                if left <= end0 - 1 < right:
                    boundaries.add(end0 - 1)
            query += length
            position += length
        elif operation == 1:
            if position in (start0, end0):
                edge_insertion = True
            elif start0 < position < end0:
                indices.extend(range(query, query + length))
            query += length
        elif operation in (2, 3):
            if operation == 3 and max(position, start0) < min(position + length, end0):
                refskip = True
            position += length
        elif operation == 4:
            query += length
        # Hard clip and pad consume neither sequence nor reference coordinates.
    if query != len(sequence) or any(index < 0 or index >= len(sequence) for index in indices):
        return reject('CIGAR_query_length_mismatch')
    if position < end0:
        return reject('does_not_span_window_end')
    if refskip:
        return reject('reference_skip_in_window')
    if edge_insertion:
        return reject('boundary_insertion')
    if boundaries != {start0, end0 - 1}:
        return reject('reference_boundary_not_aligned_to_query')
    result['sequence'] = ''.join(sequence[index] for index in indices).upper()
    result['query_indices'] = indices
    result['base_qualities'] = [int(qualities[index]) for index in indices]
    result['minimum_raw_BQ'] = min(result['base_qualities']) if indices else None
    if any(value < 0 or value > 254 for value in result['base_qualities']):
        return reject('invalid_base_quality')
    if any(base not in 'ACGT' for base in result['sequence']):
        return reject('nonconcrete_query_sequence')
    if result['minimum_raw_BQ'] is None or result['minimum_raw_BQ'] < min_bq:
        return reject('raw_base_quality_below_threshold')
    result.update(status='accepted', reason=None)
    return result


def aggregate_fragment_evidence(reads, start0, end0, min_bq, min_mq, expected_sequences=None):
    """Retain per-observation evidence and collapse only matching (RG,QNAME).

expected_sequences is a mapping from exact sequence to arbitrary JSON-safe
provenance labels, or an iterable of exact sequences. Conflicting eligible mates
are retained as conflicts and never assigned to one sequence. All raw input
observations remain distinguishable from fragments and from accepted observations.
"""
    _integer(min_mq, 'minimum mapping quality')
    if min_mq > 254:
        raise ValueError('mapping quality 255 is unavailable, not high quality')
    _integer(min_bq, 'minimum raw base quality')
    _integer(start0, 'window start')
    _integer(end0, 'window end', start0 + 1)
    expected = ({sequence: deepcopy(labels) for sequence, labels in expected_sequences.items()}
                if isinstance(expected_sequences, dict)
                else {sequence: [] for sequence in (expected_sequences or ())})
    for sequence in expected:
        _dna(sequence, 'expected haplotype sequence', allow_empty=True)
    observations, by_fragment = [], defaultdict(list)
    counts, rejection_counts = Counter(), Counter()
    observation_sequences, fragment_sequences = Counter(), Counter()
    strands, accepted_strands, read_groups = Counter(), Counter(), defaultdict(Counter)
    for ordinal, read in enumerate(reads, 1):
        counts['raw_observations'] += 1
        try:
            rg = read.get_tag('RG')
        except (KeyError, ValueError):
            rg = None
        qname = read.query_name
        strand = 'reverse' if read.is_reverse else 'forward'
        strands[strand] += 1
        observation = {'observation_id': ordinal, 'RG': rg, 'QNAME': qname,
                       'flag': read.flag, 'mapping_quality': read.mapping_quality,
                       'strand': strand, 'is_read1': bool(read.is_read1), 'is_read2': bool(read.is_read2),
                       'reference_start0': read.reference_start, 'reference_end0': read.reference_end,
                       'CIGAR': [list(pair) for pair in (read.cigartuples or ())],
                       'status': 'rejected', 'reason': None, 'sequence': None,
                       'matched_hypothesis': False, 'hypothesis_labels': []}
        if not isinstance(rg, str) or not rg:
            observation['reason'] = 'missing_or_invalid_RG'
        elif not isinstance(qname, str) or not qname:
            observation['reason'] = 'missing_or_invalid_QNAME'
        elif read.flag & EXCLUDED_FLAGS:
            observation['reason'] = 'excluded_flag_0xF04'
        elif read.mapping_quality == 255:
            observation['reason'] = 'mapping_quality_unavailable_255'
        elif read.mapping_quality < min_mq:
            observation['reason'] = 'mapping_quality_below_threshold'
        else:
            observation.update(extract_spanning_read(read, start0, end0, min_bq))
        if isinstance(rg, str) and rg:
            read_groups[rg]['raw_observations'] += 1
        if isinstance(rg, str) and rg and isinstance(qname, str) and qname:
            by_fragment[(rg, qname)].append(observation)
        if observation['status'] == 'accepted':
            sequence = observation['sequence']
            observation['matched_hypothesis'] = sequence in expected
            observation['hypothesis_labels'] = deepcopy(expected.get(sequence, []))
            counts['accepted_observations'] += 1
            counts['matched_observations' if sequence in expected else 'unmatched_observations'] += 1
            observation_sequences[sequence] += 1
            accepted_strands[strand] += 1
            read_groups[rg]['accepted_observations'] += 1
        else:
            counts['rejected_observations'] += 1
            rejection_counts[observation['reason']] += 1
        observations.append(observation)
    fragments = []
    for (rg, qname), raw in sorted(by_fragment.items()):
        accepted = [observation for observation in raw if observation['status'] == 'accepted']
        sequences = sorted({observation['sequence'] for observation in accepted})
        row = {'RG': rg, 'QNAME': qname, 'observation_ids': [observation['observation_id'] for observation in raw],
               'accepted_observation_ids': [observation['observation_id'] for observation in accepted],
               'accepted_observations': len(accepted), 'sequences': sequences,
               'sequence': sequences[0] if len(sequences) == 1 else None,
               'sequence_hypotheses': [{'sequence': sequence, 'matched_hypothesis': sequence in expected,
                                         'hypothesis_labels': deepcopy(expected.get(sequence, []))}
                                        for sequence in sequences],
               'observation_sequence_counts': dict(Counter(observation['sequence'] for observation in accepted)),
               'strand_counts': dict(Counter(observation['strand'] for observation in accepted)),
               'status': 'no_eligible_spanning_observation', 'matched_hypothesis': False,
               'hypothesis_labels': []}
        counts['fragments_with_identity'] += 1
        read_groups[rg]['fragments_with_identity'] += 1
        if not sequences:
            counts['fragments_without_eligible_span'] += 1
        elif len(sequences) > 1:
            row['status'] = 'mate_or_repeated_observation_conflict'
            counts['conflicting_fragments'] += 1
            read_groups[rg]['conflicting_fragments'] += 1
        else:
            sequence = sequences[0]
            row['status'] = 'matched' if sequence in expected else 'unmatched'
            row['matched_hypothesis'] = sequence in expected
            row['hypothesis_labels'] = deepcopy(expected.get(sequence, []))
            counts['nonconflicting_spanning_fragments'] += 1
            counts['matched_fragments' if sequence in expected else 'unmatched_fragments'] += 1
            fragment_sequences[sequence] += 1
            read_groups[rg]['nonconflicting_spanning_fragments'] += 1
        fragments.append(row)
    for name in ('raw_observations', 'accepted_observations', 'rejected_observations',
                 'matched_observations', 'unmatched_observations', 'fragments_with_identity',
                 'fragments_without_eligible_span', 'conflicting_fragments',
                 'nonconflicting_spanning_fragments', 'matched_fragments', 'unmatched_fragments'):
        counts.setdefault(name, 0)
    return {'window': {'start0': start0, 'end0': end0},
            'thresholds': {'min_raw_BQ': min_bq, 'min_MQ': min_mq, 'excluded_flags': EXCLUDED_FLAGS,
                           'reject_MQ_255': True, 'BAQ_applied': False},
            'counts': dict(counts), 'rejection_counts': dict(rejection_counts),
            'observation_sequence_counts': dict(observation_sequences),
            'fragment_sequence_counts': dict(fragment_sequences),
            'observation_strand_counts': dict(strands),
            'accepted_observation_strand_counts': dict(accepted_strands),
            'read_group_counts': {rg: dict(counter) for rg, counter in sorted(read_groups.items())},
            'observations': observations, 'fragments': fragments,
            'interpretation': 'Exact raw-quality sequence evidence in this window; observations and fragments are distinct. No spanning evidence is not a reference genotype or proof of absence.'}
