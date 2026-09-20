"""Small fictional sequence/fragment controls; no case coordinates or data."""
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from joint_haplotype_helpers import (aggregate_fragment_evidence, enumerate_haplotypes,
                                      extract_spanning_read, minimal_replacement)


def edit(caller, record, pos, ref, alt, index=1, called=True):
    return {'caller': caller, 'source_record_id': record, 'source_alt_index': index,
            'pos': pos, 'ref': ref, 'alt': alt, 'called': called}


class Read:
    def __init__(self, sequence, cigar=None, start=0, qname='fragment', rg='RG1',
                 quality=30, mq=60, flag=0):
        self.query_sequence = sequence
        self.query_qualities = [quality] * len(sequence)
        self.cigartuples = cigar or [(0, len(sequence))]
        self.reference_start = start
        self.reference_end = start + sum(size for op, size in self.cigartuples if op in (0, 2, 3, 7, 8))
        self.query_name, self.mapping_quality, self.flag, self.rg = qname, mq, flag, rg
        self.is_reverse = bool(flag & 16)
        self.is_read1, self.is_read2 = bool(flag & 64), bool(flag & 128)

    def get_tag(self, tag):
        if tag != 'RG' or self.rg is None:
            raise KeyError(tag)
        return self.rg


class JointHaplotypeTests(unittest.TestCase):
    def test_insertion_deletion_equals_four_SNPs_without_duplicate_fragment_support(self):
        reference = 'A' * 20
        edits = [edit('callerA', 'A1', 5, 'A', 'ACGTC'),
                 edit('callerA', 'A2', 8, 'AAAAA', 'A')]
        edits += [edit('callerB', 'B%d' % pos, pos, 'A', alt)
                  for pos, alt in zip((6, 7, 8, 9), 'CGTC')]
        result = enumerate_haplotypes('SYNTH', 0, reference, edits)
        sequence = 'A' * 5 + 'CGTC' + 'A' * 11
        same = next(row for row in result['haplotypes'] if row['sequence'] == sequence)
        self.assertEqual([(row['caller'], len(row['edits'])) for row in same['source_combinations']],
                         [('callerA', 2), ('callerB', 4)])
        self.assertEqual(minimal_replacement('SYNTH', 0, reference, sequence),
                         {'chrom': 'SYNTH', 'pos': 6, 'ref': 'AAAA', 'alt': 'CGTC'})
        evidence = aggregate_fragment_evidence(
            [Read(sequence, flag=64), Read(sequence, flag=128 | 16)],
            0, 20, 20, 20, {sequence: same['source_combinations']})
        self.assertEqual(evidence['counts']['accepted_observations'], 2)
        self.assertEqual(evidence['counts']['matched_fragments'], 1)
        self.assertEqual(len(evidence['fragments'][0]['hypothesis_labels']), 2)
        self.assertEqual(evidence['fragment_sequence_counts'][sequence], 1)

    def test_same_source_multiALT_exclusivity_and_uncalled_provenance(self):
        result = enumerate_haplotypes('SYNTH', 0, 'AAAA', [
            edit('A', 'source1', 2, 'A', 'C', index=1),
            edit('A', 'source1', 2, 'A', 'G', index=2, called=False)])
        self.assertEqual(result['counters']['attempted_combinations'], 3)
        self.assertEqual({row['sequence'] for row in result['haplotypes']}, {'AAAA', 'ACAA', 'AGAA'})
        uncalled = next(row for row in result['haplotypes'] if row['sequence'] == 'AGAA')
        self.assertTrue(uncalled['source_combinations'][0]['contains_uncalled_ALT'])
        self.assertFalse(uncalled['source_combinations'][0]['edits'][0]['called'])

    def test_conflicting_REF_spans_are_not_combined_and_callers_never_mix(self):
        changes = [edit('A', 'r1', 2, 'AAA', 'A'), edit('A', 'r2', 3, 'A', 'C')]
        result = enumerate_haplotypes('SYNTH', 0, 'AAAAAA', changes)
        self.assertEqual(result['counters']['conflicting_REF_span_combinations'], 1)
        self.assertEqual(result['counters']['compatible_combinations'], 3)
        self.assertEqual(result['rejected_combinations'], [{
            'caller': 'A', 'edits': changes, 'reason': 'conflicting_reference_spans'}])
        self.assertEqual([row['source_record_id'] for row in result['rejected_combinations'][0]['edits']],
                         ['r1', 'r2'])
        changes[1]['caller'] = 'B'
        result = enumerate_haplotypes('SYNTH', 0, 'AAAAAA', changes)
        self.assertEqual(result['counters']['conflicting_REF_span_combinations'], 0)
        self.assertEqual(result['rejected_combinations'], [])
        for row in result['haplotypes']:
            self.assertTrue(all(len(group['edits']) <= 1 for group in row['source_combinations']))

    def test_same_sequence_preserves_duplicate_source_labels_and_empty_caller(self):
        result = enumerate_haplotypes('SYNTH', 0, 'AAAA', [
            edit('A', 'r1', 2, 'A', 'C'), edit('A', 'r2', 2, 'A', 'C')], callers=['A', 'B'])
        alternate = next(row for row in result['haplotypes'] if row['sequence'] == 'ACAA')
        self.assertEqual({group['edits'][0]['source_record_id'] for group in alternate['source_combinations']}, {'r1', 'r2'})
        reference = next(row for row in result['haplotypes'] if row['sequence'] == 'AAAA')
        self.assertEqual([group['caller'] for group in reference['source_combinations']], ['A', 'B'])
        self.assertEqual(len(enumerate_haplotypes('SYNTH', 0, 'AAAA', [], callers=['B'])['haplotypes']), 1)

    def test_ref_bounds_empty_reference_and_state_base_guards(self):
        for changes in ([edit('A', 'r1', 2, 'C', 'G')], [edit('A', 'r1', 4, 'AA', 'A')]):
            with self.assertRaises(ValueError):
                enumerate_haplotypes('SYNTH', 0, 'AAAA', changes)
        with self.assertRaises(ValueError):
            enumerate_haplotypes('SYNTH', 0, '', [], callers=['A'])
        changes = [edit('A', 'r%d' % pos, pos, 'A', 'C') for pos in range(1, 5)]
        with self.assertRaises(ValueError):
            enumerate_haplotypes('SYNTH', 0, 'AAAA', changes, max_states=15)
        with self.assertRaises(ValueError):
            enumerate_haplotypes('SYNTH', 0, 'AAAA', [edit('A', 'r1', 2, 'A', 'ACCCCC')], max_haplotype_bases=8)
        with self.assertRaises(ValueError):
            enumerate_haplotypes('SYNTH', 0, 'AAAA', [changes[0], changes[0]])

    def test_minimal_substitution_insertion_deletion_and_identity(self):
        self.assertIsNone(minimal_replacement('SYNTH', 10, 'ACGT', 'ACGT'))
        self.assertEqual(minimal_replacement('SYNTH', 10, 'ACGT', 'ATGT'),
                         {'chrom': 'SYNTH', 'pos': 12, 'ref': 'C', 'alt': 'T'})
        self.assertEqual(minimal_replacement('SYNTH', 10, 'ACGT', 'ACTGT'),
                         {'chrom': 'SYNTH', 'pos': 12, 'ref': 'C', 'alt': 'CT'})
        self.assertEqual(minimal_replacement('SYNTH', 10, 'ACGT', 'AGT'),
                         {'chrom': 'SYNTH', 'pos': 11, 'ref': 'AC', 'alt': 'A'})

    def test_left_right_anchoring_and_unanchorable_entire_contig_deletion(self):
        self.assertEqual(minimal_replacement('SYNTH', 0, 'ACG', 'CG'),
                         {'chrom': 'SYNTH', 'pos': 1, 'ref': 'AC', 'alt': 'C'})
        self.assertEqual(minimal_replacement('SYNTH', 0, 'ACG', 'TACG'),
                         {'chrom': 'SYNTH', 'pos': 1, 'ref': 'A', 'alt': 'TA'})
        self.assertEqual(minimal_replacement('SYNTH', 10, 'ACG', '', left_anchor='T'),
                         {'chrom': 'SYNTH', 'pos': 10, 'ref': 'TACG', 'alt': 'T'})
        self.assertEqual(minimal_replacement('SYNTH', 0, 'ACG', '', right_anchor='T'),
                         {'chrom': 'SYNTH', 'pos': 1, 'ref': 'ACGT', 'alt': 'T'})
        with self.assertRaises(ValueError):
            minimal_replacement('SYNTH', 0, 'ACG', '')
        with self.assertRaises(ValueError):
            minimal_replacement('SYNTH', 0, 'ACG', 'CG', left_anchor='T')

    def test_internal_CIGAR_insertions_and_deletions_are_literal(self):
        insertion = Read('AAAACCAAAAAA', cigar=[(0, 4), (1, 2), (0, 6)])
        result = extract_spanning_read(insertion, 2, 8, 20)
        self.assertEqual(result['status'], 'accepted')
        self.assertEqual(result['sequence'], 'AACCAAAA')
        deletion = Read('AAAAAAAA', cigar=[(0, 4), (2, 2), (0, 4)])
        result = extract_spanning_read(deletion, 2, 8, 20)
        self.assertEqual(result['status'], 'accepted')
        self.assertEqual(result['sequence'], 'AAAA')
        self.assertFalse(result['BAQ_applied'])

    def test_boundary_insertions_reference_skip_and_boundary_deletion_rejected(self):
        cases = [
            (Read('A' * 12, cigar=[(0, 2), (1, 2), (0, 8)]), 'boundary_insertion'),
            (Read('A' * 11, cigar=[(0, 8), (1, 1), (0, 2)]), 'boundary_insertion'),
            (Read('A' * 8, cigar=[(0, 4), (3, 2), (0, 4)]), 'reference_skip_in_window'),
            (Read('A' * 8, cigar=[(0, 2), (2, 2), (0, 6)]), 'reference_boundary_not_aligned_to_query'),
        ]
        for read, reason in cases:
            with self.subTest(reason=reason):
                self.assertEqual(extract_spanning_read(read, 2, 8, 20)['reason'], reason)

    def test_raw_BQ_threshold_rejection_retains_sequence_and_observation(self):
        read = Read('AACCGGTTAA', quality=25)
        accepted = aggregate_fragment_evidence([read], 2, 8, 20, 20, ['CCGGTT'])
        rejected = aggregate_fragment_evidence([read], 2, 8, 30, 20, ['CCGGTT'])
        self.assertEqual(accepted['counts']['matched_fragments'], 1)
        self.assertEqual(rejected['counts']['matched_fragments'], 0)
        self.assertEqual(rejected['observations'][0]['sequence'], 'CCGGTT')
        self.assertEqual(rejected['observations'][0]['base_qualities'], [25] * 6)
        self.assertEqual(rejected['fragments'][0]['status'], 'no_eligible_spanning_observation')
        self.assertFalse(rejected['thresholds']['BAQ_applied'])

    def test_flags_MQ255_missing_RG_and_QNAME_are_rejected(self):
        reads = [Read('AAAA', flag=flag, qname='flag%d' % flag) for flag in (4, 256, 512, 1024, 2048)]
        reads += [Read('AAAA', mq=255, qname='mq255'), Read('AAAA', rg=None), Read('AAAA', qname=None)]
        result = aggregate_fragment_evidence(reads, 0, 4, 20, 0, ['AAAA'])
        self.assertEqual(result['counts']['raw_observations'], 8)
        self.assertEqual(result['counts']['accepted_observations'], 0)
        self.assertEqual(result['rejection_counts']['excluded_flag_0xF04'], 5)
        self.assertEqual(result['rejection_counts']['mapping_quality_unavailable_255'], 1)
        self.assertEqual(result['rejection_counts']['missing_or_invalid_RG'], 1)
        self.assertEqual(result['rejection_counts']['missing_or_invalid_QNAME'], 1)

    def test_MQ_views_keep_exact_threshold_and_raw_counts(self):
        reads = [Read('AAAA', qname='q%d' % mq, mq=mq) for mq in (0, 20, 60, 255)]
        for minimum, count in ((0, 3), (20, 2), (60, 1)):
            result = aggregate_fragment_evidence(reads, 0, 4, 20, minimum, ['AAAA'])
            self.assertEqual(result['thresholds']['min_MQ'], minimum)
            self.assertEqual(result['counts']['raw_observations'], 4)
            self.assertEqual(result['counts']['matched_fragments'], count)

    def test_opposite_fragments_and_conflicting_mates_are_distinct(self):
        separate = aggregate_fragment_evidence([Read('AAAA', qname='q1'), Read('ACAA', qname='q2')],
                                               0, 4, 20, 20, ['AAAA', 'ACAA'])
        self.assertEqual(separate['fragment_sequence_counts'], {'AAAA': 1, 'ACAA': 1})
        conflict = aggregate_fragment_evidence([Read('AAAA', flag=64), Read('ACAA', flag=128)],
                                               0, 4, 20, 20, ['AAAA', 'ACAA'])
        self.assertEqual(conflict['counts']['accepted_observations'], 2)
        self.assertEqual(conflict['counts']['conflicting_fragments'], 1)
        self.assertEqual(conflict['counts']['matched_fragments'], 0)
        self.assertEqual(conflict['fragment_sequence_counts'], {})
        self.assertEqual(conflict['fragments'][0]['observation_ids'], [1, 2])
        self.assertEqual(conflict['fragments'][0]['sequences'], ['AAAA', 'ACAA'])

    def test_RG_namespace_unmatched_sequence_and_strands_preserved(self):
        result = aggregate_fragment_evidence([Read('AAAA', rg='one', flag=16), Read('ACAA', rg='two')],
                                             0, 4, 20, 20, {'AAAA': ['reference-hypothesis']})
        self.assertEqual(result['counts']['fragments_with_identity'], 2)
        self.assertEqual(result['counts']['matched_fragments'], 1)
        self.assertEqual(result['counts']['unmatched_fragments'], 1)
        self.assertEqual(result['fragment_sequence_counts'], {'AAAA': 1, 'ACAA': 1})
        self.assertEqual(result['observation_strand_counts'], {'reverse': 1, 'forward': 1})
        self.assertEqual(set(result['read_group_counts']), {'one', 'two'})

    def test_missing_qualities_malformed_CIGAR_and_no_spanning_evidence(self):
        missing = Read('AAAA')
        missing.query_qualities = None
        self.assertEqual(extract_spanning_read(missing, 0, 4, 20)['reason'], 'missing_base_qualities')
        malformed = Read('AAAA', cigar=[(0, 5)])
        self.assertEqual(extract_spanning_read(malformed, 0, 4, 20)['reason'], 'CIGAR_query_length_mismatch')
        result = aggregate_fragment_evidence([Read('AAA', start=1)], 0, 4, 20, 20, ['AAAA'])
        self.assertEqual(result['counts']['fragments_without_eligible_span'], 1)
        self.assertEqual(result['counts']['matched_fragments'], 0)
        self.assertIn('not a reference genotype', result['interpretation'])

    def test_terminal_clips_allowed_and_internal_clips_rejected(self):
        terminal = Read('CCAAAAAAAAAAGG', cigar=[(5, 3), (4, 2), (0, 10), (4, 2), (5, 3)])
        self.assertEqual(extract_spanning_read(terminal, 2, 8, 20)['sequence'], 'AAAAAA')
        internal = Read('A' * 12, cigar=[(0, 5), (4, 2), (0, 5)])
        self.assertEqual(extract_spanning_read(internal, 2, 8, 20)['reason'], 'internal_soft_or_hard_clip')

    def test_native_pysam_read_fields_and_CIGAR(self):
        try:
            import pysam
        except ImportError:
            self.skipTest('native pysam control requires the installed Studio runtime')
        read = pysam.AlignedSegment()
        read.query_name = 'fictional-fragment'
        read.query_sequence = 'AAAACCAAAAAA'
        read.query_qualities = [30] * 12
        read.reference_id, read.reference_start = 0, 0
        read.cigartuples = [(0, 4), (1, 2), (0, 6)]
        read.flag, read.mapping_quality = 64, 60
        read.set_tag('RG', 'fictional-RG')
        result = aggregate_fragment_evidence([read], 2, 8, 20, 20, {'AACCAAAA': ['insertion']})
        self.assertEqual(result['counts']['matched_fragments'], 1)
        self.assertEqual(result['observations'][0]['query_indices'], list(range(2, 10)))
        self.assertEqual(result['fragments'][0]['QNAME'], 'fictional-fragment')
        read.mapping_quality = 255
        result = aggregate_fragment_evidence([read], 2, 8, 20, 0, ['AACCAAAA'])
        self.assertEqual(result['rejection_counts']['mapping_quality_unavailable_255'], 1)


if __name__ == '__main__':
    unittest.main()
