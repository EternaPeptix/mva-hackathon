from pathlib import Path
import sys
import tempfile
import unittest

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import audit_baseline_reads as audit


def write_alignment(root, rows):
    reference = root / 'ref.fa'
    reference.write_text('>1\n' + 'A' * 200 + '\n')
    pysam.faidx(str(reference))
    header = {'HD': {'VN': '1.6', 'SO': 'coordinate'}, 'SQ': [{'SN': '1', 'LN': 200}],
              'RG': [{'ID': 'lane1'}, {'ID': 'lane2'}]}
    bam = root / 'test.bam'
    with pysam.AlignmentFile(str(bam), 'wb', header=header) as stream:
        for name, sequence, cigar, flag, mq, low_position, group in rows:
            read = pysam.AlignedSegment(stream.header)
            read.query_name, read.query_sequence = name, sequence
            read.flag, read.reference_id, read.reference_start = flag, 0, 20
            read.mapping_quality, read.cigarstring = mq, cigar
            qualities = [40] * len(sequence)
            if low_position is not None:
                qualities[low_position] = 10
            read.query_qualities = qualities
            read.set_tag('RG', group)
            stream.write(read)
    pysam.index(str(bam))
    return bam, reference


class BaselineReadTests(unittest.TestCase):
    def test_snv_cooccurrence_deduplicates_mates_preserves_read_groups_and_rejects_weak_reads(self):
        def sequence(first=False, second=False, extra=False):
            text = list('A' * 40)
            if first: text[11] = 'C'
            if second: text[17] = 'G'
            if extra: text[14] = 'T'
            return ''.join(text)
        rows = [('reference', sequence(), '40M', 0, 60, None, 'lane1'),
                ('reference', sequence(), '40M', 16, 60, None, 'lane1'),
                ('cis', sequence(True, True, True), '40M', 0, 60, None, 'lane1'),
                ('first', sequence(True), '40M', 0, 60, None, 'lane1'),
                ('second', sequence(False, True), '40M', 0, 60, None, 'lane1'),
                ('sameid', sequence(True, True), '40M', 0, 60, None, 'lane1'),
                ('sameid', sequence(True, True), '40M', 0, 60, None, 'lane2'),
                ('conflict', sequence(), '40M', 0, 60, None, 'lane1'),
                ('conflict', sequence(True, True), '40M', 16, 60, None, 'lane1'),
                ('lowmap', sequence(True, True), '40M', 0, 0, None, 'lane1'),
                ('lowbase', sequence(True, True), '40M', 0, 60, 11, 'lane1'),
                ('duplicate', sequence(True, True), '40M', 1024, 60, None, 'lane1'),
                ('deletion', 'A' * 39, '11M1D28M', 0, 60, None, 'lane1')]
        with tempfile.TemporaryDirectory() as temp:
            bam, reference = write_alignment(Path(temp), rows)
            result = audit.count_snv_configuration(bam, reference, [('1', 32, 'A', 'C'), ('1', 38, 'A', 'G')])
            self.assertEqual(result['counts'], {'0|0': 1, '0|1': 1, '1|0': 1, '1|1': 3, 'other': 0, 'mate_conflict': 1})
            self.assertEqual(result['strand_counts']['0|0']['both'], 1)
            self.assertEqual(result['excluded_read_counts']['flags_or_mapping_quality'], 2)
            self.assertEqual(result['excluded_read_counts']['low_target_base_quality'], 1)
            self.assertEqual(result['excluded_read_counts']['not_covering_all_sites_or_deletion'], 1)
            with self.assertRaisesRegex(ValueError, 'reference mismatch'):
                audit.count_snv_configuration(bam, reference, [('1', 32, 'T', 'C'), ('1', 38, 'A', 'G')])

    def test_exact_indel_window_preserves_other_sequences_and_conflicting_mates(self):
        reference_seq = 'A' * 40
        alt_seq = 'A' * 12 + 'TG' + 'A' * 28
        other_seq = 'A' * 10 + 'C' + 'A' * 29
        rows = [('ref', reference_seq, '40M', 0, 60, None, 'lane1'),
                ('alt', alt_seq, '12M2I28M', 0, 60, None, 'lane1'),
                ('alt', alt_seq, '12M2I28M', 16, 60, None, 'lane1'),
                ('different', other_seq, '40M', 0, 60, None, 'lane1'),
                ('conflict', reference_seq, '40M', 0, 60, None, 'lane1'),
                ('conflict', alt_seq, '12M2I28M', 16, 60, None, 'lane1'),
                ('duplicate', alt_seq, '12M2I28M', 1024, 60, None, 'lane1')]
        with tempfile.TemporaryDirectory() as temp:
            bam, reference = write_alignment(Path(temp), rows)
            result = audit.indel_window_evidence(bam, reference, ('1', 32, 'A', 'ATG'))
            self.assertEqual(result['counts'], {'REF': 1, 'ALT': 1, 'other': 1, 'mate_conflict': 1})
            self.assertEqual(result['strand_counts']['ALT']['both'], 1)


if __name__ == '__main__':
    unittest.main()
