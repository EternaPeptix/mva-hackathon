from pathlib import Path
import sys
import unittest

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from diagnose_indel_windows import window_state, target_cycle


class IndelWindowDiagnosticTests(unittest.TestCase):
    def read(self, cigar='30M', length=30, flag=0, mq=60):
        read = pysam.AlignedSegment()
        read.query_name = 'synthetic'
        read.query_sequence = 'A' * length
        read.query_qualities = [40] * length
        read.reference_id, read.reference_start = 0, 20
        read.flag, read.mapping_quality, read.cigarstring = flag, mq, cigar
        return read

    def test_quality_and_missing_boundaries_remain_distinct(self):
        read = self.read()
        self.assertEqual(window_state(read, 25, 40), 'eligible')
        qualities = [40] * 30
        qualities[10] = 10
        read.query_qualities = qualities
        self.assertEqual(window_state(read, 25, 40), 'low_window_quality')
        self.assertEqual(window_state(read, 10, 40), 'missing_left')
        self.assertEqual(window_state(read, 25, 60), 'missing_right')
        self.assertEqual(window_state(read, 10, 60), 'missing_both')
        self.assertEqual(window_state(self.read('5M1D25M'), 25, 40), 'missing_left')
        self.assertEqual(window_state(self.read('10M1N20M'), 25, 40), 'reference_skip')
        self.assertEqual(window_state(self.read(mq=255), 25, 40), 'low_or_unknown_mapping_quality')

    def test_original_orientation_cycle_with_softclip_and_deleted_base(self):
        forward, reverse = self.read('5S25M'), self.read('5S25M', flag=16)
        self.assertEqual(target_cycle(forward, 25), 11)
        self.assertEqual(target_cycle(reverse, 25), 20)
        self.assertIsNone(target_cycle(self.read('5M1D25M'), 25))


if __name__ == '__main__':
    unittest.main()
