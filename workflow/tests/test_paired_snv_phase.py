import copy
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from audit_paired_snv_phase import merge_mate_calls


class PairedSnpPhaseTests(unittest.TestCase):
    def test_partial_mates_combine_but_missing_bases_are_not_reference(self):
        pair = [{'mate_number': 1, 'start': 10, 'mate_start': 90, 'template_length': 150, 'calls': {0: '1'}},
                {'mate_number': 2, 'start': 90, 'mate_start': 10, 'template_length': -150, 'calls': {1: '1'}}]
        self.assertEqual(merge_mate_calls(pair, 2), ('1|1', None))
        self.assertEqual(merge_mate_calls(pair, 3), (None, 'incomplete_target_coverage'))
        pair[1]['calls'][0] = '0'
        self.assertEqual(merge_mate_calls(pair, 2), ('mate_conflict', None))
        pair[1]['calls'][0] = '1'
        pair[1]['calls'][1] = '?'
        self.assertEqual(merge_mate_calls(pair, 2), ('other', None))

    def test_unrelated_or_duplicate_mates_cannot_create_phase(self):
        pair = [{'mate_number': 1, 'start': 10, 'mate_start': 90, 'template_length': 150, 'calls': {0: '1'}},
                {'mate_number': 2, 'start': 90, 'mate_start': 10, 'template_length': -150, 'calls': {1: '0'}}]
        self.assertEqual(merge_mate_calls(pair, 2), ('1|0', None))
        self.assertEqual(merge_mate_calls(pair[:1], 2), (None, 'missing_or_ambiguous_mate'))
        bad = copy.deepcopy(pair);bad[1]['mate_start'] = 11
        self.assertEqual(merge_mate_calls(bad, 2), (None, 'inconsistent_mate_geometry'))
        bad = copy.deepcopy(pair);bad[1]['template_length'] = -900
        self.assertEqual(merge_mate_calls(bad, 2), (None, 'inconsistent_mate_geometry'))
        self.assertEqual(merge_mate_calls(pair + pair[:1], 2), (None, 'missing_or_ambiguous_mate'))


if __name__ == '__main__':
    unittest.main()
