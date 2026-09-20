import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from review_four_lane_copy_number import overlap, union, reasons


class ReviewControls(unittest.TestCase):
    def test_inclusive_boundary_and_union(self):
        self.assertEqual(overlap(1, 100, 100, 200), 1)
        self.assertEqual(overlap(1, 100, 101, 200), 0)
        self.assertEqual(union([(11,20),(1,10),(7,15),(30,30)]), [[1,20],[30,30]])

    def test_missing_q0_and_reference_overhang_never_pass(self):
        policy = dict(maximum_q0=0.1, maximum_gap_fraction=0.01, maximum_e1_exclusive=0.05,
                      maximum_deletion_ratio=0.75, minimum_duplication_ratio=1.25)
        row = dict(within_reference_bounds=True, MQ_zero_fraction=0.02, reference_gap_fraction=0,
                   e1=0, relative_depth=0.5, type='deletion')
        self.assertEqual(reasons(row, policy), [])
        self.assertEqual(reasons(dict(row, MQ_zero_fraction=-1), policy), ['q0_missing_or_invalid'])
        self.assertEqual(reasons(dict(row, within_reference_bounds=False), policy), ['outside_reference_bounds'])
        self.assertIn('e1_outside_limit', reasons(dict(row, e1={'nonfinite':'nan'}), policy))

    def test_neutral_depth_call_label_and_zero_statistic_do_not_override_effect(self):
        policy = dict(maximum_q0=0.1, maximum_gap_fraction=0.01, maximum_e1_exclusive=0.05,
                      maximum_deletion_ratio=0.75, minimum_duplication_ratio=1.25)
        row = dict(within_reference_bounds=True, MQ_zero_fraction=0, reference_gap_fraction=0,
                   e1=0, relative_depth=1.01, type='duplication')
        self.assertEqual(reasons(row, policy), ['relative_depth_effect_below_limit'])


if __name__ == '__main__':
    unittest.main()
