from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from call_ploidy_regions import validate_regions, ploidy_at


class PloidyRegionTests(unittest.TestCase):
    def test_half_open_transition_matches_both_sides(self):
        config = {'regions': [{'chrom': 'X', 'length': 100, 'intervals_0based_half_open': [[0, 20, 1], [20, 40, 2], [40, 100, 1]]}]}
        regions = validate_regions(config, {'X': 100})
        self.assertEqual([ploidy_at(regions, 'X', p)[0] for p in (0, 19, 20, 39, 40, 99)], [1, 1, 2, 2, 1, 1])

    def test_gaps_overlaps_invalid_ploidy_and_length_fail(self):
        for intervals in [[[0, 20, 1], [21, 100, 2]], [[0, 21, 1], [20, 100, 2]], [[0, 100, 0]], [[0, 99, 1]]]:
            with self.assertRaises(ValueError):
                validate_regions({'regions': [{'chrom': 'X', 'length': 100, 'intervals_0based_half_open': intervals}]}, {'X': 100})

    def test_outside_contig_is_not_defaulted_to_diploid(self):
        regions = [{'chrom': 'X', 'length': 100, 'intervals_0based_half_open': [[0, 100, 1]]}]
        for chrom, pos in [('Y', 1), ('X', -1), ('X', 100)]:
            with self.assertRaises(ValueError):
                ploidy_at(regions, chrom, pos)


if __name__ == '__main__':
    unittest.main()
