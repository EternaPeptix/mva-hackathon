"""Literal half-open BED edge cases; native matching may still cross these edges."""
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from summarize_public_representation import overlaps, read_regions


class ScopeGeometry(unittest.TestCase):
    def test_exclusive_edges_and_inside_bases(self):
        regions, starts = [(100, 110), (200, 210)], [100, 200]
        for left, right in [(99, 100), (110, 111), (199, 200), (210, 211)]:
            self.assertFalse(overlaps(regions, starts, left, right))
        for left, right in [(100, 101), (109, 110), (200, 201), (209, 210)]:
            self.assertTrue(overlaps(regions, starts, left, right))

    def test_deletion_anchor_outside_but_REF_span_inside(self):
        self.assertFalse(overlaps([(100, 110)], [100], 99, 100))
        self.assertTrue(overlaps([(100, 110)], [100], 99, 102))

    def test_overlapping_BED_is_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'regions.bed'
            path.write_text('chr20\t100\t110\nchr20\t109\t120\n')
            with self.assertRaisesRegex(ValueError, 'nonoverlapping'):
                read_regions(path)


if __name__ == '__main__':
    unittest.main()
