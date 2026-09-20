import sys
from pathlib import Path
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from project_cluster_target_bases import reoriented_query_position, complement


class ProjectionControls(unittest.TestCase):
    def test_all_four_strand_combinations_and_boundaries(self):
        for old in [False,True]:
            for new in [False,True]:
                expected=[0,10,148] if old==new else [148,138,0]
                self.assertEqual([reoriented_query_position(p,149,old,new) for p in [0,10,148]],expected)
        with self.assertRaises(ValueError):reoriented_query_position(149,149,False,False)

    def test_base_complement_involution(self):
        self.assertEqual(complement('ACGTN'),'TGCAN')
        self.assertEqual(complement(complement('ACGTN')),'ACGTN')


if __name__=='__main__':unittest.main()
