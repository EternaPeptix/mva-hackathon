import sys
from pathlib import Path
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from compare_cluster_public_sequence import parse_embl, complete_exact, plausible_pair, score_difference


class CloneControls(unittest.TestCase):
    def test_accession_version_length_and_sequence(self):
        text='ID   X123; SV 2; linear; genomic DNA; STD; HUM; 6 BP.\nSQ   Sequence 6 BP;\n     acg tnn 6\n//\n'
        self.assertEqual(parse_embl(text),('X123.2','ACGTNN'))
        with self.assertRaises(ValueError):parse_embl(text.replace('6 BP.','7 BP.'))
        with self.assertRaises(ValueError):parse_embl(text.replace('//',''))

    def row(self,start=100,flag=0,nm=0,aligned=149):
        return dict(flag=flag,chrom='clone',start_0based=start,end_exclusive=start+aligned,
                    query_alignment_length=aligned,query_length=149,tags={'NM':nm})

    def test_exact_requires_full_query_and_no_edits(self):
        self.assertTrue(complete_exact(self.row()))
        self.assertFalse(complete_exact(self.row(nm=1)))
        self.assertFalse(complete_exact(self.row(aligned=100)))
        self.assertFalse(complete_exact(self.row(flag=4)))

    def test_pair_orientation_span_and_clipping(self):
        a,b=self.row(),self.row(start=350,flag=16)
        self.assertTrue(plausible_pair(a,b,2000))
        self.assertTrue(plausible_pair(b,a,2000))
        self.assertFalse(plausible_pair(a,self.row(start=350),2000))
        self.assertFalse(plausible_pair(a,self.row(start=1,flag=16),2000))
        self.assertFalse(plausible_pair(a,self.row(start=3000,flag=16),2000))
        self.assertFalse(plausible_pair(a,self.row(start=350,flag=16,aligned=100),2000))

    def test_unmapped_score_tag_is_not_a_valid_comparison(self):
        a,b=self.row(),self.row()
        a['tags']['AS']=149;b['tags']['AS']=145
        self.assertEqual(score_difference(a,b),4)
        b['flag']=4;b['tags']['AS']=0
        self.assertIsNone(score_difference(a,b))
        b['flag']=0;del b['tags']['AS']
        self.assertIsNone(score_difference(a,b))


if __name__=='__main__':unittest.main()
