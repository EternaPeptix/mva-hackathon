from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from review_four_lane_structural import usable_support


class StructuralContextControls(unittest.TestCase):
    def base(self):
        return dict(filter=['PASS'],genotype={'GT':[0,1],'FT':'PASS'},
            endpoint_reference_bases_in_query_intervals=[True,True],native_diploid_sex_genotype_uninterpreted=False)

    def test_site_PASS_is_not_nonreference_genotype_support(self):
        row=self.base(); self.assertTrue(usable_support(row))
        for gt in ([0,0],[1,None],None):
            row['genotype']['GT']=gt; self.assertFalse(usable_support(row))
        row=self.base();row['genotype']['FT']='LowQual';self.assertFalse(usable_support(row))

    def test_sex_query_boundaries_and_native_filters_retained_outside_subset(self):
        for key,value in [('filter',['LowQual']),('endpoint_reference_bases_in_query_intervals',[True,False]),
                          ('native_diploid_sex_genotype_uninterpreted',True)]:
            row=self.base();row[key]=value;self.assertFalse(usable_support(row))


if __name__=='__main__':unittest.main()
