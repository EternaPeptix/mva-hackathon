"""Controls for phase counting and conservative reported-alignment comparisons."""
import sys
from pathlib import Path
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import pysam
from review_paralog_candidate_pair import (compare_reported_mappings, haplotypes,
    merge_intervals, phase_observation, summarize_phase)
from summarize_paralog_mates import summarize_mates


def row(name, label='1|1', mq=40, flag=0, reverse=False, rg='lane1'):
    return {'id': name + str(flag) + str(reverse) + rg, 'query_name':name, 'RG':rg,
            'phase_label':label, 'mapping_quality':mq, 'flag':flag, 'reverse':reverse,
            'tags':{'NM':4}, 'cigar':'2M1D2M3D2M'}


class ParalogPairTests(unittest.TestCase):
    def test_fragment_dedup_conflicts_and_separate_read_groups(self):
        rows = [row('cis'), row('cis',reverse=True), row('conflict'), row('conflict','0|0',reverse=True),
                row('cis',rg='lane2'), row('dup',flag=1024)]
        result = summarize_phase(rows,20)
        self.assertEqual(result['counts']['1|1'],2)
        self.assertEqual(result['counts']['mate_conflict'],1)
        self.assertEqual(result['details']['1|1']['strand'],{'both':1,'forward':1})
        self.assertEqual(result['excluded_record_reasons_nonexclusive'],{'excluded_flag':1})

    def test_mapping_threshold_and_unknown_quality(self):
        rows = [row('x',mq=19),row('y',mq=20),row('z',mq=60),row('unknown',mq=255),row('not_span',label=None)]
        self.assertEqual(summarize_phase(rows,20)['counts']['1|1'],2)
        self.assertEqual(summarize_phase(rows,60)['counts']['1|1'],1)

    def test_exact_phase_and_base_quality_controls(self):
        expected=haplotypes('CAGCTAAGTT',20,[(22,'AG','A'),(25,'TAAG','T')])
        read=pysam.AlignedSegment()
        read.query_name='cis';read.reference_id=0;read.reference_start=20
        read.query_sequence='CACTTT';read.cigarstring='2M1D2M3D2M';read.query_qualities=[40]*6
        self.assertEqual(phase_observation(read,20,30,expected,20),'1|1')
        read.query_qualities=[40,40,19,40,40,40]
        self.assertIsNone(phase_observation(read,20,30,expected,20))
        read.query_qualities=[40]*6;read.cigarstring='2M1N2M3D2M'
        self.assertIsNone(phase_observation(read,20,30,expected,20))

    def test_reported_tie_missing_alternative_and_supplementary(self):
        locus={'chrom':'1','start_0based':100,'end_exclusive':500}
        def a(start,score,flag):
            return {'chrom':'1','start_0based':start,'end_exclusive':start+100,
                    'flag':flag,'mapping_quality':0,'tags':{'AS':score}}
        one=a(100,90,0)
        result=compare_reported_mappings([one,a(900,90,256),a(1500,110,2048)],locus)
        self.assertEqual(result['gene_minus_other_AS'],0)
        self.assertEqual(len(result['best_reported_other_locations']),1)
        self.assertIsNone(compare_reported_mappings([one],locus)['gene_minus_other_AS'])
        with self.assertRaises(ValueError):compare_reported_mappings([one,a(900,90,0)],locus)

    def test_regions_merge_without_duplicate_overlaps(self):
        self.assertEqual(merge_intervals([[100,200],[90,110],[200,250],[400,500]]),[[90,250],[400,500]])
        with self.assertRaises(ValueError):merge_intervals([[10,10]])

    def test_mate_support_requires_reciprocity_and_keeps_clipping(self):
        locus={'chrom':'1','start_0based':100,'end_exclusive':500}
        observation={'id':'r','RG':'L','query_name':'q','chrom':'1','start_0based':120,
                     'mate_chrom':'1','mate_start_0based':300,'mate_id':'m'}
        mate={'id':'m','RG':'L','query_name':'q','read_number':2,'chrom':'1','start_0based':300,
              'mate_chrom':'1','mate_start_0based':120,'mapping_quality':60,'flag':147,'tags':{'NM':2}}
        review={'observations':[observation],'observed_mates':[mate],
                'phase_thresholds':{'20':{'fragments':[{'RG':'L','query_name':'q','label':'1|1','observation_ids':['r']}]}}}
        remap={'comparisons':{'m':{'gene_minus_other_AS':None}}}
        primary={'m':{'chrom':'1','start_0based':300,'end_exclusive':350,'MQ':60,'flag':16,
                      'query_length':150,'aligned_query_length':50}}
        result=summarize_mates(review,remap,primary,locus)['groups']['1|1']['thresholds']['20']
        self.assertEqual(result['fragments_with_any_eligible_mate_primary_overlapping_gene'],1)
        self.assertEqual(result['fragments_with_gene_mate_and_at_least_90_percent_available_query_aligned'],0)
        mate['mate_start_0based']=121
        result=summarize_mates(review,remap,primary,locus)['groups']['1|1']
        self.assertEqual(result['fragments_with_nonreciprocal_mate'],1)
        self.assertEqual(result['thresholds']['20']['fragments_with_any_eligible_mate_primary_overlapping_gene'],0)


if __name__ == '__main__':
    unittest.main()
