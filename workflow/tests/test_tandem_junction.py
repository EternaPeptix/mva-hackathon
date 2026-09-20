from pathlib import Path
import random
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from review_tandem_junction import decompose_consensus, exact_high_quality, inspect_reads, query_span, reverse_complement, transcript_geometry


class TandemJunctionControls(unittest.TestCase):
    def test_native_BAM_evidence_counts_fragments_and_preserves_excluded_reads(self):
        try:import pysam
        except ImportError:self.skipTest('Native pysam environment required')
        needle='ACCGATATTCGAGTCGCTAGTTCAGGATACCTAGGCATGA'
        header={'HD':{'VN':'1.6'},'SQ':[{'SN':'1','LN':5000}],'RG':[{'ID':'lane','SM':'fictional'}]}
        config={'excluded_flags':1796,'minimum_mapping_quality':20,'minimum_base_quality':20,
                'split_coordinate_tolerance':3,'primary_read_anchor_length':20,'paired_endpoint_distance':1000}
        with tempfile.TemporaryDirectory() as tmp:
            bam=Path(tmp)/'fixture.bam'
            with pysam.AlignmentFile(str(bam),'wb',header=header) as out:
                for name,flag,mq in [('one_fragment',65,60),('one_fragment',129,60),('duplicate',1089,60),('unknown_MQ',65,255),('low_MQ',65,10)]:
                    r=pysam.AlignedSegment(out.header);r.query_name=name;r.flag=flag;r.reference_id=0;r.reference_start=100
                    r.mapping_quality=mq;r.cigarstring=str(len(needle))+'M';r.query_sequence=needle;r.query_qualities=[35]*len(needle);r.set_tag('RG','lane');out.write(r)
            result=inspect_reads(bam,'1',1001,2000,{'junction_20':needle},config)
            self.assertEqual(result['summary']['extracted_alignment_records'],5)
            self.assertEqual(result['summary']['evidence']['junction_20']['fragments'],1)
            self.assertEqual(result['summary']['excluded_record_reasons'],{'excluded_flag':1,'mapping_quality':2})

    def test_known_tandem_junction_and_reverse_complement_resolve_same_coordinates(self):
        rng=random.Random(7319);ref=''.join(rng.choice('ACGT') for _ in range(400))
        bases=list(ref)
        bases[79]=next(c for c in 'ACGT' if c!=ref[239])
        bases[240]=next(c for c in 'ACGT' if c!=ref[80])
        ref=''.join(bases)
        windows={'left':(50,ref[50:150]),'right':(180,ref[180:280])}
        junction=ref[200:240]+ref[80:120]
        for seq,direction in [(junction,'native'),(reverse_complement(junction),'reverse_complement')]:
            self.assertEqual(decompose_consensus(seq,windows,20),[{'orientation':direction,'order':['right','left'],'split_offset':40,
                'prefix_last_reference_base_1based':240,'suffix_first_reference_base_1based':81}])
        mutant=junction[:10]+next(c for c in 'ACGT' if c!=junction[10])+junction[11:]
        self.assertEqual(decompose_consensus(mutant,windows,20),[])

    def test_microhomology_retains_all_equivalent_solutions(self):
        rng=random.Random(7319);ref=''.join(rng.choice('ACGT') for _ in range(400))
        junction=ref[200:240]+ref[80:120]
        matches=decompose_consensus(junction,{'left':(50,ref[50:150]),'right':(180,ref[180:280])},20)
        self.assertEqual([(m['prefix_last_reference_base_1based'],m['suffix_first_reference_base_1based']) for m in matches],
                         [(237,78),(238,79),(239,80),(240,81)])

    def test_reference_sequence_low_quality_and_missing_quality_do_not_support_junction(self):
        needle='ACCGATATTCGAGTCGCTAG'+'TTCAGGATACCTAGGCATGA'
        read='GGA'+needle+'TCT';qualities=[35]*len(read)
        self.assertTrue(exact_high_quality(read,qualities,needle,20))
        self.assertTrue(exact_high_quality(reverse_complement(read),qualities[::-1],needle,20))
        self.assertFalse(exact_high_quality('ACCGATATTCGAGTCGCTAG'*3,[35]*60,needle,20))
        qualities[20]=10;self.assertFalse(exact_high_quality(read,qualities,needle,20))
        self.assertFalse(exact_high_quality(read,None,needle,20))

    def test_hard_and_soft_clips_share_original_query_coordinates(self):
        self.assertEqual(query_span([(0,50),(4,100)]),(0,50,150))
        self.assertEqual(query_span([(5,50),(0,100)]),(50,150,150))
        self.assertEqual(query_span([(4,50),(0,100)]),(50,150,150))
        self.assertEqual(query_span([(5,10),(4,15),(0,30),(1,5),(0,50),(4,40)]),(25,110,150))

    def test_extra_partial_copy_is_distinct_from_insertion_inside_original_transcript(self):
        exons=[{'start':a,'end':a+9,'attributes':{'rank':str(i)}} for i,a in enumerate([100,200,300],1)]
        transcript={'start':100,'end':310,'strand':'+'}
        result=transcript_geometry(transcript,exons,exons,195,500)
        self.assertTrue(result['original_transcript_span_uninterrupted_in_simple_tandem_model'])
        self.assertFalse(result['whole_transcript_span_copied_under_literal_tandem_model'])
        self.assertEqual(result['complete_exon_ranks_inside_literal_copy_interval'],[2,3])
        self.assertEqual(result['copied_CDS_bp_under_literal_tandem_model'],20)
        intragenic=transcript_geometry(transcript,exons,exons,195,305)
        self.assertFalse(intragenic['original_transcript_span_uninterrupted_in_simple_tandem_model'])
        self.assertEqual(intragenic['partial_exon_ranks_inside_literal_copy_interval'],[3])
        for exon in exons:exon['attributes']['rank']=str(4-int(exon['attributes']['rank']))
        negative=transcript_geometry(dict(transcript,strand='-'),exons,exons,195,500)
        self.assertEqual(negative['complete_exon_ranks_inside_literal_copy_interval'],[1,2])


if __name__=='__main__':unittest.main()
