import sys
from pathlib import Path
import unittest

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from review_snv_cluster_mapping import snv_states, fragment_configuration, location


class ClusterControls(unittest.TestCase):
    def test_clips_deletion_and_low_quality_are_not_reference(self):
        read = pysam.AlignedSegment()
        read.query_sequence = 'AACCGG'
        read.query_qualities = [35,35,35,5,35,35]
        read.reference_start = 100
        read.cigartuples = [(4,1),(0,2),(2,1),(0,3)]
        self.assertEqual(snv_states(read,[(100,'A','T'),(101,'A','T'),(102,'C','T'),
                                          (103,'A','T'),(104,'C','T'),(105,'G','T')],20),
                         [None,'0','0',None,None,'0'])

    def row(self,identifier,states,group='g',name='q',mq=60,flag=0):
        return dict(id=identifier,states=states,RG=group,query_name=name,mapping_quality=mq,flag=flag)

    def test_template_collapse_and_read_group_identity(self):
        rows=[self.row('a',['1','1']),self.row('b',['1','1']),
              self.row('c',['1','1'],group='h'),self.row('d',['0','0'],name='r')]
        counts=fragment_configuration(rows,[0,1],20)['counts']
        self.assertEqual((counts['1|1'],counts['0|0'],counts['mate_conflict']),(2,1,0))

    def test_partial_mates_do_not_create_a_spanning_observation(self):
        rows=[self.row('a',['1',None]),self.row('b',[None,'1']),
              self.row('c',['0','0'],name='bad_mq',mq=255),
              self.row('d',['1','1'],name='dup',flag=1024)]
        self.assertEqual(sum(fragment_configuration(rows,[0,1],20)['counts'].values()),0)

    def test_conflicting_observations_and_other_base(self):
        rows=[self.row('a',['1','1']),self.row('b',['0','1']),
              self.row('c',['?','0'],name='other')]
        counts=fragment_configuration(rows,[0,1],20)['counts']
        self.assertEqual(counts['mate_conflict'],1)
        self.assertEqual(counts['other'],1)
        self.assertEqual(counts['1|1'],0)

    def test_half_open_annotation_boundaries_and_unmapped(self):
        loci={'G':dict(chrom='1',start_0based=100,end_exclusive=200)}
        row=dict(flag=0,chrom='1',start_0based=200,end_exclusive=210)
        self.assertEqual(location(row,loci),'outside_annotated_loci')
        row['start_0based']=199
        self.assertEqual(location(row,loci),'G')
        row['flag']=4
        self.assertEqual(location(row,loci),'unmapped')


if __name__ == '__main__':
    unittest.main()
