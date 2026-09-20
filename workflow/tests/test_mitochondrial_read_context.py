from copy import deepcopy
from pathlib import Path
import sys
import tempfile
import unittest

import pysam

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from review_mitochondrial_read_context import base_context, select_fragments, project_anchor, identity


class ContextControls(unittest.TestCase):
    def read(self):
        read=pysam.AlignedSegment()
        read.query_name='fictional';read.flag=65;read.reference_start=100
        read.cigartuples=[(4,2),(0,3),(1,1),(0,2),(2,1),(0,2)]
        read.query_sequence='TTACGTACGT';read.query_qualities=[35]*10;read.set_tag('RG','synthetic')
        return read

    def test_cigar_geometry_and_other_edits(self):
        read=self.read()
        row=base_context(read,'ACGACAGT',100,104,[(104,'A','C'),(106,'A','C')])
        self.assertEqual(row['anchor']['base'],'A')
        self.assertEqual(row['anchor']['query_position_0based'],6)
        self.assertEqual(row['sites']['106'],None)
        self.assertEqual((row['other_mismatches'],row['insertion_deletion_bases'],row['other_edit_bases']),(0,2,2))
        self.assertEqual(row['soft_clipped_bases'],2)

    def test_reverse_cycle_and_deleted_anchor_are_explicit(self):
        read=self.read();read.flag|=16
        row=base_context(read,'ACGACAGT',100,104,[(104,'A','C')])
        self.assertEqual((row['anchor']['cycle_1based'],row['anchor']['distance_to_nearest_read_end']),(4,3))
        self.assertIsNone(base_context(read,'ACGACAGT',100,106,[(106,'A','C')])['anchor'])
        self.assertEqual(identity(read),('synthetic','fictional',1))

    def row(self, identifier, base='G', rg='g', name='q', reverse=False):
        return {'id':identifier,'RG':rg,'query_name':name,'flag':99,'mapping_quality':60,'reverse':reverse,
            'context':{'anchor':{'base':base,'raw_BQ':35,'distance_to_nearest_read_end':20},'other_edit_bases':0},
            'BAQ':{'100':{'base':base,'BAQ_BQ':35,'deletion':False,'refskip':False}}}

    def policy(self, **changes):
        return {'MQ':20,'BQ':20,'quality':'BAQ',**changes}

    def test_overlap_conflict_and_read_group_units(self):
        rows=[self.row('a'),self.row('b',reverse=True),self.row('c',base='T',name='conflict'),
              self.row('d',name='conflict'),self.row('e',rg='h')]
        result=select_fragments(rows,100,'T','G',self.policy())
        self.assertEqual((result['counts']['ALT'],result['counts']['mate_conflict']),(2,1))
        self.assertEqual(next(f for f in result['fragments'] if f['query_name']=='q' and f['RG']=='g')['strand'],'both')

    def test_absent_baq_low_quality_and_excluded_flags_are_not_ref(self):
        rows=[self.row(str(i),name=str(i)) for i in range(5)]
        rows[0]['BAQ']={};rows[1]['BAQ']['100']['BAQ_BQ']=5
        rows[2]['flag']|=1024;rows[3]['mapping_quality']=255
        result=select_fragments(rows,100,'T','G',self.policy())
        self.assertEqual(sum(result['counts'].values()),1)
        self.assertEqual(result['counts']['ALT'],1)
        self.assertEqual(select_fragments(rows[:1],100,'T','G',self.policy(quality='raw'))['counts']['ALT'],1)

    def test_context_filters_apply_to_ref_and_alt_with_inclusive_boundary(self):
        rows=[self.row('r','T',name='r'),self.row('a','G',name='a')]
        self.assertEqual(select_fragments(rows,100,'T','G',self.policy(trim=20))['counts']['ALT'],1)
        self.assertEqual(sum(select_fragments(rows,100,'T','G',self.policy(trim=21))['counts'].values()),0)
        for row in rows:row['context']['other_edit_bases']=1
        self.assertEqual(sum(select_fragments(rows,100,'T','G',self.policy(max_other_edits=0))['counts'].values()),0)
        counts=select_fragments(rows,100,'T','G',self.policy(max_other_edits=1))['counts']
        self.assertEqual((counts['REF'],counts['ALT']),(1,1))

    def test_deletion_does_not_become_reference(self):
        row=self.row('d');row['context']['anchor']=None
        row['BAQ']['100'].update(base=None,BAQ_BQ=None,deletion=True)
        counts=select_fragments([row],100,'T','G',self.policy())['counts']
        self.assertEqual((counts['deletion'],counts['REF']),(1,0))

    def test_projection_reverses_base_and_detects_clipping(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'ref.fa';path.write_text('>nuclear\nAACGTTAC\n');pysam.faidx(str(path))
            header=pysam.AlignmentHeader.from_references(['nuclear'],[8])
            read=pysam.AlignedSegment(header);read.query_name='q';read.reference_id=0;read.reference_start=0
            read.query_sequence='AACG';read.cigarstring='4M';read.flag=16
            old={'reverse':False,'context':{'read_length':4,'anchor':{'base':'G','query_position_0based':1}}}
            with pysam.FastaFile(str(path)) as fasta:
                projection=project_anchor(old,read,100,fasta)
                self.assertEqual((projection['status'],projection['position_1based'],projection['oriented_base']),('matched',3,'C'))
                read.cigarstring='3S1M'
                self.assertEqual(project_anchor(old,read,100,fasta)['status'],'anchor_clipped_or_inserted')
                read.cigarstring='1H3M'
                self.assertEqual(project_anchor(old,read,100,fasta)['status'],'hard_clip_not_projected')


if __name__=='__main__':unittest.main()
