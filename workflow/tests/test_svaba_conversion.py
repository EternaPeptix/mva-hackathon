"""Independent haplotype and malformed-input checks for SvABA conversion."""
import copy
import sys
from pathlib import Path
import unittest
from types import SimpleNamespace
import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from convert_svaba_pilot import (decode_alt, linear_allele, paired_events,
                                reference_status, index_assembly_records, UnresolvedJunction)


class Reference:
    sequence = 'AACCTTGGTTCCAA'
    def fetch(self, chrom, start, end):
        if chrom != '1': raise ValueError('unknown chromosome')
        return self.sequence[start:end]
    def get_reference_length(self, chrom): return len(self.sequence)


def event(kind, p, q, insertion='', homology=''):
    return {'type':kind,'connection_type':'3to5' if kind=='DEL' else '5to3',
            'a':{'chrom':'1','pos':p,'info':{'SCTG':'assembled','INSERTION':insertion,'HOMSEQ':homology}},
            'b':{'chrom':'1','pos':q}}


def resulting_haplotype(allele):
    seq=Reference.sequence;start=allele['pos']-1
    assert seq[start:start+len(allele['ref'])] == allele['ref']
    return seq[:start]+allele['alt']+seq[start+len(allele['ref']):]


def assembly(e, sequence, left_qend, right_qstart, reverse=False):
    """Independent literal contig and SAM alignments for the tiny test chromosome."""
    header=pysam.AlignmentHeader.from_references(['1'],[14]); reads=[]
    left_end=e['a']['pos'] if e['type']=='DEL' else e['b']['pos']
    right_start=(e['b']['pos'] if e['type']=='DEL' else e['a']['pos'])-1
    for start,cigar in [(left_end-left_qend,f'{left_qend}M{len(sequence)-left_qend}S'),
                        (right_start,f'{right_qstart}S{len(sequence)-right_qstart}M')]:
        r=pysam.AlignedSegment(header);r.query_name='assembled';r.reference_id=0;r.reference_start=start
        r.flag=16 if reverse else 0;r.query_sequence=sequence;r.cigarstring=cigar;reads.append(r)
    return index_assembly_records(reads)


class ConversionTests(unittest.TestCase):
    def test_left_virtual_boundary_has_no_reference_base(self):
        row=SimpleNamespace(contig='1',start=-1,pos=0,ref='N',info={'SVTYPE':'BND'})
        self.assertEqual(reference_status(row,Reference()),'virtual_boundary_no_reference_base')

    def test_right_virtual_boundary_has_no_reference_base(self):
        row=SimpleNamespace(contig='1',start=14,pos=15,ref='N',info={'SVTYPE':'BND'})
        self.assertEqual(reference_status(row,Reference()),'virtual_boundary_no_reference_base')

    def test_ordinary_position_zero_is_rejected(self):
        row=SimpleNamespace(contig='1',start=-1,pos=0,ref='N',info={'SVTYPE':'DEL'})
        with self.assertRaisesRegex(ValueError,'bounds'):reference_status(row,Reference())

    def test_interior_reference_mismatch_still_rejected(self):
        row=SimpleNamespace(contig='1',start=1,pos=2,ref='C',info={'SVTYPE':'BND'})
        with self.assertRaisesRegex(ValueError,'REF mismatch'):reference_status(row,Reference())

    def test_deletion_haplotype(self):
        e=event('DEL',4,9);a=assembly(e,'AACCTTCCAA',4,4)
        self.assertEqual(resulting_haplotype(linear_allele(e,Reference(),a)), 'AACCTTCCAA')

    def test_deletion_shared_overlap_counted_once(self):
        e=event('DEL',6,9,homology='TT');a=assembly(e,'AACCTTCCAA',6,4)
        self.assertEqual(resulting_haplotype(linear_allele(e,Reference(),a)), 'AACCTTCCAA')

    def test_duplication_haplotype(self):
        e=event('DUP',5,8);a=assembly(e,'AACCTTGGTTGGTTCCAA',8,8)
        self.assertEqual(resulting_haplotype(linear_allele(e,Reference(),a)), 'AACCTTGGTTGGTTCCAA')

    def test_duplication_shared_overlap_counted_once(self):
        e=event('DUP',5,10,homology='TT');a=assembly(e,'AACCTTGGTTGGTTCCAA',10,8)
        self.assertEqual(resulting_haplotype(linear_allele(e,Reference(),a)), 'AACCTTGGTTGGTTCCAA')

    def test_inserted_sequence_forward(self):
        e=event('DEL',4,5,'GAT');a=assembly(e,'AACCGATTTGGTTCCAA',4,7)
        self.assertEqual(resulting_haplotype(linear_allele(e,Reference(),a)), 'AACCGATTTGGTTCCAA')

    def test_inserted_sequence_reverse(self):
        e=event('DEL',4,5,'ATC');a=assembly(e,'AACCGATTTGGTTCCAA',4,7,True)
        self.assertEqual(resulting_haplotype(linear_allele(e,Reference(),a)), 'AACCGATTTGGTTCCAA')

    def test_deletion_with_insertion(self):
        e=event('DEL',4,9,'GAT');a=assembly(e,'AACCGATTTCCAA',4,7)
        self.assertEqual(resulting_haplotype(linear_allele(e,Reference(),a)), 'AACCGATTTCCAA')

    def test_unknown_insertion_orientation_rejected(self):
        with self.assertRaisesRegex(ValueError,'orientation'):linear_allele(event('DEL',4,5,'GAT'),Reference(),{})

    def test_deletion_overlap_mismatch_preserves_assembled_bases(self):
        e=event('DEL',6,8,homology='TT');a=assembly(e,'AACCTTTCCAA',6,4)
        allele=linear_allele(e,Reference(),a)
        self.assertEqual(resulting_haplotype(allele),'AACCTTTCCAA')
        self.assertFalse(allele['overlap_matches_both_reference_flanks'])

    def test_duplication_overlap_mismatch_preserves_assembled_bases(self):
        e=event('DUP',5,10,homology='AC');a=assembly(e,'AACCTTGGACGGTTCCAA',10,8)
        self.assertEqual(resulting_haplotype(linear_allele(e,Reference(),a)),'AACCTTGGACGGTTCCAA')

    def test_reverse_overlap_sequence(self):
        e=event('DUP',5,10,homology='GT');a=assembly(e,'AACCTTGGACGGTTCCAA',10,8,True)
        self.assertEqual(resulting_haplotype(linear_allele(e,Reference(),a)),'AACCTTGGACGGTTCCAA')

    def test_stale_homology_is_not_added_to_insertion(self):
        e=event('DEL',4,9,'GAT','CC');a=assembly(e,'AACCGATTTCCAA',4,7)
        allele=linear_allele(e,Reference(),a)
        self.assertEqual(resulting_haplotype(allele),'AACCGATTTCCAA')
        self.assertEqual(allele['assembly_junction']['inactive_annotation_fields'],['HOMSEQ'])

    def test_stale_insertion_is_not_added_to_overlap(self):
        e=event('DEL',6,9,'GAT','TT');a=assembly(e,'AACCTTCCAA',6,4)
        allele=linear_allele(e,Reference(),a)
        self.assertEqual(resulting_haplotype(allele),'AACCTTCCAA')
        self.assertEqual(allele['assembly_junction']['inactive_annotation_fields'],['INSERTION'])

    def test_annotation_must_match_assembled_sequence(self):
        e=event('DEL',4,9,'CAT');a=assembly(e,'AACCGATTTCCAA',4,7)
        with self.assertRaisesRegex(UnresolvedJunction,'sequence'):linear_allele(e,Reference(),a)

    def test_imprecise_coordinates_are_an_explicit_assumption(self):
        e=event('DEL',4,9);e['a']['info']['IMPRECISE']=True
        allele=linear_allele(e,Reference(),{})
        self.assertEqual(resulting_haplotype(allele),'AACCTTCCAA')
        self.assertEqual(allele['assembly_junction']['source'],'approximate_discordant_coordinates')

    def test_imprecise_sequence_is_not_invented(self):
        e=event('DEL',4,9,'GAT');e['a']['info']['IMPRECISE']=True
        with self.assertRaisesRegex(UnresolvedJunction,'unverified'):linear_allele(e,Reference(),{})

    def test_gapped_overlap_is_retained_as_unresolved(self):
        e=event('DEL',6,7,homology='CTTG');a=assembly(e,'AACCTTGCCAA',7,3)
        for v in a.values():
            for r in v:
                if r['ordinal']==0:r['cigar']=[(0,4),(1,1),(0,2),(4,4)];r['cigar_string']='4M1I2M4S'
        with self.assertRaisesRegex(UnresolvedJunction,'gapped'):linear_allele(e,Reference(),a)

    def test_a_single_alignment_cannot_pair_with_itself(self):
        e=event('DUP',5,10,homology='TTGGTT');a=assembly(e,'TTGGTT',6,0)
        for key in a:a[key]=[r for r in a[key] if r['ordinal']==0]
        with self.assertRaisesRegex(UnresolvedJunction,'sequence'):linear_allele(e,Reference(),a)

    def test_all_four_breakend_orientations(self):
        for alt,local,remote in [('C[1:9[','3','5'),('C]1:9]','3','3'),(']1:9]C','5','3'),('[1:9[C','5','5')]:
            d=decode_alt(alt,'C');self.assertEqual((d['local_side'],d['remote_side']),(local,remote))

    def pair(self):
        a={'id':'a','chrom':'1','pos':4,'ref':'C','alt':'C[1:9[','info':{'MATEID':'b'},'gt':[0,1],'filters':['PASS']}
        b={'id':'b','chrom':'1','pos':9,'ref':'T','alt':']1:4]T','info':{'MATEID':'a'},'gt':[0,1],'filters':['PASS']}
        return a,b

    def test_mates_count_once_independent_of_input_order(self):
        a,b=self.pair();self.assertEqual(paired_events([a,b],{'1':0}),paired_events([b,a],{'1':0}))
        self.assertEqual(len(paired_events([a,b],{'1':0})),1)

    def test_missing_mate_rejected(self):
        with self.assertRaisesRegex(ValueError,'mate'):paired_events([self.pair()[0]],{'1':0})

    def test_wrong_reciprocal_coordinate_rejected(self):
        a,b=self.pair();b['alt']=']1:3]T'
        with self.assertRaisesRegex(ValueError,'coordinate'):paired_events([a,b],{'1':0})

    def test_discordant_mate_genotype_rejected(self):
        a,b=self.pair();b['gt']=[1,1]
        with self.assertRaisesRegex(ValueError,'genotype'):paired_events([a,b],{'1':0})


if __name__=='__main__': unittest.main()
