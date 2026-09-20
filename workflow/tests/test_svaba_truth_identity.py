"""Do not merge distinct truth alleles that share position, ID, type and length."""
import sys
import tempfile
import unittest
from pathlib import Path

import pysam

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from compare_svaba_public import allele_records, key


class TruthIdentityTests(unittest.TestCase):
    def read_rows(self, genotypes=((0,1),(1,0)), alleles=('AAC','AGT')):
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/'truth.vcf';h=pysam.VariantHeader();h.contigs.add('chr20',length=1000)
            h.formats.add('GT',1,'String','Genotype');h.add_sample('CONTROL')
            with pysam.VariantFile(str(path),'w',header=h) as out:
                for alt,gt in zip(alleles,genotypes):
                    r=out.new_record(contig='chr20',start=100,alleles=('A',alt));r.samples['CONTROL']['GT']=gt;out.write(r)
            return allele_records(path)

    def test_same_length_distinct_alternates_remain_separate(self):
        rows=self.read_rows()
        self.assertEqual((rows[0]['pos'],rows[0]['id'],rows[0]['type'],rows[0]['size']),
                         (rows[1]['pos'],rows[1]['id'],rows[1]['type'],rows[1]['size']))
        self.assertEqual(len({key(r) for r in rows}),2)

    def test_genotype_order_does_not_change_allele_identity(self):
        rows=self.read_rows(alleles=('AAC','AAC'))
        self.assertEqual(key(rows[0]),key(rows[1]))

    def test_repeat_readback_identity_is_stable(self):
        self.assertEqual([key(r) for r in self.read_rows()],[key(r) for r in self.read_rows()])


if __name__=='__main__':unittest.main()
