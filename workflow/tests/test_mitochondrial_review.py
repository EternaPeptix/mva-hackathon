from pathlib import Path
import os
import sys
import tempfile
import unittest

import pysam

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from review_mitochondrial_calls import normalized_records


class MitochondrialReviewTests(unittest.TestCase):
    def test_native_alias_and_multiallelic_identity_preserve_original_context(self):
        with tempfile.TemporaryDirectory() as d:
            r=Path(d);ref=r/'ref.fa';ref.write_text('>chrM\nACGTACGTACGT\n');pysam.faidx(str(ref))
            source=r/'input.vcf'
            source.write_text('##fileformat=VCFv4.2\n##contig=<ID=M,length=12>\n'
                '##INFO=<ID=AS_FilterStatus,Number=A,Type=String,Description="Allele filters">\n'
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Fraction">\n'
                '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tfictional\n'
                'M\t2\t.\tC\tA,G\t.\tPASS\tAS_FilterStatus=weak_evidence,PASS\tGT:AF\t1/2:0.01,0.90\n')
            gz=Path(pysam.tabix_index(str(source),preset='vcf',force=True))
            result=normalized_records(gz,ref,Path(os.environ['MVA_BCFTOOLS']),r/'out')
            self.assertEqual([x['key'] for x in result['normalized_identities']],['chrM:2:C:A','chrM:2:C:G'])
            self.assertEqual([x['original_ALT_index'] for x in result['normalized_identities']],[1,2])
            self.assertEqual(result['original_records'][0]['chrom'],'M')
            self.assertEqual(result['original_records'][0]['samples']['fictional']['GT'],(1,2))
            self.assertEqual(result['original_records'][0]['info']['AS_FilterStatus'],('weak_evidence','PASS'))

    def test_wrong_reference_is_not_silently_relabelled(self):
        with tempfile.TemporaryDirectory() as d:
            r=Path(d);ref=r/'ref.fa';ref.write_text('>chrM\nACGT\n');pysam.faidx(str(ref))
            source=r/'input.vcf';source.write_text('##fileformat=VCFv4.2\n##contig=<ID=MT,length=4>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nMT\t2\t.\tA\tG\t.\tPASS\t.\n')
            gz=Path(pysam.tabix_index(str(source),preset='vcf',force=True))
            with self.assertRaises(ValueError):
                normalized_records(gz,ref,Path(os.environ['MVA_BCFTOOLS']),r/'out')


if __name__ == '__main__':
    unittest.main()
