"""Check aggregate counting and reference mismatches with a tiny synthetic genome."""
import importlib.util
from pathlib import Path
import tempfile
import unittest


class VendorAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if importlib.util.find_spec('pysam') is None:
            raise unittest.SkipTest('pysam is installed in the Studio runtime')
        import pysam
        cls.pysam = pysam
        spec = importlib.util.spec_from_file_location('vendor_audit', Path(__file__).resolve().parents[1] / 'scripts/audit_vendor_vcf.py')
        cls.audit = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.audit)

    def test_counts_phase_and_reference_mismatch(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            reference = root/'synthetic.fa'
            reference.write_text('>15\nAAAAAAAAAA\n')
            self.pysam.faidx(str(reference))
            vcf = root/'synthetic.vcf'
            vcf.write_text('''##fileformat=VCFv4.2
##contig=<ID=15,length=10>
##FILTER=<ID=q10,Description="Synthetic failed record">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PGT,Number=1,Type=String,Description="Physical genotype">
##FORMAT=<ID=PID,Number=1,Type=String,Description="Physical phase block">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC
15\t1\t.\tA\tT\t50\tPASS\t.\tGT:PGT:PID\t0/1:0|1:block
15\t2\t.\tA\tG\t50\tPASS\t.\tGT:PGT:PID\t0/1:0|1:block
15\t3\t.\tC\tT\t50\tPASS\t.\tGT\t0/1
15\t4\t.\tA\tT\t5\tq10\t.\tGT\t0/1
''')
            result = self.audit.inspect(vcf, reference, sample_size=10)
            self.assertEqual(result['record_counts']['records'], 4)
            self.assertEqual(result['record_counts']['pass_records'], 3)
            self.assertEqual(result['reference_check']['matched'], 2)
            self.assertEqual(result['reference_check']['mismatched'], 1)
            self.assertEqual(result['vendor_phase']['pass_phase_blocks'], 1)
            self.assertEqual(result['vendor_phase']['maximum_observed_variant_span_bp'], 1)


if __name__ == '__main__':
    unittest.main()
