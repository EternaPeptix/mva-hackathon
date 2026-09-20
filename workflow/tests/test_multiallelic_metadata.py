"""Ensure representation checks allow serialization rounding, not biological edits."""
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from prepare_public_multiallelic_sensitivity import record_fingerprint, verify_metadata


class MetadataChecks(unittest.TestCase):
    def compare(self, left, right):
        header = ('##fileformat=VCFv4.2\n##contig=<ID=chr20,length=1000>\n'
                  '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                  '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n')
        with tempfile.TemporaryDirectory() as folder:
            paths = [Path(folder) / name for name in ('raw.vcf', 'retained.vcf')]
            for path, rows in zip(paths, (left, right)):
                path.write_text(header + rows)
            return verify_metadata(*record_fingerprint(paths[0]), *record_fingerprint(paths[1]))

    def test_native_vcf_quality_rounding_is_reported(self):
        result = self.compare('chr20\t101\t.\tGAAAA\tGAA,G\t10289.64\t.\t.\tGT\t1/2\n',
                              'chr20\t101\t.\tGAAAA\tGAA,G\t10289.6\t.\t.\tGT\t1/2\n')
        self.assertEqual(result['parsed_QUAL_numeric_multiset_records_changed'], 1)
        self.assertTrue(result['FORMAT_serialized_QUAL_FILTER_multiset_preserved'])

    def test_changed_reported_quality_rejected(self):
        with self.assertRaisesRegex(ValueError, 'serialized QUAL'):
            self.compare('chr20\t101\t.\tG\tA\t10289.64\t.\t.\tGT\t0/1\n',
                         'chr20\t101\t.\tG\tA\t10290.6\t.\t.\tGT\t0/1\n')

    def test_changed_multiallelic_genotype_rejected(self):
        with self.assertRaises(ValueError):
            self.compare('chr20\t101\t.\tGAAAA\tGAA,G\t99\t.\t.\tGT\t1/2\n',
                         'chr20\t101\t.\tGAAAA\tGAA,G\t99\t.\t.\tGT\t0/2\n')

    def test_changed_phase_rejected(self):
        with self.assertRaises(ValueError):
            self.compare('chr20\t101\t.\tGAAAA\tGAA,G\t99\t.\t.\tGT\t1/2\n',
                         'chr20\t101\t.\tGAAAA\tGAA,G\t99\t.\t.\tGT\t1|2\n')

    def test_lost_record_rejected(self):
        row = 'chr20\t101\t.\tG\tA\t99\t.\t.\tGT\t0/1\n'
        with self.assertRaises(ValueError):
            self.compare(row + row, row)


if __name__ == '__main__':
    unittest.main()
