import gzip
from pathlib import Path
import tempfile
import unittest
from verify_completion import compare_streams, records

HEADER = '##fileformat=VCFv4.2\n##contig=<ID=1,length=100>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tfictional\n'
A = b'1\t10\t.\tA\tC\t.\t.\t.\tGT\t0/1\n'
B = b'1\t20\t.\tG\tT\t.\t.\t.\tGT\t1/1\n'


class Controls(unittest.TestCase):
    def test_full_payload_and_duplicate_multiplicity(self):
        result = compare_streams([[A, A], [B]], [A, A, B], [('1', 100)])
        self.assertEqual(result['records'], 3)
        self.assertEqual(result['records_by_contig'], {'1': 3})

    def test_loss_or_extra_record_refused(self):
        for final in ([A], [A, B, B]):
            with self.assertRaisesRegex(ValueError, 'count differs'):
                compare_streams([[A, B]], final, [('1', 100)])

    def test_field_change_refused(self):
        with self.assertRaisesRegex(ValueError, 'payload'):
            compare_streams([[A]], [A.replace(b'0/1', b'1/1')], [('1', 100)])

    def test_order_and_coordinate_refused(self):
        for values, dictionary in (([B, A], [('1', 100)]), ([A], [('1', 5)]), ([A], [('2', 100)])):
            with self.assertRaises(ValueError):
                compare_streams([values], values, dictionary)

    def test_gzip_header_and_literal_payload(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory) / 'fictional.vcf.gz'
            with gzip.open(p, 'wb') as stream:
                stream.write(HEADER.encode() + A)
            self.assertEqual(list(records(p, 'fictional', [('1', 100)])), [A])
            for sample, dictionary in (('other', [('1', 100)]), ('fictional', [('1', 99)])):
                with self.assertRaises(ValueError):
                    list(records(p, sample, dictionary))

    def test_unterminated_or_missing_header_refused(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory) / 'fictional.vcf.gz'
            for content in (HEADER.encode() + A.rstrip(b'\n'), A):
                with gzip.open(p, 'wb') as stream:
                    stream.write(content)
                with self.assertRaises(ValueError):
                    list(records(p, 'fictional', [('1', 100)]))


if __name__ == '__main__':
    unittest.main()
