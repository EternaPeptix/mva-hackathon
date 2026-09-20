import hashlib
from pathlib import Path
import sys
import tempfile
import unittest

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from run_mitochondrial_discovery import verify_rotation, original_position, reference_partition, load_depth, primary_identity


class MitochondrialMethodTests(unittest.TestCase):
    def test_rotation_is_exact_and_wrong_shift_fails(self):
        bases = 'A' * 8000 + 'C' * 8569
        self.assertTrue(verify_rotation(bases, bases[8000:] + bases[:8000]))
        with self.assertRaises(ValueError):
            verify_rotation(bases, bases[7999:] + bases[:7999])

    def test_partition_covers_every_reference_base_exactly_once(self):
        a, b = reference_partition()
        self.assertEqual(len(a), 15449)
        self.assertEqual(len(b), 1120)
        self.assertFalse(a & b)
        self.assertEqual(a | b, set(range(1, 16570)))
        self.assertEqual([original_position(p) for p in [1, 8569, 8570, 16569]], [8001, 16569, 1, 8000])

    def test_bad_shifted_positions_rejected(self):
        for p in [0, -1, 16570]:
            with self.assertRaises(ValueError):
                original_position(p)

    def test_coverage_requires_explicit_zero_and_unique_positions(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'depth.tsv'
            p.write_text('chrM\t1\t0\nchrM\t2\t8\n')
            self.assertEqual(load_depth(p, 2), {1: 0, 2: 8})
            for text in ['chrM\t2\t8\n', 'chrM\t1\t0\nchrM\t1\t8\n', 'MT\t1\t0\nMT\t2\t8\n']:
                p.write_text(text)
                with self.assertRaises(ValueError):
                    load_depth(p, 2)

    def write_pair(self, path, orphan=False, repeated=False):
        with pysam.AlignmentFile(str(path), 'wb', header={'HD': {'VN': '1.6'}, 'SQ': [{'SN': 'MT', 'LN': 100}], 'RG': [{'ID': 'a', 'SM': 'fictional'}]}) as out:
            for i in ([0] if orphan else [0, 1]) + ([0] if repeated else []):
                read = pysam.AlignedSegment(out.header)
                read.query_name, read.query_sequence = 'pair', ('ACGTA' if i == 0 else 'TACGT')
                read.query_qualities = ([31, 32, 33, 34, 35] if i == 0 else [35, 34, 33, 32, 31])
                read.flag, read.reference_id, read.reference_start = (65 if i == 0 else 145), 0, 20
                read.cigarstring = '5M'
                read.set_tag('RG', 'a')
                out.write(read)

    def test_forward_identity_retains_reverse_qualities(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / 'pair.bam'
            self.write_pair(path)
            rows, counts = primary_identity(path)
            self.assertEqual(counts['complete_pairs'], 1)
            self.assertEqual(len(set(rows.values())), 1)
            self.assertEqual(next(iter(rows.values())), hashlib.sha256(b'ACGTA\x00' + bytes([31,32,33,34,35])).hexdigest())

    def test_orphan_and_duplicate_primary_reads_fail(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / 'pair.bam'
            for kwargs in [{'orphan': True}, {'repeated': True}]:
                self.write_pair(path, **kwargs)
                with self.assertRaises(ValueError):
                    primary_identity(path)

    def test_reverted_unmapped_bam_has_no_reference_dictionary(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / 'unmapped.bam'
            with pysam.AlignmentFile(str(path), 'wb', header={'HD': {'VN': '1.6'}, 'RG': [{'ID':'a', 'SM':'fictional'}]}) as out:
                for flag in (77, 141):
                    read = pysam.AlignedSegment(out.header)
                    read.query_name, read.query_sequence, read.query_qualities = 'pair', 'ACGTA', [31,32,33,34,35]
                    read.flag = flag
                    read.set_tag('RG','a')
                    out.write(read)
            rows, counts = primary_identity(path)
            self.assertEqual(counts['primary_records'], 2)
            self.assertEqual(len(set(rows.values())), 1)


if __name__ == '__main__':
    unittest.main()
