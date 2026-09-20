import io
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from survey_chromosome_depth import collect_depth, grid_windows, summarize


class ChromosomeDepthSurveyTests(unittest.TestCase):
    def test_grid_is_prespecified_nonoverlapping_and_bounded(self):
        rows = grid_windows('X', 10, 1010, 10, 50, 'nonPAR_X')
        self.assertEqual(len(rows), 10)
        self.assertEqual((rows[0]['start'], rows[-1]['end']), (35, 985))
        self.assertTrue(all(a['end'] <= b['start'] for a, b in zip(rows, rows[1:])))
        with self.assertRaises(ValueError):
            grid_windows('X', 10, 110, 10, 20, 'nonPAR_X')

    def test_explicit_zeros_reference_ambiguity_and_missing_positions(self):
        windows = [{'chrom': '1', 'start': 0, 'end': 3}]
        data = '1\t1\t0\n1\t2\t8\n1\t3\t20\n'
        rows, count = collect_depth(io.StringIO(data), windows, ['ANC'])
        self.assertEqual(count, 3)
        self.assertEqual(rows[0]['span_depth']['sum_depth'], 28)
        self.assertEqual(rows[0]['ACGT_depth']['histogram'], {0: 1, 20: 1})
        self.assertEqual(rows[0]['ACGT_depth']['bases_at_least']['20'], 1)
        for bad in ['1\t1\t0\n1\t3\t20\n', data + '1\t4\t0\n', data[:-7]]:
            with self.assertRaises(ValueError):
                collect_depth(io.StringIO(bad), windows, ['ANC'])
        with self.assertRaises(ValueError):
            collect_depth(io.StringIO(data), windows, [])

    def test_zero_depth_does_not_change_reference_eligibility(self):
        rows = []
        for category, mean, eligible in [('autosome', 0, True), ('autosome', 20, True),
                                          ('autosome', 100, False), ('nonPAR_X', 5, True)]:
            rows.append({'category': category, 'reference_eligible': eligible, 'ACGT_depth': {'mean': mean},
                         'span': 100, 'ACGT_bases': 100})
        result = summarize(rows)
        self.assertEqual(result['autosome']['median_reference_eligible_window_mean'], 10)
        self.assertEqual(result['nonPAR_X']['median_over_autosomal_median'], 0.5)
        self.assertEqual(result['autosome']['reference_eligible_windows'], 2)

    def test_native_flags_quality_overlap_deletions_and_empty_contig(self):
        samtools = os.environ['MVA_SAMTOOLS']
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            header = {'HD': {'VN': '1.6', 'SO': 'coordinate'},
                      'SQ': [{'SN': '1', 'LN': 60}, {'SN': '2', 'LN': 20}]}
            records = [('pair', 8, '12M', 12, 99, 60, 40, 12),
                       ('lowmap', 8, '10M', 10, 0, 0, 40, -1),
                       ('lowbase', 8, '10M', 10, 0, 60, 10, -1),
                       ('duplicate', 8, '10M', 10, 1024, 60, 40, -1),
                       ('supplementary', 8, '10M', 10, 2048, 60, 40, -1),
                       ('deletion', 8, '3M2D5M', 8, 0, 60, 40, -1),
                       ('skip', 8, '3M2N5M', 8, 0, 60, 40, -1),
                       ('pair', 12, '12M', 12, 147, 60, 40, 8)]
            bam = root / 'synthetic.bam'
            with pysam.AlignmentFile(str(bam), 'wb', header=header) as out:
                for name, start, cigar, length, flag, mq, bq, mate in records:
                    read = pysam.AlignedSegment(out.header)
                    read.query_name, read.query_sequence = name, 'A' * length
                    read.query_qualities = [bq] * length
                    read.reference_id, read.reference_start = 0, start
                    read.cigarstring, read.flag, read.mapping_quality = cigar, flag, mq
                    if mate >= 0:
                        read.next_reference_id, read.next_reference_start = 0, mate
                        read.template_length = 16 if start == 8 else -16
                    out.write(read)
            pysam.index(str(bam))
            bed = root / 'regions.bed'
            bed.write_text('1\t8\t18\n1\t25\t29\n2\t0\t4\n')
            for mq, bq, expected in [(20, 20, [3, 3, 3, 1, 1, 3, 3, 3, 3, 3]),
                                      (0, 1, [5, 5, 5, 3, 3, 5, 5, 5, 5, 5])]:
                command = [samtools, 'depth', '-aa', '-b', str(bed), '-q', str(bq), '-Q', str(mq), '-G', '0xF04', '-s', str(bam)]
                p = subprocess.run(command, text=True, capture_output=True)
                self.assertEqual(p.returncode, 0, p.stderr)
                rows = [line.split() for line in p.stdout.splitlines()]
                self.assertEqual(len(rows), 18)
                self.assertEqual([int(r[2]) for r in rows], expected + [0] * 8)
                self.assertEqual(rows[-1], ['2', '4', '0'])


if __name__ == '__main__':
    unittest.main()
