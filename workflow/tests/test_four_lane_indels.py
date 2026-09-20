from pathlib import Path
import sys
import tempfile
import unittest

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from audit_final_alignment_indels import window_evidence, call_context
from audit_baseline_reads import exact_calls


def fixture(root, rows):
    reference = root / 'reference.fa'
    reference.write_text('>1\n' + 'G' * 20 + 'A' * 20 + 'C' * 160 + '\n')
    pysam.faidx(str(reference))
    bam = root / 'reads.bam'
    header = {'HD': {'VN': '1.6', 'SO': 'coordinate'}, 'SQ': [{'SN': '1', 'LN': 200}],
              'RG': [{'ID': 'lane1', 'SM': 'test'}, {'ID': 'lane2', 'SM': 'test'}]}
    with pysam.AlignmentFile(str(bam), 'wb', header=header) as out:
        for name, sequence, cigar, flag, mq, low, lane in rows:
            read = pysam.AlignedSegment(out.header)
            read.query_name, read.query_sequence = name, sequence
            read.flag, read.reference_id, read.reference_start = flag, 0, 20
            read.cigarstring, read.mapping_quality = cigar, mq
            qualities = [40] * len(sequence)
            if low is not None:
                qualities[low] = 10
            read.query_qualities = qualities
            read.set_tag('RG', lane)
            out.write(read)
    pysam.index(str(bam))
    target = {'key': '1:31:A:AA', 'chrom': '1', 'pos': 31, 'ref': 'A', 'alt': 'AA',
              'other_hypotheses': {'ALT2': 'AAA'}}
    return bam, reference, target


class FourLaneIndelTests(unittest.TestCase):
    def test_shifted_repeat_insertions_and_distinct_alternatives(self):
        rows = [('ref', 'A' * 20 + 'C' * 40, '60M', 0, 60, None, 'lane1'),
                ('alt', 'A' * 21 + 'C' * 40, '11M1I49M', 0, 60, None, 'lane1'),
                ('shifted', 'A' * 21 + 'C' * 40, '20M1I40M', 16, 60, None, 'lane1'),
                ('different', 'A' * 22 + 'C' * 40, '11M2I49M', 0, 60, None, 'lane2')]
        with tempfile.TemporaryDirectory() as directory:
            bam, reference, target = fixture(Path(directory), rows)
            result = window_evidence(bam, reference, target, 10, 10)
            self.assertEqual(result['counts'], {'REF': 1, 'ALT': 2, 'ALT2': 1, 'ambiguous': 0, 'other': 0, 'mate_conflict': 0})
            self.assertEqual(result['strand_counts']['ALT'], {'forward': 1, 'reverse': 1})

    def test_fragment_conflicts_read_groups_and_missing_quality(self):
        ref, alt = 'A' * 20 + 'C' * 40, 'A' * 21 + 'C' * 40
        rows = [('conflict', ref, '60M', 0, 60, None, 'lane1'),
                ('conflict', alt, '11M1I49M', 16, 60, None, 'lane1'),
                ('same', alt, '11M1I49M', 0, 60, None, 'lane1'),
                ('same', alt, '11M1I49M', 16, 60, None, 'lane1'),
                ('same', alt, '11M1I49M', 0, 60, None, 'lane2'),
                ('unknownmq', alt, '11M1I49M', 0, 255, None, 'lane1'),
                ('lowmq', alt, '11M1I49M', 0, 0, None, 'lane1'),
                ('lowbq', alt, '11M1I49M', 0, 60, 11, 'lane1'),
                ('dup', alt, '11M1I49M', 1024, 60, None, 'lane1')]
        with tempfile.TemporaryDirectory() as directory:
            bam, reference, target = fixture(Path(directory), rows)
            result = window_evidence(bam, reference, target, 1, 2)
            self.assertEqual(result['eligible_fragments'], 3)
            self.assertEqual(result['counts']['ALT'], 2)
            self.assertEqual(result['counts']['mate_conflict'], 1)
            self.assertEqual(result['excluded_read_counts'], {'low_or_unknown_mapping_quality': 2,
                'not_spanning_or_low_quality_or_reference_skip': 1, 'excluded_flags': 1})

    def test_longer_window_nearby_change_and_missing_boundaries_not_reference(self):
        sequence = list('A' * 21 + 'C' * 40)
        sequence[26] = 'G'
        rows = [('neighbor', ''.join(sequence), '11M1I49M', 0, 60, None, 'lane1'),
                ('deleted_boundary', 'A' * 19 + 'C' * 40, '9M1D50M', 0, 60, None, 'lane1'),
                ('skip', 'A' * 19 + 'C' * 40, '9M1N50M', 0, 60, None, 'lane1')]
        with tempfile.TemporaryDirectory() as directory:
            bam, reference, target = fixture(Path(directory), rows)
            short = window_evidence(bam, reference, target, 1, 2)
            long = window_evidence(bam, reference, target, 10, 25)
            self.assertEqual(short['counts']['ALT'], 1)
            self.assertEqual(short['counts']['REF'], 0)
            self.assertEqual(short['excluded_read_counts']['not_spanning_or_low_quality_or_reference_skip'], 2)
            self.assertEqual(long['counts']['ALT'], 0)
            self.assertEqual(long['counts']['REF'], 0)
            self.assertEqual(long['counts']['other'], 2)

    def test_original_multiallelic_genotype_survives_identity_normalization(self):
        header = ('##fileformat=VCFv4.2\n##contig=<ID=1,length=200>\n'
                  '##INFO=<ID=SF_ORIG,Number=1,Type=String,Description="Original record">\n'
                  '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                  '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
                  '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depth">\n'
                  '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttest\n')
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw, normalized = root / 'raw.vcf', root / 'norm.vcf'
            raw.write_text(header + '1\t31\t.\tA\tAA,AAA\t50\tPASS\t.\tGT:DP:AD\t1/2:12:0,7,5\n')
            normalized.write_text(header + '1\t31\t.\tA\tAA\t50\tPASS\tSF_ORIG=1|31|A|AA,AAA|1\tGT:DP:AD\t1/0:12:0,7\n')
            paths = []
            for path in (raw, normalized):
                pysam.tabix_compress(str(path), str(path) + '.gz')
                pysam.tabix_index(str(path) + '.gz', preset='vcf')
                paths.append(Path(str(path) + '.gz'))
            target = {'key': '1:31:A:AA', 'chrom': '1', 'pos': 31, 'ref': 'A', 'alt': 'AA'}
            call = exact_calls(paths[0], paths[1], [target], 'test')[0]['exact_normalized_calls'][0]
            self.assertEqual(call['GT_values'], [1, 2])
            self.assertEqual(call['alt_copies'], 1)
            self.assertEqual(call['source_all_alts'], ['AA', 'AAA'])
            self.assertEqual(call['AD'], [0, 7, 5])
            self.assertEqual(call_context(paths[0], [target], 'test', 0)[target['key']][0]['GT'], (1, 2))


if __name__ == '__main__':
    unittest.main()
