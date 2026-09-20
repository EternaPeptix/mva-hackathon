from pathlib import Path
import random
import sys
import tempfile
import unittest

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from review_changed_call_reads import (call_context, fragment_evidence, merged_regions,
                                       require_known_mapping_quality, validate_selection)


def targets():
    return [{'key': '1:100:T:'+alt, 'chrom': '1', 'pos': 100, 'ref': 'T', 'alt': alt}
            for alt in ('A', 'C')]


def reference(root):
    rng = random.Random(520731)
    seq = list(''.join(rng.choice('ACGT') for _ in range(300)))
    seq[99] = 'T'
    seq = ''.join(seq)
    p = root/'reference.fa'
    p.write_text('>1\n'+seq+'\n')
    pysam.faidx(str(p))
    return p, seq


class ChangedCallReadTests(unittest.TestCase):
    def test_competing_alternates_share_one_extraction_interval(self):
        with tempfile.TemporaryDirectory() as tmp:
            fa, _ = reference(Path(tmp))
            rows = targets()
            validate_selection(rows, fa)
            self.assertEqual(merged_regions(rows, 50, fa), [['1', 49, 150]])
            with self.assertRaisesRegex(ValueError, 'duplicate'):
                validate_selection(rows+[rows[0]], fa)
            changed = dict(rows[0], ref='G', key='1:100:G:A')
            with self.assertRaisesRegex(ValueError, 'reference'):
                validate_selection([changed], fa)

    def test_haploid_context_preserves_uncalled_alternate_and_likelihoods(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'calls.vcf'
            path.write_text('##fileformat=VCFv4.2\n##contig=<ID=1,length=300>\n'
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">\n'
                '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Likelihoods">\n'
                '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTH\n'
                '1\t100\t.\tT\tA,C\t50\t.\t.\tGT:AD:PL\t2:0,3,9:40,50,0\n')
            gz = pysam.tabix_index(str(path), preset='vcf')
            result = call_context(gz, targets(), 'SYNTH', 0)
            for target in targets():
                row = result[target['key']][0]
                self.assertEqual(row['GT'], (2,))
                self.assertEqual(row['alts'], ['A', 'C'])
                self.assertEqual(row['AD'], (0, 3, 9))
                self.assertEqual(row['PL'], (40, 50, 0))
            empty = dict(targets()[0], pos=200, key='1:200:T:A')
            self.assertEqual(call_context(gz, [empty], 'SYNTH', 0)[empty['key']], [])

    def test_three_bases_mate_conflicts_lanes_and_unknown_mapping_quality(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fa, seq = reference(root)
            header = {'HD': {'VN': '1.6', 'SO': 'coordinate'},
                      'SQ': [{'SN': '1', 'LN': 300}],
                      'RG': [{'ID': x, 'SM': 'SYNTH'} for x in ('L1', 'L2')]}
            rows = [('same', 'A', 'L1', 0, 60, 40), ('same', 'A', 'L1', 16, 60, 40),
                    ('same', 'A', 'L2', 0, 60, 40), ('reference', 'T', 'L1', 0, 60, 40),
                    ('third', 'C', 'L1', 16, 60, 40), ('conflict', 'A', 'L1', 0, 60, 40),
                    ('conflict', 'C', 'L1', 16, 60, 40), ('duplicate', 'A', 'L1', 1024, 60, 40),
                    ('lowmap', 'A', 'L1', 0, 0, 40), ('lowbase', 'A', 'L1', 0, 60, 10),
                    ('orphan', 'A', 'L1', 1, 60, 40)]

            def write(path, specifications):
                with pysam.AlignmentFile(str(path), 'wb', header=header) as out:
                    for name, base, lane, flag, mq, bq in specifications:
                        read = pysam.AlignedSegment(out.header)
                        bases = list(seq[60:140]); bases[39] = base
                        read.query_name, read.query_sequence = name, ''.join(bases)
                        qualities = [40]*80; qualities[39] = bq
                        read.query_qualities = qualities
                        read.reference_id, read.reference_start = 0, 60
                        read.cigarstring, read.flag, read.mapping_quality = '80M', flag, mq
                        read.set_tag('RG', lane)
                        out.write(read)
                pysam.index(str(path))

            bam = root/'reads.bam'; write(bam, rows)
            require_known_mapping_quality(bam)
            filters = {'minimum_mapping_quality': 20, 'minimum_base_quality': 20,
                       'exclude_flags': '0xF04'}
            a, c = [fragment_evidence(bam, fa, t, filters) for t in targets()]
            self.assertEqual(a['counts'], {'REF': 1, 'ALT': 2, 'deletion': 0,
                                          'refskip': 0, 'other': 1, 'mate_conflict': 1})
            self.assertEqual(c['counts'], {'REF': 1, 'ALT': 1, 'deletion': 0,
                                          'refskip': 0, 'other': 2, 'mate_conflict': 1})
            self.assertEqual(a['eligible_fragments'], c['eligible_fragments'])
            unknown = root/'unknown.bam'
            write(unknown, rows+[('unknown', 'A', 'L1', 0, 255, 40)])
            with self.assertRaisesRegex(ValueError, 'MQ255'):
                require_known_mapping_quality(unknown)


if __name__ == '__main__':
    unittest.main()
