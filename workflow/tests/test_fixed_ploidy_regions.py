import json
from pathlib import Path
import sys
import tempfile
import unittest

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from call_fixed_ploidy_regions import merge_fixed_calls, record_contract


class FixedPloidyContractTests(unittest.TestCase):
    def setUp(self):
        self.header = pysam.VariantHeader()
        self.header.contigs.add('control', length=6000)
        for name, number, kind in [('GT', 1, 'String'), ('AD', 'R', 'Integer'), ('PL', 'G', 'Integer')]:
            self.header.formats.add(name, number, kind, name)
        self.header.add_sample('fictional')
        self.regions = [{'chrom': 'control', 'length': 6000,
            'intervals_0based_half_open': [[0, 3000, 2], [3000, 6000, 1]]}]

    def record(self, pos, alleles, gt, pl, ad):
        record = self.header.new_record(contig='control', start=pos - 1, alleles=alleles)
        record.samples['fictional']['GT'] = gt
        record.samples['fictional']['PL'] = pl
        record.samples['fictional']['AD'] = ad
        return record

    def test_original_multialleles_and_boundary_span_preserved(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            records = {
                1: [self.record(3000, ('A', 'C'), (1,), (90, 0), (0, 12)),
                    self.record(3001, ('A', 'C', 'G'), (2,), (90, 50, 0), (0, 0, 12))],
                2: [self.record(3000, ('AA', 'A'), (0, 1), (90, 0, 90), (6, 6)),
                    self.record(3001, ('A', 'C', 'G'), (2, 2), (90, 90, 90, 90, 90, 0), (0, 0, 12))]}
            paths = {}
            for p, values in records.items():
                paths[p] = root / f'p{p}.vcf'
                with pysam.VariantFile(str(paths[p]), 'w', header=self.header) as stream:
                    for row in values:
                        stream.write(row)
            evidence = merge_fixed_calls(paths, root / 'owned.vcf.gz', self.regions, 'fictional', root / 'excluded.jsonl')
            with pysam.VariantFile(str(root / 'owned.vcf.gz')) as stream:
                owned = list(stream)
            self.assertEqual([r.pos for r in owned], [3000, 3001])
            self.assertEqual(owned[0].samples['fictional']['GT'], (0, 1))
            self.assertEqual(owned[1].alleles, ('A', 'C', 'G'))
            self.assertEqual(owned[1].samples['fictional']['GT'], (2,))
            self.assertEqual(owned[1].samples['fictional']['PL'], (90, 50, 0))
            self.assertEqual(owned[1].samples['fictional']['AD'], (0, 0, 12))
            self.assertEqual(len(evidence['reference_span_crosses_ploidy_boundary']), 1)
            self.assertEqual(evidence['reference_span_crosses_ploidy_boundary'][0]['pos'], 3000)
            excluded = [json.loads(line) for line in (root / 'excluded.jsonl').read_text().splitlines()]
            self.assertEqual(len(excluded), 2)
            self.assertEqual({r['pos'] for r in excluded}, {3000, 3001})

    def test_diploid_shape_cannot_pass_as_haploid(self):
        row = self.record(3001, ('A', 'C'), (1, 1), (90, 40, 0), (0, 12))
        with self.assertRaisesRegex(ValueError, 'genotype contract'):
            record_contract(row, 'fictional', 1)

    def test_partial_genotype_is_not_reference_evidence(self):
        row = self.record(3000, ('A', 'C'), (0, None), (0, 40, 90), (12, 0))
        with self.assertRaisesRegex(ValueError, 'genotype contract'):
            record_contract(row, 'fictional', 2)

    def test_missing_likelihoods_fail_contract(self):
        row = self.header.new_record(contig='control', start=3000, alleles=('A', 'C'))
        row.samples['fictional']['GT'] = (1,)
        with self.assertRaisesRegex(ValueError, 'likelihood cardinality'):
            record_contract(row, 'fictional', 1)


if __name__ == '__main__':
    unittest.main()
