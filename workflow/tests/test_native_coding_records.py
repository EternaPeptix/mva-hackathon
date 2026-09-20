"""Fictional native-call representation controls, not sensitivity measurements."""
from copy import deepcopy
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from native_coding_records import (allele_call, build_models, original_alt_index,
                                   review_axes, source_context)


class Sample(dict):
    phased = False


class Record:
    def __init__(self, pos=10, ref='A', alts=('C',), gt=(0, 1), fields=None,
                 phased=False, filters=(), rendered=None, info=None):
        self.contig, self.pos, self.ref, self.alts = '1', pos, ref, alts
        sample = Sample({'GT': gt, 'DP': 20, 'GQ': 35, 'AD': (10, 10), 'PL': (100, 0, 100)})
        if fields is not None:
            sample = Sample(fields)
        sample.phased = phased
        self.samples, self.filter, self.qual = {'SYNTH': sample}, filters, 100.0
        self.info = info or {}
        self.rendered = rendered or 'fictional original record at %d\n' % pos

    def __str__(self):
        return self.rendered


def context(record=None, caller='gatk', source='a', ordinal=0):
    return source_context(record or Record(), 'SYNTH', caller, source * 64, ordinal)


def variant(record=None, index=1, caller='gatk', source='a', ordinal=0,
            key=None, frequency=0.0, group=2, transcript='TR1', duplicate=False):
    call = allele_call(context(record, caller, source, ordinal), index)
    allele = deepcopy(call['source_allele'])
    return {'key': key or '%s:%d:%s:%s' % (allele['chrom'], allele['pos'], allele['ref'], allele['alt']),
            'allele': allele, 'source_call': call,
            'genes': {'GENE': {transcript: {'group': group, 'annotations': ['fictional']}}},
            'frequency': {'status': 'unresolved' if frequency is None else 'observed',
                          'max_AF_fraction': frequency, 'counts': []},
            'duplicate_normalized_key': duplicate, 'review_axes': review_axes(call, duplicate)}


class NativeCodingRecordTests(unittest.TestCase):
    def test_full_multiallelic_native_GT_PL_AD_and_extra_fields_survive(self):
        record = Record(alts=('C', 'G'), fields={'GT': (1, 2), 'DP': 22, 'GQ': None,
                       'AD': (1, 8, 13), 'PL': (91, 82, 73, 64, 0, 55), 'EXTRA': ('retained', None)})
        original = context(record)
        first, second = allele_call(original, 1), allele_call(original, 2)
        for call in (first, second):
            self.assertEqual(call['GT_values'], [1, 2])
            self.assertEqual(call['AD'], [1, 8, 13])
            self.assertEqual(call['PL'], [91, 82, 73, 64, 0, 55])
            self.assertEqual(call['source_all_alts'], ['C', 'G'])
            self.assertEqual(call['native_context']['FORMAT']['EXTRA'], ['retained', None])
            self.assertEqual(call['native_context']['rendered_record'], str(record))
        self.assertEqual((first['alt_depth'], second['alt_depth']), (8, 13))
        first['PL'][0] = -1
        self.assertEqual(original['FORMAT']['PL'][0], 91)
        self.assertEqual(second['PL'][0], 91)

    def test_original_ALT_mapping_checks_multiallelic_and_left_shift_tags(self):
        original = context(Record(alts=('C', 'G'), gt=(1, 2)))
        split = Record(alts=('G',), info={'NC_ORIG': '1|10|A|C,G|2'})
        self.assertEqual(original_alt_index(split, original), 2)
        split.info['NC_ORIG'] = ('1|10|A|C', 'G|2')
        self.assertEqual(original_alt_index(split, original), 2)
        shifted_context = context(Record(pos=12, ref='AA', alts=('A',)))
        shifted = Record(pos=8, ref='AA', alts=('A',), info={'NC_ORIG': '1|12|AA|A'})
        self.assertEqual(original_alt_index(shifted, shifted_context), 1)
        self.assertEqual(original_alt_index(Record(), context()), 1)
        for bad in ('1|10|A|G,C|2', '1|10|A|C,G|3', '1|9|A|C,G|2'):
            with self.subTest(tag=bad), self.assertRaises(ValueError):
                original_alt_index(Record(info={'NC_ORIG': bad}), original)
        with self.assertRaises(ValueError):
            original_alt_index(Record(), original)
        with self.assertRaises(ValueError):
            original_alt_index(Record(pos=8, ref='AA', alts=('A',)), shifted_context)
        for malformed in ((), ('1|10|A|C', None), ('1|10|A|C', 2), '', ('',)):
            with self.subTest(malformed=malformed), self.assertRaises(ValueError):
                original_alt_index(Record(info={'NC_ORIG': malformed}), original)

    def test_duplicate_native_records_have_distinct_identity_and_review_axes(self):
        first, second = context(ordinal=3), context(ordinal=4)
        self.assertNotEqual(first['source_record_id'], second['source_record_id'])
        self.assertEqual(first['source_namespace'], second['source_namespace'])
        self.assertNotEqual(first['source_namespace'], context(caller='bcftools')['source_namespace'])
        self.assertNotEqual(first['source_namespace'], context(source='b')['source_namespace'])
        rows = [variant(ordinal=3, duplicate=True), variant(ordinal=4, duplicate=True)]
        models, reviews, counts = build_models(rows, 0.001)
        self.assertEqual((models, reviews), ([], []))
        self.assertEqual(counts['structural_reason_duplicate_normalized_key'], 2)
        self.assertTrue(rows[0]['review_axes']['duplicate_normalized_key'])

    def test_missing_explicit_zero_and_below_threshold_are_distinct_and_retained(self):
        states = []
        for fields in ({'GT': (1, 1)}, {'GT': (1, 1), 'DP': None, 'GQ': None},
                       {'GT': (1, 1), 'DP': 0, 'GQ': 0}, {'GT': (1, 1), 'DP': 9, 'GQ': 19}):
            row = variant(Record(fields=fields))
            axes = row['review_axes']
            states.append((axes['DP']['state'], axes['DP']['FORMAT_presence']))
            self.assertEqual(axes['filter_state'], 'unfiltered')
            self.assertTrue(axes['structural_model_eligible'])
            self.assertEqual(len(build_models([row], 0.001)[0]), 1)
        self.assertEqual(states, [('missing', 'absent'), ('missing', 'present_missing'),
                                  ('explicit_zero', 'present_value'),
                                  ('below_historical_threshold', 'present_value')])
        named = variant(Record(gt=(1, 1), filters=('LowQual',)))
        self.assertEqual(named['review_axes']['filter_state'], 'named_filters')
        self.assertEqual(len(build_models([named], 0.001)[0]), 1)

    def test_partial_haploid_nonconcrete_and_uncalled_ALT_stay_out_of_models(self):
        rows = [variant(Record(gt=(None, 1)), ordinal=0),
                variant(Record(gt=(1,)), ordinal=1),
                variant(Record(alts=('C', '*'), gt=(1, 2)), ordinal=2),
                variant(Record(alts=('C', 'G'), gt=(0, 1)), index=2, ordinal=3)]
        for row in rows:
            self.assertFalse(row['review_axes']['structural_model_eligible'])
        self.assertEqual(rows[0]['review_axes']['genotype_state'], 'partial_or_missing')
        self.assertEqual(rows[1]['review_axes']['genotype_state'], 'haploid')
        self.assertTrue(rows[2]['review_axes']['nonconcrete_in_GT'])
        self.assertEqual(rows[3]['review_axes']['allele_participation'], 'ALT_absent_from_GT')
        self.assertEqual(build_models(rows, 0.001)[0], [])

    def test_concrete_1_2_retains_two_full_calls_and_source_opposite_relationship(self):
        record = Record(alts=('C', 'G'), gt=(1, 2))
        rows = [variant(record, index=1), variant(record, index=2, group=1)]
        models, reviews, _ = build_models(rows, 0.001)
        self.assertEqual(len(models), 1)
        self.assertEqual(reviews, [])
        self.assertEqual(models[0]['phase'], 'opposite_alleles_in_source_GT')
        self.assertEqual(models[0]['tier'], 2)
        self.assertEqual(len(models[0]['allele_evidence']), 2)
        self.assertFalse(models[0]['phase_evidence']['parental_phase_established'])
        self.assertEqual(models[0]['allele_evidence'][0]['source_call']['GT_values'], [1, 2])

    def test_callers_and_source_streams_do_not_share_phase(self):
        phased = {'GT': (0, 1), 'DP': 20, 'GQ': 30, 'PS': 7}
        first = variant(Record(pos=10, fields=phased, phased=True))
        other_caller = variant(Record(pos=20, fields=phased, phased=True), caller='bcftools', ordinal=1)
        self.assertEqual(build_models([first, other_caller], 0.001)[0], [])
        other_stream = variant(Record(pos=20, fields=phased, phased=True), source='b', ordinal=1)
        models = build_models([first, other_stream], 0.001)[0]
        self.assertEqual(len(models), 1)
        self.assertEqual(models[0]['phase'], 'unknown')

    def test_conflicting_PS_and_PID_remains_unknown_with_reason(self):
        first_fields = {'GT': (0, 1), 'PS': 7, 'PGT': '0|1', 'PID': 'BLOCK'}
        second_fields = {'GT': (0, 1), 'PS': 7, 'PGT': '1|0', 'PID': 'BLOCK'}
        rows = [variant(Record(pos=10, fields=first_fields, phased=True)),
                variant(Record(pos=20, fields=second_fields, phased=True), ordinal=1)]
        model = build_models(rows, 0.001)[0][0]
        self.assertEqual(model['phase'], 'unknown')
        self.assertIn('conflicting_PS_PID_pair_relations', model['phase_evidence']['issues'])

    def test_explicit_cis_and_overlap_pairs_retained_for_review(self):
        fields = {'GT': (0, 1), 'PS': 7}
        first = variant(Record(pos=10, fields=fields, phased=True))
        second = variant(Record(pos=20, fields=fields, phased=True), ordinal=1)
        models, reviews, _ = build_models([first, second], 0.001)
        self.assertEqual(models, [])
        self.assertEqual(reviews[0]['review_reasons'], ['explicit_source_cis'])
        first = variant(Record(pos=10, ref='AAA', alts=('A',)))
        second = variant(Record(pos=11), ordinal=1)
        models, reviews, _ = build_models([first, second], 0.001)
        self.assertEqual(models, [])
        self.assertEqual(reviews[0]['review_reasons'], ['overlapping_representations'])

    def test_frequency_common_unknown_and_observed_zero_remain_distinct(self):
        record = Record(gt=(1, 1))
        common = variant(record, frequency=0.02)
        self.assertEqual(build_models([common], 0.001)[0], [])
        observed = build_models([variant(record, frequency=0.0)], 0.001)[0][0]
        unknown = build_models([variant(record, frequency=None)], 0.001)[0][0]
        self.assertEqual(observed['frequency_group'], 'observed_rare')
        self.assertEqual(unknown['frequency_group'], 'frequency_unresolved')
        self.assertIsNone(unknown['allele_evidence'][0]['frequency']['max_AF_fraction'])
        invalid = variant(record, frequency=None)
        invalid['frequency']['max_AF_fraction'] = 0
        with self.assertRaises(ValueError):
            build_models([invalid], 0.001)

    def test_transcript_sharing_is_required_and_duplicate_rows_fail(self):
        first = variant(Record(pos=10), transcript='TR1')
        second = variant(Record(pos=20), ordinal=1, transcript='TR2')
        models, _, counts = build_models([first, second], 0.001)
        self.assertEqual(models, [])
        self.assertEqual(counts['heterozygous_pairs_without_shared_transcript'], 1)
        with self.assertRaises(ValueError):
            build_models([first, deepcopy(first)], 0.001)

    def test_native_pysam_rendering_and_fields(self):
        try:
            import pysam
        except ImportError:
            self.skipTest('native pysam control requires the installed Studio runtime')
        text = ('##fileformat=VCFv4.2\n##contig=<ID=1,length=100>\n'
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
                '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">\n'
                '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Likelihoods">\n'
                '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTH\n'
                '1\t10\t.\tA\tC,G\t100\t.\t.\tGT:DP:AD:PL\t1/2:22:1,8,13:91,82,73,64,0,55\n')
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'fictional.vcf'
            path.write_text(text)
            with pysam.VariantFile(path) as reader:
                record = next(reader)
                saved = source_context(record, 'SYNTH', 'gatk', 'a' * 64, 0)
                self.assertEqual(saved['rendered_record'], text.split('#CHROM')[1].split('\n', 1)[1])
                self.assertEqual(saved['FORMAT']['PL'], [91, 82, 73, 64, 0, 55])
                self.assertEqual(saved['FORMAT']['GT'], [1, 2])
                self.assertNotIn('GQ', saved['FORMAT'])
                self.assertEqual(allele_call(saved, 2)['alt_depth'], 13)
            self.assertEqual(path.read_text(), text)

    def test_native_pysam_comma_origin_in_VCF_and_BCF(self):
        try:
            import pysam
        except ImportError:
            self.skipTest('native pysam control requires the installed Studio runtime')
        text = ('##fileformat=VCFv4.2\n##contig=<ID=1,length=100>\n'
                '##INFO=<ID=NC_ORIG,Number=1,Type=String,Description="Original record">\n'
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTH\n'
                '1\t10\t.\tA\tG\t100\t.\tNC_ORIG=1|10|A|C,G|2\tGT\t0/1\n')
        original = context(Record(alts=('C', 'G'), gt=(1, 2)))
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'fictional-normalized.vcf'
            bcf = Path(folder) / 'fictional-normalized.bcf'
            path.write_text(text)
            with pysam.VariantFile(path) as reader:
                record = next(reader)
                self.assertEqual(record.info['NC_ORIG'], ('1|10|A|C', 'G|2'))
                self.assertEqual(original_alt_index(record, original), 2)
                with pysam.VariantFile(bcf, 'wb', header=reader.header) as writer:
                    writer.write(record)
            with pysam.VariantFile(bcf) as reader:
                record = next(reader)
                self.assertEqual(record.info['NC_ORIG'], ('1|10|A|C', 'G|2'))
                self.assertEqual(original_alt_index(record, original), 2)
            self.assertEqual(path.read_text(), text)


if __name__ == '__main__':
    unittest.main()
