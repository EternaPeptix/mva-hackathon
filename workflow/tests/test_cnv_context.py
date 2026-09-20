import base64
from email.message import Message
import hashlib
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from fetch_cnv_context_resources import head_metadata, verify_publisher_md5
from review_cnv_context import (IntervalIndex, allele_fraction, marker_summary, population_relation,
                               biallelic_frequency, read_clingen, annotate_dosage, summarize_population)
from summarize_cnv_context import copy_state_distribution


class ContextControls(unittest.TestCase):
    def test_nested_interval_and_inclusive_boundaries(self):
        rows = [dict(start_1based=1, end_inclusive=1000, id='long'),
                dict(start_1based=100, end_inclusive=150, id='short'),
                dict(start_1based=300, end_inclusive=300, id='point')]
        index = IntervalIndex(rows)
        self.assertEqual([r['id'] for r in index.fetch(300, 300)], ['long', 'point'])
        self.assertEqual([r['id'] for r in index.fetch(1000, 1001)], ['long'])
        self.assertEqual(index.fetch(1001, 1002), [])
        self.assertEqual(index.fetch(5, 4), [])

    def record(self, **changes):
        return dict(ID='public1', affected_start_1based=101, affected_end_inclusive=200, SVTYPE='DEL',
                    ALT=['<DEL>'], FILTER=['PASS'], INFO={'AN': 100, 'AC': (2,), 'AF': (0.02,)}) | changes

    def call(self, **changes):
        return dict(start_1based=101, end_inclusive=200, size_bp=100, type='deletion', bin_size=10) | changes

    def test_nominal_interval_endpoints_and_tolerance_not_equivalence(self):
        self.assertEqual(population_relation(self.call(), self.record())['minimum_reciprocal_overlap'], 1)
        with self.assertRaises(ValueError):
            population_relation(self.call(start_1based=1, end_inclusive=100), self.record())
        relation = population_relation(self.call(start_1based=111, end_inclusive=210), self.record())
        self.assertEqual(relation['minimum_reciprocal_overlap'], .9)
        self.assertTrue(relation['both_endpoints_within_one_bin'])
        self.assertNotIn('equivalent', relation)

    def test_multiallelic_copy_states_are_not_biallelic_AF(self):
        record = self.record(SVTYPE='CNV', ALT=['<CN0>', '<CN2>', '<CN3>'],
                             INFO={'CN_NUMBER': 100, 'CN_COUNT': (0, 2, 90, 8), 'CN_FREQ': (0, .02, .9, .08)})
        self.assertFalse(population_relation(self.call(), record)['same_simple_type'])
        self.assertIsNone(biallelic_frequency(record))
        self.assertIsNone(biallelic_frequency(self.record(INFO={'MULTIALLELIC': True})))

    def test_missing_frequency_is_unknown_and_variant_denominator_used(self):
        self.assertIsNone(biallelic_frequency(self.record(INFO={'AN': 0, 'AC': (0,), 'AF': (None,)})))
        self.assertIsNone(biallelic_frequency(self.record(INFO={})))
        self.assertEqual(biallelic_frequency(self.record())['AC_divided_by_AN'], .02)
        with self.assertRaises(ValueError):
            biallelic_frequency(self.record(INFO={'AN': 10, 'AC': (2,), 'AF': (.02,)}))

    def test_nonPASS_overlap_not_promoted_and_empty_frequency_not_zero(self):
        record = self.record(FILTER=['HIGH_NCR'])
        relation = population_relation(self.call(), record)
        record['biallelic_frequency'] = biallelic_frequency(record)
        result = summarize_population([relation], {'public1': record},
            {'population_reciprocal_overlap_thresholds': [.5, .8], 'maximum_population_examples_per_call': 5})
        self.assertEqual(result['any_nominal_span_overlap'], 1)
        self.assertEqual(result['same_simple_type_PASS'], 0)
        self.assertIsNone(result['same_simple_type_PASS_reciprocal_overlap']['0.5']['maximum_population_record_AF'])

    def test_excluded_high_depth_and_missing_markers_remain(self):
        base = dict(pos=100, AD=(40, 80), DP=120, PL=(20, 0, 30), GQ=20, category='heterozygous',
                    usable=False, excluded=['depth_outside_declared_range'])
        missing = dict(pos=200, AD=None, DP=None, PL=None, GQ=None, category='missing', usable=False,
                       excluded=['no_usable_depth'])
        result = marker_summary([base, missing], 10)
        self.assertEqual(result['markers'], 2)
        self.assertEqual(result['native_diploid_genotypes']['missing'], 1)
        self.assertEqual(result['usable_for_previous_ROH_model'], 0)
        self.assertAlmostEqual(result['native_heterozygous_AD_at_least_10']['ALT_fraction']['median'], 2/3)
        self.assertIsNone(result['native_heterozygous_ROH_usable_AD_at_least_10']['ALT_fraction']['median'])
        self.assertIsNone(allele_fraction((0, 0)))
        self.assertIsNone(allele_fraction((None, 4)))
        self.assertEqual(allele_fraction((0, 4)), 1)

    def test_dosage_raw_scores_unknown_locations_and_partial_overlap_retained(self):
        header = '#ISCA ID\tISCA Region Name\tGenomic Location\tHaploinsufficiency Score\tHaploinsufficiency Description\tTriplosensitivity Score\tTriplosensitivity Description\n'
        text = '#GRCh38\n' + header + 'one\tregion\tchr1:1-1000\t40\tunlikely\tNot yet evaluated\tunknown\n' + 'two\tpending\ttbd\t\t\t0\tno evidence\n'
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'regions.tsv'; p.write_text(text)
            rows, unmapped = read_clingen(p, 'region')
        self.assertEqual(len(unmapped), 1)
        annotated = annotate_dosage(self.call(type='duplication'), rows)[0]
        self.assertEqual(annotated['relevant_score'], 'Not yet evaluated')
        self.assertEqual(annotated['curation']['Haploinsufficiency Score'], '40')
        self.assertEqual(annotated['fraction_of_curated_span'], .1)

    def test_repeated_gcs_hash_headers_and_publisher_checksum_mismatch(self):
        body = b'public source fixture'
        message = Message(); message['x-goog-hash'] = 'crc32c=example'
        message['x-goog-hash'] = 'md5=' + base64.b64encode(hashlib.md5(body).digest()).decode()
        metadata = head_metadata(message)
        self.assertIn('crc32c=example', metadata['x-goog-hash'])
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'fixture'; p.write_bytes(body)
            self.assertEqual(verify_publisher_md5(p, metadata), hashlib.md5(body).hexdigest())
            p.write_bytes(body + b'corruption')
            with self.assertRaises(ValueError):
                verify_publisher_md5(p, metadata)

    def test_copy_state_denominator_is_samples_and_malformed_counts_fail(self):
        record = self.record(SVTYPE='CNV', INFO={'CN_NUMBER': 100, 'CN_COUNT': (0, 2, 90, 8), 'CN_FREQ': (0, .02, .9, .08)})
        self.assertEqual(copy_state_distribution(record)['count_divided_by_sample_denominator'][3], .08)
        record['INFO']['CN_COUNT'] = (0, 2, 90, 9)
        with self.assertRaises(ValueError):
            copy_state_distribution(record)


if __name__ == '__main__':
    unittest.main()
