"""Guard interpretation errors in synthetic structural-event comparisons."""
import copy
import json
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).parents[1] / 'scripts'))
from probe_structural_variants import compare_truth
from probe_sv_insertion import local_haplotype_matches


def event(svtype='BND'):
    value = {'id': 'event', 'type': svtype, 'chrom': '1', 'left_cut_0based': 1000,
             'right_cut_0based': 2000, 'length_bp': 1000}
    if svtype == 'BND':
        value.update(partner_chrom='2', partner_cut_0based=2000,
                     expected_connection_types=['3to5'])
    elif svtype == 'INV':
        value['expected_connection_types'] = ['3to3', '5to5']
    elif svtype == 'INS':
        value.update(right_cut_0based=1000, length_bp=100)
    return value


def call(svtype='BND'):
    return {'id': 'call', 'type': svtype, 'chrom': '1', 'pos_1based': 1000,
            'end_chrom': '2' if svtype == 'BND' else '1', 'end_position': 2000,
            'connection_type': '3to5', 'length': 100 if svtype == 'INS' else 1000,
            'genotype': {'GT': (0, 1)}, 'filter': ['PASS']}


def compare(value, calls, tolerance=150):
    return compare_truth({'events_on_one_haplotype': [value],
                          'unaltered_control_chromosome': '3'}, calls, tolerance)


class StructuralComparison(unittest.TestCase):
    def test_reversed_breakend_order_requires_reversed_orientation(self):
        original = call()
        reversed_call = copy.deepcopy(original)
        reversed_call.update(chrom='2', pos_1based=2000, end_chrom='1', end_position=1000,
                             connection_type='5to3')
        retained = copy.deepcopy(reversed_call)
        a = compare(event(), [original])['events'][0]
        b = compare(event(), [reversed_call])['events'][0]
        self.assertEqual(a['matched_endpoint_views'], b['matched_endpoint_views'])
        self.assertTrue(b['required_connection_types_detected'])
        self.assertEqual(reversed_call, retained)
        reversed_call['connection_type'] = '3to5'
        wrong = compare(event(), [reversed_call])['events'][0]
        self.assertTrue(wrong['site_detected_within_tolerance'])
        self.assertFalse(wrong['required_connection_types_detected'])

    def test_wrong_partner_chromosome_and_wrong_partner_position_do_not_match(self):
        for changes in ({'end_chrom': '4'}, {'end_position': 5000}):
            with self.subTest(changes=changes):
                wrong = call(); wrong.update(changes)
                result = compare(event(), [wrong])
                self.assertFalse(result['events'][0]['site_detected_within_tolerance'])
                self.assertEqual(result['unmatched_calls'], [wrong])

    def test_single_inversion_connection_is_incomplete(self):
        a = call('INV'); a['connection_type'] = '3to3'
        result = compare(event('INV'), [a])['events'][0]
        self.assertTrue(result['site_detected_within_tolerance'])
        self.assertFalse(result['required_connection_types_detected'])
        b = copy.deepcopy(a); b.update(id='second', connection_type='5to5')
        self.assertTrue(compare(event('INV'), [a, b])['events'][0]['required_connection_types_detected'])

    def test_insertion_needs_length_evidence_as_well_as_position(self):
        for length in (None, 1, 89, 111, 1000):
            with self.subTest(length=length):
                row = call('INS'); row['length'] = length
                self.assertFalse(compare(event('INS'), [row])['events'][0]['site_detected_within_tolerance'])
        for length in (90, 100, 110):
            row = call('INS'); row['length'] = length
            self.assertTrue(compare(event('INS'), [row])['events'][0]['site_detected_within_tolerance'])
        row['pos_1based'] = 1151
        self.assertFalse(compare(event('INS'), [row])['events'][0]['site_detected_within_tolerance'])

    def test_empty_callset_retains_every_missing_event(self):
        values = [event(kind) for kind in ('DEL', 'DUP', 'INV', 'INS', 'BND')]
        result = compare_truth({'events_on_one_haplotype': values, 'unaltered_control_chromosome': '3'}, [])
        self.assertEqual(len(result['events']), 5)
        for row in result['events']:
            self.assertFalse(row['site_detected_within_tolerance'])
            self.assertFalse(row['all_matched_PASS'])
            self.assertFalse(row['all_matched_genotypes_heterozygous'])

    def test_json_genotypes_match_native_tuples_but_missing_and_homozygous_do_not(self):
        native = call()
        loaded = json.loads(json.dumps(native))
        for row in (native, loaded):
            self.assertTrue(compare(event(), [row])['events'][0]['all_matched_genotypes_heterozygous'])
        for gt in (None, (None, None), (0, 0), [1, 1], [0], [0, 1, 1]):
            row = copy.deepcopy(native); row['genotype']['GT'] = gt
            self.assertFalse(compare(event(), [row])['events'][0]['all_matched_genotypes_heterozygous'])

    def test_filtered_or_control_chromosome_call_remains_visible(self):
        matched = call(); matched['filter'] = ['LowQual']
        control = copy.deepcopy(matched); control.update(id='control', chrom='3', end_chrom='3')
        result = compare(event(), [matched, control])
        self.assertTrue(result['events'][0]['site_detected_within_tolerance'])
        self.assertFalse(result['events'][0]['all_matched_PASS'])
        self.assertEqual(result['unmatched_calls'], [control])
        self.assertEqual(result['control_chromosome_calls'], [control])

    def test_equivalent_insertion_shift_is_sequence_checked(self):
        reference = 'GGGACACACTTT'
        # AC inserted at either repeat boundary gives the same resulting haplotype.
        self.assertTrue(local_haplotype_matches(reference, 0, 9, 'AC', 9, 'C', 'CAC'))
        self.assertTrue(local_haplotype_matches(reference, 0, 9, 'AC', 3, 'G', 'GAC'))
        self.assertFalse(local_haplotype_matches(reference, 0, 9, 'AC', 3, 'G', 'GAG'))
        self.assertFalse(local_haplotype_matches(reference, 0, 9, 'AC', 3, 'A', 'AAC'))
        self.assertFalse(local_haplotype_matches(reference, 0, 9, 'AC', 3, 'G', '<INS>'))
        self.assertFalse(local_haplotype_matches(reference, 0, 9, 'AC', 0, 'G', 'GAC'))


if __name__ == '__main__':
    unittest.main()
