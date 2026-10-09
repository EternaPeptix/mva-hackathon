import copy
import hashlib
import json
import unittest

from verify_cpdat_bub1b_screen import decision_snapshot


def ledger_bytes(ledger, **kwargs):
    return (json.dumps(ledger, indent=2, ensure_ascii=False, **kwargs) + '\n').encode('utf-8')


class DecisionSnapshotTests(unittest.TestCase):
    def setUp(self):
        self.ledger = {
            'date': '2026-09-28',
            'status': 'Dated screen',
            'decisions': [
                {'decision_id': 'one', 'name': 'β compound', 'status': 'review'},
                {'decision_id': 'two', 'name': 'Other compound', 'status': 'reject'},
            ],
        }
        self.recorded_count = len(self.ledger['decisions'])
        self.recorded_sha256 = hashlib.sha256(ledger_bytes(self.ledger)).hexdigest()

    def snapshot(self, ledger):
        return decision_snapshot(ledger_bytes(ledger), self.recorded_count, self.recorded_sha256)

    def appended_ledger(self):
        ledger = copy.deepcopy(self.ledger)
        ledger['decisions'].append({'decision_id': 'later', 'status': 'review'})
        return ledger

    def test_exact_snapshot_and_later_appends_pass(self):
        self.assertEqual(self.snapshot(self.ledger), self.ledger['decisions'])
        ledger = self.appended_ledger()
        ledger['decisions'].append({'decision_id': 'still-later', 'status': 'reject'})
        self.assertEqual(self.snapshot(ledger), self.ledger['decisions'])

    def test_changed_historical_decision_fails_even_with_append(self):
        ledger = self.appended_ledger()
        ledger['decisions'][0]['status'] = 'reject'
        with self.assertRaisesRegex(AssertionError, 'snapshot changed'):
            self.snapshot(ledger)

    def test_reordered_historical_decisions_fail(self):
        ledger = self.appended_ledger()
        ledger['decisions'][:2] = reversed(ledger['decisions'][:2])
        with self.assertRaisesRegex(AssertionError, 'snapshot changed'):
            self.snapshot(ledger)

    def test_changed_ledger_metadata_fails_even_with_append(self):
        ledger = self.appended_ledger()
        ledger['date'] = '2026-09-29'
        with self.assertRaisesRegex(AssertionError, 'snapshot changed'):
            self.snapshot(ledger)

    def test_truncated_snapshot_fails(self):
        ledger = copy.deepcopy(self.ledger)
        ledger['decisions'].pop()
        with self.assertRaisesRegex(AssertionError, 'snapshot is truncated'):
            self.snapshot(ledger)

    def test_unreproducible_original_serialization_is_not_normalized(self):
        original_bytes = json.dumps(self.ledger, indent=4, ensure_ascii=False).encode('utf-8')
        recorded_hash = hashlib.sha256(original_bytes).hexdigest()
        self.assertEqual(decision_snapshot(original_bytes, self.recorded_count, recorded_hash), self.ledger['decisions'])
        with self.assertRaisesRegex(AssertionError, 'snapshot changed'):
            decision_snapshot(ledger_bytes(self.appended_ledger()), self.recorded_count, recorded_hash)

    def audited_refinement(self):
        ledger = self.appended_ledger()
        ledger['date'] = '2026-10-04'
        ledger['decisions'][0]['status'] = 'source-qualified revision'
        current_hash = hashlib.sha256(json.dumps(ledger['decisions'][0], ensure_ascii=False,
                                                sort_keys=True, separators=(',', ':')).encode()).hexdigest()
        audit = {
            'recorded_decision_count': self.recorded_count,
            'recorded_ledger_sha256': self.recorded_sha256,
            'metadata_refinements': {'date': {'original': self.ledger['date'], 'accepted_current': ledger['date']}},
            'decision_refinements': [{'index': 0, 'decision_id': 'one',
                                     'original_object': copy.deepcopy(self.ledger['decisions'][0]),
                                     'accepted_current_object_sha256': current_hash}],
        }
        return ledger, audit

    def test_exact_audited_refinement_and_later_append_pass(self):
        ledger, audit = self.audited_refinement()
        result = decision_snapshot(ledger_bytes(ledger), self.recorded_count, self.recorded_sha256, audit)
        self.assertEqual(result, self.ledger['decisions'])
        self.assertEqual(ledger['decisions'][0]['status'], 'source-qualified revision')
        self.assertEqual(self.snapshot(self.ledger), self.ledger['decisions'])

    def test_audited_refinement_rejects_unrecorded_edits_and_wrong_original(self):
        for change in ('metadata', 'admitted_object', 'other_object', 'original_object', 'order'):
            with self.subTest(change=change):
                ledger, audit = self.audited_refinement()
                if change == 'metadata':
                    ledger['date'] = '2026-10-09'
                elif change == 'admitted_object':
                    ledger['decisions'][0]['status'] = 'unreviewed'
                elif change == 'other_object':
                    ledger['decisions'][1]['status'] = 'unreviewed'
                elif change == 'original_object':
                    audit['decision_refinements'][0]['original_object']['status'] = 'invented baseline'
                else:
                    ledger['decisions'][:2] = reversed(ledger['decisions'][:2])
                with self.assertRaises(AssertionError):
                    decision_snapshot(ledger_bytes(ledger), self.recorded_count, self.recorded_sha256, audit)


if __name__ == '__main__':
    unittest.main()
