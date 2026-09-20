import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import pysam

SCRIPT = Path(__file__).resolve().parents[1] / 'scripts/build_submission.py'
spec = importlib.util.spec_from_file_location('submission', SCRIPT)
submission = importlib.util.module_from_spec(spec)
spec.loader.exec_module(submission)


class SubmissionTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.fa = self.root / 'tiny.fa'
        self.fa.write_text('>15\nACGTACGTACGT\n')
        pysam.faidx(str(self.fa))
        self.reference = pysam.FastaFile(str(self.fa))
        self.addCleanup(self.reference.close)
        self.review = {
            'schema_version': 1, 'assembly': 'GRCh38', 'proband_id': 'PROBAND01',
            'team': 'fixture', 'approach': 'validation', 'review_status': 'reviewed',
            'reviewed_by': 'synthetic fixture', 'reviewed_at': '2026-09-19',
            'reference_sha256': hashlib.sha256(self.fa.read_bytes()).hexdigest(),
            'candidates': [{'epcr': .5, 'epcr_rationale': 'Synthetic test value, not a patient estimate.',
                            'finding_type': 'primary', 'notes': 'Synthetic alleles.', 'phase': 'unknown',
                            'variants': [{'chrom': 'chr15', 'pos': 2, 'ref': 'C', 'alt': 'T'},
                                         {'chrom': 'chr15', 'pos': 5, 'ref': 'A', 'alt': 'G'}]}]}

    def test_phase_and_reference_fail_closed(self):
        rows = submission.build_rows(self.review, self.reference)
        self.assertIn('phase unknown', rows[0]['notes'])
        for value in ('cis', 'trans', None):
            changed = copy.deepcopy(self.review)
            changed['candidates'][0]['phase'] = value
            with self.assertRaises(ValueError):
                submission.build_rows(changed, self.reference)
        changed = copy.deepcopy(self.review)
        changed['candidates'][0]['variants'][0]['ref'] = 'G'
        with self.assertRaisesRegex(ValueError, 'REF mismatch'):
            submission.build_rows(changed, self.reference)

    def test_probability_is_explicit_and_pair_identity_unordered(self):
        for probability in (None, True, float('nan'), float('inf'), 0, 1.1):
            changed = copy.deepcopy(self.review)
            changed['candidates'][0]['epcr'] = probability
            with self.assertRaises(ValueError):
                submission.build_rows(changed, self.reference)
        changed = copy.deepcopy(self.review)
        second = copy.deepcopy(changed['candidates'][0])
        second['variants'].reverse()
        changed['candidates'].append(second)
        with self.assertRaisesRegex(ValueError, 'duplicate candidate'):
            submission.build_rows(changed, self.reference)
        changed = copy.deepcopy(self.review)
        changed['candidates'][0]['epcr_rationale'] = ''
        with self.assertRaisesRegex(ValueError, 'epcr_rationale'):
            submission.build_rows(changed, self.reference)

    def test_honest_ties_and_secondary_priority_are_preserved(self):
        second = copy.deepcopy(self.review['candidates'][0])
        second.update(epcr=.8, finding_type='secondary', phase='not_applicable')
        second['variants'] = [{'chrom': 'chr15', 'pos': 7, 'ref': 'G', 'alt': 'A'}]
        self.review['candidates'].append(second)
        rows = submission.build_rows(self.review, self.reference)
        self.assertEqual(rows[0]['finding_type'], 'secondary')
        self.assertTrue(all(rows[0][k] == '' for k in ('chrom_2', 'pos_2', 'ref_2', 'alt_2')))
        second['epcr'] = .5
        rows = submission.build_rows(self.review, self.reference)
        self.assertEqual([r['finding_type'] for r in rows], ['primary', 'secondary'])

    def test_cli_provenance_and_preservation(self):
        review_file = self.root / 'review.json'
        review_file.write_text(json.dumps(self.review))
        output = self.root / 'fixture_validation.csv'
        cmd = [sys.executable, str(SCRIPT), '--review', str(review_file), '--reference', str(self.fa), '--out', str(output)]
        result = subprocess.run(cmd, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        before = output.read_bytes()
        receipt = json.loads(Path(str(output) + '.provenance.json').read_text())
        self.assertEqual(receipt['csv_sha256'], hashlib.sha256(before).hexdigest())
        self.assertNotEqual(subprocess.run(cmd, capture_output=True).returncode, 0)
        self.assertEqual(output.read_bytes(), before)
        self.review['review_status'] = 'draft'
        review_file.write_text(json.dumps(self.review))
        output.unlink()
        Path(str(output) + '.provenance.json').unlink()
        self.assertNotEqual(subprocess.run(cmd, capture_output=True).returncode, 0)
        self.assertFalse(output.exists())


if __name__ == '__main__':
    unittest.main()
