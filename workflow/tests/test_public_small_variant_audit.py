"""Check overlap and GT diagnostics using already-executed native RTG controls.

The two caller labels here select fictional comparator outputs, not real caller
performance. Literal one-site outcomes are independent of the audit implementation.
"""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


class NativeOutputAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root = Path(os.environ.get('MVA_PROJECT_ROOT', '/Users/jeweled/Documents/mva-hackathon'))
        cls.controls = cls.root / 'results/germline-public-comparator-controls-20260919'
        if not (cls.controls / 'summary.json').is_file():
            raise unittest.SkipTest('Completed native RTG control outputs required')
        recorded = json.loads((cls.controls / 'summary.json').read_text())
        if recorded['status'] != 'complete' or recorded['cases_passed'] != 12:
            raise RuntimeError('Native control study is not complete')
        cls.cases = {row['name']: row['comparison'] for row in recorded['cases']}
        cls.script = Path(__file__).parents[1] / 'scripts/audit_public_small_variants.py'

    def exercise(self, first, second, expected, mismatch_gt=None, expect_failure=False):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            comparison = base / 'comparison'
            comparison.mkdir()
            def pin(path):
                raw = path.read_bytes()
                return {'path': str(path), 'bytes': len(raw), 'sha256': hashlib.sha256(raw).hexdigest()}
            summary = {'status': 'complete', 'identical_scored_truth_universe': True,
                       'scope': 'Fictional native-output audit control; caller labels are bookkeeping only',
                       'callers': {'gatk': self.cases[first], 'bcftools': self.cases[second]},
                       'inputs': {'truth': pin(self.controls / first / 'baseline.vcf.gz'),
                                  'gatk': pin(self.controls / first / 'calls.vcf.gz'),
                                  'bcftools': pin(self.controls / second / 'calls.vcf.gz')}}
            (comparison / 'summary.json').write_text(json.dumps(summary))
            command = [sys.executable, str(self.script), '--comparison', str(comparison),
                       '--output', str(base / 'audit'), '--sites-per-class', '2']
            process = subprocess.run(command, capture_output=True, text=True)
            if expect_failure:
                self.assertNotEqual(process.returncode, 0)
                self.assertIn('Caller truth universes differ', process.stderr)
                self.assertFalse((base / 'audit/summary.json').exists())
                return
            self.assertEqual(process.returncode, 0, process.stdout + process.stderr)
            result = json.loads((base / 'audit/summary.json').read_text())
            self.assertEqual(result['scored_truth_records'], 1)
            self.assertEqual({k: v['records'] for k, v in result['truth_partitions'].items()}, expected)
            self.assertIsNone(result['matched_by_at_least_one_caller']['merged_precision'])
            self.assertFalse(result['matched_by_at_least_one_caller']['is_validated_merged_callset'])
            if mismatch_gt is not None:
                self.assertEqual({k: v['fp_exact_raw_alleles_with_discordant_GT']
                                  for k, v in result['output_diagnostics'].items()}, mismatch_gt)
            if first == 'mnp_as_two_snvs':
                self.assertEqual(result['output_diagnostics']['gatk']['tp_calls']['records'], 2)
                self.assertEqual(result['truth_partitions']['matched_both']['record_classes'], {'mnp': 1})
            # A repeated output path must not overwrite the completed audit.
            before = (base / 'audit/summary.json').read_bytes()
            repeated = subprocess.run(command, capture_output=True, text=True)
            self.assertNotEqual(repeated.returncode, 0)
            self.assertEqual((base / 'audit/summary.json').read_bytes(), before)

    def test_match_by_both_preserves_one_truth_record(self):
        self.exercise('zero_qual_is_included', 'zero_qual_is_included',
                      {'matched_both': 1, 'gatk_only': 0, 'bcftools_only': 0, 'missed_both': 0})

    def test_first_only_and_same_allele_wrong_GT(self):
        self.exercise('zero_qual_is_included', 'correct_allele_wrong_genotype',
                      {'matched_both': 0, 'gatk_only': 1, 'bcftools_only': 0, 'missed_both': 0},
                      {'gatk': 0, 'bcftools': 1})

    def test_second_only_is_not_assigned_to_first(self):
        self.exercise('correct_allele_wrong_genotype', 'zero_qual_is_included',
                      {'matched_both': 0, 'gatk_only': 0, 'bcftools_only': 1, 'missed_both': 0},
                      {'gatk': 1, 'bcftools': 0})

    def test_both_misses_distinguish_GT_from_ALT_diagnostic(self):
        self.exercise('wrong_alternate_allele', 'correct_allele_wrong_genotype',
                      {'matched_both': 0, 'gatk_only': 0, 'bcftools_only': 0, 'missed_both': 1},
                      {'gatk': 0, 'bcftools': 1})

    def test_two_call_records_do_not_inflate_truth_overlap(self):
        self.exercise('mnp_as_two_snvs', 'mnp_as_two_snvs',
                      {'matched_both': 1, 'gatk_only': 0, 'bcftools_only': 0, 'missed_both': 0})

    def test_different_genotype_truth_universe_is_rejected(self):
        self.exercise('zero_qual_is_included', 'correct_homozygous_genotype', {}, expect_failure=True)


if __name__ == '__main__':
    unittest.main()
