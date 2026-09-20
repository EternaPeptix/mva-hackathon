import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / 'scripts/summarize_lirical.py'
spec = importlib.util.spec_from_file_location('summary', SCRIPT)
summary = importlib.util.module_from_spec(spec)
spec.loader.exec_module(summary)


class LiricalSummaryTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.hpo = ['HP:0002859']
        self.data = {
            'analysisData': {'sampleId': 'fixture', 'observedPhenotypicFeatures': self.hpo,
                             'excludedPhenotypicFeatures': []},
            'analysisMetadata': {'isGlobalAnalysisMode': False},
            'analysisResults': []}
        # Intentionally put the best diagnosis last and hide it in the TSV,
        # matching the real distinction between JSON and a filtered report.
        for ident, lr in [('OMIM:257300', 2.), ('OMIM:241200', 9.)]:
            self.data['analysisResults'].append({
                'diseaseId': ident, 'pretestProbability': .1, 'compositeLR': lr,
                'posttestProbability': lr / (9 + lr),
                'observedPhenotypicFeatures': [{'query': self.hpo[0], 'lr': lr}],
                'excludedPhenotypicFeatures': [],
                'genotypeLR': {'geneId': {'id': 'NCBIGene:1', 'symbol': 'FIXTURE'}, 'lr': 1.}})
        self.state = {'status': 'complete', 'mode': 'genotype_aware', 'observed': self.hpo,
                      'negated': [], 'command': ['--sample-id', 'fixture'],
                      'started_at': 'synthetic', 'completed_at': 'synthetic', 'elapsed_seconds': 1.}
        (self.root / 'lirical.html').write_text('synthetic fixture')
        (self.root / 'lirical.log').write_text('WARN source - Unknown contig M\n')
        with (self.root / 'lirical.tsv').open('w') as stream:
            writer = csv.DictWriter(stream, fieldnames=['rank', 'diseaseName', 'diseaseCurie',
                                   'pretestprob', 'posttestprob', 'compositeLR'], delimiter='\t')
            writer.writeheader()
            writer.writerow({'rank': '1', 'diseaseName': 'Synthetic disease',
                             'diseaseCurie': 'OMIM:257300', 'pretestprob': '1/10',
                             'posttestprob': '18.18%', 'compositeLR': '0.301'})
        self.save()

    def save(self):
        (self.root / 'lirical.json').write_text(json.dumps(self.data))
        self.state['output_hashes'] = {name: summary.digest(self.root / name)
            for name in ('lirical.json', 'lirical.tsv', 'lirical.html')}
        (self.root / 'run-status.json').write_text(json.dumps(self.state))

    def test_full_json_order_and_hidden_disease_are_preserved(self):
        result = summary.summarize_run(self.root)
        self.assertEqual(result['top_json_results'][0]['diseaseId'], 'OMIM:241200')
        self.assertIsNone(result['top_json_results'][0]['native_tsv_rank'])
        self.assertEqual(result['json_diseases_absent_from_tsv'], 1)
        self.assertEqual(result['warnings']['unknown_primary_contig_records'], {'M': 1})

    def test_tampering_and_wrong_phenotype_fail_validation(self):
        (self.root / 'lirical.tsv').write_text('tampered')
        with self.assertRaisesRegex(ValueError, 'checksum mismatch'):
            summary.summarize_run(self.root)
        self.save()
        self.data['analysisData']['observedPhenotypicFeatures'] = ['HP:0000121']
        self.save()
        with self.assertRaisesRegex(ValueError, 'HPO features'):
            summary.summarize_run(self.root)

    def test_inconsistent_model_probability_is_rejected(self):
        self.data['analysisResults'][0]['posttestProbability'] = .95
        self.save()
        with self.assertRaisesRegex(ValueError, 'post-test probability'):
            summary.summarize_run(self.root)


if __name__ == '__main__':
    unittest.main()
