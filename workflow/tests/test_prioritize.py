"""Synthetic checks for command/config consistency and honest completion state."""
import argparse
import importlib.util
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import yaml

SPEC = importlib.util.spec_from_file_location(
    'prioritize', Path(__file__).resolve().parents[1] / 'scripts/prioritize.py')
RUN = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RUN)


class PrioritizeTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name).resolve()
        (self.root/'input.vcf').write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC\n')
        (self.root/'app').mkdir()
        (self.root/'app/exomiser-cli-test.jar').write_bytes(b'not a real jar')
        (self.root/'app/application.properties').write_text('exomiser.hg19.data-version=2512\n')
        for name in ('2512_hg38', '2512_phenotype'):
            (self.root/'data'/name).mkdir(parents=True)
        config = {'analysis': {'genomeAssembly':'hg38', 'proband':'SYNTHETIC',
                              'vcf':'/stale/path.vcf', 'hpoIds':['HP:0000001']},
                  'outputOptions': {'outputDirectory':'/stale/output'}}
        (self.root/'config.yml').write_text(yaml.safe_dump(config))
        self.args = argparse.Namespace(vcf=str(self.root/'input.vcf'), output=str(self.root/'output'),
            exomiser=str(self.root/'app'), data=str(self.root/'data'), analysis=str(self.root/'config.yml'),
            data_version='2512', threads=8, heap='16g', java='synthetic-java', add_hpo=[])

    def status(self):
        return json.loads((self.root/'output/analysis-status.json').read_text())

    def test_input_arguments_replace_stale_config_and_pin_runtime(self):
        def execute(command, **kwargs):
            actual = yaml.safe_load((self.root/'output/analysis.yml').read_text())
            self.assertEqual(actual['analysis']['vcf'], self.args.vcf)
            self.assertEqual(actual['outputOptions']['outputDirectory'], self.args.output)
            self.assertIn('-XX:ActiveProcessorCount=8', command)
            self.assertEqual(command[command.index('-jar') + 2], 'analyse')
            self.assertIn('-Dexomiser.hg38.data-version=2512', command[:command.index('-jar')])
            self.assertFalse(any(arg.startswith('--exomiser.') or arg.startswith('--spring.') for arg in command))
            self.assertNotIn('--exomiser.hg19.data-version=', command)
            self.assertNotIn('exomiser.hg19.data-version=', (self.root/'output/application.properties').read_text())
            (self.root/'output/synthetic.genes.tsv').write_text('#RANK\tGENE_SYMBOL\tMOI\tEXOMISER_GENE_COMBINED_SCORE\n1\tSYNTHETIC\tAR\t0.5\n')
        with patch.object(RUN.subprocess, 'run', side_effect=execute):
            result = RUN.run(self.args)
        self.assertEqual(self.status()['status'], 'complete')
        self.assertEqual(result['top_genes'][0]['GENE_SYMBOL'], 'SYNTHETIC')
        self.assertEqual(result['top_genes'][0]['MOI'], 'AR')
        with self.assertRaisesRegex(ValueError, 'existing run preserved'):
            RUN.run(self.args)
        self.assertEqual(self.status()['status'], 'complete')

    def test_zero_exit_without_results_is_a_failed_run(self):
        with patch.object(RUN.subprocess, 'run'), self.assertRaisesRegex(ValueError, 'gene ranking table'):
            RUN.run(self.args)
        self.assertEqual(self.status()['status'], 'failed')
        self.assertFalse((self.root/'output/gene-summary.json').exists())

    def test_java_failure_cannot_be_reported_as_complete(self):
        with patch.object(RUN.subprocess, 'run', side_effect=subprocess.CalledProcessError(1, ['synthetic-java'])), self.assertRaises(subprocess.CalledProcessError):
            RUN.run(self.args)
        self.assertEqual(self.status()['status'], 'failed')

    def test_sample_mismatch_fails_before_execution(self):
        config = yaml.safe_load((self.root/'config.yml').read_text())
        config['analysis']['proband'] = 'ABSENT'
        (self.root/'config.yml').write_text(yaml.safe_dump(config))
        with patch.object(RUN.subprocess, 'run') as execute, self.assertRaisesRegex(ValueError, 'absent'):
            RUN.run(self.args)
        execute.assert_not_called()

    def test_unlabelled_existing_outputs_are_preserved(self):
        output = self.root/'output'
        output.mkdir()
        (output/'old.genes.tsv').write_text('retained result')
        with patch.object(RUN.subprocess, 'run') as execute, self.assertRaisesRegex(ValueError, 'existing run preserved'):
            RUN.run(self.args)
        execute.assert_not_called()
        self.assertEqual((output/'old.genes.tsv').read_text(), 'retained result')

    def test_valid_header_without_candidates_records_an_empty_result(self):
        def execute(*args, **kwargs):
            (self.root/'output/synthetic.genes.tsv').write_text('#RANK\tGENE_SYMBOL\n')
        with patch.object(RUN.subprocess, 'run', side_effect=execute):
            result = RUN.run(self.args)
        self.assertEqual(self.status()['status'], 'complete')
        self.assertEqual(result['reported_gene_rows'], 0)
        self.assertEqual(result['top_genes'], [])


if __name__ == '__main__':
    unittest.main()
