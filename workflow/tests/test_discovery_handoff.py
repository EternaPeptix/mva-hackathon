"""Native, fictional handoff controls; no caller execution or case data required.

The GATK/bcftools labels describe file containers in a miniature completed shard.
They are not measurements of either caller's sensitivity or biological accuracy.
Every normalization and replay uses the installed Studio bcftools executable.
"""
import copy
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import pysam


SCRIPT = Path(__file__).resolve().parents[1] / 'scripts/audit_discovery_handoff.py'
AUTOSOMES = [str(number) for number in range(1, 23)]


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


class FictionalShard:
    """One completed whole-autosome shard and 21 explicit pending autosomes."""

    def __init__(self, folder, bcftools, empty=False):
        self.folder, self.bcftools = folder.resolve(), bcftools.resolve()
        self.run = self.folder / 'calling'
        self.attempt = self.run / 'shard-0000/attempt-fictional'
        self.attempt.mkdir(parents=True)
        self.reference = self.folder / 'reference.fa'
        self.reference.write_text(''.join('>' + name + '\n' + 'A' * 200 + '\n'
                                          for name in AUTOSOMES))
        pysam.faidx(str(self.reference))
        self.definition = {
            'reference_sha256': sha256(self.reference), 'sample': 'SYNTH',
            'contigs': 22, 'reference_bases': 4400, 'ploidy': 2,
            'shards': [[[name, 200]] for name in AUTOSOMES],
        }
        self.definition_path = self.run / 'run-definition.json'
        write_json(self.definition_path, self.definition)
        self.quality_path = self.folder / 'quality.json'
        write_json(self.quality_path, {'sample': 'SYNTH', 'quality': {
            'explicit_PASS': True, 'fully_called_diploid': True,
            'minimum_DP': 10, 'minimum_GQ': 20,
        }})
        self.header = (
            '##fileformat=VCFv4.2\n'
            + ''.join('##contig=<ID=' + name + ',length=200>\n' for name in AUTOSOMES)
            + '##FILTER=<ID=PASS,Description="All filters passed">\n'
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
            '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n'
            '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">\n'
            '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Likelihoods">\n'
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTH\n'
        )
        self.rows = [] if empty else [
            '1\t10\t.\tA\tC,G\t100\tPASS\t.\tGT:DP:GQ:AD:PL\t1/2:21:35:1,10,10:100,100,100,100,0,100\n',
            '1\t20\t.\tA\tC\t100\tPASS\t.\tGT:DP:GQ:AD:PL\t0/1:20:40:10,10:100,0,100\n',
            # GQ is absent from FORMAT, rather than inferred from the sharp PL.
            '1\t30\t.\tA\tG\t100\t.\t.\tGT:DP:AD:PL\t0/1:20:10,10:100,0,100\n',
            # Explicit zero GQ is present and remains distinguishable from missing.
            '1\t40\t.\tA\tT\t100\tPASS\t.\tGT:DP:GQ:AD:PL\t0/1:20:0:10,10:100,0,100\n',
        ]
        self.seed = self.folder / 'fictional.vcf'
        self.seed.write_text(self.header + ''.join(self.rows))
        self.raw = {'gatk': self.attempt / 'gatk.raw.vcf.gz',
                    'bcftools': self.attempt / 'bcftools.raw.bcf'}
        self.rewrite_gatk_raw(self.rows)
        self.native(['view', '-Ob', '-o', str(self.raw['bcftools']), str(self.seed)])
        self.native(['index', str(self.raw['bcftools'])])
        commands, self.receipt, self.normalized = [], {}, {}
        for method in ('gatk', 'bcftools'):
            normalized = self.attempt / (method + '.normalized.raw.vcf.gz')
            self.normalized[method] = normalized
            normalize = ['bcftools', 'norm', '--check-ref', 'e', '-f', str(self.reference),
                         '-m', '-any', '-Oz', '-o', str(normalized), str(self.raw[method])]
            commands.append(normalize)
            self.native(normalize[1:])
            index = ['bcftools', 'index', '-t', str(normalized)]
            commands.append(index)
            self.native(index[1:])
            with pysam.VariantFile(str(normalized)) as reader:
                records = sum(1 for _ in reader.fetch())
            self.receipt[method] = {
                'path': str(normalized.relative_to(self.attempt.parent)),
                'sha256': sha256(normalized), 'counts': {'records': records},
            }
        self.commands_path = self.attempt / 'commands.json'
        write_json(self.commands_path, commands)
        self.receipt['commands_sha256'] = sha256(self.commands_path)
        self.receipt_path = self.attempt.parent / 'complete.json'
        write_json(self.receipt_path, self.receipt)

    def native(self, arguments):
        subprocess.run([str(self.bcftools), *arguments], capture_output=True,
                       text=True, check=True)

    def rewrite_gatk_raw(self, rows):
        self.seed.write_text(self.header + ''.join(rows))
        pysam.tabix_compress(str(self.seed), str(self.raw['gatk']), force=True)
        pysam.tabix_index(str(self.raw['gatk']), preset='vcf', force=True)

    def execute(self, name, require_complete=False):
        output = self.folder / name
        command = [sys.executable, str(SCRIPT), '--run', str(self.run),
                   '--reference', str(self.reference), '--bcftools', str(self.bcftools),
                   '--quality-config', str(self.quality_path), '--output', str(output),
                   '--expected-definition-sha256', sha256(self.definition_path)]
        if require_complete:
            command.append('--require-complete')
        result = subprocess.run(command, capture_output=True, text=True)
        return result, output

    def snapshot(self):
        """Detect any byte or modification-time change to the small source fixture."""
        paths = [path for path in self.run.rglob('*') if path.is_file()]
        paths.extend([self.reference, Path(str(self.reference) + '.fai'), self.quality_path])
        return {str(path): (path.read_bytes(), path.stat().st_mtime_ns) for path in paths}


class NativeDiscoveryHandoffTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root = Path(os.environ.get('MVA_PROJECT_ROOT', '/Users/jeweled/Documents/mva-hackathon'))
        cls.bcftools = root / 'tools/bio/bin/bcftools'
        if not cls.bcftools.is_file():
            raise unittest.SkipTest('installed Studio pysam and bcftools runtime required')

    def assert_failed(self, result, output, reason):
        self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        state = json.loads((output / 'status.json').read_text())
        self.assertEqual(state['status'], 'failed')
        self.assertIn(reason, state['error'])
        self.assertFalse((output / 'summary.json').exists())

    def test_native_multiallelic_quality_missingness_and_input_preservation(self):
        with tempfile.TemporaryDirectory() as temporary:
            fixture = FictionalShard(Path(temporary), self.bcftools)
            before = fixture.snapshot()
            process, output = fixture.execute('audit')
            self.assertEqual(process.returncode, 0, process.stdout + process.stderr)
            summary = json.loads((output / 'summary.json').read_text())
            self.assertEqual(summary['selected_autosomes'], ['1'])
            self.assertEqual(summary['pending_autosomes'], AUTOSOMES[1:])
            self.assertFalse(summary['all_autosomes_complete'])
            self.assertEqual(len(summary['shards']), 1)
            self.assertEqual(summary['shards'][0]['contigs'], [['1', 200]])
            self.assertTrue(any('absence of a shard/record is not reference evidence' in item
                                for item in summary['limitations']))
            for method, reviewed in summary['shards'][0]['methods'].items():
                with self.subTest(method=method):
                    raw = reviewed['original_raw']
                    self.assertEqual(raw['counts']['records'], 4)
                    self.assertEqual(raw['counts']['retained_ALT_alleles'], 5)
                    self.assertEqual(raw['counts']['GT_positive_ALT_alleles'], 5)
                    self.assertEqual(raw['counts']['two_distinct_nonreference_diploid_GT_records'], 1)
                    self.assertEqual(raw['counts']['multiallelic_records'], 1)
                    self.assertEqual(raw['filters'], {'PASS': 3, '.': 1})
                    self.assertEqual(raw['FORMAT_availability']['GQ'], {'present': 3, 'missing': 1})
                    self.assertEqual(raw['counts']['records_passing_existing_preannotation_quality_gate'], 2)
                    self.assertEqual(raw['counts']['GT_positive_concrete_ALT_alleles_passing_existing_preannotation_quality_gate'], 3)
                    failures = raw['existing_gate_nonexclusive_record_failure_counts']
                    self.assertEqual(failures['record_not_explicit_PASS'], 1)
                    self.assertEqual(failures['genotype_quality_missing_or_below_threshold'], 2)
                    self.assertEqual(reviewed['normalized_split']['counts']['records'], 5)
                    self.assertTrue(reviewed['replayed_split_payload_identical'])
                    # Read native source directly: split summaries cannot replace GT/AD/PL.
                    with pysam.VariantFile(str(fixture.raw[method])) as reader:
                        records = list(reader)
                    self.assertEqual(records[0].alts, ('C', 'G'))
                    self.assertEqual(records[0].samples['SYNTH']['GT'], (1, 2))
                    self.assertEqual(records[0].samples['SYNTH']['AD'], (1, 10, 10))
                    self.assertEqual(records[0].samples['SYNTH']['PL'], (100, 100, 100, 100, 0, 100))
                    self.assertIsNone(records[2].samples['SYNTH'].get('GQ'))
                    self.assertEqual(records[3].samples['SYNTH']['GQ'], 0)
            self.assertEqual(fixture.snapshot(), before)
            # A second invocation must preserve the completed output directory too.
            summary_before = (output / 'summary.json').read_bytes()
            repeated, _ = fixture.execute('audit')
            self.assertNotEqual(repeated.returncode, 0)
            self.assertEqual((output / 'summary.json').read_bytes(), summary_before)
            self.assertEqual(fixture.snapshot(), before)

    def test_changed_raw_genotype_cannot_acquire_trust_from_a_new_fingerprint(self):
        with tempfile.TemporaryDirectory() as temporary:
            fixture = FictionalShard(Path(temporary), self.bcftools)
            frozen_receipt = fixture.receipt_path.read_bytes()
            frozen_normalized = {method: path.read_bytes() for method, path in fixture.normalized.items()}
            original_raw_hash = sha256(fixture.raw['gatk'])
            changed = list(fixture.rows)
            changed[0] = changed[0].replace('1/2:21:', '0/1:21:')
            fixture.rewrite_gatk_raw(changed)
            self.assertNotEqual(sha256(fixture.raw['gatk']), original_raw_hash)
            before_audit = fixture.snapshot()
            process, output = fixture.execute('changed-raw-audit')
            self.assert_failed(process, output, 'does not reproduce frozen split record payloads')
            self.assertEqual(fixture.snapshot(), before_audit)
            self.assertEqual(fixture.receipt_path.read_bytes(), frozen_receipt)
            for method, path in fixture.normalized.items():
                self.assertEqual(path.read_bytes(), frozen_normalized[method])

    def test_pending_duplicate_mixed_and_escaping_inputs_are_refused(self):
        with tempfile.TemporaryDirectory() as temporary:
            fixture = FictionalShard(Path(temporary), self.bcftools)
            before = fixture.snapshot()
            process, output = fixture.execute('require-complete', require_complete=True)
            self.assert_failed(process, output, 'autosomal calling is incomplete')
            self.assertEqual(fixture.snapshot(), before)

            duplicate = copy.deepcopy(fixture.definition)
            duplicate['shards'][1] = [['1', 200]]
            mixed = copy.deepcopy(fixture.definition)
            mixed['shards'][0].append(['X', 200])
            mixed['contigs'], mixed['reference_bases'] = 23, 4600
            for name, definition, reason in (
                    ('duplicate', duplicate, 'inconsistent discovery definition'),
                    ('mixed', mixed, 'autosome/nonautosome mixed shard')):
                with self.subTest(definition=name):
                    write_json(fixture.definition_path, definition)
                    changed_before = fixture.snapshot()
                    process, output = fixture.execute(name)
                    self.assert_failed(process, output, reason)
                    self.assertEqual(fixture.snapshot(), changed_before)
            write_json(fixture.definition_path, fixture.definition)

            for name, escaped in (('relative-escape', '../outside.vcf.gz'),
                                  ('absolute-escape', str(fixture.folder / 'outside.vcf.gz'))):
                with self.subTest(path=name):
                    receipt = copy.deepcopy(fixture.receipt)
                    receipt['gatk']['path'] = escaped
                    write_json(fixture.receipt_path, receipt)
                    changed_before = fixture.snapshot()
                    process, output = fixture.execute(name)
                    self.assert_failed(process, output, 'shard output path escapes its directory')
                    self.assertEqual(fixture.snapshot(), changed_before)

    def test_empty_indexed_shard_is_valid_but_does_not_complete_pending_chromosomes(self):
        with tempfile.TemporaryDirectory() as temporary:
            fixture = FictionalShard(Path(temporary), self.bcftools, empty=True)
            before = fixture.snapshot()
            process, output = fixture.execute('empty-audit')
            self.assertEqual(process.returncode, 0, process.stdout + process.stderr)
            summary = json.loads((output / 'summary.json').read_text())
            self.assertEqual(summary['selected_autosomes'], ['1'])
            self.assertEqual(summary['pending_autosomes'], AUTOSOMES[1:])
            self.assertFalse(summary['all_autosomes_complete'])
            for reviewed in summary['shards'][0]['methods'].values():
                self.assertEqual(reviewed['original_raw']['counts']['records'], 0)
                self.assertEqual(reviewed['normalized_split']['counts']['records'], 0)
                self.assertEqual(reviewed['normalized_split']['rendered_record_payload_sha256'],
                                 hashlib.sha256(b'').hexdigest())
                self.assertTrue(reviewed['replayed_split_payload_identical'])
            self.assertEqual(fixture.snapshot(), before)

    def test_readable_empty_index_cannot_validate_a_nonempty_normalized_shard(self):
        with tempfile.TemporaryDirectory() as temporary:
            fixture = FictionalShard(Path(temporary), self.bcftools)
            empty = fixture.folder / 'empty-index-source.vcf'
            empty.write_text(fixture.header)
            empty_gz = Path(str(empty) + '.gz')
            pysam.tabix_compress(str(empty), str(empty_gz), force=True)
            pysam.tabix_index(str(empty_gz), preset='vcf', force=True)
            # The substitute is a readable real tabix index, not arbitrary damage.
            with pysam.VariantFile(str(empty_gz)) as reader:
                self.assertIsNotNone(reader.index)
                self.assertEqual(list(reader.fetch('1', 0, 200)), [])
            normalized = fixture.normalized['gatk']
            normalized_bytes = normalized.read_bytes()
            target_index = Path(str(normalized) + '.tbi')
            original_index = target_index.read_bytes()
            target_index.write_bytes(Path(str(empty_gz) + '.tbi').read_bytes())
            self.assertNotEqual(target_index.read_bytes(), original_index)
            before_audit = fixture.snapshot()
            process, output = fixture.execute('stale-index-audit')
            self.assertNotEqual(process.returncode, 0, process.stdout + process.stderr)
            status = json.loads((output / 'status.json').read_text())
            self.assertEqual(status['status'], 'failed')
            self.assertRegex(status['error'].lower(),
                             'record count mismatch|index|iterator|fetch|truncated|read|seek|invalid|bgzf')
            self.assertFalse((output / 'summary.json').exists())
            self.assertEqual(normalized.read_bytes(), normalized_bytes)
            self.assertEqual(fixture.snapshot(), before_audit)


if __name__ == '__main__':
    unittest.main()
