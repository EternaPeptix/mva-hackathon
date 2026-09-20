"""Native integration controls using a fictional gene and miniature handoff.

These checks exercise real bcftools normalization and local consequence annotation.
The model combinations are descriptive test objects, never diagnoses or estimates
of caller accuracy. No case genome, actual calling job or public database is used.
"""
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_discovery_handoff import AUTOSOMES, FictionalShard, sha256, write_json


WORKFLOW = Path(__file__).resolve().parents[1]
SCRIPT = WORKFLOW / 'scripts/prepare_native_coding_review.py'
EXPECTED_KEYS = {'1:10:A:C', '1:10:A:G', '1:20:A:C', '1:30:A:G', '1:40:A:T'}


def read_jsonl(path):
    with gzip.open(path, 'rt') as stream:
        return [json.loads(line) for line in stream]


def pinned_bytes(paths):
    return {str(path): (path.read_bytes(), path.stat().st_mtime_ns) for path in paths}


class FictionalCodingSnapshot:
    def __init__(self, folder, runtime, empty=False):
        self.runtime = runtime.resolve()
        self.fixture = FictionalShard(folder, runtime / 'tools/bio/bin/bcftools', empty=empty)
        self.folder = self.fixture.folder
        # A 99-base complete CDS: start ATG, terminal TAA and mostly lysines.
        # Codon 28..30 is ATA: the GQ-absent 30:A>G call becomes missense,
        # rather than disappearing as the synonymous AAA>AAG consequence.
        coding = list('ATG' + 'A' * 93 + 'TAA')
        coding[28] = 'T'
        sequence = ''.join(coding) + 'A' * 101
        self.fixture.reference.write_text(''.join(
            '>' + chrom + '\n' + (sequence if chrom == '1' else 'A' * 200) + '\n'
            for chrom in AUTOSOMES))
        pysam.faidx(str(self.fixture.reference))
        self.fixture.definition['reference_sha256'] = sha256(self.fixture.reference)
        write_json(self.fixture.definition_path, self.fixture.definition)
        # The original SNP REF/ALT and normalized payloads remain unchanged.
        # The actual handoff audit must replay them against this pinned FASTA.
        process, handoff_output = self.fixture.execute('completed-handoff')
        if process.returncode:
            raise AssertionError('fictional handoff failed:\n' + process.stdout + process.stderr)
        self.handoff = handoff_output / 'summary.json'
        self.gff = self.folder / 'fictional-coding.gff3'
        self.gff.write_text(
            '##gff-version 3\n'
            '##sequence-region 1 1 200\n'
            '1\tensembl\tgene\t1\t99\t.\t+\t.\tID=gene:ENSG00000000001;Name=FICTIONAL;biotype=protein_coding\n'
            '1\tensembl\tmRNA\t1\t99\t.\t+\t.\tID=transcript:ENST00000000001;Parent=gene:ENSG00000000001;Name=FICTIONAL-001;biotype=protein_coding\n'
            '1\tensembl\texon\t1\t99\t.\t+\t.\tID=exon:ENSE00000000001;Parent=transcript:ENST00000000001\n'
            '1\tensembl\tCDS\t1\t99\t.\t+\t0\tID=CDS:ENSP00000000001;Parent=transcript:ENST00000000001;protein_id=ENSP00000000001\n')
        self.config = self.folder / 'coding-config.json'
        configuration = json.loads((WORKFLOW / 'configs/native-coding-review.json').read_text())
        configuration['sample'] = 'SYNTH'
        configuration['frequency']['source'] = 'Fictional test AC/AN values; no population database used'
        write_json(self.config, configuration)
        self.frequencies = self.folder / 'fictional-frequency-responses.jsonl'
        self.frequency_rows = []
        for key in ([] if empty else sorted(EXPECTED_KEYS)):
            row = {'key': key, 'allele_found': False, 'frequencies': []}
            if key == '1:20:A:C':
                row.update(allele_found=True,
                           frequencies=[{'source': 'FICTIONAL', 'ac': 1, 'an': 10000}])
            self.frequency_rows.append(row)
        self.frequencies.write_text(''.join(json.dumps(row) + '\n' for row in self.frequency_rows))

    def execute(self, name, caller='gatk', allow_partial=True, expected_gff=None):
        output = self.folder / name
        command = [sys.executable, str(SCRIPT), '--root', str(self.runtime),
                   '--handoff', str(self.handoff), '--expected-handoff-sha256', sha256(self.handoff),
                   '--reference', str(self.fixture.reference), '--gff', str(self.gff),
                   '--expected-gff-sha256', expected_gff or sha256(self.gff),
                   '--config', str(self.config), '--caller', caller,
                   '--frequency-responses', str(self.frequencies),
                   '--expected-frequency-sha256', sha256(self.frequencies), '--output', str(output)]
        if allow_partial:
            command.append('--allow-partial-snapshot')
        return subprocess.run(command, capture_output=True, text=True), output

    def source_paths(self):
        paths = [Path(name) for name in self.fixture.snapshot()]
        paths.extend(path for path in self.handoff.parent.rglob('*') if path.is_file())
        paths.extend([self.gff, self.config, self.frequencies])
        return paths


class NativeCodingPreparationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.runtime = Path(os.environ.get('MVA_PROJECT_ROOT', '/Users/jeweled/Documents/mva-hackathon'))
        if not (cls.runtime / 'tools/bio/bin/bcftools').is_file():
            raise unittest.SkipTest('installed Studio pysam and bcftools runtime required')

    def assert_complete(self, process, output):
        logs = '\n'.join(path.read_text(errors='replace')[-4000:] for path in output.rglob('*.log'))
        self.assertEqual(process.returncode, 0, process.stdout + process.stderr + logs)
        summary = json.loads((output / 'summary.json').read_text())
        self.assertEqual(summary['status'], 'complete')
        return summary

    def assert_failed(self, process, output, reason):
        self.assertNotEqual(process.returncode, 0, process.stdout + process.stderr)
        status = json.loads((output / 'status.json').read_text())
        self.assertEqual(status['status'], 'failed')
        self.assertIn(reason, status['error'])
        self.assertFalse((output / 'summary.json').exists())

    def test_native_annotation_recovers_both_original_alleles_and_keeps_quality_missingness(self):
        with tempfile.TemporaryDirectory() as temporary:
            fixture = FictionalCodingSnapshot(Path(temporary), self.runtime)
            paths = fixture.source_paths()
            before = pinned_bytes(paths)
            source_ids = {}
            for caller in ('gatk', 'bcftools'):
                with self.subTest(caller=caller):
                    process, output = fixture.execute(caller + '-coding', caller=caller)
                    summary = self.assert_complete(process, output)
                    self.assertEqual(summary['caller'], caller)
                    self.assertEqual(summary['selected_autosomes'], ['1'])
                    self.assertEqual(summary['pending_autosomes'], AUTOSOMES[1:])
                    self.assertFalse(summary['all_autosomes_complete'])
                    self.assertEqual(summary['coding_allele_rows'], 5)
                    self.assertEqual(summary['source_records_with_coding_rows'], 4)
                    self.assertEqual(summary['frequency_responses'], 5)
                    self.assertEqual(summary['unique_normalized_keys'], 5)
                    self.assertEqual(summary['duplicate_normalized_keys'], 0)
                    self.assertTrue(any('not callable-reference evidence' in value
                                        for value in summary['limitations']))
                    shard = summary['shards'][0]
                    self.assertTrue(shard['original_copy_check']['derived_copy_payload_preserved'])
                    self.assertEqual(shard['original_copy_check']['records'], 4)
                    self.assertEqual(shard['annotation_preservation_check']['records'], 5)
                    self.assertEqual(shard['counters']['protein_altering_allele_records'], 5)
                    self.assertIn('--old-rec-tag', shard['commands'][0])
                    self.assertIn('--local-csq', shard['commands'][1])
                    rows = {row['key']: row for row in read_jsonl(output / 'coding-alleles.jsonl.gz')}
                    self.assertEqual(set(rows), EXPECTED_KEYS)
                    first, second = rows['1:10:A:C'], rows['1:10:A:G']
                    for index, row in enumerate((first, second), 1):
                        source = row['source_call']
                        self.assertEqual(source['caller'], caller)
                        self.assertEqual(source['GT'], '1/2')
                        self.assertEqual(source['GT_values'], [1, 2])
                        self.assertEqual(source['source_all_alts'], ['C', 'G'])
                        self.assertEqual(source['source_alt_index'], index)
                        self.assertEqual(source['alt_copies'], 1)
                        self.assertEqual(source['AD'], [1, 10, 10])
                        self.assertEqual(source['PL'], [100, 100, 100, 100, 0, 100])
                        self.assertEqual(source['native_context']['source_ordinal'], 1)
                        self.assertIn('NC_ORIG=', row['normalized_record'])
                        self.assertIn('FICTIONAL', row['genes'])
                        self.assertEqual(row['frequency']['status'], 'unresolved')
                        self.assertIsNone(row['frequency']['max_AF_fraction'])
                    self.assertEqual(first['source_call']['source_record_id'],
                                     second['source_call']['source_record_id'])
                    source_ids[caller] = first['source_call']['source_record_id']
                    missing, zero = rows['1:30:A:G'], rows['1:40:A:T']
                    self.assertEqual(missing['source_call']['FILTER'], [])
                    self.assertIsNone(missing['source_call']['GQ'])
                    self.assertEqual(missing['review_axes']['filter_state'], 'unfiltered')
                    self.assertEqual(missing['review_axes']['GQ']['state'], 'missing')
                    self.assertEqual(missing['review_axes']['GQ']['FORMAT_presence'], 'absent')
                    self.assertTrue(missing['review_axes']['structural_model_eligible'])
                    self.assertEqual(zero['source_call']['GQ'], 0)
                    self.assertEqual(zero['review_axes']['GQ']['state'], 'explicit_zero')
                    self.assertTrue(zero['review_axes']['structural_model_eligible'])
                    observed = rows['1:20:A:C']['frequency']
                    self.assertEqual(observed['status'], 'observed')
                    self.assertEqual(observed['max_AF_fraction'], 0.0001)
                    self.assertEqual(observed['counts'][0]['an'], 10000)
                    models = read_jsonl(output / 'descriptive-models-AF0.001.jsonl.gz')
                    same_record = [model for model in models if set(model['allele_keys']) == {'1:10:A:C', '1:10:A:G'}]
                    self.assertEqual(len(same_record), 1)
                    self.assertEqual(same_record[0]['phase'], 'opposite_alleles_in_source_GT')
                    self.assertFalse(same_record[0]['phase_evidence']['parental_phase_established'])
                    self.assertTrue(any('1:30:A:G' in model['allele_keys'] for model in models))
                    self.assertTrue(any('1:40:A:T' in model['allele_keys'] for model in models))
                    self.assertEqual(pinned_bytes(paths), before)
            self.assertNotEqual(source_ids['gatk'], source_ids['bcftools'])

    def test_partial_snapshot_requires_explicit_permission_to_use_that_scope(self):
        with tempfile.TemporaryDirectory() as temporary:
            fixture = FictionalCodingSnapshot(Path(temporary), self.runtime)
            paths, before = fixture.source_paths(), pinned_bytes(fixture.source_paths())
            process, output = fixture.execute('partial-refused', allow_partial=False)
            self.assert_failed(process, output, 'incomplete autosomal snapshot requires explicit')
            self.assertEqual(pinned_bytes(paths), before)

    def test_raw_source_changes_and_unpinned_annotation_are_refused(self):
        with tempfile.TemporaryDirectory() as temporary:
            fixture = FictionalCodingSnapshot(Path(temporary), self.runtime)
            paths = fixture.source_paths()
            before = pinned_bytes(paths)
            process, output = fixture.execute('bad-gff-pin', expected_gff='0' * 64)
            self.assert_failed(process, output, 'input checksum mismatch')
            self.assertIn('fictional-coding.gff3', json.loads((output / 'status.json').read_text())['error'])
            self.assertEqual(pinned_bytes(paths), before)

            with self.subTest(rejection='forged normalized count'):
                original_handoff = fixture.handoff.read_bytes()
                forged = json.loads(original_handoff)
                forged['shards'][0]['methods']['gatk']['normalized_split']['counts']['records'] += 1
                write_json(fixture.handoff, forged)
                forged_before = pinned_bytes(paths)
                # execute() supplies this altered summary's exact checksum; its
                # inconsistent record count must still fail native reconciliation.
                process, output = fixture.execute('forged-normalized-count')
                self.assert_failed(process, output, 'normalized record count differs from completed handoff')
                self.assertEqual(pinned_bytes(paths), forged_before)
                fixture.handoff.write_bytes(original_handoff)

            with self.subTest(rejection='changed reference index'):
                reference_index = Path(str(fixture.fixture.reference) + '.fai')
                original_index = reference_index.read_bytes()
                reference_index.write_bytes(original_index + b'\n')
                changed_index_before = pinned_bytes(paths)
                process, output = fixture.execute('changed-reference-index')
                self.assert_failed(process, output, 'input checksum mismatch')
                self.assertIn('reference.fa.fai', json.loads((output / 'status.json').read_text())['error'])
                self.assertEqual(pinned_bytes(paths), changed_index_before)
                reference_index.write_bytes(original_index)

            original_hash = sha256(fixture.fixture.raw['gatk'])
            changed = list(fixture.fixture.rows)
            changed[0] = changed[0].replace('1/2:21:', '0/1:21:')
            fixture.fixture.rewrite_gatk_raw(changed)
            self.assertNotEqual(sha256(fixture.fixture.raw['gatk']), original_hash)
            changed_before = pinned_bytes(paths)
            process, output = fixture.execute('changed-raw-refused')
            self.assert_failed(process, output, 'input checksum mismatch')
            self.assertIn('gatk.raw.vcf.gz', json.loads((output / 'status.json').read_text())['error'])
            self.assertEqual(pinned_bytes(paths), changed_before)

    def test_empty_native_snapshot_has_zero_coding_rows_with_explicit_pending_scope(self):
        with tempfile.TemporaryDirectory() as temporary:
            fixture = FictionalCodingSnapshot(Path(temporary), self.runtime, empty=True)
            paths = fixture.source_paths()
            before = pinned_bytes(paths)
            self.assertEqual(fixture.frequencies.read_bytes(), b'')
            for caller in ('gatk', 'bcftools'):
                with self.subTest(caller=caller):
                    process, output = fixture.execute('empty-' + caller, caller=caller)
                    summary = self.assert_complete(process, output)
                    self.assertEqual(summary['coding_allele_rows'], 0)
                    self.assertEqual(summary['source_records_with_coding_rows'], 0)
                    self.assertEqual(summary['frequency_responses'], 0)
                    self.assertEqual(summary['pending_autosomes'], AUTOSOMES[1:])
                    self.assertFalse(summary['all_autosomes_complete'])
                    self.assertEqual(read_jsonl(output / 'coding-alleles.jsonl.gz'), [])
                    self.assertEqual(read_jsonl(output / 'descriptive-models-AF0.001.jsonl.gz'), [])
                    self.assertTrue(summary['shards'][0]['original_copy_check']['derived_copy_payload_preserved'])
                    self.assertEqual(pinned_bytes(paths), before)

    def test_repeat_shifted_deletion_matches_sequence_but_wrong_origin_alt_is_rejected(self):
        sys.path.insert(0, str(WORKFLOW / 'scripts'))
        from native_coding_records import allele_call, original_alt_index, source_context
        from prepare_native_coding_review import verify_allele_haplotype

        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary).resolve()
            reference = folder / 'repeat.fa'
            reference.write_text('>1\nCAAAAAGTC\n')
            pysam.faidx(str(reference))
            source = folder / 'repeat.vcf'
            source.write_text(
                '##fileformat=VCFv4.2\n'
                '##contig=<ID=1,length=9>\n'
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTH\n'
                '1\t4\t.\tAA\tA,G\t100\tPASS\t.\tGT\t1/2\n')
            normalized = folder / 'repeat.normalized.vcf'
            result = subprocess.run([
                str(self.runtime / 'tools/bio/bin/bcftools'), 'norm', '--check-ref', 'e',
                '-f', str(reference), '-m', '-any', '--old-rec-tag', 'NC_ORIG',
                '-Ov', '-o', str(normalized), str(source)], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            with pysam.VariantFile(str(source)) as reader:
                context = source_context(next(reader), 'SYNTH', 'gatk', sha256(source), 1)
            with pysam.VariantFile(str(normalized)) as reader:
                deletions = [record.copy() for record in reader
                             if record.pos == 1 and record.ref == 'CA' and record.alts == ('C',)]
            self.assertEqual(len(deletions), 1)
            deletion = deletions[0]
            self.assertEqual(original_alt_index(deletion, context), 1)
            paths = [reference, Path(str(reference) + '.fai'), source, normalized]
            before = pinned_bytes(paths)
            with pysam.FastaFile(str(reference)) as fasta:
                verify_allele_haplotype(fasta, deletion, allele_call(context, 1))
                # An in-range origin index can be syntactically valid but wrong.
                # Force it to the replacement ALT G, which is not this deletion.
                deletion.info['NC_ORIG'] = '1|4|AA|A,G|2'
                wrong_index = original_alt_index(deletion, context)
                self.assertEqual(wrong_index, 2)
                with self.assertRaisesRegex(ValueError, 'not sequence-equivalent'):
                    verify_allele_haplotype(fasta, deletion, allele_call(context, wrong_index))
            self.assertEqual(pinned_bytes(paths), before)


if __name__ == '__main__':
    unittest.main()
