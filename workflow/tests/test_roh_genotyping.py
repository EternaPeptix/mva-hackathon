"""Native integration and provenance failures affecting genotype/ROH interpretation."""
import copy
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import pysam

sys.path.insert(0, str(Path(__file__).parents[1] / 'scripts'))
from genotype_roh_panel import (PATIENT_SCOPE, fingerprint, qualify_calls, run_analysis,
                                save, validate_config, verify_completion)
from probe_roh_genotyping import create_fixture, check_result
from segmented_roh import POPULATION_TAGS, run_segmented


class IntegratedROH(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        default = '/Users/jeweled/Documents/mva-hackathon/tools/bio/bin/'
        self.bcf = Path(os.environ.get('MVA_BCFTOOLS', default + 'bcftools'))
        self.sam = Path(os.environ.get('MVA_SAMTOOLS', default + 'samtools'))

    def tearDown(self):
        self.tmp.cleanup()

    def fixture(self, mode='controls'):
        return create_fixture(self.root / 'fixture', mode, self.bcf, self.sam)

    def run_fixture(self, mode='controls'):
        rows, alignment, reference, panel, config = self.fixture(mode)
        result = run_analysis(alignment, reference, panel, config, self.bcf, self.sam, self.root / 'analysis')
        check_result(rows, mode, result, self.root / 'analysis')
        return rows, alignment, reference, panel, config, result

    def test_native_controls_preserve_genotypes_and_source_population_fields(self):
        rows, _, _, panel, _, result = self.run_fixture()
        self.assertEqual(result['marker_qualification']['counts']['qualified_markers'], 3)
        self.assertEqual(result['read_extraction']['native_counts']['records_filter_accepted'], 183)
        self.assertEqual(result['read_extraction']['native_counts']['records_filter_rejected'], 180)
        with pysam.VariantFile(str(self.root / 'analysis/calls.bcf')) as original, \
                pysam.VariantFile(str(self.root / 'analysis/annotated-all-markers.vcf.gz')) as annotated, \
                pysam.VariantFile(str(panel / 'panel.vcf.gz')) as population:
            for raw, observed, source in zip(original, annotated, population):
                for key in ('GT', 'DP', 'AD', 'GQ', 'PL'):
                    self.assertEqual(raw.samples['SYNTH'].get(key), observed.samples['SYNTH'].get(key))
                for tag in POPULATION_TAGS:
                    self.assertEqual(observed.info[tag], source.info[tag])
        self.assertEqual(len(result['homozygosity_models']), 8)
        for model in result['homozygosity_models']:
            summary = json.loads(Path(model['summary']['path']).read_text())
            for segment in summary['segments']:
                args = segment['command']
                self.assertEqual(args[args.index('--AF-tag') + 1], model['population_af_tag'])

    def test_empty_evidence_is_missing_and_unanalyzed(self):
        _, _, _, _, _, result = self.run_fixture('empty')
        self.assertEqual(result['marker_qualification']['counts']['raw_genotype_missing'], 14)
        self.assertEqual(result['analysis_status'], 'no_qualified_markers_no_models')
        self.assertFalse(result['homozygosity_models'])

    def test_wrong_pins_and_wrong_sample_fail_before_read_extraction(self):
        _, alignment, reference, panel, config = self.fixture()
        for i, (field, value) in enumerate([('reference_sha256', '0' * 64),
                ('panel_summary_sha256', '0' * 64), ('bcftools_sha256', '0' * 64),
                ('samtools_sha256', '0' * 64), ('sample', 'WRONG')]):
            with self.subTest(field=field):
                bad = copy.deepcopy(config); bad[field] = value
                folder = self.root / ('failed-' + str(i))
                with self.assertRaises(ValueError):
                    run_analysis(alignment, reference, panel, bad, self.bcf, self.sam, folder)
                self.assertFalse((folder / 'marker-reads.bam').exists())
                self.assertEqual(json.loads((folder / 'status.json').read_text())['status'], 'failed')

    def test_panel_sidecars_require_exact_marker_identity_even_with_updated_hash(self):
        _, alignment, reference, panel, config = self.fixture()
        p = panel / 'sites.bed'; p.write_text(p.read_text().replace('1\t1000\t1001', '1\t1001\t1002'))
        summary = json.loads((panel / 'summary.json').read_text())
        summary['outputs']['sites.bed'] = fingerprint(p); save(panel / 'summary.json', summary)
        config['panel_summary_sha256'] = fingerprint(panel / 'summary.json')['sha256']
        with self.assertRaisesRegex(ValueError, 'representations disagree'):
            run_analysis(alignment, reference, panel, config, self.bcf, self.sam, self.root / 'bad')
        self.assertFalse((self.root / 'bad/marker-reads.bam').exists())

    def test_partial_or_mismatched_patient_completion_cannot_pass(self):
        alignment = self.root / 'final.cram'
        identity = {'bytes': 100, 'sha256': '1' * 64}
        config = {'scope': PATIENT_SCOPE}
        receipt = self.root / 'completion.json'
        with self.assertRaisesRegex(ValueError, 'receipt required'):
            verify_completion(config, alignment, identity, None)
        valid = {'status': 'complete', 'expected_lanes': 4, 'cram_path': str(alignment),
                 'cram_bytes': 100, 'cram_sha256': '1' * 64}
        save(receipt, valid)
        self.assertEqual(verify_completion(config, alignment, identity, receipt)['sha256'], fingerprint(receipt)['sha256'])
        for field, bad_value in [('status', 'running'), ('expected_lanes', 2), ('cram_bytes', 99),
                ('cram_sha256', '2' * 64), ('cram_path', str(self.root / 'retained-two-lane.cram'))]:
            with self.subTest(field=field):
                invalid = dict(valid); invalid[field] = bad_value; save(receipt, invalid)
                with self.assertRaisesRegex(ValueError, 'does not match'):
                    verify_completion(config, alignment, identity, receipt)
        prototype = json.loads(Path(__file__).parents[1].joinpath('configs/roh-genotyping.json').read_text())
        prototype['chromosomes'] = ['1', '2']
        with self.assertRaisesRegex(ValueError, 'all 22'):
            validate_config(prototype)

    def test_dropped_caller_marker_is_an_error_not_an_imputed_reference(self):
        rows, _, _, panel, config, _ = self.run_fixture()
        path = self.root / 'dropped.bcf'
        with pysam.VariantFile(str(self.root / 'analysis/calls.bcf')) as source, \
                pysam.VariantFile(str(path), 'wb', header=source.header) as writer:
            for i, row in enumerate(source):
                if i != 3: writer.write(row)
        with pysam.VariantFile(str(panel / 'panel.vcf.gz')) as source:
            marker_rows = list(source); header = source.header.copy()
        folder = self.root / 'qualification'; folder.mkdir()
        with self.assertRaisesRegex(ValueError, 'omitted, reordered or added'):
            qualify_calls(marker_rows, header, self.root / 'analysis/pileup.bcf', path, 'SYNTH', config['evidence'], folder)

    def test_multi_region_extraction_does_not_repeat_a_read_overlapping_two_markers(self):
        _, alignment, reference, _, _ = self.fixture()
        bed = self.root / 'adjacent.bed'; bed.write_text('1\t1000\t1001\n1\t1010\t1011\n')
        out = self.root / 'adjacent.bam'
        subprocess.run([str(self.sam), 'view', '-b', '-M', '-L', str(bed), '-T', str(reference),
                        '-q', '20', '-F', '0xF04', '-e', 'mapq != 255', '-o', str(out), str(alignment)], check=True)
        with pysam.AlignmentFile(str(out)) as reader:
            self.assertEqual(sum(1 for _ in reader), 30)

    def test_external_frequency_profile_changes_native_support_with_fixed_genotypes(self):
        path = self.root / 'likelihoods.vcf'
        header = pysam.VariantHeader(); header.contigs.add('1', length=20000000)
        for tag in ('POP_AF', 'POP_AFR_AF'): header.info.add(tag, 1, 'Float', 'Synthetic frequency sensitivity control')
        header.formats.add('GT', 1, 'String', 'Synthetic genotype'); header.formats.add('PL', 'G', 'Integer', 'Synthetic likelihoods')
        header.add_sample('SYNTH')
        with pysam.VariantFile(str(path), 'w', header=header) as writer:
            for pos in range(5001, 20000000, 10000):
                r = writer.new_record(contig='1', start=pos - 1, alleles=('A', 'G'))
                r.info['POP_AF'] = 0.5; r.info['POP_AFR_AF'] = 0.001
                r.samples['SYNTH']['GT'] = (0, 0); r.samples['SYNTH']['PL'] = (0, 90, 219); writer.write(r)
        primary = run_segmented(path, 'SYNTH', self.bcf, self.root / 'primary')
        alternate = run_segmented(path, 'SYNTH', self.bcf, self.root / 'alternate', af_tag='POP_AFR_AF')
        # For homozygous reference evidence, the AZ/HW emission ratio is
        # 1/(1-AF): AF=.5 supplies stronger support than AF=.001. Neither
        # frequency guarantees no interval in this all-homozygous fixture.
        a = primary['segments'][0]['regions']; b = alternate['segments'][0]['regions']
        self.assertEqual(len(a), 1); self.assertEqual(len(b), 1)
        self.assertEqual((a[0]['start_1based'], a[0]['end_inclusive']), (b[0]['start_1based'], b[0]['end_inclusive']))
        self.assertGreater(a[0]['mean_model_phred'], b[0]['mean_model_phred'])
        self.assertEqual(primary['counts'], alternate['counts'])


if __name__ == '__main__': unittest.main()
