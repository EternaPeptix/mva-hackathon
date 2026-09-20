import importlib.util
import json
from pathlib import Path
import shutil
import tempfile
from types import SimpleNamespace
import unittest

import pysam

SCRIPT = Path(__file__).resolve().parents[1] / 'scripts/screen_secondary_findings.py'
spec = importlib.util.spec_from_file_location('secondary', SCRIPT)
secondary = importlib.util.module_from_spec(spec)
spec.loader.exec_module(secondary)


class SecondaryFindingsTests(unittest.TestCase):
    def test_review_strata_and_non_germline_assertions(self):
        for status, expected in [
            ('criteria_provided,_multiple_submitters,_no_conflicts', 'aggregate_plp_2plus_stars'),
            ('criteria_provided,_single_submitter', 'aggregate_plp_1star'),
            ('no_assertion_criteria_provided', 'aggregate_plp_0star_review')]:
            self.assertEqual(secondary.germline_stratum(
                {'CLNSIG': ('Pathogenic/Likely_pathogenic',), 'CLNREVSTAT': status})[0], expected)
        self.assertIsNone(secondary.germline_stratum(
            {'CLNSIG': 'Uncertain_significance', 'ONC': 'Oncogenic',
             'ONCREVSTAT': 'reviewed_by_expert_panel'}))
        self.assertIsNone(secondary.germline_stratum({'CLNSIG': 'Likely_benign'}))

    def test_conflicting_and_haplotype_assertions_never_become_clean_plp(self):
        info = {'CLNSIG': 'Conflicting_classifications_of_pathogenicity',
                'CLNSIGCONF': 'Pathogenic(3)|Uncertain_significance(2)',
                'CLNREVSTAT': 'criteria_provided,_conflicting_classifications'}
        self.assertEqual(secondary.germline_stratum(info)[0], 'conflicting_or_mixed_classification_review')
        for term in ('Likely_pathogenic_(1)|Uncertain_significance_(1)',
                     'Benign_(2)|Pathogenic_(1)|Likely_pathogenic_(1)',
                     'Likely_benign_(2),_Pathogenic_(1)'):
            info['CLNSIGCONF'] = term
            self.assertEqual(secondary.germline_stratum(info)[0], 'conflicting_or_mixed_classification_review')
        self.assertFalse(secondary.has_plp_signal('Conflicting_classifications_of_pathogenicity|Likely_benign_(1)'))
        info = {'CLNSIG': 'Pathogenic|risk_factor', 'CLNREVSTAT': 'reviewed_by_expert_panel'}
        self.assertEqual(secondary.germline_stratum(info)[0], 'conflicting_or_mixed_classification_review')
        info = {'CLNSIGINCL': '123:Pathogenic', 'CLNREVSTAT': 'no_classification_for_the_individual_variant'}
        self.assertEqual(secondary.germline_stratum(info)[0], 'haplotype_or_genotype_only_review')
        info['CLNREVSTAT'] = 'no_classification_for_the_single_variant'
        self.assertEqual(secondary.germline_stratum(info)[0], 'haplotype_or_genotype_only_review')
        with self.assertRaisesRegex(ValueError, 'unrecognized'):
            secondary.germline_stratum({'CLNSIG': 'Pathogenic', 'CLNREVSTAT': 'new_unknown_status'})

    def test_quality_missingness_and_homref_are_not_positive_calls(self):
        call = {'alt_copies': 1, 'fully_called': True, 'FILTER': ['PASS'], 'DP': 20, 'GQ': 99}
        self.assertEqual(secondary.quality_failures(call), [])
        for key, value, expected in [('alt_copies', 0, 'no_called_alternate'),
            ('fully_called', False, 'partial_or_missing_genotype'),
            ('FILTER', [], 'record_not_explicit_PASS'),
            ('GQ', None, 'genotype_quality_missing_or_below_threshold'),
            ('DP', 9, 'depth_missing_or_below_threshold')]:
            changed = dict(call, **{key: value})
            self.assertIn(expected, secondary.quality_failures(changed))

    def test_panel_count_guard_and_exact_gene_names(self):
        self.assertEqual(secondary.gene_symbols({'GENEINFO': 'GENE1:10|GENE10:100'}), {'GENE1', 'GENE10'})
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'panel.json'
            path.write_text(json.dumps({'total': 2, 'rows': [{'symbol': 'A'}, {'symbol': 'A'}]}))
            with self.assertRaisesRegex(ValueError, 'membership'):
                secondary.load_panel(path, 2)

    def test_real_normalization_preserves_multiallelic_source_and_phase(self):
        bcftools = shutil.which('bcftools')
        if not bcftools:
            self.skipTest('requires bcftools on PATH for the integration fixture')
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reference = root / 'reference.fa'
            reference.write_text('>1\nCAAAAAGTACGTACGTACGT\n')
            pysam.faidx(str(reference))
            vcf = root / 'sample.vcf'
            vcf.write_text('''##fileformat=VCFv4.2
##contig=<ID=1,length=19>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Quality">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depth">
##FORMAT=<ID=PGT,Number=1,Type=String,Description="Physical phase">
##FORMAT=<ID=PID,Number=1,Type=String,Description="Physical block">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S
1	4	.	AA	A	60	PASS	.	GT:DP:GQ:AD:PGT:PID	0/1:30:99:15,15:0|1:block
1	7	.	G	A,T	60	PASS	.	GT:DP:GQ:AD:PGT:PID	1/2:30:99:0,12,18:1|2:block
''')
            compressed = str(vcf) + '.gz'
            pysam.tabix_compress(str(vcf), compressed, force=True)
            pysam.tabix_index(compressed, preset='vcf', force=True)
            calls = []
            with pysam.VariantFile(compressed) as original:
                with secondary.normalized_stream(bcftools, reference, compressed, root / 'norm.log') as (stream, commands):
                    for record in stream:
                        calls.append((secondary.allele_key(record), secondary.original_call(record, original, 'S')))
            self.assertEqual(calls[0][0], ('1', 1, 'CA', 'C'))
            self.assertEqual(calls[0][1]['source_allele']['pos'], 4)
            self.assertEqual(len(calls), 3)
            self.assertEqual([c[1]['source_alt_index'] for c in calls[1:]], [1, 2])
            self.assertEqual([c[1]['alt_depth'] for c in calls[1:]], [12, 18])
            for _, call in calls[1:]:
                self.assertEqual(call['GT'], '1/2')
                self.assertEqual(call['PGT'], '1|2')
                self.assertEqual(call['PID'], 'block')
                self.assertEqual(call['alt_copies'], 1)
                self.assertEqual(secondary.quality_failures(call), [])

    def test_complete_screen_with_missing_public_contig_headers(self):
        bcftools = shutil.which('bcftools')
        if not bcftools:
            self.skipTest('requires bcftools on PATH for the integration fixture')
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reference = root / 'reference.fa'
            reference.write_text(''.join('>' + c + '\nCAAAAAGTACGTACGTACGT\n' for c in secondary.PRIMARY))
            pysam.faidx(str(reference))
            panel = root / 'panel.json'
            panel.write_text(json.dumps({'total': 2, 'rows': [{'symbol': 'G1'}, {'symbol': 'G2'}]}))
            public = root / 'clinvar.vcf'
            public.write_text('''##fileformat=VCFv4.2
##reference=GRCh38
##fileDate=synthetic
##INFO=<ID=GENEINFO,Number=1,Type=String,Description="Genes">
##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Classification">
##INFO=<ID=CLNSIGCONF,Number=.,Type=String,Description="Conflict">
##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description="Review">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	1	1	C	A	.	.	GENEINFO=G1:1;CLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts
1	7	2	G	T	.	.	GENEINFO=G1:1;CLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_single_submitter
1	8	3	T	G	.	.	GENEINFO=G2:2;CLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_expert_panel
1	9	4	A	C	.	.	GENEINFO=G2:2;CLNSIG=Conflicting_classifications_of_pathogenicity;CLNSIGCONF=Pathogenic(1)|Uncertain_significance(2);CLNREVSTAT=criteria_provided,_conflicting_classifications
''')
            subject = root / 'subject.vcf'
            subject.write_text('##fileformat=VCFv4.2\n' + ''.join(
                '##contig=<ID=' + c + ',length=19>\n' for c in secondary.PRIMARY) + '''##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="GQ">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S
1	1	.	C	A	60	PASS	.	GT:DP:GQ	0/1:30:99
1	7	.	G	T	60	PASS	.	GT:DP:GQ	0/1:5:99
1	8	.	T	G	60	PASS	.	GT:DP:GQ	0/0:30:99
1	9	.	A	C	60	PASS	.	GT:DP:GQ	0/1:30:99
''')
            for path in (public, subject):
                pysam.tabix_compress(str(path), str(path) + '.gz', force=True)
                pysam.tabix_index(str(path) + '.gz', preset='vcf', force=True)
            secondary.run(SimpleNamespace(output=root / 'out', sample='S', min_dp=10, min_gq=20,
                panel=panel, expected_genes=2, vcf=Path(str(subject) + '.gz'),
                clinvar=Path(str(public) + '.gz'), reference=reference, bcftools=Path(bcftools)))
            result = json.loads((root / 'out/triage.json').read_text())
            self.assertEqual(len(result['matches']), 3)
            self.assertEqual(result['sample_counts']['matching_records_no_called_alternate'], 1)
            self.assertEqual(result['quality_pass_unique_alleles_by_stratum'],
                {'aggregate_plp_2plus_stars': 1, 'conflicting_or_mixed_classification_review': 1})
            self.assertTrue(all(h['reportability'].startswith('pending_') for h in result['matches']))


if __name__ == '__main__':
    unittest.main()
