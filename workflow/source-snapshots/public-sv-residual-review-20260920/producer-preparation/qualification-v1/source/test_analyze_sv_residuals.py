"""Small wholly fictional saved-VCF controls; no real source path is opened."""
import copy
import gzip
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

SPEC = importlib.util.spec_from_file_location('residual', Path(__file__).with_name('analyze_sv_residuals.py'))
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)


def write_json(path, obj):
    path.write_text(json.dumps(obj, sort_keys=True) + '\n')


class Fixture:
    def __init__(self, root):
        self.root = Path(root)
        self.paths = {name: self.root/(name + ('.vcf.gz' if 'truth' in name else '.json')) for name in M.ROLES}
        self.paths['official_chr20_benchmark'] = self.root/'benchmark.bed'
        self.header = ['##fileformat=VCFv4.3', '##contig=<ID=chr20,length=10000>',
                       '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                       '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
                       '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SVTYPE">',
                       '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SVLEN">']
        for name, kind in M.ANNOTATIONS.items():
            self.header.append(f'##INFO=<ID={name},Number={0 if kind == "Flag" else 1},Type={kind},Description="Fictional annotation">')
        self.columns = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002'
        def row(pos, ref, alt, info, gt='0/1'):
            return '\t'.join(['chr20', str(pos), '.', ref, alt, '.', 'PASS', info, 'GT:AD', gt+':.'])
        self.original = [row(20, 'A', 'T', '.', '0/1'),
                         row(100, 'A', 'A'+'C'*50, 'SVTYPE=INS;SVLEN=50;TRF;TRFperiod=1'),
                         row(100, 'A', 'A'+'G'*50, 'SVTYPE=INS;SVLEN=50;RM_clsfam=LINE/L1;LCR=0;REMAP=partial'),
                         row(2000, 'A'*51, 'A', 'SVTYPE=DEL;SVLEN=50;TRF;RM_clsfam=Simple_repeat;LCR=0.5'),
                         row(4000, 'A'*101, 'A', 'SVTYPE=DEL;SVLEN=100;LCR=.;RM_score=.')]
        def final(values):
            result = []
            for raw, bd in zip(self.original[1:], values):
                f = raw.split('\t');f[8] += ':BD';f[9] += ':'+bd;result.append('\t'.join(f))
            return result
        self.delly = final(['TP', 'TP', 'FN', 'FN'])
        self.svaba = final(['TP', 'FN', 'TP', 'FN'])
        self.expected = {'truth_records': 4, 'original_full_context_records': 5,
                         'benchmark_intervals': 1, 'benchmark_bases': 10000,
                         'classes': {'DEL': 2, 'INS': 2},
                         'outcomes': dict.fromkeys(M.OUTCOMES, 1),
                         'outcomes_by_class': {
                             'INS': dict(zip(M.OUTCOMES, [1, 1, 0, 0])),
                             'DEL': dict(zip(M.OUTCOMES, [0, 0, 1, 1]))},
                         'per_caller_by_size': {
                             'DELLY': {'INS_like:50-99:TP': 2, 'DEL:50-99:FN': 1, 'DEL:100-499:FN': 1},
                             'SvABA': {'INS_like:50-99:TP': 1, 'INS_like:50-99:FN': 1, 'DEL:50-99:TP': 1, 'DEL:100-499:FN': 1}}}
        self.paths['official_chr20_benchmark'].write_text('chr20\t0\t10000\n')
        self.write_vcfs()
        self.make_receipts()

    def write_vcf(self, role, rows, final=False):
        header = self.header + (['##FORMAT=<ID=BD,Number=1,Type=String,Description="Decision">'] if final else [])
        with gzip.open(self.paths[role], 'wt') as writer:
            writer.write('\n'.join(header+[self.columns]+rows)+'\n')

    def write_vcfs(self):
        self.write_vcf('original_truth_chr20_no_star', self.original)
        self.write_vcf('delly_final_truth', self.delly, True)
        self.write_vcf('svaba_final_truth', self.svaba, True)

    def make_receipts(self):
        write_json(self.paths['comparison_config'], {'truvari_version': '5.5.0'})
        cp = M.fingerprint(self.paths['comparison_config'])
        for caller in ['delly', 'svaba']:
            write_json(self.paths[caller+'_status'], {'status': 'complete'})
            write_json(self.paths[caller+'_variant_summary'], {'TP-base': 2, 'FN': 2})
        ds = {'comparison': {'refined': {'TP-base': 2, 'FN': 2}}, 'config': cp,
              'inputs': {'baseline': M.fingerprint(self.paths['original_truth_chr20_no_star']),
                         'benchmark_bed': M.fingerprint(self.paths['official_chr20_benchmark'])}}
        ss = {'comparison': {'refined': {'TP-base': 2, 'FN': 2}}, 'configuration': cp,
              'by_size': {'baseline': self.expected['per_caller_by_size']['SvABA']}}
        write_json(self.paths['delly_summary'], ds);write_json(self.paths['svaba_summary'], ss)
        audit = {key: True for key in ['truth_alleles_equal_between_callers', 'truth_alleles_unique',
                 'truth_genotypes_unchanged_between_callers', 'every_final_truth_allele_occurs_exactly_in_original_input',
                 'counts_reconcile', 'overlap_totals_reconcile']}
        audit.update(comparison_status='complete', artifacts={
            'benchmark/refine.base.vcf.gz': M.fingerprint(self.paths['svaba_final_truth']),
            'summary.json': M.fingerprint(self.paths['svaba_summary'])})
        write_json(self.paths['prior_final_audit'], audit)

    def config(self):
        return {'ready': True, 'scope': 'fictional only', 'expected': self.expected,
                'inputs': {k: M.fingerprint(v) for k, v in self.paths.items()}}

    def build(self):
        return M.build_ledger(self.paths, self.expected)


class ResidualControls(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.f = Fixture(self.root)

    def tearDown(self):
        self.temp.cleanup()

    def test_four_outcomes_and_original_ordinals(self):
        rows, summary, headers = self.f.build()
        self.assertEqual([r['original_ordinal'] for r in rows], [2, 3, 4, 5])
        self.assertEqual([r['outcome'] for r in rows], list(M.OUTCOMES))
        self.assertEqual(summary['total']['outcomes'], dict.fromkeys(M.OUTCOMES, 1))
        self.assertEqual(len(summary['class_size']), 10)
        self.assertEqual(len(summary['TRF_RM_by_class_size']), 60)
        self.assertEqual(len(summary['raw_LCR_by_class_size_and_outcome']), 40)
        self.assertEqual(set(headers), {'original', 'DELLY', 'SvABA'})

    def test_same_position_same_size_distinct_alt_not_collapsed(self):
        rows, _, _ = self.f.build()
        self.assertEqual(rows[0]['pos_1based'], rows[1]['pos_1based'])
        self.assertNotEqual(rows[0]['allele_sha256'], rows[1]['allele_sha256'])
        self.assertEqual(rows[0]['outcome'], 'matched_by_both')
        self.assertEqual(rows[1]['outcome'], 'DELLY_only')

    def test_absence_zero_and_explicit_missing_are_distinct(self):
        rows, summary, _ = self.f.build()
        self.assertEqual([M.value_state(r['annotations']['LCR']) for r in rows],
                         [('not_encoded', None), ('value', 0.0), ('value', 0.5), ('explicit_missing', None)])
        self.assertEqual(M.repeat_state(rows[1]), ('flag_not_encoded', 'value_present'))
        self.assertEqual(M.repeat_state(rows[3]), ('flag_not_encoded', 'fields_explicit_missing'))
        for table in ('class_size', 'TRF_RM_by_class_size', 'RM_class_family_by_class_size', 'REMAP_by_class_size'):
            self.assertEqual(sum(r['truth_records'] for r in summary[table]), 4)

    def test_empty_stratum_fraction_is_null(self):
        _, summary, _ = self.f.build()
        row = next(r for r in summary['class_size'] if r['truth_records'] == 0)
        self.assertIsNone(row['shared_miss_fraction'])
        self.assertEqual(set(row['outcomes'].values()), {0})

    def test_genotype_order_and_phase_encoding_retained_without_phase_claim(self):
        self.f.delly[0] = self.f.delly[0].replace('0/1:', '1|0:')
        self.f.write_vcfs()
        rows, _, _ = self.f.build()
        self.assertEqual(rows[0]['callers']['DELLY']['GT'], '1|0')
        self.assertEqual(rows[0]['original_GT'], '0/1')
        self.assertFalse(rows[0]['GT_encoding_matches_original']['DELLY'])

    def test_changed_genotype_allele_content_refused(self):
        self.f.delly[0] = self.f.delly[0].replace('0/1:', '1/1:');self.f.write_vcfs()
        with self.assertRaisesRegex(ValueError, 'genotype differs'):self.f.build()

    def test_duplicate_final_identity_refused(self):
        self.f.delly.append(self.f.delly[0]);self.f.write_vcfs()
        with self.assertRaisesRegex(ValueError, 'duplicate final'):self.f.build()

    def test_duplicate_original_identity_refused(self):
        self.f.original.append(self.f.original[1]);self.f.write_vcfs()
        with self.assertRaisesRegex(ValueError, 'duplicate original'):self.f.build()

    def test_missing_final_allele_refused(self):
        self.f.svaba.pop();self.f.write_vcfs()
        with self.assertRaisesRegex(ValueError, 'universes differ'):self.f.build()

    def test_missing_original_allele_refused(self):
        self.f.original.pop();self.f.expected['original_full_context_records'] = 4;self.f.write_vcfs()
        with self.assertRaisesRegex(ValueError, 'original truth missing'):self.f.build()

    def test_changed_same_length_allele_refused(self):
        self.f.svaba[0] = self.f.svaba[0].replace('A'+'C'*50, 'A'+'T'*50);self.f.write_vcfs()
        with self.assertRaisesRegex(ValueError, 'universes differ'):self.f.build()

    def test_changed_annotation_refused(self):
        self.f.svaba[1] = self.f.svaba[1].replace('LCR=0;', 'LCR=1;');self.f.write_vcfs()
        with self.assertRaisesRegex(ValueError, 'annotation value'):self.f.build()

    def test_missing_flag_not_converted_to_false_value(self):
        rows, _, _ = self.f.build()
        self.assertEqual(rows[0]['annotations']['TRF'], {'present': True, 'value': True, 'raw': None})
        self.assertEqual(rows[1]['annotations']['TRF'], {'present': False, 'value': None, 'raw': None})

    def test_flag_equals_zero_refused(self):
        self.f.delly[0] = self.f.delly[0].replace(';TRF;', ';TRF=0;');self.f.write_vcfs()
        with self.assertRaisesRegex(ValueError, 'flag must'):self.f.build()

    def test_nan_numeric_annotation_refused(self):
        self.f.delly[1] = self.f.delly[1].replace('LCR=0;', 'LCR=NaN;');self.f.write_vcfs()
        with self.assertRaisesRegex(ValueError, 'nonfinite'):self.f.build()

    def test_duplicate_info_refused(self):
        self.f.delly[0] = self.f.delly[0].replace(';TRF;', ';TRF;TRF;');self.f.write_vcfs()
        with self.assertRaisesRegex(ValueError, 'duplicate INFO'):self.f.build()

    def test_annotation_header_type_change_refused(self):
        text = gzip.open(self.f.paths['svaba_final_truth'], 'rt').read().replace('ID=LCR,Number=1,Type=Float', 'ID=LCR,Number=1,Type=String')
        with gzip.open(self.f.paths['svaba_final_truth'], 'wt') as out:out.write(text)
        with self.assertRaisesRegex(ValueError, 'annotation schema'):self.f.build()

    def test_unexpected_chromosome_refused(self):
        self.f.delly[0] = self.f.delly[0].replace('chr20\t', 'chr21\t');self.f.write_vcfs()
        with self.assertRaisesRegex(ValueError, 'unexpected chromosome'):self.f.build()

    def test_label_change_not_reclassified(self):
        self.f.delly[0] = self.f.delly[0][:-2]+'FP';self.f.write_vcfs()
        with self.assertRaisesRegex(ValueError, 'unsupported final'):self.f.build()

    def test_published_delly_size_margin_enforced(self):
        self.f.expected['per_caller_by_size']['DELLY']['INS_like:50-99:TP'] = 1
        with self.assertRaisesRegex(ValueError, 'class/size margins'):self.f.build()

    def test_expected_bool_is_not_count(self):
        self.f.expected['outcomes']['matched_by_both'] = True
        with self.assertRaisesRegex(ValueError, 'nonnegative integer'):self.f.build()

    def test_all_size_boundaries(self):
        values = [50, 99, 100, 499, 500, 999, 1000, 9999, 10000, 50000]
        self.assertEqual([M.size_bin(v) for v in values], [x for x in M.BINS for _ in range(2)])
        for value in [49, 50001, True]:
            with self.subTest(value=value), self.assertRaises(ValueError):M.size_bin(value)

    def test_bed_boundary_is_descriptive_not_new_exclusion(self):
        self.f.paths['official_chr20_benchmark'].write_text('chr20\t0\t100\n')
        self.f.expected['benchmark_bases'] = 100
        rows, summary, _ = self.f.build()
        self.assertEqual(len(rows), 4)
        self.assertFalse(rows[2]['REF_span_contained_in_one_benchmark_interval'])
        self.assertEqual(summary['total']['truth_records'], 4)

    def test_integration_output_and_no_input_mutation(self):
        cfg = self.root/'config.json';write_json(cfg, self.f.config())
        before = {k: M.fingerprint(p) for k, p in self.f.paths.items()}
        output = self.root/'output';status = M.run(cfg, M.fingerprint(cfg)['sha256'], output)
        self.assertEqual(status['status'], 'complete')
        self.assertEqual(status['truth_records'], 4)
        self.assertEqual(len((output/'ledger.jsonl').read_text().splitlines()), 4)
        self.assertEqual(before, {k: M.fingerprint(p) for k, p in self.f.paths.items()})
        self.assertEqual(set(status['outputs']), {'ledger.jsonl', 'summary.json', 'headers.json'})

    def test_disabled_config_refuses_before_input_access(self):
        config = self.f.config();config['ready'] = False;config['inputs']['delly_final_truth']['path'] = '/does/not/exist'
        cfg = self.root/'config.json';write_json(cfg, config)
        with self.assertRaisesRegex(ValueError, 'not approved'):M.run(cfg, M.fingerprint(cfg)['sha256'], self.root/'output')

    def test_wrong_input_pin_no_output(self):
        config = self.f.config();config['inputs']['delly_final_truth']['sha256'] = '0'*64
        cfg = self.root/'config.json';write_json(cfg, config)
        with self.assertRaisesRegex(ValueError, 'input SHA'):M.run(cfg, M.fingerprint(cfg)['sha256'], self.root/'output')
        self.assertFalse((self.root/'output').exists())

    def test_wrong_config_pin_and_existing_output_preserved(self):
        cfg = self.root/'config.json';write_json(cfg, self.f.config())
        with self.assertRaisesRegex(ValueError, 'config SHA'):M.run(cfg, '0'*64, self.root/'output')
        out = self.root/'output';out.mkdir();(out/'keep').write_text('keep')
        with self.assertRaisesRegex(ValueError, 'existing output'):M.run(cfg, M.fingerprint(cfg)['sha256'], out)
        self.assertEqual((out/'keep').read_text(), 'keep')

    def test_stale_historical_audit_pin_refused(self):
        path = self.f.paths['prior_final_audit'];value = json.loads(path.read_text())
        value['artifacts']['benchmark/refine.base.vcf.gz']['sha256'] = '0'*64;write_json(path, value)
        cfg = self.root/'config.json';write_json(cfg, self.f.config())
        with self.assertRaisesRegex(ValueError, 'historical audit artifact pin'):M.run(cfg, M.fingerprint(cfg)['sha256'], self.root/'output')
        self.assertFalse((self.root/'output').exists())

    def test_nonfinite_and_duplicate_json_refused(self):
        for raw in ['{"x":NaN}', '{"x":Infinity}', '{"x":1,"x":2}']:
            with self.subTest(raw=raw), self.assertRaises(ValueError):M.parse_json(raw)


if __name__ == '__main__':
    unittest.main()
