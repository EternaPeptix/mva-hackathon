#!/usr/bin/env python3
"""Construction/decision controls, plus fictional CRAM3.0 roundtrip; no GATK or case data."""
import collections,copy,importlib.util,json,os,tempfile,unittest
from fractions import Fraction
from pathlib import Path
import fictional_controls as f
s=importlib.util.spec_from_file_location('pilot',Path(__file__).with_name('run-fictional-pilot.py'));r=importlib.util.module_from_spec(s);s.loader.exec_module(r)
class Controls(unittest.TestCase):
    def test_genotype_joint_grid_marginals(self):
        pairs=f.genotype_grid(20260920);self.assertEqual(len(pairs),10000)
        for side in [0,1]:self.assertEqual(collections.Counter(p[side] for p in pairs),{0:8100,1:1800,2:100})
        self.assertEqual(collections.Counter(pairs)[(2,2)],1);self.assertEqual(collections.Counter(pairs)[(0,0)],6561)
    def test_seed_only_permutes_order(self):
        a=f.genotype_grid(20260920);b=f.genotype_grid(20260921);self.assertNotEqual(a,b);self.assertEqual(collections.Counter(a),collections.Counter(b))
    def test_known_mixture_exact_integer_composition(self):
        spec=next(s for s in f.SCENARIOS if s['name']=='mixture_0.02_d100');rows,truth=f.table_rows(spec)
        self.assertEqual(len(rows),10000);self.assertTrue(all(not t['rounded'] for t in truth))
        for row,t in zip(rows,truth):self.assertEqual(row['alt_count'],49*t['recipient_ALT_copies']+t['distinct_donor_ALT_copies']);self.assertEqual(row['ref_count']+row['alt_count'],100)
    def test_no_mixture_and_same_donor_identical(self):
        clean,_=f.table_rows(f.SCENARIOS[0]);same,_=f.table_rows(next(s for s in f.SCENARIOS if s['kind']=='same_donor'));self.assertEqual(clean,same)
    def test_cn_depth_and_allelic_imbalance_separate(self):
        rows,truth=f.table_rows(next(s for s in f.SCENARIOS if s['kind']=='cn_region'))
        self.assertTrue(all(r['ref_count']+r['alt_count']==450 for r in rows[:5000]));self.assertTrue(all(r['ref_count']+r['alt_count']==300 for r in rows[5000:]));self.assertTrue(all(t['copy_number_pattern']==3 for t in truth[:5000]))
    def test_empty_uninformative_and_rounding(self):
        self.assertEqual(f.table_rows(next(s for s in f.SCENARIOS if s['kind']=='empty')),([],[]))
        rows,_=f.table_rows(next(s for s in f.SCENARIOS if s['kind']=='uninformative'));self.assertTrue(all(r['alt_count']==50 for r in rows))
        self.assertEqual(f.round_half_up(Fraction(1,2)),1);self.assertEqual(f.round_half_up(Fraction(1,3)),0)
    def test_site_REF_orientation_and_bounds_refusal(self):
        f.validate_site('ACGT',2,'C','G')
        for vals in [(0,'A','C'),(2,'G','C'),(2,'C','C'),(True,'A','C')]:
            with self.assertRaises(ValueError):f.validate_site('ACGT',*vals)
    def result(self,d,c,e):
        p=Path(d)/'result.tsv';p.write_text('sample\tcontamination\terror\nFICTIONAL\t'+str(c)+'\t'+str(e)+'\n');return p
    def test_numeric_tolerances_fixed_and_fallback_never_clean(self):
        with tempfile.TemporaryDirectory() as d:
            spec=dict(kind='diploid',alpha='0.1')
            self.assertEqual(r.assess_estimator(self.result(d,.1,.002),spec)['classification'],'component_numeric_pass')
            for c,e in [(.14,.002),(0,1),(float('nan'),.1)]:
                with self.assertRaises(ValueError):r.assess_estimator(self.result(d,c,e),spec)
    def test_uninformative_and_CN_never_numeric_pass(self):
        with tempfile.TemporaryDirectory() as d:
            for kind,classification in [('uninformative','ineligible_low_information'),('cn_region','descriptive_confounding_control')]:
                got=r.assess_estimator(self.result(d,0,1),dict(kind=kind,alpha='0'));self.assertEqual(got['classification'],classification);self.assertIsNone(got['numerical_pass'])
    def test_exact_count_membership_and_zero_omission(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'counts.tsv';p.write_text('contig\tposition\tref_count\talt_count\tother_alt_count\tallele_frequency\nfictional1\t100\t1\t2\t0\t0.1\n')
            expected=[dict(contig='fictional1',position=100,ref_count=1,alt_count=2,other_alt_count=0,allele_frequency='0.1')];r.compare_counts(p,expected)
            for patch in [dict(position=101),dict(ref_count=0)]:
                with self.assertRaises(ValueError):r.compare_counts(p,[expected[0]|patch])
            with self.assertRaises(ValueError):r.compare_counts(p,expected+[expected[0]])
    def test_native_fictional_CRAM30_full_original_payload_roundtrip(self):
        try:import pysam
        except ImportError:self.skipTest('native pysam required')
        if os.environ.get('MVA_PILOT_CONSTRUCTOR_OUTPUT'):out=Path(os.environ['MVA_PILOT_CONSTRUCTOR_OUTPUT'])
        else:
            temporary=tempfile.TemporaryDirectory();self.addCleanup(temporary.cleanup);out=Path(temporary.name)/'fixture'
        fixture=f.count_fixture(out)
        with Path(fixture['cram']).open('rb') as stream:self.assertEqual(stream.read(6),b'CRAM\x03\x00')
        self.assertEqual(fixture['cram_version'],'3.0');self.assertTrue(fixture['decoded_payload_equal'])
        payload=json.loads((out/'decoded-payloads.json').read_text());original=[x['decoded_record'] for x in json.loads((out/'read-construction.json').read_text())]
        canonical=lambda rows:collections.Counter(json.dumps(x,sort_keys=True) for x in rows)
        self.assertEqual(canonical(original),canonical(payload['bam']));self.assertEqual(canonical(original),canonical(payload['cram']))
        self.assertEqual(fixture['records'],len(original));self.assertEqual(fixture['records'],len(payload['bam']))
        missing=next(x for x in original if x['QNAME']=='FICTIONAL_missing_BQ');partial=next(x for x in original if x['QNAME']=='FICTIONAL_partial_BQ255')
        self.assertIsNone(missing['qualities']);self.assertEqual(partial['qualities'][74],255)
        unmapped=next(x for x in original if x['FLAG']&4);self.assertEqual(unmapped['MQ'],0);self.assertIsNone(unmapped['CIGAR'])
        self.assertTrue(any(any(c[0]==3 for c in (x['CIGAR'] or [])) for x in original))
if __name__=='__main__':unittest.main(verbosity=2)
