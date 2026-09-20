"""Fictional adversarial audit controls. The producer fixture supplies inputs only.

The audit module never imports the producer. Integration fixture execution is
explicitly confined to the test and preserves generated files when a control
artifact directory is supplied.
"""
import copy
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
from types import ModuleType, SimpleNamespace
import unittest

try:
    import pysam
    NATIVE=hasattr(pysam,'VariantFile')
except ImportError:
    sys.modules['pysam']=ModuleType('pysam');NATIVE=False

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('independent_secondary_audit',HERE/'verify-secondary.py')
a=importlib.util.module_from_spec(spec);spec.loader.exec_module(a)
QUALIFICATION=os.environ.get('MVA_SECONDARY_QUALIFICATION')
NATIVE=NATIVE and bool(QUALIFICATION) and shutil.which('bcftools') is not None


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path,value):Path(path).write_text(json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n')


def json_rows(path):
    with gzip.open(path,'rt') as f:return [json.loads(line) for line in f]


def save_rows(path,rows):
    with gzip.open(path,'wt') as f:
        for row in rows:f.write(json.dumps(row,sort_keys=True,allow_nan=False)+'\n')


class PureControls(unittest.TestCase):
    def test_typed_equality(self):
        for v in [True,1.0]:
            with self.assertRaises(ValueError):a.equal(v,1,'typed')
    def test_nonfinite_and_duplicate_json(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'v.json'
            for value in ['{"v":NaN}','{"v":Infinity}','{"v":1,"v":2}']:
                p.write_text(value)
                with self.assertRaises(ValueError):a.read_json(p)
    def test_missing_GT_slot_is_not_ploidy(self):
        for f in [{},{'GT':[None]},{'GT':[None,None]}]:
            self.assertIsNone(a.shape({'FORMAT':f},2)['observed_GT_ploidy'])
        with self.assertRaises(ValueError):a.shape({'FORMAT':{'GT':[1]}},2)
    def test_quality_missing_absent_zero_distinct(self):
        got=[a.axis({'FORMAT':f},'GQ',20) for f in [{},{'GQ':None},{'GQ':0}]]
        self.assertEqual([v['FORMAT_presence'] for v in got],['absent','present_missing','present_value'])
        self.assertEqual([v['state'] for v in got],['missing','missing','explicit_zero'])
    def test_stratum_not_somatic_or_benign(self):
        self.assertIsNone(a.classification({'ONC':'Oncogenic','CLNSIG':'Benign'}))
        self.assertEqual(a.classification({'CLNSIG':'Pathogenic','CLNREVSTAT':'criteria_provided,_multiple_submitters,_no_conflicts'}),('aggregate_plp_2plus_stars',2))
        self.assertEqual(a.classification({'CLNSIG':['Pathogenic','Benign'],'CLNREVSTAT':'criteria_provided,_single_submitter'})[0],'conflicting_or_mixed_classification_review')
    def test_unknown_review_status_refused(self):
        with self.assertRaises(ValueError):a.classification({'CLNSIG':'Pathogenic','CLNREVSTAT':'invented'})
    def test_distinct_literal_edits_refused(self):
        fasta=SimpleNamespace(fetch=lambda c,s,e:'ACGT'[s:e])
        normalized=SimpleNamespace(contig='1',pos=2,ref='C',alts=('A',))
        with self.assertRaises(ValueError):a.literal(fasta,normalized,dict(chrom='1',pos=2,ref='C',alts=['G']),1)
    def test_stale_hash_refused(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'x';p.write_text('a');expected={'sha256':sha(p)};p.write_text('b')
            with self.assertRaises(ValueError):a.Files().hash(p,expected)
    def test_same_bytes_wrong_path_refused(self):
        with tempfile.TemporaryDirectory() as d:
            p,q=Path(d)/'x',Path(d)/'y';p.write_text('a');q.write_text('a')
            with self.assertRaises(ValueError):a.Files().hash(q,{'path':str(p),'sha256':sha(p)})
    def metadata_contract(self):
        return dict(expected_genes=84,public_snapshot=dict(clinvar_fileDate='2026-09-13',clinvar_snapshot_label='20260913',clinvar_resource_filename='clinvar_20260913.vcf.gz'),inputs={k:dict(sha256=v,path='/public/clinvar_20260913.vcf.gz' if k=='clinvar' else '/public/'+k) for k,v in a.FROZEN_PUBLIC.items()})
    def test_exact_header_date_distinct_from_snapshot_label(self):
        a.public_contract(self.metadata_contract())
    def test_compact_header_date_and_wrong_filename_refused(self):
        for key,value in [('clinvar_fileDate','20260913'),('clinvar_snapshot_label','2026-09-13'),('clinvar_resource_filename','clinvar_20260914.vcf.gz')]:
            config=self.metadata_contract();config['public_snapshot'][key]=value
            with self.subTest(key=key),self.assertRaises(ValueError):a.public_contract(config)
    def test_zero_counts_keep_callers(self):
        groups,counts=a.grouped([])
        self.assertEqual(groups,[])
        self.assertEqual(set(counts),{'gatk','bcftools'})
        self.assertTrue(all(v==0 for c in counts.values() for v in c.values()))


@unittest.skipUnless(NATIVE,'native pysam/bcftools and explicit qualified fictional fixture required')
class IntegrationControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.qualification=Path(QUALIFICATION).resolve()
        cls.qualification_sha=sha(cls.qualification)
        workflow=cls.qualification.parent/'workflow'
        fixture_script=workflow/'tests/test_native_secondary.py'
        spec=importlib.util.spec_from_file_location('secondary_fixture_input_only',fixture_script)
        fixture=importlib.util.module_from_spec(spec);spec.loader.exec_module(fixture);cls.fixture=fixture
        declared=os.environ.get('MVA_SECONDARY_CONTROL_OUTPUT')
        cls.temporary=None
        if declared:cls.root=Path(declared).resolve();cls.root.mkdir(parents=True,exist_ok=False)
        else:cls.temporary=tempfile.TemporaryDirectory();cls.root=Path(cls.temporary.name).resolve()
        cls.fixture_root=cls.root/'fictional-inputs';cls.fixture_root.mkdir()
        cls.config,cls.config_path=fixture.fixture(cls.fixture_root)
        cls.producer_run=cls.root/'producer-fixture'
        fixture.a.run(cls.config_path,sha(cls.config_path),cls.producer_run)
        cls.input_pins={str(p):sha(p) for p in cls.fixture_root.rglob('*') if p.is_file()}
    @classmethod
    def tearDownClass(cls):
        if cls.temporary:cls.temporary.cleanup()
    def setUp(self):
        self.folder=self.root/self._testMethodName;self.folder.mkdir()
        self.run=self.folder/'run';shutil.copytree(self.producer_run,self.run)
        self.repin()
    def repin(self):
        status=json.loads((self.run/'status.json').read_text())
        files=a.Files()
        status['outputs']=[files.hash(p) for p in sorted(self.run.rglob('*')) if p.is_file() and p.name!='status.json']
        status['summary_sha256']=sha(self.run/'summary.json');save(self.run/'status.json',status)
    def run_audit(self):
        return a.audit(self.config_path,sha(self.config_path),self.run,sha(self.run/'status.json'),self.qualification,self.qualification_sha,self.folder/'audit')
    def rows(self):return json_rows(self.run/'exact-matching-occurrences.jsonl.gz')
    def change_rows(self,mutate):
        rows=self.rows();mutate(rows);save_rows(self.run/'exact-matching-occurrences.jsonl.gz',rows);self.repin()
    def refusal(self):
        with self.assertRaises((ValueError,KeyError,TypeError)):self.run_audit()
    def test_success_complete_fixture_preserves_inputs(self):
        result=self.run_audit();self.assertEqual(result['source_files'],46)
        self.assertEqual(result['counts_by_caller']['gatk']['exact_matching_occurrences'],12)
        self.assertEqual(self.input_pins,{p:sha(p) for p in self.input_pins})
        self.assertIn('same pinned bcftools',result['normalization_engine'])
    def test_omitted_occurrence_even_rehashed_refused(self):
        self.change_rows(lambda rows:rows.pop());self.refusal()
    def test_extra_duplicate_even_rehashed_refused(self):
        self.change_rows(lambda rows:rows.append(copy.deepcopy(rows[0])));self.refusal()
    def test_source_ordinal_mutation_refused(self):
        self.change_rows(lambda rows:rows[0]['source_call']['native_context'].update(source_ordinal=999));self.refusal()
    def test_source_full_PL_mutation_refused(self):
        self.change_rows(lambda rows:rows[0]['source_call']['native_context']['FORMAT']['PL'].__setitem__(0,0));self.refusal()
    def test_extra_FORMAT_payload_mutation_refused(self):
        self.change_rows(lambda rows:rows[0]['source_call']['native_context']['FORMAT'].update(XX='changed'));self.refusal()
    def test_missing_GQ_imputation_refused(self):
        self.change_rows(lambda rows:rows[0]['quality']['GQ'].update(state='explicit_zero',value=0));self.refusal()
    def test_missing_GT_reclassified_as_called_refused(self):
        def mutation(rows):
            row=next(r for r in rows if r['quality']['genotype']['state']=='fully_missing')
            row['quality']['called_ALT_present']=True
        self.change_rows(mutation);self.refusal()
    def test_public_original_INFO_loss_refused(self):
        self.change_rows(lambda rows:rows[0]['clinvar'][0]['source']['INFO'].pop('EXTRA'));self.refusal()
    def test_public_multiallelic_overclaim_refused(self):
        def mutation(rows):
            row=next(r for r in rows if any(c['public_original_multiallelic'] for c in r['clinvar']))
            row['clinvar'][0]['stratum']='aggregate_plp_1star'
        self.change_rows(mutation);self.refusal()
    def test_wrong_ALT_lineage_refused(self):
        self.change_rows(lambda rows:rows[0]['source_call'].update(source_alt_index=2));self.refusal()
    def test_group_duplicate_count_mutation_refused(self):
        p=self.run/'exact-identity-groups.json';groups=json.loads(p.read_text());groups[0]['source_occurrences']+=1;save(p,groups);self.repin();self.refusal()
    def test_summary_called_count_mutation_refused(self):
        p=self.run/'summary.json';v=json.loads(p.read_text());v['counts_by_caller']['gatk']['unique_called_ALT_identities']+=1;save(p,v);self.repin();self.refusal()
    def test_summary_quality_count_mutation_refused(self):
        p=self.run/'summary.json';v=json.loads(p.read_text());v['exact_occurrence_quality_axes']['gatk']['GQ']['missing']=999;save(p,v);self.repin();self.refusal()
    def test_native_tagged_omission_refused(self):
        p=self.run/'gatk-1/ordinal-tagged-original.bcf'
        with pysam.VariantFile(p) as reader:header=reader.header.copy();records=[r.copy() for r in reader]
        with pysam.VariantFile(p,'wb',header=header) as writer:
            for r in records[:-1]:writer.write(r)
        self.repin();self.refusal()
    def test_public_selected_omission_refused(self):
        p=self.run/'public-selected.vcf.gz'
        with pysam.VariantFile(p) as reader:header=reader.header.copy();records=[r.copy() for r in reader]
        with pysam.VariantFile(p,'wz',header=header) as writer:
            for r in records[1:]:writer.write(r)
        self.repin();self.refusal()
    def test_unpinned_output_refused(self):
        (self.run/'extra.txt').write_text('unexpected');self.refusal()
    def test_stale_output_hash_refused(self):
        with (self.run/'exact-identity-groups.json').open('a') as f:f.write(' ')
        self.refusal()
    def test_wrong_summary_link_refused(self):
        p=self.run/'status.json';v=json.loads(p.read_text());v['summary_sha256']='0'*64;save(p,v);self.refusal()
    def test_false_scope_completion_refused(self):
        handoff=copy.deepcopy(json.loads(Path(self.config['inputs']['handoff']['path']).read_text()))
        handoff['all_autosomes_complete']=1
        original=a.Files.get
        def substituted(obj,pin):
            result=original(obj,pin)
            return handoff if pin['path']==self.config['inputs']['handoff']['path'] else result
        from unittest.mock import patch
        with patch.object(a.Files,'get',substituted):self.refusal()
    def test_duplicate_json_occurrence_key_refused(self):
        path=self.run/'exact-matching-occurrences.jsonl.gz'
        with gzip.open(path,'rt') as f:lines=f.readlines()
        lines[0]='{"caller":"gatk",'+lines[0].lstrip()[1:]
        with gzip.open(path,'wt') as f:f.writelines(lines)
        self.repin();self.refusal()
    def test_native_source_omitted_from_input_inventory_refused(self):
        p=self.run/'summary.json';v=json.loads(p.read_text())
        target=next(k for k in v['inputs'] if k.endswith('/gatk-22.vcf.gz'))
        del v['inputs'][target];save(p,v);self.repin();self.refusal()
    def test_missing_public_CLNALLE_header_end_to_end(self):
        root=self.folder/'without-public-header';root.mkdir()
        config,config_path=self.fixture.fixture(root)
        path=Path(config['inputs']['clinvar']['path'])
        with pysam.VariantFile(str(path)) as reader:records=[str(r).rstrip('\n').replace(';CLNALLE=1','') for r in reader]
        config['inputs']['clinvar']=self.fixture.write_vcf(path,records,[],public=True,public_clnalle_header=False)
        config['inputs']['clinvar_index']=self.fixture.pin(str(path)+'.tbi');save(config_path,config)
        with pysam.VariantFile(str(path)) as reader:self.assertNotIn('CLNALLE',reader.header.info)
        run=root/'producer';self.fixture.a.run(config_path,sha(config_path),run)
        result=a.audit(config_path,sha(config_path),run,sha(run/'status.json'),self.qualification,self.qualification_sha,root/'independent-audit')
        self.assertEqual(result['status'],'complete')
        contexts=[c for row in json_rows(run/'exact-matching-occurrences.jsonl.gz') for c in row['clinvar']]
        self.assertTrue(contexts)
        self.assertTrue(all(c['assertion_mapping']['normalized_CLNALLE'] is None and c['assertion_mapping']['original_CLNALLE'] is None for c in contexts))
        multi=[c for c in contexts if c['public_original_multiallelic']]
        self.assertTrue(multi);self.assertTrue(all(c['stratum']=='multiallelic_assertion_mapping_review' for c in multi))
    def test_existing_audit_directory_preserved(self):
        out=self.folder/'audit';out.mkdir();(out/'sentinel').write_text('keep')
        self.refusal();self.assertEqual((out/'sentinel').read_text(),'keep')


if __name__=='__main__':unittest.main(verbosity=2)
