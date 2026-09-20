import copy
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

SOURCE=Path(__file__).resolve().parents[1]/'scripts/preflight_public_targets.py'
spec=importlib.util.spec_from_file_location('public_target_tested',SOURCE)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
try:
    import pysam
except ImportError:
    pysam=None


def save(path,value):path.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
def pin(path):return dict(path=str(path.resolve()),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),bytes=path.stat().st_size)


def fixture(root, native=False):
    root=Path(root);ref=root/'reference.fa';ref.write_text('>X\nAAAAANAAAAAAAAAAAAAAAAAA\n')
    if native:pysam.faidx(str(ref))
    else:Path(str(ref)+'.fai').write_text('X\t24\t3\t24\t25\n')
    panel=root/'panel.json';save(panel,dict(total=2,rows=[dict(symbol='ALPHA',hgnc_id='HGNC:1'),dict(symbol='BETA',hgnc_id='HGNC:2')]))
    gff=root/'annotation.gff3.gz'
    with gzip.open(gff,'wt') as stream:
        stream.write('##gff-version 3\nX\tx\tgene\t1\t20\t.\t+\t.\tID=gene:A;Name=ALPHA;biotype=protein_coding;description=fictional [Source:HGNC Symbol%3BAcc:HGNC:1];version=2\nX\tx\tmRNA\t1\t20\t.\t+\t.\tID=transcript:A;Parent=gene:A\nX\tx\tCDS\t3\t8\t.\t+\t0\tParent=transcript:A\nX\tx\tgene\t1\t20\t.\t-\t.\tID=gene:B;Name=BETA;biotype=protein_coding;description=fictional [Source:HGNC Symbol%3BAcc:HGNC:2]\nX\tx\tmRNA\t1\t20\t.\t-\t.\tID=transcript:B;Parent=gene:B\nX\tx\tCDS\t7\t10\t.\t-\t0\tParent=transcript:B\n')
    sources={}
    for name in m.QUALIFIED:
        p=root/name;p.write_bytes((m.VENDOR/name).read_bytes() if name=='coding_depth_targets.py' else b'# fictional unused qualification dependency\n')
        sources[name]=pin(p)
    q=root/'controls.json';save(q,dict(status='complete',exit_code=0,source_sha256={'workflow/scripts/'+k:v['sha256'] for k,v in sources.items()}))
    cfg=dict(schema_version=1,public_preflight_authorized=True,fixture_only=True,
        inputs=dict(panel=pin(panel),gff=pin(gff),reference=pin(ref),reference_index=pin(Path(str(ref)+'.fai'))),
        genes=['ALPHA','BETA'],expected_genes=2,producer_qualification=pin(q),qualified_sources=sources,
        expected_gene_identity={'ALPHA':dict(gene_id='gene:A',hgnc_id='HGNC:1'),'BETA':dict(gene_id='gene:B',hgnc_id='HGNC:2')},
        output=str(root/'targets'),approval_note='Fictional fixture only.')
    cp=root/'config.json';save(cp,cfg)
    return cfg,cp


class Contract(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name);self.cfg,self.cp=fixture(self.root)
    def config(self,cfg=None):
        save(self.cp,cfg or self.cfg)
        return m.configuration(self.cp,pin(self.cp)['sha256'],Path(self.cfg['output']))
    def test_public_only_inputs_exact_and_case_path_not_touched(self):
        self.cfg['inputs']['alignment']=dict(path='/forbidden/patient.cram',sha256='0'*64)
        with patch.object(m.guard,'fingerprint',wraps=m.guard.fingerprint) as seen:
            with self.assertRaisesRegex(ValueError,'public-only input'):self.config()
        self.assertEqual([str(call.args[0]) for call in seen.call_args_list],[str(self.cp)])
    def test_disabled_public_authorization_and_typed_flags(self):
        for key,value in [('public_preflight_authorized',False),('public_preflight_authorized',1),('fixture_only',0),('expected_genes',True),('schema_version',True)]:
            bad=copy.deepcopy(self.cfg);bad[key]=value
            with self.subTest(key=key,value=value),self.assertRaises(ValueError):self.config(bad)
    def test_strict_json_nan_and_duplicate_keys(self):
        for text in ['{"a":NaN}','{"a":Infinity}','{"a":1,"a":2}']:
            with self.assertRaises(ValueError):m.guard.parse(text)
    def test_wrong_config_hash_fails_before_interpretation(self):
        with self.assertRaisesRegex(ValueError,'sha256'):m.configuration(self.cp,'0'*64,Path(self.cfg['output']))
    def test_source_pin_changed(self):
        self.cfg['qualified_sources']['coding_depth_targets.py']['sha256']='0'*64
        cfg,p=self.config()
        with self.assertRaisesRegex(ValueError,'sha256'):m.bind_inputs(cfg,p)
    def test_public_pin_changed(self):
        self.cfg['inputs']['gff']['sha256']='0'*64
        cfg,p=self.config()
        with self.assertRaisesRegex(ValueError,'sha256'):m.bind_inputs(cfg,p)
    def test_panel_omission_duplicate_and_hgnc_refused(self):
        p=Path(self.cfg['inputs']['panel']['path'])
        for rows in [[dict(symbol='ALPHA',hgnc_id='HGNC:1')],[dict(symbol='ALPHA',hgnc_id='HGNC:1')]*2,[dict(symbol='ALPHA',hgnc_id='HGNC:999'),dict(symbol='BETA',hgnc_id='HGNC:2')]]:
            save(p,dict(total=2,rows=rows));self.cfg['inputs']['panel']=pin(p);cfg,cp=self.config()
            with self.subTest(rows=rows),self.assertRaises(ValueError):m.bind_inputs(cfg,cp)
    def test_qualification_link_or_status_refused(self):
        p=Path(self.cfg['producer_qualification']['path']);original=json.loads(p.read_text())
        for bad in [dict(original,status='failed'),dict(original,exit_code=True),dict(original,source_sha256={})]:
            save(p,bad);self.cfg['producer_qualification']=pin(p);cfg,cp=self.config()
            with self.assertRaises((ValueError,KeyError)):m.bind_inputs(cfg,cp)
    def test_existing_output_wrong_output_and_index_refused(self):
        out=Path(self.cfg['output']);out.mkdir()
        with self.assertRaisesRegex(ValueError,'existing output'):self.config()
        out.rmdir();bad=copy.deepcopy(self.cfg);bad['output']=str(self.root/'wrong')
        with self.assertRaisesRegex(ValueError,'output path'):self.config(bad)
        wrong=self.root/'wrong.fai';wrong.write_bytes(Path(self.cfg['inputs']['reference_index']['path']).read_bytes());self.cfg['inputs']['reference_index']=pin(wrong)
        cfg,cp=self.config()
        with self.assertRaisesRegex(ValueError,'index binding'):m.bind_inputs(cfg,cp)
    def test_omitted_duplicated_and_unsorted_genes_refused(self):
        for genes in [['ALPHA'],['ALPHA','ALPHA'],['BETA','ALPHA']]:
            bad=copy.deepcopy(self.cfg);bad['genes']=genes
            with self.assertRaises(ValueError):self.config(bad)
    def test_immutable_vendor_and_proposed_limits(self):
        for name,digest in m.PINS.items():self.assertEqual(pin(m.VENDOR/name)['sha256'],digest)
        self.assertEqual(m.LIMITS,dict(wall_seconds=1800,group_RSS_bytes=16*1024**3,output_bytes=2*1024**3,minimum_free_bytes=40*1024**3))
        for key in ['wall_seconds','group_RSS_bytes','output_bytes']:
            values=dict(elapsed=0,rss=0,size=0,free=40*1024**3);values[{'wall_seconds':'elapsed','group_RSS_bytes':'rss','output_bytes':'size'}[key]]=m.LIMITS[key]
            self.assertEqual(m.guard.limit_reason(**values,limits=m.LIMITS),key)
    def test_worker_cannot_be_directly_invoked(self):
        r=subprocess.run([sys.executable,'-B',str(SOURCE),'--config',str(self.cp),'--expected-config-sha256',pin(self.cp)['sha256'],'--expected-source-sha256',pin(SOURCE)['sha256'],'--output',self.cfg['output'],'--internal-worker'],capture_output=True,text=True)
        self.assertNotEqual(r.returncode,0);self.assertIn('owned parent/session',r.stderr);self.assertFalse(Path(self.cfg['output']).exists())


@unittest.skipIf(pysam is None,'native pysam integration unavailable')
class Native(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name);self.cfg,self.cp=fixture(self.root,True)
    def prepare(self):
        save(self.cp,self.cfg);cfg,p=m.configuration(self.cp,pin(self.cp)['sha256'],Path(self.cfg['output']));m.prepare(cfg,p,Path(self.cfg['output']))
    def test_reference_only_x_overlap_n_and_exact_outputs(self):
        self.prepare();out=Path(self.cfg['output']);s=json.loads((out/'target-summary.json').read_text());d=json.loads((out/'reference-design.json').read_text())
        self.assertEqual((s['genes'],s['windows'],s['panel_union_bases'],s['gene_union_bases_sum']),(2,3,8,10))
        self.assertEqual((s['panel_ACGT_bases'],s['panel_non_ACGT_bases']),(7,1));self.assertIs(s['reads_inspected'],False)
        self.assertEqual(d['windows'][1]['gene_ids'],['gene:A','gene:B']);self.assertTrue(all(g['original_gene_feature'] for g in d['genes']))
        st=json.loads((out/'status.json').read_text());self.assertEqual(set(st['outputs']),{'reference-design.json','regions.bed','target-summary.json'})
        for name,expected in st['outputs'].items():self.assertEqual(pin(out/name),expected)
    def test_missing_cds_ambiguous_gene_and_bad_parent(self):
        p=Path(self.cfg['inputs']['gff']['path']);text=gzip.open(p,'rt').read()
        changes=[text.replace('Parent=transcript:A','Parent=transcript:missing'), '\n'.join(line for line in text.splitlines() if 'Parent=transcript:A' not in line)+'\n',text+'X\tx\tgene\t1\t20\t.\t+\t.\tID=gene:extra;Name=ALPHA;biotype=protein_coding\n']
        for value in changes:
            with gzip.open(p,'wt') as f:f.write(value)
            self.cfg['inputs']['gff']=pin(p)
            with self.assertRaises(ValueError):self.prepare()
            self.assertFalse(Path(self.cfg['output']).exists())
    def test_gene_identity_or_gff_hgnc_mismatch(self):
        for key,value in [('gene_id','gene:wrong'),('hgnc_id','HGNC:999')]:
            old=self.cfg['expected_gene_identity']['ALPHA'][key];self.cfg['expected_gene_identity']['ALPHA'][key]=value
            if key=='hgnc_id':
                pp=Path(self.cfg['inputs']['panel']['path']);panel=json.loads(pp.read_text());panel['rows'][0]['hgnc_id']=value;save(pp,panel);self.cfg['inputs']['panel']=pin(pp)
            with self.assertRaises(ValueError):self.prepare()
            self.cfg['expected_gene_identity']['ALPHA'][key]=old
    def test_guarded_fictional_end_to_end(self):
        if m.shutil.disk_usage(self.root).free<m.LIMITS['minimum_free_bytes']:self.skipTest('guard requires 40GiB free')
        r=subprocess.run([sys.executable,'-B',str(SOURCE),'--config',str(self.cp),'--expected-config-sha256',pin(self.cp)['sha256'],'--expected-source-sha256',pin(SOURCE)['sha256'],'--output',self.cfg['output']],capture_output=True,text=True)
        self.assertEqual(r.returncode,0,r.stderr);launch=Path(self.cfg['output']+'-launch');state=json.loads((launch/'status.json').read_text())
        self.assertEqual(state['status'],'complete');self.assertFalse(state['reads_inspected'])
        self.assertEqual(state['result_status'],pin(Path(self.cfg['output'])/'status.json'))


if __name__=='__main__':unittest.main()
