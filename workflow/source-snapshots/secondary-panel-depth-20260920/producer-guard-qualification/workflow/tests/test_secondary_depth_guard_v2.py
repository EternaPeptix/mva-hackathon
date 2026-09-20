"""Fictional metadata and fail-closed wrapper controls; no clinical inputs."""
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

SOURCE=Path(__file__).resolve().parents[1]/'scripts/launch_secondary_depth_v2.py'
spec=importlib.util.spec_from_file_location('secondary_depth_tested',SOURCE)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def pin(path):return dict(path=str(path.resolve()),bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest())
def save(path,obj):path.write_bytes(m.encoded(obj));return pin(path)


class Fixture:
    def __init__(self,root):
        self.root=root; self.public=root/'public';self.public.mkdir();self.src=root/'source';self.src.mkdir()
        self.output=root/'depth';self.inputs={}
        for name in ['alignment','alignment_index','reference','reference_index','gff','samtools','panel']:
            p=root/name;p.write_text('FICTIONAL '+name+'\n');self.inputs[name]=pin(p)
        sources={}
        for name in m.QUALIFIED:
            p=self.src/name;p.write_text('# fictional '+name+'\n');sources[name]=pin(p)
        q=save(root/'qualification.json',dict(status='complete',exit_code=0,source_sha256={'workflow/scripts/'+n:p['sha256'] for n,p in sources.items()}))
        self.design=dict(status='prepared',reads_inspected=False,genes=[dict(name='FAKE',id='gene:F',chrom='1',strand='+',CDS_union_bases=3,transcript_ids=['transcript:T'])],windows=[dict(chrom='1',start=3,end=6,gene_ids=['gene:F'],reference_ACGT_bases=3,reference_non_ACGT_bases=0,reference_sequence_sha256='f'*64)],panel_union_bases=3,gene_union_bases_sum=3)
        self.target=dict(genes=1,windows=1,panel_union_bases=3,gene_union_bases_sum=3,panel_ACGT_bases=3,panel_non_ACGT_bases=0)
        config=save(self.public/'config.json',dict(genes=['FAKE']))
        d=save(self.public/'reference-design.json',self.design)
        summary=save(self.public/'target-summary.json',dict(status='complete',reads_inspected=False,**self.target))
        bed=self.public/'regions.bed';bed.write_text('1\t3\t6\n')
        state=save(self.public/'status.json',dict(status='complete',reads_inspected=False,inputs=dict(config=config,**{n:self.inputs[n] for n in ['reference','reference_index','gff','panel']}),outputs={'reference-design.json':d,'target-summary.json':summary,'regions.bed':pin(bed)}))
        launch=save(self.public/'launch.json',dict(status='complete',exit_code=0,result_status=state))
        inner=dict(schema_version=1,assembly='GRCh38',sample='FICTIONAL',expected_read_groups=1,genes=['FAKE'],inputs={n:p for n,p in self.inputs.items() if n!='panel'},gene_list_source=self.inputs['panel'],source_sha256={n:p['sha256'] for n,p in sources.items()},depth_modes=copy.deepcopy(m.MODES))
        definition=save(root/'definition.json',dict(shards=[[['FICTIONAL'+str(i),100]] for i in range(64)],sample='FICTIONAL',alignment_sha256=self.inputs['alignment']['sha256'],reference_sha256=self.inputs['reference']['sha256']))
        caller=save(root/'caller.json',dict(status='complete',completed_shards=list(range(64)),expected_shards=64,definition_sha256=definition['sha256']))
        self.cfg=dict(schema_version=1,ready=True,preparation_note='Fictional control',fixture_only=True,output=str(self.output),limits=copy.deepcopy(m.LIMITS),producer_configuration=inner,producer_config_sha256=hashlib.sha256(m.encoded(inner)).hexdigest(),qualified_sources=sources,producer_qualification=q,public_preflight={'config':config,'status':state,'launch':launch,'reference-design.json':d,'target-summary.json':summary,'regions.bed':pin(bed)},runtime=pin(Path(sys.executable).resolve()),completed_caller_status=caller,caller_definition=definition,expected_target=self.target)

    def mutate_public(self,name,value,bind_status=True):
        p=Path(self.cfg['public_preflight'][name]['path']);self.cfg['public_preflight'][name]=save(p,value)
        if name in ('reference-design.json','target-summary.json','regions.bed') and bind_status:
            state=json.loads(Path(self.cfg['public_preflight']['status']['path']).read_text());state['outputs'][name]=self.cfg['public_preflight'][name];self.mutate_public('status',state)
        if name=='status' and bind_status:
            launch=json.loads(Path(self.cfg['public_preflight']['launch']['path']).read_text());launch['result_status']=self.cfg['public_preflight']['status'];self.mutate_public('launch',launch)

    def inner(self):self.cfg['producer_config_sha256']=hashlib.sha256(m.encoded(self.cfg['producer_configuration'])).hexdigest()

    def terminal(self):
        out=self.output;out.mkdir();inner=out.parent/'producer-config.json';innerpin=save(inner,self.cfg['producer_configuration'])
        selection=m.core_design(self.design);selection['config']=self.cfg['producer_configuration'];save(out/'selection.json',selection)
        (out/'regions.bed').write_text('1\t3\t6\n')
        inputs={'config':innerpin}
        summary=dict(status='complete',inputs=inputs,panel_union_bases=3,gene_union_bases_sum=3,depth_modes={v['name']:dict(mode=v,position_count=3) for v in m.MODES})
        save(out/'summary.json',summary)
        for n in ['selected.bam','selected.bam.bai','commands.json','extraction-counts.json','extract-preserving-records.log','quickcheck.log','index.log']:
            (out/n).write_text('fictional placeholder, not alignment data\n')
        for mode in m.MODES:
            for suffix in ['.json','.depth.tsv.gz','.log']:(out/(mode['name']+suffix)).write_text('fictional\n')
        save(out/'status.json',dict(status='complete',summary_sha256=pin(out/'summary.json')['sha256'],inputs=inputs,outputs={p.name:pin(p) for p in out.iterdir()}))
        return innerpin


class Controls(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup);self.f=Fixture(Path(self.tmp.name));self.cfg=self.f.cfg
    def test_complete_fictional_bind_without_opening_alignment(self):
        m.shape(self.cfg,self.f.output)
        real=m.guard.fingerprint
        def guarded(path,*a,**k):
            self.assertNotIn(str(path),[self.f.inputs[n]['path'] for n in ['alignment','alignment_index']]);return real(path,*a,**k)
        with mock.patch.object(m.guard,'fingerprint',side_effect=guarded):
            d,p=m.bind(self.cfg)
        self.assertEqual(d,self.f.design);self.assertIn('completed_caller',p)
    def test_disabled_and_nonboolean_authorization(self):
        for value in [False,1,'true']:
            self.cfg['ready']=value
            with self.assertRaises(ValueError):m.shape(self.cfg,self.f.output)
    def test_existing_output_or_launch_preserved(self):
        for p in [self.f.output,self.f.output.with_name('depth-launch')]:
            p.mkdir()
            with self.assertRaises(ValueError):m.shape(self.cfg,self.f.output)
            p.rmdir()
    def test_inner_config_pin_and_modes(self):
        self.cfg['producer_configuration']['depth_modes'][0]['minimum_BQ']=0
        with self.assertRaises(ValueError):m.shape(self.cfg,self.f.output)
        self.f.inner()
        with self.assertRaises(ValueError):m.shape(self.cfg,self.f.output)
    def test_approval_delta_rejects_scope_and_boolean_change(self):
        draft=copy.deepcopy(self.cfg);draft['ready']=False;m.guard.approved_delta(draft,self.cfg)
        changed=copy.deepcopy(self.cfg);changed['producer_configuration']['genes']=['OTHER']
        with self.assertRaises(ValueError):m.guard.approved_delta(draft,changed)
        changed=copy.deepcopy(self.cfg);changed['ready']=1
        with self.assertRaises(ValueError):m.guard.approved_delta(draft,changed)
    def test_stale_public_output_pin(self):
        Path(self.cfg['public_preflight']['target-summary.json']['path']).write_text('{}')
        with self.assertRaises(ValueError):m.bind(self.cfg)
    def test_public_status_and_guard_must_be_complete(self):
        for name in ['status','launch']:
            old=json.loads(Path(self.cfg['public_preflight'][name]['path']).read_text());v=copy.deepcopy(old);v['status']='failed';self.f.mutate_public(name,v)
            with self.assertRaises(ValueError):m.bind(self.cfg)
            self.f.mutate_public(name,old)
    def test_public_guard_wrong_path_even_same_bytes(self):
        s=Path(self.cfg['public_preflight']['status']['path']);other=s.with_name('copy.json');other.write_bytes(s.read_bytes())
        launch=json.loads(Path(self.cfg['public_preflight']['launch']['path']).read_text());launch['result_status']=pin(other);self.f.mutate_public('launch',launch)
        with self.assertRaises(ValueError):m.bind(self.cfg)
    def test_gene_scope_and_input_refusal(self):
        self.cfg['producer_configuration']['genes']=['OTHER'];self.f.inner()
        with self.assertRaises(ValueError):m.bind(self.cfg)
        self.cfg['producer_configuration']['genes']=['FAKE'];self.cfg['producer_configuration']['inputs']['gff']['sha256']='0'*64
        with self.assertRaises(ValueError):m.bind(self.cfg)
    def test_denominator_type_and_target_drift(self):
        for field,value in [('genes',True),('panel_union_bases',4)]:
            summary=json.loads(Path(self.cfg['public_preflight']['target-summary.json']['path']).read_text());summary[field]=value;self.f.mutate_public('target-summary.json',summary)
            with self.assertRaises(ValueError):m.bind(self.cfg)
            summary[field]=self.target_value(field);self.f.mutate_public('target-summary.json',summary)
    def target_value(self,field):return self.f.target[field]
    def test_bed_drift(self):
        p=Path(self.cfg['public_preflight']['regions.bed']['path']);p.write_text('1\t3\t7\n');self.cfg['public_preflight']['regions.bed']=pin(p)
        s=json.loads(Path(self.cfg['public_preflight']['status']['path']).read_text());s['outputs']['regions.bed']=pin(p);self.f.mutate_public('status',s)
        with self.assertRaises(ValueError):m.bind(self.cfg)
    def test_qualified_source_and_receipt_link(self):
        p=Path(self.cfg['qualified_sources']['measure_coding_depth.py']['path']);p.write_text('# mutation\n')
        with self.assertRaises(ValueError):m.bind(self.cfg)
        self.cfg['qualified_sources']['measure_coding_depth.py']=pin(p)
        with self.assertRaises(ValueError):m.bind(self.cfg)
    def test_original_caller_completed_scope(self):
        for n,value in [('status','running'),('completed_shards',63),('expected_shards',True)]:
            path=Path(self.cfg['completed_caller_status']['path']);original=json.loads(path.read_text());v=copy.deepcopy(original);v[n]=value;self.cfg['completed_caller_status']=save(path,v)
            with self.assertRaises(ValueError):m.bind(self.cfg)
            self.cfg['completed_caller_status']=save(path,original)
    def test_completed_shard_duplicate_refused(self):
        path=Path(self.cfg['completed_caller_status']['path']);value=json.loads(path.read_text());value['completed_shards'][-1]=62;self.cfg['completed_caller_status']=save(path,value)
        with self.assertRaises(ValueError):m.bind(self.cfg)
    def test_completed_shard_missing_refused(self):
        path=Path(self.cfg['completed_caller_status']['path']);value=json.loads(path.read_text());value['completed_shards'].pop();self.cfg['completed_caller_status']=save(path,value)
        with self.assertRaises(ValueError):m.bind(self.cfg)
    def test_completed_shard_noninteger_refused(self):
        path=Path(self.cfg['completed_caller_status']['path']);original=json.loads(path.read_text())
        for item in [False,0.0,'0']:
            value=copy.deepcopy(original);value['completed_shards'][0]=item;self.cfg['completed_caller_status']=save(path,value)
            with self.assertRaises(ValueError):m.bind(self.cfg)
    def test_runtime_pin_and_path(self):
        p=Path(self.cfg['runtime']['path']);other=self.f.root/'python-copy';other.write_bytes(p.read_bytes());self.cfg['runtime']=pin(other)
        with self.assertRaises(ValueError):m.bind(self.cfg)
    def test_strict_json_nonfinite_and_duplicate(self):
        for value in ['{"x":NaN}','{"x":Infinity}','{"x":1,"x":2}']:
            with self.assertRaises(ValueError):m.guard.parse(value)
    def test_complete_terminal_then_target_mutation_refused(self):
        inner=self.f.terminal();m.terminal(self.cfg,self.f.output,inner,self.f.design)
        p=self.f.output/'selection.json';v=json.loads(p.read_text());v['windows'][0]['end']=7;save(p,v)
        state=json.loads((self.f.output/'status.json').read_text());state['outputs']['selection.json']=pin(p);save(self.f.output/'status.json',state)
        with self.assertRaises(ValueError):m.terminal(self.cfg,self.f.output,inner,self.f.design)
    def test_limits_and_owned_rss(self):
        self.assertIsNone(m.guard.limit_reason(1,1,1,100*1024**3,m.LIMITS))
        for key in m.LIMITS:
            values=dict(wall_seconds=0,group_RSS_bytes=0,output_bytes=0,minimum_free_bytes=100*1024**3);values[key]=m.LIMITS[key]-(1 if key=='minimum_free_bytes' else 0)
            self.assertEqual(m.guard.limit_reason(*values.values(),m.LIMITS),key)
        self.assertEqual(m.guard.owned_rss([dict(pgid=1,RSS_bytes=7),dict(pgid=2,RSS_bytes=100)],1),7)
    def test_worker_direct_cli_refused_without_parent(self):
        p=subprocess.run([sys.executable,'-B',str(SOURCE),'--config','missing','--proposal','missing','--output',str(self.f.output),'--qualification','missing','--expected-config-sha256','0'*64,'--expected-proposal-sha256','0'*64,'--expected-source-sha256','0'*64,'--expected-qualification-sha256','0'*64,'--internal-worker'],capture_output=True,text=True)
        self.assertNotEqual(p.returncode,0);self.assertIn('worker requires owned parent/session',p.stderr);self.assertFalse(self.f.output.exists())


if __name__=='__main__':unittest.main()
