#!/usr/bin/env python3
"""Fictional process controls for whole-pilot guard; no native GATK/patient work."""
import copy,hashlib,importlib.util,json,os,subprocess,sys,tempfile,time,unittest
from pathlib import Path
from unittest.mock import patch
s=importlib.util.spec_from_file_location('guard',Path(__file__).with_name('guard-fictional-pilot.py'));g=importlib.util.module_from_spec(s);s.loader.exec_module(g)
class Controls(unittest.TestCase):
    def limits(self,**patch):return dict(g.LIMITS,minimum_free_bytes=0,maximum_fixture_total_bytes=1024**2,maximum_group_RSS_bytes=256*1024**2,total_wall_seconds=5,**patch) if not patch else (dict(g.LIMITS)|dict(minimum_free_bytes=0,maximum_fixture_total_bytes=1024**2,maximum_group_RSS_bytes=256*1024**2,total_wall_seconds=5)|patch)
    def command(self,out,code):return [sys.executable,'-B','-c',code,str(out)]
    def test_exact_caller_command_and_start_identity(self):
        line='79898 79800 Sun Sep 20 01:02:03 2026 '+' '.join(g.EXPECTED_CALLER_ARGV);before=g.parse_caller(line);g.exact_identity(before,dict(before))
        for source in [line.replace('--workers 4','--workers 5'),line.replace('from-fastq-cram30','other'),line.replace('79898 ','79899 ',1)]:
            with self.assertRaises(ValueError):g.parse_caller(source)
        for patch in [dict(started='new'),dict(pgid=2),dict(pid=1)]:
            with self.assertRaises(ValueError):g.exact_identity(before,before|patch)
    def test_all_resource_boundaries(self):
        a=g.LIMITS;self.assertIsNone(g.limit_reason(0,0,0,a['minimum_free_bytes'],a))
        for i,key in enumerate(['total_wall_seconds','maximum_group_RSS_bytes','maximum_fixture_total_bytes']):
            vals=[0,0,0,a['minimum_free_bytes'],a];vals[i]=a[key];self.assertEqual(g.limit_reason(*vals),key)
        self.assertEqual(g.limit_reason(0,0,0,a['minimum_free_bytes']-1,a),'minimum_free_bytes')
    def test_full_python_and_inherited_child_success(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);out=root/'out';launch=root/'launch'
            code="import pathlib,subprocess,sys; p=pathlib.Path(sys.argv[1]);p.mkdir();(p/'fixture').write_bytes(b'x'*100);subprocess.run([sys.executable,'-c','import time;time.sleep(.1)'],check=True)"
            v=g.supervise(self.command(out,code),out,launch,self.limits(),root);self.assertEqual(v['status'],'complete');self.assertGreater(v['peak_sampled_RSS_bytes'],0);self.assertFalse(g.group(v['pgid']))
    def test_python_output_growth_guard_before_any_external_tool(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);out=root/'out';launch=root/'launch';code="import pathlib,sys,time;p=pathlib.Path(sys.argv[1]);p.mkdir();(p/'fixture').write_bytes(b'x'*65536);time.sleep(3)"
            with self.assertRaisesRegex(ValueError,'maximum_fixture_total_bytes'):g.supervise(self.command(out,code),out,launch,self.limits(maximum_fixture_total_bytes=32768),root)
            v=json.loads((launch/'status.json').read_text());self.assertFalse(g.group(v['pgid']));self.assertEqual((out/'fixture').stat().st_size,65536)
    def test_python_RSS_guard_before_any_external_tool(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);out=root/'out';launch=root/'launch';code="import time;v=bytearray(80*1024**2);time.sleep(3)"
            with self.assertRaisesRegex(ValueError,'maximum_group_RSS_bytes'):g.supervise(self.command(out,code),out,launch,self.limits(maximum_group_RSS_bytes=64*1024**2),root)
            self.assertFalse(g.group(json.loads((launch/'status.json').read_text())['pgid']))
    def test_timeout_kills_inherited_descendant_preserves_unrelated(self):
        other=subprocess.Popen([sys.executable,'-c','import time;time.sleep(30)'],start_new_session=True)
        try:
            with tempfile.TemporaryDirectory() as d:
                root=Path(d);out=root/'out';launch=root/'launch';code="import subprocess,sys,time;subprocess.Popen([sys.executable,'-c','import time;time.sleep(30)']);time.sleep(30)"
                with self.assertRaisesRegex(ValueError,'total_wall_seconds'):g.supervise(self.command(out,code),out,launch,self.limits(total_wall_seconds=.3),root)
                self.assertIsNone(other.poll());self.assertFalse(g.group(json.loads((launch/'status.json').read_text())['pgid']))
        finally:other.terminate();other.wait(timeout=5)
    def test_no_retry_no_output_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);out=root/'out';out.mkdir();(out/'sentinel').write_text('keep')
            with self.assertRaises(ValueError):g.supervise(['false'],out,root/'launch',self.limits(),root)
            self.assertEqual((out/'sentinel').read_text(),'keep')
    def test_child_nonzero_failure_preserved(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);out=root/'out';launch=root/'launch'
            with self.assertRaisesRegex(ValueError,'runner exit 7'):g.supervise(self.command(out,'raise SystemExit(7)'),out,launch,self.limits(),root)
            self.assertEqual(json.loads((launch/'status.json').read_text())['exit_code'],7)
    def test_strict_json_hash_and_symlink_refusal(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);p=root/'data.json'
            for text in ['{"a":1,"a":2}','{"a":NaN}']:
                p.write_text(text)
                with self.assertRaises(ValueError):g.read_json(p)
            p.write_text('{}');self.assertEqual(len(g.sha(p)),64);(root/'link').symlink_to(p)
            with self.assertRaises(ValueError):g.size(root)
    def terminal_fixture(self):
        definition=dict(shards=[['fictional']]*64,sample='WGS_EX2312012')
        shards=[dict(completed_at='fictional',commands_sha256='a'*64,**{c:dict(counts=dict(records=2)) for c in ['gatk','bcftools']}) for _ in range(64)]
        methods={c:dict(path=str(g.CALLER_RUN/('WGS_EX2312012.'+c+'.normalized.raw.vcf.gz')),sha256='b'*64,counts=dict(records=128)) for c in ['gatk','bcftools']}
        state=dict(status='complete',pid=79898,definition_sha256=g.DEFINITION_SHA,expected_shards=64,completed_shards=list(range(64)),completed_at='later',started_at='earlier',outputs=copy.deepcopy(methods))
        return state,definition,shards,methods
    def test_natural_completion_exact_metadata_gate(self):
        values=self.terminal_fixture();self.assertEqual(g.validate_terminal_metadata(*values,g.CALLER_RUN),dict(gatk=128,bcftools=128))
    def test_incomplete_or_changed_completion_refused(self):
        for kind in ['missing_shard','missing_member','wrong_definition','wrong_PID','wrong_path','wrong_total','receipt_disagreement','failed','bool_count']:
            state,definition,shards,methods=self.terminal_fixture()
            if kind=='missing_shard':shards.pop()
            elif kind=='missing_member':state['completed_shards'].pop()
            elif kind=='wrong_definition':state['definition_sha256']='c'*64
            elif kind=='wrong_PID':state['pid']=79999
            elif kind=='wrong_path':state['outputs']['gatk']['path']='/another/run/file.vcf.gz';methods['gatk']=copy.deepcopy(state['outputs']['gatk'])
            elif kind=='wrong_total':state['outputs']['gatk']['counts']['records']=127;methods['gatk']=copy.deepcopy(state['outputs']['gatk'])
            elif kind=='receipt_disagreement':methods['gatk']['sha256']='c'*64
            elif kind=='failed':state['status']='failed'
            elif kind=='bool_count':shards[0]['gatk']['counts']['records']=True
            with self.subTest(kind=kind),self.assertRaises(ValueError):g.validate_terminal_metadata(state,definition,shards,methods,g.CALLER_RUN)
    def test_running_to_complete_transition_and_no_restart(self):
        identity=dict(pid=79898,started_at='original',definition_sha256=g.DEFINITION_SHA,run=str(g.CALLER_RUN),expected_shards=64)
        live=dict(mode='running',identity=identity,process=dict(pid=79898,started='original_process'),metadata_pins={})
        complete=dict(mode='complete_metadata_only',identity=identity,process=None,metadata_pins={'terminal':'fixed'})
        g.exact_identity(live,complete);g.exact_identity(complete,copy.deepcopy(complete))
        for before,after in [(complete,live),(live,complete|dict(identity=identity|dict(started_at='restart'))),(complete,complete|dict(metadata_pins={'terminal':'changed'}))]:
            with self.assertRaises(ValueError):g.exact_identity(before,after)
    def test_completed_terminal_SHA_mismatch_refused_without_case_files(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);definition=root/'run-definition.json';definition.write_text('{}')
            digest=hashlib.sha256(definition.read_bytes()).hexdigest();state=root/'status.json';state.write_text(json.dumps(dict(status='complete',pid=79898,definition_sha256=digest,expected_shards=64,started_at='fictional')))
            with patch.multiple(g,CALLER_RUN=root,DEFINITION_SHA=digest,TERMINAL_STATUS_SHA='0'*64),patch.object(g.subprocess,'run',return_value=type('Result',(),dict(returncode=1,stdout=''))()),self.assertRaisesRegex(ValueError,'terminal caller status SHA'):
                g.caller_identity()
    def test_complete_metadata_process_cannot_restart_or_change(self):
        done=dict(mode='complete_metadata_only',identity={'run':'fixed'},process=None,metadata_pins={'terminal':'fixed'})
        draining=done|dict(process=dict(pid=79898,started='original',pgid=79800,argv=['fixed']))
        g.exact_identity(draining,copy.deepcopy(draining));g.exact_identity(draining,done);g.exact_identity(done,copy.deepcopy(done))
        for before,after in [(done,draining),(draining,draining|dict(process=draining['process']|dict(started='restarted'))),(draining,draining|dict(process=draining['process']|dict(pgid=1)))]:
            with self.assertRaises(ValueError):g.exact_identity(before,after)
if __name__=='__main__':unittest.main(verbosity=2)
