#!/usr/bin/env python3
"""One explicitly authorized public transformation and independent reconciliation."""
import argparse,hashlib,importlib.util,json,os,platform,subprocess,sys
from pathlib import Path

HERE=Path(__file__).resolve().parent
GUARD_SHA='61eb3dffa51c23f30e203036d43d9fbdee76d314b5e869daf9fefbbd7aaacb18'
TRANSFORM_SHA='8ac5a7397284252a2ad259ea58ded0f887809e26b900ced97961bd844000bdab'
PLAN_SHA='be9ed02e121ab5b276aa3418917d11e2e358ec30e0b0fa97376e3bc0d49e3bbd'
ROOT=Path('/Users/jeweled/Documents/mva-hackathon')
OUT=ROOT/'results/nuclear-mixture-pilot-20260920/public-panel-v1'
LAUNCH=ROOT/'methods/nuclear-mixture-pilot-20260920/public-panel-launch-v1'
guard_path=HERE/'guard-fictional-pilot.py'
if hashlib.sha256(guard_path.read_bytes()).hexdigest()!=GUARD_SHA:raise ValueError('qualified guard differs')
spec=importlib.util.spec_from_file_location('qualified_public_guard',guard_path);g=importlib.util.module_from_spec(spec);spec.loader.exec_module(g)

def bind(path,expected):
    g.require(g.sha(path)==expected,'authorization SHA differs');a=g.read_json(path)
    g.require(set(a)=={'schema_version','public_transform_authorized','patient_assay_authorized','plan','source_sha256','output','launch','limits'},'authorization fields differ')
    g.require(type(a['schema_version']) is int and a['schema_version']==1,'authorization schema')
    g.require(a['public_transform_authorized'] is True and a['patient_assay_authorized'] is False,'public-only authorization required')
    g.require(a['output']==str(OUT) and a['launch']==str(LAUNCH) and a['limits']==g.LIMITS,'scope or caps differ')
    g.require(a['plan']['sha256']==PLAN_SHA and g.sha(a['plan']['path'])==PLAN_SHA,'unchanged reviewed resource plan required')
    plan=g.read_json(a['plan']['path']);g.require(plan['case_execution_authorized'] is False and plan['ready_for_native_pilot'] is False,'resource proposal changed')
    g.require(set(a['source_sha256'])=={'guard-public-panel.py','guard-fictional-pilot.py','prepare-public-panel.py','verify-public-panel.py'},'source inventory differs')
    g.require(a['source_sha256']['prepare-public-panel.py']==TRANSFORM_SHA and a['source_sha256']['guard-fictional-pilot.py']==GUARD_SHA,'unchanged transformer/guard required')
    for name,digest in a['source_sha256'].items():g.require(g.sha(HERE/name)==digest,'source changed '+name)
    return a

def worker(a):
    g.require(os.environ.get('MVA_PUBLIC_PANEL_PARENT')==str(os.getppid()) and os.getpgrp()==os.getpid(),'owned guard parent required')
    OUT.mkdir(parents=True,exist_ok=False);logs=OUT/'commands';logs.mkdir();commands=[]
    try:
        common=['--plan',a['plan']['path'],'--expected-plan-sha256',PLAN_SHA]
        for name,arguments in [('prepare-public-panel.py',common+['--output',str(OUT/'derivative')]),('verify-public-panel.py',common+['--run',str(OUT/'derivative'),'--output',str(OUT/'verification')])]:
            command=[sys.executable,'-B',str(HERE/name)]+arguments
            with (logs/(name+'.stdout.log')).open('xb') as stdout,(logs/(name+'.stderr.log')).open('xb') as stderr:code=subprocess.run(command,stdout=stdout,stderr=stderr,timeout=a['limits']['per_command_wall_seconds']).returncode
            commands.append(dict(command=command,exit_code=code));g.save(OUT/'commands.json',commands);g.require(code==0,'public step failed '+name)
        result=dict(status='complete',completed_at=g.stamp(),summary_sha256=g.sha(OUT/'derivative/summary.json'),verification_sha256=g.sha(OUT/'verification/verification.json'),patient_data_accessed=False,estimator_executed=False)
        g.save(OUT/'status.json',result)
    except BaseException as e:g.save(OUT/'status.json',dict(status='failed',failed_at=g.stamp(),error=str(e),commands=commands));raise

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--authorization',type=Path,required=True);p.add_argument('--expected-authorization-sha256',required=True);p.add_argument('--worker',action='store_true');a=p.parse_args()
    g.require(platform.node().split('.')[0]=='512S2','Mac Studio512S2 required');authorization=bind(a.authorization,a.expected_authorization_sha256)
    if a.worker:worker(authorization);return
    os.nice(10);before=g.caller_identity();os.environ['MVA_PUBLIC_PANEL_PARENT']=str(os.getpid());state=None
    try:
        state=g.supervise([sys.executable,'-B',str(Path(__file__).resolve()),'--authorization',str(a.authorization.resolve()),'--expected-authorization-sha256',a.expected_authorization_sha256,'--worker'],OUT,LAUNCH,g.LIMITS,ROOT)
        terminal=g.read_json(OUT/'status.json');g.require(terminal['status']=='complete','public terminal incomplete')
        g.require(g.sha(OUT/'derivative/summary.json')==terminal['summary_sha256'] and g.sha(OUT/'verification/verification.json')==terminal['verification_sha256'],'public terminal links differ')
        after=g.caller_identity();g.exact_identity(before,after);state.update(original_caller_before=before,original_caller_after=after,authorization_sha256=a.expected_authorization_sha256,source_sha256=g.sha(__file__),result_status_sha256=g.sha(OUT/'status.json'))
    except BaseException as e:
        if LAUNCH.exists():state=g.read_json(LAUNCH/'status.json');state.update(status='failed',error=str(e),original_caller_before=before)
        raise
    finally:
        if state is not None:g.save(LAUNCH/'status.json',state)
if __name__=='__main__':main()
