#!/usr/bin/env python3
"""One separately authorized read-only independent depth audit; bounded owned group."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time

SOURCE=Path(__file__).resolve()
HELPER=SOURCE.with_name('guard_helpers.py')
HELPER_SHA='e30e5ad4819d488e2509feffff1027bab6102e473bb44ed33bcf4dbfc734ffb6'
if hashlib.sha256(HELPER.read_bytes()).hexdigest()!=HELPER_SHA:raise ValueError('guard helper changed')
spec=importlib.util.spec_from_file_location('depth_audit_guard',HELPER)
g=importlib.util.module_from_spec(spec);spec.loader.exec_module(g)
ROOT=Path('/Users/jeweled/Documents/mva-hackathon')
METHOD=ROOT/'methods/secondary-panel-depth-preparation-20260920/independent-verification-preparation'
QUAL=METHOD/'qualification-v3'
QUAL_SHA='a0603fb9cf8fa064ad384127d610c9c088ee62f7abc77d4c2e853f4d5bb3c568'
VERIFIER_SHA='dac844cd69fb16c4431a6ac907f3dcef72b5a24e6d31f2c29156b6a43e804ec0'
PRODUCER=ROOT/'results/secondary-panel-cds-depth-20260920/review-v2'
OUTPUT=ROOT/'results/secondary-panel-cds-depth-20260920/independent-verification/actual-v1'
LIMITS=dict(wall_seconds=5400,group_RSS_bytes=16*1024**3,output_bytes=2*1024**3,minimum_free_bytes=40*1024**3)


def same(a,b,message):g.require(json.dumps(a,sort_keys=True,allow_nan=False)==json.dumps(b,sort_keys=True,allow_nan=False),message)


def complete_control(q,count):
    same(q['status'],'complete','qualification incomplete')
    same(q['tests_run'],count,'qualification count differs')
    for k in ('errors','failures','skipped'):same(q[k],0,'qualification '+k+' differs')


def fresh_output(output=OUTPUT):
    same(str(output.resolve()),str(OUTPUT),'wrong owned output')
    g.require(not output.exists() and not output.with_name(output.name+'-launch').exists(),'existing audit preserved')


def binding(args):
    same(platform.node().split('.')[0],'512S2','Mac Studio512S2 only')
    source=g.fingerprint(SOURCE,{'sha256':args.expected_source_sha256})
    controlpin=g.fingerprint(args.controls,{'sha256':args.expected_controls_sha256})
    controls=g.parse(args.controls.read_text());complete_control(controls,8)
    same(sorted(controls['source_sha256']),['guard-verification.py','guard_helpers.py','test_guard_verification.py'],'guard source membership differs')
    for name,sha in controls['source_sha256'].items():g.fingerprint(SOURCE.parent/name,{'sha256':sha})
    qualpin=g.fingerprint(QUAL/'controls.json',{'sha256':QUAL_SHA})
    qual=g.parse((QUAL/'controls.json').read_text());complete_control(qual,40)
    same(sorted(qual['source_sha256']),['qualify.py','test_verify_depth.py','verify-depth.py'],'verifier source membership differs')
    for name,sha in qual['source_sha256'].items():g.fingerprint(QUAL/name,{'sha256':sha})
    verifier=g.fingerprint(QUAL/'verify-depth.py',{'sha256':VERIFIER_SHA})
    runtime=g.fingerprint(Path(sys.executable).resolve(),qual['runtime'])
    launch=PRODUCER.with_name(PRODUCER.name+'-launch')
    producerpin=g.fingerprint(launch/'status.json',{'sha256':args.expected_launch_status_sha256})
    st=g.parse((launch/'status.json').read_text());same(st['status'],'complete','producer launch incomplete');same(st['exit_code'],0,'producer launch failed')
    workerpin=g.fingerprint(launch/'worker-binding.json',st['worker_binding']);worker=g.parse((launch/'worker-binding.json').read_text());same(worker['status'],'complete','producer binding incomplete')
    g.fingerprint(PRODUCER/'status.json',st['producer_status'])
    command=['/usr/bin/nice','-n','10',runtime['path'],'-B',verifier['path'],'--root',str(ROOT),'--output',str(OUTPUT),'--expected-launch-status-sha256',producerpin['sha256']]
    return command,dict(launcher=source,guard_controls=controlpin,verifier=verifier,verifier_controls=qualpin,runtime=runtime,producer_launch=producerpin,producer_worker=workerpin)


def run(args):
    fresh_output();command,pins=binding(args)
    g.require(shutil.disk_usage(ROOT).free>=LIMITS['minimum_free_bytes'],'insufficient free space')
    OUTPUT.parent.mkdir(parents=True,exist_ok=True)
    launch=OUTPUT.with_name(OUTPUT.name+'-launch');launch.mkdir()
    state=dict(status='launching',started_at=g.stamp(),inputs=pins,command=command,limits=LIMITS,guard_sampling_seconds=2,before_processes=g.snapshot())
    g.save(launch/'status.json',state)
    process=None;start=time.monotonic();peaks=dict(group_RSS_bytes=0,output_bytes=0)
    try:
        with (launch/'stdout.txt').open('xb') as stdout,(launch/'stderr.txt').open('xb') as stderr:
            process=subprocess.Popen(command,stdin=subprocess.DEVNULL,stdout=stdout,stderr=stderr,start_new_session=True)
            state.update(status='running',pid=process.pid,pgid=process.pid)
            while True:
                rows=g.snapshot();rss=g.owned_rss(rows,process.pid);size=g.tree_bytes(OUTPUT)+g.tree_bytes(launch);free=shutil.disk_usage(ROOT).free;elapsed=time.monotonic()-start
                peaks=dict(group_RSS_bytes=max(peaks['group_RSS_bytes'],rss),output_bytes=max(peaks['output_bytes'],size))
                state.update(elapsed_seconds=elapsed,group_RSS_bytes=rss,output_bytes=size,free_bytes=free,peaks=peaks)
                reason=g.limit_reason(elapsed,rss,size,free,LIMITS);g.require(reason is None,'guard limit '+str(reason))
                code=process.poll()
                if code is not None:
                    g.require(not any(r['pgid']==process.pid and r['pid']!=process.pid for r in rows),'owned child remains')
                    g.require(code==0,'verifier exit '+str(code));break
                g.save(launch/'status.json',state);time.sleep(2)
        result=g.parse((OUTPUT/'status.json').read_text());same(result['status'],'complete','audit incomplete')
        receipt=g.fingerprint(OUTPUT/'verification.json',{'sha256':result['verification_sha256']})
        for pin in pins.values():g.fingerprint(pin['path'],pin)
        state.update(status='complete',completed_at=g.stamp(),exit_code=0,result_status=g.fingerprint(OUTPUT/'status.json'),verification=receipt)
    except BaseException as error:
        if process is not None:g.stop_owned(process)
        state.update(status='failed',failed_at=g.stamp(),error=str(error),exit_code=process.returncode if process else None)
        raise
    finally:
        state['after_processes']=g.snapshot();g.save(launch/'status.json',state)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--controls',type=Path,required=True)
    for n in ('source','controls','launch-status'):p.add_argument('--expected-'+n+'-sha256',required=True)
    run(p.parse_args())
