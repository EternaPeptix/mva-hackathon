"""Run one reviewed saved-public-VCF analysis in an isolated, bounded group."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--expected-plan-sha256', required=True)
    parser.add_argument('--expected-guard-sha256', required=True)
    args = parser.parse_args()
    import hashlib
    digest = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
    if digest(__file__) != args.expected_guard_sha256 or digest(args.plan) != args.expected_plan_sha256:
        raise ValueError('Launch source or plan differs')
    plan = json.loads(args.plan.read_text())
    if plan['ready'] is not True:
        raise ValueError('Execution plan is disabled')
    if platform.node().split('.')[0] != '512S2' or platform.machine() != 'arm64':
        raise ValueError('Wrong host')
    if subprocess.check_output(['sysctl','-n','hw.model'],text=True).strip() != 'Mac15,14':
        raise ValueError('Wrong Mac Studio model')
    helper = plan['guard_helper']
    if digest(helper['path']) != helper['sha256']:
        raise ValueError('Guard helper changed')
    spec = importlib.util.spec_from_file_location('saved_review_helpers',helper['path'])
    h = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(h)
    for pin in plan['pins']:
        h.fingerprint(pin['path'],pin)
    root = Path(plan['root']).resolve()
    output,launch = (Path(plan[k]).resolve() for k in ('output','launch'))
    h.require(output.is_relative_to(root/'results') and launch.is_relative_to(root/'methods'),'Outputs outside project')
    h.require(not output.exists() and not launch.exists(),'Existing result or launch preserved; no automatic retry')
    h.require(plan['command'].count('--output') == 1,'Exactly one producer output required')
    oi=plan['command'].index('--output')
    h.require(oi+1<len(plan['command']) and Path(plan['command'][oi+1]).resolve()==output,'Producer output differs from guarded output')
    limits = plan['limits']
    h.require(limits == dict(wall_seconds=300,group_RSS_bytes=2*1024**3,output_bytes=64*1024**2,minimum_free_bytes=40*1024**3),'Unexpected resource limits')
    h.require(shutil.disk_usage(root).free >= limits['minimum_free_bytes'],'Insufficient free space')
    command = ['/usr/bin/nice','-n','10'] + plan['command']
    launch.mkdir(parents=True,exist_ok=False)
    state = dict(status='launching',started_at=h.stamp(),plan=h.fingerprint(args.plan),guard=h.fingerprint(__file__),
                 command=command,limits=limits,sampling_seconds=0.25,scope=plan['scope'],
                 busy_processes_before=[r for r in h.snapshot() if r['cpu_percent']>=50])
    h.save(launch/'status.json',state)
    process=None
    tick=time.monotonic()
    try:
        with (launch/'stdout.log').open('xb') as out,(launch/'stderr.log').open('xb') as err:
            process=subprocess.Popen(command,stdout=out,stderr=err,start_new_session=True)
            state.update(status='running',producer_pid=process.pid,producer_pgid=process.pid)
            peaks=dict(group_RSS_bytes=0,output_bytes=0)
            while True:
                rows=h.snapshot();rss=h.owned_rss(rows,process.pid)
                size=h.tree_bytes(output)+h.tree_bytes(launch)
                elapsed=time.monotonic()-tick;free=shutil.disk_usage(root).free
                peaks['group_RSS_bytes']=max(peaks['group_RSS_bytes'],rss)
                peaks['output_bytes']=max(peaks['output_bytes'],size)
                state.update(checked_at=h.stamp(),elapsed_seconds=elapsed,group_RSS_bytes=rss,output_bytes=size,free_disk_bytes=free,peaks=peaks)
                reason=h.limit_reason(elapsed,rss,size,free,limits)
                h.require(reason is None,'Guard limit '+str(reason))
                code=process.poll()
                if code is not None:
                    h.require(code==0,'Producer exit '+str(code))
                    h.require(not any(r['pgid']==process.pid and r['pid']!=process.pid for r in rows),'Owned children remain')
                    break
                h.save(launch/'status.json',state);time.sleep(0.25)
        terminal=h.parse((output/'status.json').read_text())
        h.require(terminal['status']=='complete','Producer did not complete')
        for pin in plan['pins']:
            h.fingerprint(pin['path'],pin)
        h.fingerprint(args.plan,{'sha256':args.expected_plan_sha256})
        h.fingerprint(__file__,{'sha256':args.expected_guard_sha256})
        h.fingerprint(helper['path'],helper)
        state.update(status='complete',completed_at=h.stamp(),exit_code=0,producer_status=h.fingerprint(output/'status.json'),producer_summary=h.fingerprint(output/'summary.json'),all_pins_stable=True)
    except BaseException as error:
        if process is not None:h.stop_owned(process)
        state.update(status='failed',failed_at=h.stamp(),error_type=type(error).__name__,error=str(error),exit_code=process.returncode if process else None)
        raise
    finally:
        state['busy_processes_after']=[r for r in h.snapshot() if r['cpu_percent']>=50]
        state['free_disk_after']=shutil.disk_usage(root).free
        h.save(launch/'status.json',state)
    print(json.dumps({k:state[k] for k in ('status','exit_code','elapsed_seconds','peaks')}))


if __name__=='__main__':
    main()
