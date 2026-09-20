#!/usr/bin/env python3
"""One authorized immutable producer invocation with a dedicated guarded group."""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import signal
import subprocess
import time

LIMITS = dict(wall_seconds=5400, group_RSS_bytes=16*1024**3,
              output_bytes=20*1024**3, minimum_free_bytes=40*1024**3)


def require(ok, message):
    if not ok: raise ValueError(message)


def parse(text):
    def pairs(items):
        result = {}
        for key, value in items:
            require(key not in result, 'duplicate JSON key')
            result[key] = value
        return result
    def bad(value): raise ValueError('nonfinite JSON '+value)
    return json.loads(text, object_pairs_hook=pairs, parse_constant=bad)


def stamp(): return datetime.datetime.now(datetime.timezone.utc).isoformat()


def save(path, value):
    temporary=path.with_name(path.name+'.tmp')
    temporary.write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n')
    temporary.replace(path)


def fingerprint(path, expected=None):
    path=Path(path).resolve();before=path.stat();digest=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):digest.update(block)
    after=path.stat()
    require((before.st_size,before.st_mtime_ns)==(after.st_size,after.st_mtime_ns),'file changed during hash')
    result=dict(path=str(path),bytes=after.st_size,sha256=digest.hexdigest())
    if expected:
        if 'path' in expected: require(path==Path(expected['path']).resolve(),'pin path differs')
        for key in ('bytes','sha256'):
            if key in expected:require(result[key]==expected[key],'pin '+key+' differs')
    return result


def approved_delta(draft, approved):
    require(draft.get('ready') is False and approved.get('ready') is True,'approval ready transition')
    require(isinstance(approved.get('preparation_note'),str) and bool(approved['preparation_note']),'approval note missing')
    def remaining(value):return {k:v for k,v in value.items() if k not in ('ready','preparation_note')}
    require(json.dumps(remaining(draft),sort_keys=True,allow_nan=False)==json.dumps(remaining(approved),sort_keys=True,allow_nan=False),'unapproved configuration change')


def qualification_ok(controls, proposal):
    require(controls['status']=='complete' and type(controls['tests_run']) is int and controls['tests_run']==32,'wrong source qualification')
    require(controls['skipped']==[] and all(type(controls[k]) is int and controls[k]==0 for k in ('errors','failures','exit_code')),'qualification not wholly passing')
    require(set(controls['source_sha256'])==set(proposal['qualified_sources']),'qualification source membership differs')
    for name,pin in proposal['qualified_sources'].items():require(pin['sha256']==controls['source_sha256'][name],'qualified source differs')


def snapshot():
    rows=[]
    for line in subprocess.check_output(['ps','-axo','pid=,ppid=,pgid=,stat=,rss=,%cpu=,comm='],text=True).splitlines():
        part=line.split(None,6)
        if len(part)==7:
            rows.append(dict(pid=int(part[0]),ppid=int(part[1]),pgid=int(part[2]),state=part[3],RSS_bytes=1024*int(part[4]),cpu_percent=float(part[5]),executable=part[6]))
    return rows


def owned_rss(rows, group):return sum(r['RSS_bytes'] for r in rows if r['pgid']==group)


def tree_bytes(path):
    if not path.exists():return 0
    total=0
    for member in path.rglob('*'):
        require(not member.is_symlink(),'unexpected output symlink')
        if member.is_file():total+=member.stat().st_size
    return total


def limit_reason(elapsed, rss, size, free, limits=LIMITS):
    for field,value in [('wall_seconds',elapsed),('group_RSS_bytes',rss),('output_bytes',size)]:
        if value>=limits[field]:return field
    if free<limits['minimum_free_bytes']:return 'minimum_free_bytes'
    return None


def stop_owned(process):
    require(process.pid>1 and process.pid!=os.getpgrp(),'refuse unsafe process group')
    # Only the dedicated session created by this launcher is ever signalled.
    try:os.killpg(process.pid,signal.SIGTERM)
    except ProcessLookupError:pass
    deadline=time.monotonic()+5
    while time.monotonic()<deadline:
        process.poll()  # Reap our exited child before checking group survivors.
        if not any(r['pgid']==process.pid for r in snapshot()):break
        time.sleep(0.1)
    if any(r['pgid']==process.pid for r in snapshot()):
        try:os.killpg(process.pid,signal.SIGKILL)
        except ProcessLookupError:pass
    try:process.wait(timeout=5)
    except subprocess.TimeoutExpired:raise RuntimeError('owned process group did not exit')


def bind_plan(proposal_path, proposal_sha, approved_path, approved_sha):
    pins={}
    def pin(path, expected):
        value=fingerprint(path,expected);pins[value['path']]=value;return value
    pin(proposal_path,{'sha256':proposal_sha});proposal=parse(proposal_path.read_text())
    require(proposal['ready'] is False and proposal['case_screen_executed'] is False,'proposal already used')
    for name in ('config_draft','source_qualification','qualified_entrypoint','python','input_preflight'):
        value=proposal[name];pin(value['path'],value)
    for value in proposal['qualified_sources'].values():pin(value['path'],value)
    controls=parse(Path(proposal['source_qualification']['path']).read_text());qualification_ok(controls,proposal)
    require(proposal['qualified_entrypoint']==proposal['qualified_sources']['scripts/review_native_secondary.py'],'entrypoint not qualified')
    pin(approved_path,{'sha256':approved_sha})
    approved=parse(approved_path.read_text());draft=parse(Path(proposal['config_draft']['path']).read_text());approved_delta(draft,approved)
    require(approved['fixture_only'] is False,'real authorized configuration required')
    root=Path(approved['root']).resolve()
    output=Path(proposal['proposed_output']).resolve()
    require(output==root/'results/native-secondary-review-20260920/review-v1','wrong approved output')
    require(not output.exists(),'existing result preserved')
    command=[proposal['python']['path'],'-B',proposal['qualified_entrypoint']['path'],'--config',str(approved_path),'--expected-config-sha256',approved_sha,'--output',str(output)]
    return proposal,approved,root,output,command,pins


def run(args):
    require(platform.node().split('.')[0]=='512S2','512S2 only')
    source=fingerprint(__file__,{'sha256':args.expected_launcher_sha256})
    proposal,config,root,output,command,pins=bind_plan(args.proposal,args.expected_proposal_sha256,args.config,args.expected_config_sha256)
    method=root/'methods/native-secondary-review-20260920'
    launch=method/'launch-v1';require(not launch.exists(),'existing launch preserved; no retry')
    require(Path(__file__).resolve().is_relative_to(method),'launcher outside owned methods')
    free=shutil.disk_usage(root).free;require(free>=LIMITS['minimum_free_bytes'],'insufficient disk before start')
    before=snapshot();caller=[r for r in before if r['pid']==79898]
    require(len(caller)==1,'original caller is not live; review required')
    caller_command=subprocess.check_output(['ps','-p','79898','-o','command='],text=True).strip()
    require('/workflow/scripts/call_genome.py' in caller_command and str(root/'calls/from-fastq-cram30-20260920') in caller_command,'original caller identity differs')
    launch.mkdir();state=dict(status='launching',started_at=stamp(),source=source,bindings=pins,command=command,
        limits=LIMITS,guard_sampling_seconds=2,caller_before=caller[0],caller_command_before=caller_command,
        caller_children_before=[r for r in before if r['ppid']==79898],free_disk_before=free,
        other_busy_processes_before=[r for r in before if r['cpu_percent']>=50 and r['pid']!=79898 and r['ppid']!=79898],
        scope='One authorized secondary inventory; no genotyper, alignment read, download or publication.')
    save(launch/'status.json',state);process=None;started=time.monotonic()
    try:
        with (launch/'producer.stdout.log').open('xb') as stdout,(launch/'producer.stderr.log').open('xb') as stderr:
            process=subprocess.Popen(command,stdout=stdout,stderr=stderr,start_new_session=True)
            state.update(status='running',producer_pid=process.pid,producer_pgid=process.pid,launched_at=stamp())
            save(launch/'status.json',state)
            peaks=dict(group_RSS_bytes=0,output_bytes=0)
            while True:
                rows=snapshot();rss=owned_rss(rows,process.pid);size=tree_bytes(output)+tree_bytes(launch)
                free=shutil.disk_usage(root).free;elapsed=time.monotonic()-started
                peaks={k:max(peaks[k],v) for k,v in [('group_RSS_bytes',rss),('output_bytes',size)]}
                reason=limit_reason(elapsed,rss,size,free)
                state.update(checked_at=stamp(),elapsed_seconds=elapsed,group_RSS_bytes=rss,output_bytes=size,free_disk_bytes=free,peaks=peaks)
                if reason:state['limit_triggered']=reason;raise RuntimeError('guard limit '+reason)
                code=process.poll()
                if code is not None:
                    require(not any(r['pgid']==process.pid for r in rows if r['pid']!=process.pid),'producer exited with owned children still live')
                    require(code==0,'producer exit '+str(code));break
                save(launch/'status.json',state);time.sleep(2)
        terminal=parse((output/'status.json').read_text());require(terminal['status']=='complete','producer terminal status incomplete')
        fingerprint(output/'summary.json',{'sha256':terminal['summary_sha256']})
        state.update(status='complete',completed_at=stamp(),exit_code=0,producer_status=fingerprint(output/'status.json'),producer_summary=fingerprint(output/'summary.json'))
    except BaseException as error:
        if process is not None:stop_owned(process)
        state.update(status='failed',failed_at=stamp(),error_type=type(error).__name__,error=str(error),exit_code=process.returncode if process else None)
        raise
    finally:
        after=snapshot();state['caller_after']=[r for r in after if r['pid']==79898];state['caller_children_after']=[r for r in after if r['ppid']==79898]
        state['free_disk_after']=shutil.disk_usage(root).free
        save(launch/'status.json',state)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('proposal','config'):parser.add_argument('--'+name,type=Path,required=True)
    for name in ('proposal','config','launcher'):parser.add_argument('--expected-'+name+'-sha256',required=True)
    run(parser.parse_args())
