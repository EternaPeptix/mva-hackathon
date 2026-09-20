#!/usr/bin/env python3
"""Guard only an owned saved-depth mapping or fictional qualification process."""
import argparse
import importlib.util
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parent/'vendor'))
from guard_helpers import (fingerprint, limit_reason, owned_rss, parse, require, save,
                           snapshot, stamp, stop_owned, tree_bytes)

LIMITS=dict(wall_seconds=900,group_RSS_bytes=4*1024**3,output_bytes=512*1024**2,minimum_free_bytes=40*1024**3)
SOURCES=('scripts/map_depth_gaps.py','scripts/guard_gap_map.py','scripts/vendor/coding_depth_targets.py',
         'scripts/vendor/guard_helpers.py','tests/test_depth_gap_map.py','tests/test_guard_gap_map.py')


def sources():
    base=Path(__file__).resolve().parents[1]
    return {name:fingerprint(base/name) for name in SOURCES}


def qualification_ok(receipt, current):
    require(receipt.get('status')=='complete' and receipt.get('fictional_only') is True,'qualification incomplete')
    require(type(receipt.get('tests_run')) is int and receipt['tests_run']==39,'qualification count')
    require(receipt.get('skipped')==[] and all(type(receipt.get(k)) is int and receipt[k]==0 for k in ('errors','failures','exit_code')),
            'qualification not wholly passing')
    require(set(receipt['source_files'])==set(SOURCES),'qualified source membership')
    for name,pin in current.items():
        require(receipt['source_files'][name]['sha256']==pin['sha256'] and receipt['source_files'][name]['bytes']==pin['bytes'],
                'qualified source differs')


def worker(args):
    os.nice(10)
    pins=sources()
    if args.mode=='qualification':
        output=args.output/'controls.json'
        suite=unittest.defaultTestLoader.discover(str(Path(__file__).resolve().parents[1]/'tests'),pattern='test_*gap_map.py')
        result=unittest.TextTestRunner(verbosity=2).run(suite)
        for pin in pins.values():fingerprint(pin['path'],pin)
        save(output,dict(status='complete' if result.wasSuccessful() else 'failed',completed_at=stamp(),
             fictional_only=True,tests_run=result.testsRun,errors=len(result.errors),failures=len(result.failures),
             skipped=[x[0].id() for x in result.skipped],exit_code=0 if result.wasSuccessful() else 1,
             source_files=pins,runtime=fingerprint(sys.executable),host=platform.node(),system=platform.platform()))
        return 0 if result.wasSuccessful() else 1
    require(args.config and args.expected_config_sha256 and args.qualification and args.expected_qualification_sha256,
            'case pins absent')
    qpin=fingerprint(args.qualification,{'sha256':args.expected_qualification_sha256})
    qualification_ok(parse(Path(args.qualification).read_text()),pins)
    cpin=fingerprint(args.config,{'sha256':args.expected_config_sha256})
    config=parse(Path(args.config).read_text())
    require(config.get('ready') is True,'case not approved')
    require(Path(config['output_directory']).resolve()==args.output.with_name(args.output.name.removesuffix('-launch')),
            'case output differs')
    fingerprint(sys.executable,config['python'])
    save(args.output/'worker-binding.json',dict(config=cpin,qualification=qpin,sources=pins,runtime=fingerprint(sys.executable)))
    source=Path(__file__).with_name('map_depth_gaps.py')
    cmd=[sys.executable,str(source),'--config',str(args.config),'--expected-config-sha256',args.expected_config_sha256,
         '--output',config['output_directory']]
    result=subprocess.run(cmd)
    require(result.returncode==0,'mapper exited unsuccessfully')
    status=Path(config['output_directory'])/'status.json'
    require(parse(status.read_text())['status']=='complete','mapper terminal state')
    save(args.output/'worker-result.json',dict(status='complete',mapper_status=fingerprint(status)))
    return 0


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--mode',choices=('qualification','case'),required=True)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--config',type=Path);ap.add_argument('--expected-config-sha256')
    ap.add_argument('--qualification',type=Path);ap.add_argument('--expected-qualification-sha256')
    ap.add_argument('--worker',action='store_true');args=ap.parse_args()
    args.output=args.output.resolve()
    if args.worker:return worker(args)
    require(platform.system()=='Darwin' and platform.machine()=='arm64','Mac Studio native runtime required')
    require(subprocess.check_output(['sysctl','-n','hw.model'],text=True).strip()=='Mac15,14','512S2 hardware identity')
    require(not args.output.exists(),'existing launch preserved')
    require(args.mode=='qualification' or args.output.name.endswith('-launch'),'case launch suffix')
    if args.mode=='case':require(not args.output.with_name(args.output.name.removesuffix('-launch')).exists(),'existing result preserved')
    require(shutil.disk_usage(args.output.parent).free>=LIMITS['minimum_free_bytes'],'insufficient free disk')
    args.output.mkdir(parents=True,exist_ok=False)
    before=snapshot();state=dict(status='starting',mode=args.mode,started_at=stamp(),limits=LIMITS,
                                host=platform.node(),runtime=fingerprint(sys.executable),source_files=sources(),
                                sampling_seconds=1,limits_are_sampled=True,process_snapshot_before=before)
    save(args.output/'status.json',state)
    command=[sys.executable,str(Path(__file__).resolve()),*sys.argv[1:],'--worker']
    started=time.monotonic();peak_rss=0;peak_output=0;reason=None
    with (args.output/'stdout.txt').open('w') as stdout,(args.output/'stderr.txt').open('w') as stderr:
        process=subprocess.Popen(command,start_new_session=True,stdout=stdout,stderr=stderr)
        state.update(status='running',owned_process_group=process.pid,command=command);save(args.output/'status.json',state)
        try:
            while True:
                rss=owned_rss(snapshot(),process.pid);size=tree_bytes(args.output)
                if args.mode=='case':size+=tree_bytes(args.output.with_name(args.output.name.removesuffix('-launch')))
                peak_rss=max(peak_rss,rss);peak_output=max(peak_output,size)
                reason=limit_reason(time.monotonic()-started,rss,size,shutil.disk_usage(args.output).free,LIMITS)
                if reason:stop_owned(process);break
                if process.poll() is not None:break
                time.sleep(1)
        except BaseException:
            stop_owned(process);raise
    state.update(status='complete' if process.returncode==0 and reason is None else 'failed',
                 completed_at=stamp(),exit_code=process.returncode,limit_reason=reason,elapsed_seconds=time.monotonic()-started,
                 peak_sampled_RSS_bytes=peak_rss,peak_sampled_output_bytes=peak_output,process_snapshot_after=snapshot(),
                 outputs={p.name:fingerprint(p) for p in args.output.iterdir() if p.is_file() and p.name!='status.json'})
    save(args.output/'status.json',state)
    return 0 if state['status']=='complete' else 1


if __name__=='__main__':sys.exit(main())
