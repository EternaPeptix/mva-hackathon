#!/usr/bin/env python3
"""Finish one public benchmark after the existing native caller completes.

Observe, but never start or restart, the existing public or patient callers.
Dependencies and this sequence's source/configuration are frozen at launch.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import time

from build_structural_exclusions import fingerprint, save, stamp


def run(root, output, wait_seconds):
    if output.resolve().is_relative_to(Path(__file__).resolve().parents[2]): raise ValueError('logs stay outside repository')
    pilot=root/'results/germline-public-pilot-20260919'
    prepared=root/'ref/giab-small-variants-20260919/prepared-v2'
    controls=root/'results/germline-public-comparator-controls-20260919'
    comparison=root/'results/germline-public-comparison-20260919'
    for folder in (prepared,controls):
        if json.loads((folder/'summary.json').read_text())['status']!='complete':raise ValueError('dependency incomplete')
    initial=json.loads((pilot/'status.json').read_text())
    if initial['status'] not in ('running','complete'):raise ValueError('existing public caller failed or absent')
    if comparison.exists():raise ValueError('comparison output already exists; preserve it')
    output.mkdir(parents=True,exist_ok=False);state={'status':'waiting_public_caller','pid':os.getpid(),'started_at':stamp(),'dependency_pid':initial['pid']}
    save(output/'status.json',state)
    scripts=root/'repo/workflow/scripts';config=root/'repo/workflow/configs/public-small-variant-comparison.json'
    paths=[Path(__file__),scripts/'compare_public_small_variants.py',scripts/'rtg_small_variant_helpers.py',
           scripts/'run_public_small_variant_pilot.py',scripts/'call_genome.py',config,
           prepared/'summary.json',controls/'summary.json']
    sources={str(p):fingerprint(p) for p in paths};save(output/'sources-at-launch.json',sources)
    try:
        deadline=time.monotonic()+wait_seconds
        while True:
            try:dependency=json.loads((pilot/'status.json').read_text())
            except json.JSONDecodeError:
                if time.monotonic()>deadline:raise RuntimeError('public caller status did not settle')
                time.sleep(1);continue
            if dependency.get('pid')!=initial['pid'] or dependency.get('started_at')!=initial['started_at']:
                raise ValueError('public dependency identity changed')
            if dependency['status']=='complete':break
            if dependency['status']!='running':raise RuntimeError('public caller did not finish successfully')
            os.kill(dependency['pid'],0)
            command=subprocess.check_output(['ps','-p',str(dependency['pid']),'-o','command='],text=True)
            if str(scripts/'run_public_small_variant_pilot.py') not in command:raise ValueError('dependency PID no longer owns the expected caller')
            if time.monotonic()>deadline:raise RuntimeError('observation deadline reached; existing caller left running')
            time.sleep(15)
        for path,entry in sources.items():
            if fingerprint(Path(path))['sha256']!=entry['sha256']:raise ValueError('sequence dependency changed while waiting: '+path)
        args=list(map(str,[root/'tools/bio/bin/python',scripts/'compare_public_small_variants.py',
                          '--root',root,'--output',comparison,'--config',config,
                          '--prepared',prepared,'--controls',controls,'--pilot',pilot]))
        start=time.monotonic()
        with (output/'comparison.log').open('w') as stream:
            child=subprocess.Popen(args,stdout=stream,stderr=subprocess.STDOUT)
            state.update(status='running',stage='fixed genotype-aware public comparison',child_pid=child.pid);save(output/'status.json',state)
            exit_code=child.wait()
        command={'command':args,'exit_code':exit_code,'elapsed_seconds':time.monotonic()-start,'log':fingerprint(output/'comparison.log')}
        save(output/'command.json',command)
        if exit_code:raise RuntimeError('public comparison failed; original outputs retained')
        if json.loads((comparison/'summary.json').read_text())['status']!='complete':raise ValueError('comparison receipt incomplete')
        state.update(status='complete',completed_at=stamp());save(output/'status.json',state)
        save(output/'summary.json',{'status':'complete','state':state,'command':command,'sources':sources,'comparison':fingerprint(comparison/'summary.json')})
        return state
    except BaseException as error:
        state.update(status='failed',failed_at=stamp(),error_type=type(error).__name__,error=str(error));save(output/'status.json',state);raise


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('root','output'):parser.add_argument('--'+name,required=True,type=Path)
    parser.add_argument('--wait-seconds',type=int,default=21600)
    args=parser.parse_args()
    if args.wait_seconds<=0:parser.error('--wait-seconds must be positive')
    print(json.dumps(run(args.root,args.output,args.wait_seconds)))
