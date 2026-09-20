#!/usr/bin/env python3
"""Wait for the already-running public reference preparation, then run one pilot.

This finite sequence never starts, edits or restarts the patient controllers or
reference preparation. Any failed stage stops dependent work and keeps evidence.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import time

from build_structural_exclusions import fingerprint, save, stamp


def run(root, output, wait_seconds):
    if output.resolve().is_relative_to(Path(__file__).resolve().parents[2]): raise ValueError('sequence logs stay outside Git')
    output.mkdir(parents=True,exist_ok=False)
    state={'status':'waiting_reference','started_at':stamp(),'pid':os.getpid()};save(output/'status.json',state)
    scripts=root/'repo/workflow/scripts';config=root/'repo/workflow/configs'
    code={name:fingerprint(scripts/name) for name in ('run_svaba_public_sequence.py','run_svaba_pilot.py','convert_svaba_pilot.py','compare_svaba_public.py')}
    sources={'code':code,'svaba_config':fingerprint(config/'svaba-pilot.json'),
             'comparison_config':fingerprint(config/'public-sv-pilot-comparison.json')}
    save(output/'sources-at-launch.json',sources);commands=[]
    try:
        deadline=time.monotonic()+wait_seconds;preparation=root/'ref/giab-sv-20260919/assembly-reference/status.json'
        while True:
            try: reference_state=json.loads(preparation.read_text())
            except json.JSONDecodeError:
                if time.monotonic()>deadline: raise RuntimeError('reference status did not settle')
                time.sleep(1);continue
            if reference_state['status']=='complete':break
            if reference_state['status']!='running':raise RuntimeError('reference preparation did not complete successfully')
            os.kill(reference_state['pid'],0)
            if time.monotonic()>deadline:raise RuntimeError('reference observation deadline reached; existing preparation left running')
            time.sleep(10)
        save(output/'reference-preparation.json',json.loads((preparation.parent/'summary.json').read_text()))
        for name,identity in code.items():
            if fingerprint(scripts/name)['sha256']!=identity['sha256']:raise ValueError('sequence source changed while waiting')
        for key,name in [('svaba_config','svaba-pilot.json'),('comparison_config','public-sv-pilot-comparison.json')]:
            if fingerprint(config/name)['sha256']!=sources[key]['sha256']:raise ValueError('sequence configuration changed while waiting')
        native=root/'results/svaba-public-20260919';conversion=root/'results/svaba-public-conversion-20260919'
        comparison=root/'results/svaba-public-comparison-20260919'
        stages=[('native public assembly',[root/'tools/bio/bin/python',scripts/'run_svaba_pilot.py','--root',root,'--mode','public','--config',config/'svaba-pilot.json','--output',native]),
                ('paired-breakend conversion',[root/'tools/bio/bin/python',scripts/'convert_svaba_pilot.py','--pilot',native,'--output',conversion]),
                ('unchanged public benchmark',[root/'tools/truvari-5.5.0/bin/python',scripts/'compare_svaba_public.py','--root',root,'--output',comparison,'--config',config/'public-sv-pilot-comparison.json'])]
        for i,(label,args) in enumerate(stages):
            args=list(map(str,args));tick=time.monotonic();log=output/f'step-{i}.log'
            with log.open('w') as handle:
                child=subprocess.Popen(args,stdout=handle,stderr=subprocess.STDOUT)
                state.update(status='running',stage=label,child_pid=child.pid,command=args);save(output/'status.json',state)
                exit_code=child.wait()
            commands.append({'stage':label,'command':args,'exit_code':exit_code,'elapsed_seconds':time.monotonic()-tick,'log':fingerprint(log)})
            save(output/'commands.json',commands)
            if exit_code:raise RuntimeError(label+' failed; subsequent stages not run')
        state.update(status='complete',completed_at=stamp());save(output/'status.json',state)
        save(output/'summary.json',{'state':state,'commands':commands,'sources':sources,
                                  'public_comparison':fingerprint(comparison/'summary.json')})
        return state
    except BaseException as error:
        state.update(status='failed',failed_at=stamp(),error_type=type(error).__name__,error=str(error));save(output/'status.json',state);raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('root','output'):p.add_argument('--'+name,required=True,type=Path)
    p.add_argument('--wait-seconds',type=int,default=3600)
    a=p.parse_args()
    if a.wait_seconds<=0: p.error('--wait-seconds must be positive')
    print(json.dumps(run(a.root,a.output,a.wait_seconds)))
