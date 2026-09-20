#!/usr/bin/env python3
"""Wait for the existing four-lane job, verify its CRAM, then run native discovery callers."""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import time


def stamp():
    return time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())


def save(path, value):
    temporary=Path(str(path)+'.tmp')
    temporary.write_text(json.dumps(value,indent=2)+'\n')
    os.replace(temporary,path)


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''): h.update(block)
    return h.hexdigest()


def dependency_action(state, verified_process):
    if state.get('status')=='complete':
        if not all(state.get(k) for k in ('cram_path','cram_sha256','cram_bytes')):
            raise ValueError('completed alignment lacks its output provenance')
        return 'verify_completed_alignment'
    if state.get('status')=='running' and verified_process:
        return 'wait_for_existing_alignment'
    if state.get('status')=='running':
        raise RuntimeError('alignment handle missing; inspect existing job, do not restart it automatically')
    raise RuntimeError('alignment is terminal or unrecognised; inspect its retained outputs')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    args=parser.parse_args()
    root=args.root.resolve()
    status_path=root/'state/genome-calling.json'
    os.nice(10)
    with (root/'state/genome-calling.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        if status_path.exists():
            save(root/('state/genome-calling-history-%d.json'%time.time_ns()),json.loads(status_path.read_text()))
        state={'status':'running','pid':os.getpid(),'started_at':stamp(),'controller_sha256':digest(Path(__file__)),
               'scope':'Depend on the existing four-lane CRAM, then run initial diploid small-variant discovery across all reference contigs.',
               'current_stage':'waiting for existing alignment'}
        save(status_path,state)
        try:
            while True:
                dependency=json.loads((root/'state/raw-reads.json').read_text())
                process=False
                if dependency.get('status')=='running':
                    pid=int(dependency['pid'])
                    inspected=subprocess.run(['ps','-p',str(pid),'-o','args='],capture_output=True,text=True)
                    arguments=shlex.split(inspected.stdout)
                    expected=str(root/'repo/workflow/scripts/resume_reads.py')
                    process=inspected.returncode==0 and expected in arguments and '--root' in arguments and arguments[arguments.index('--root')+1]==str(root)
                    if not process:
                        # A job can finish between the first state read and process inspection.
                        dependency=json.loads((root/'state/raw-reads.json').read_text())
                action=dependency_action(dependency,process)
                state.update(last_checked_at=stamp(),dependency_pid=dependency.get('pid'),dependency_status=dependency.get('status'))
                save(status_path,state)
                if action=='verify_completed_alignment': break
                time.sleep(30)
            cram=Path(dependency['cram_path']).resolve()
            expected_cram=root/'align/from-fastq-20260919/WGS_EX2312012.cram'
            if cram!=expected_cram or cram.stat().st_size!=dependency['cram_bytes']:
                raise ValueError('unexpected completed alignment path or size')
            state['current_stage']='verify completed CRAM checksum'
            save(status_path,state)
            if digest(cram)!=dependency['cram_sha256']:
                raise ValueError('completed CRAM checksum mismatch')
            state['current_stage']='native whole-reference diploid discovery calls'
            save(status_path,state)
            output=root/'calls/from-fastq-20260919'
            command=[str(root/'tools/bio/bin/python'),str(root/'repo/workflow/scripts/call_genome.py'),
                     '--root',str(root),'--reference',str(root/'ref/realignment/GRCh38_nochr.fa'),
                     '--alignment',str(cram),'--sample','WGS_EX2312012','--output',str(output),'--workers','4']
            state['command']=command
            save(status_path,state)
            with (root/'logs/genome-calling.log').open('a') as log:
                subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True,pass_fds=(lock.fileno(),))
            result=json.loads((output/'status.json').read_text())
            if result.get('status')!='complete': raise RuntimeError('caller returned without verified callsets')
            state.update(status='complete',completed_at=stamp(),current_stage='raw discovery complete',
                         callset_status=str(output/'status.json'),remaining_gates=result['next_gates'])
            save(status_path,state)
        except BaseException as error:
            state.update(status='failed',failed_at=stamp(),error_type=type(error).__name__)
            save(status_path,state)
            raise


if __name__=='__main__':
    main()
