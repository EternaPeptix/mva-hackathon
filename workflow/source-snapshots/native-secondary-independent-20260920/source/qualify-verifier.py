#!/usr/bin/env python3
"""Run only the independent verifier's fictional controls in a fresh directory."""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--producer-qualification',type=Path,required=True)
    p.add_argument('--expected-producer-qualification-sha256',required=True)
    a=p.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    q=a.producer_qualification.resolve()
    if sha(q)!=a.expected_producer_qualification_sha256:raise ValueError('producer qualification hash')
    qdata=json.loads(q.read_text())
    if qdata['status']!='complete' or any(qdata[k]!=0 for k in ['errors','failures','exit_code']) or qdata['skipped']:raise ValueError('producer qualification unsuccessful')
    source=out/'source';source.mkdir()
    pins={}
    for name in ['verify-secondary.py','test_verify_secondary.py','qualify-verifier.py']:
        original=Path(__file__).resolve().parent/name
        shutil.copyfile(original,source/name);pins[name]=dict(bytes=(source/name).stat().st_size,sha256=sha(source/name))
    env=os.environ.copy();env['MVA_SECONDARY_QUALIFICATION']=str(q);env['MVA_SECONDARY_CONTROL_OUTPUT']=str(out/'fictional-cases')
    started=time.monotonic()
    command=[sys.executable,'-u','-m','unittest','discover','-s',str(source),'-p','test_verify_secondary.py','-v']
    with (out/'stdout.log').open('xb') as stdout,(out/'stderr.log').open('xb') as stderr:
        result=subprocess.run(command,env=env,stdout=stdout,stderr=stderr)
    log=(out/'stderr.log').read_text()
    import re
    count=re.findall(r'^Ran (\d+) tests? in ',log,re.M)
    skips=re.findall(r'^.*\.\.\. skipped .*$',log,re.M)
    errors=len(re.findall(r'^ERROR: ',log,re.M));failures=len(re.findall(r'^FAIL: ',log,re.M))
    complete=result.returncode==0 and len(count)==1 and int(count[0])>=36 and not skips and not errors and not failures
    data=dict(status='complete' if complete else 'failed',completed_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),elapsed_seconds=time.monotonic()-started,tests_run=int(count[0]) if len(count)==1 else None,skipped=skips,errors=errors,failures=failures,exit_code=result.returncode,source_pins=pins,producer_qualification=dict(path=str(q),sha256=sha(q)),command=command,runtime=dict(host=platform.node(),machine=platform.machine(),system=platform.system(),python=sys.version,executable=sys.executable),stdout_sha256=sha(out/'stdout.log'),stderr_sha256=sha(out/'stderr.log'),scope='Fictional fixture and intentional corruptions only. No patient alignments, genotypes, case screen or audit.')
    (out/'controls.json').write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:data[k] for k in ['status','tests_run','errors','failures','skipped','elapsed_seconds']}))
    return 0 if complete else 1


if __name__=='__main__':sys.exit(main())
