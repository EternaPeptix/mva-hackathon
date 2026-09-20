#!/usr/bin/env python3
"""Finite fictional public-transform and guard controls; no actual public/case input."""
import argparse,datetime,hashlib,json,os,platform,re,shutil,subprocess,sys,time
from pathlib import Path
NAMES=['prepare-public-panel.py','verify-public-panel.py','guard-public-panel.py','guard-fictional-pilot.py','test_public_panel.py','test_pilot_guard.py','qualify-public-preparation.py']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);src=out/'source';src.mkdir();here=Path(__file__).resolve().parent
    for name in NAMES:shutil.copyfile(here/name,src/name)
    environment=os.environ.copy();environment['MVA_PUBLIC_FIXTURE_ROOT']=str(out/'fictional-transform')
    command=[sys.executable,'-B','-m','unittest','discover','-s',str(src),'-p','test*.py','-v'];start=time.monotonic()
    with (out/'stdout.log').open('xb') as stdout,(out/'stderr.log').open('xb') as stderr:r=subprocess.run(command,stdout=stdout,stderr=stderr,env=environment,timeout=120)
    log=(out/'stderr.log').read_text();count=re.findall(r'^Ran (\d+) tests? in ',log,re.M);errors=len(re.findall(r'^ERROR:',log,re.M));failures=len(re.findall(r'^FAIL:',log,re.M));skips=re.findall(r'^.*\.\.\. skipped .*$',log,re.M)
    good=r.returncode==0 and count==['25'] and errors==failures==0 and not skips
    receipt=dict(status='complete' if good else 'failed',completed_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),tests_run=int(count[0]) if len(count)==1 else None,errors=errors,failures=failures,skipped=skips,exit_code=r.returncode,elapsed_seconds=time.monotonic()-start,source_sha256={n:sha(src/n) for n in NAMES},stdout_sha256=sha(out/'stdout.log'),stderr_sha256=sha(out/'stderr.log'),runtime=dict(host=platform.node(),machine=platform.machine(),python=sys.version),command=command,scope='Generated fictional 224487-record public transform, independent per-original-row reconciliation and owned-process guard controls only. No actual public resource, GATK or patient input accessed.',fictional_files=[dict(path=str(x.relative_to(out)),bytes=x.stat().st_size,sha256=sha(x)) for x in sorted((out/'fictional-transform').rglob('*')) if x.is_file()])
    with (out/'controls.json').open('x') as f:f.write(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:receipt[k] for k in ['status','tests_run','errors','failures','skipped','elapsed_seconds']}));return 0 if good else 1
if __name__=='__main__':raise SystemExit(main())
