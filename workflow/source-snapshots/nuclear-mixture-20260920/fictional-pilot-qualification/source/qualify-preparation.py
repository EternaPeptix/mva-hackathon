#!/usr/bin/env python3
"""Run only finite pure/owned-process tests; never execute the GATK pilot."""
import argparse,datetime,hashlib,json,os,platform,re,shutil,subprocess,sys,time
from pathlib import Path
NAMES=['guard-fictional-pilot.py','run-fictional-pilot.py','fictional_controls.py','test_pilot_controls.py','test_pilot_guard.py']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);src=out/'source';src.mkdir();here=Path(__file__).resolve().parent
 for name in NAMES+['qualify-preparation.py']:shutil.copyfile(here/name,src/name)
 command=[sys.executable,'-B','-m','unittest','discover','-s',str(src),'-p','test_pilot*.py','-v'];started=time.monotonic()
 environment=os.environ.copy();environment['MVA_PILOT_CONSTRUCTOR_OUTPUT']=str(out/'fictional-constructor')
 with (out/'stdout.log').open('xb') as stdout,(out/'stderr.log').open('xb') as stderr:r=subprocess.run(command,stdout=stdout,stderr=stderr,timeout=120,env=environment)
 log=(out/'stderr.log').read_text();count=re.findall(r'^Ran (\d+) tests? in ',log,re.M);errors=len(re.findall(r'^ERROR:',log,re.M));failures=len(re.findall(r'^FAIL:',log,re.M));skips=re.findall(r'^.*\.\.\. skipped .*$',log,re.M)
 good=r.returncode==0 and count==['25'] and errors==failures==0 and not skips
 receipt=dict(status='complete' if good else 'failed',completed_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),tests_run=int(count[0]) if len(count)==1 else None,errors=errors,failures=failures,skipped=skips,exit_code=r.returncode,elapsed_seconds=time.monotonic()-started,source_sha256={n:sha(src/n) for n in NAMES+['qualify-preparation.py']},stdout_sha256=sha(out/'stdout.log'),stderr_sha256=sha(out/'stderr.log'),runtime=dict(host=platform.node(),machine=platform.machine(),python=sys.version),command=command,scope='Constructor/decision and owned fictional process-group controls plus generated fictional BAM/CRAM3.0 full-payload roundtrip. No GATK invocation, public derivative or patient input access.',fictional_constructor_files=[dict(path=str(x.relative_to(out)),bytes=x.stat().st_size,sha256=sha(x)) for x in sorted((out/'fictional-constructor').rglob('*')) if x.is_file()])
 with (out/'controls.json').open('x') as f:f.write(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
 print(json.dumps({k:receipt[k] for k in ['status','tests_run','errors','failures','skipped','elapsed_seconds']}));return 0 if good else 1
if __name__=='__main__':raise SystemExit(main())
