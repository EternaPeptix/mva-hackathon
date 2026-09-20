#!/usr/bin/env python3
"""Bounded native fictional qualification only; no patient inputs."""
from pathlib import Path
import datetime,hashlib,io,json,os,platform,sys,time,unittest
ROOT=Path('/Users/jeweled/Documents/mva-hackathon')
DIRECTORY=Path(__file__).resolve().parent

def pin(p):
 s=p.stat();return dict(path=str(p),bytes=s.st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())

if __name__=='__main__':
 import signal,pysam
 signal.alarm(300)
 if any((DIRECTORY/n).exists() for n in ('controls.json','stdout.txt','stderr.txt')):raise ValueError('qualification output already exists')
 os.environ['MVA_PROJECT_ROOT']=str(ROOT)
 os.environ['MVA_SAMTOOLS']=str(ROOT/'tools/bio/bin/samtools')
 sources={n:pin(DIRECTORY/n) for n in ['verify-depth.py','test_verify_depth.py','qualify.py']}
 old=ROOT/'methods/coding-depth-20260920/qualification-v2'
 original_sources={n:pin(old/n) for n in ['controls.json','workflow/scripts/measure_coding_depth.py','workflow/scripts/coding_depth_targets.py','workflow/scripts/survey_chromosome_depth.py','workflow/scripts/prepare_gatk_compatible_cram.py']}
 start=time.monotonic()
 with (DIRECTORY/'stdout.txt').open('x') as stdout,(DIRECTORY/'stderr.txt').open('x') as stderr:
  from contextlib import redirect_stdout,redirect_stderr
  with redirect_stdout(stdout),redirect_stderr(stderr):
   suite=unittest.defaultTestLoader.discover(str(DIRECTORY),pattern='test_verify_depth.py')
   result=unittest.TextTestRunner(stream=stderr,verbosity=2).run(suite)
 success=result.wasSuccessful() and not result.skipped
 for n,expected in sources.items():
  if pin(DIRECTORY/n)!=expected:raise ValueError('qualification source changed')
 report=dict(status='complete' if success else 'failed',completed_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),host=platform.node(),python=sys.version,pysam=pysam.__version__,runtime=pin(Path(sys.executable).resolve()),samtools=pin(ROOT/'tools/bio/bin/samtools'),scope='Fictional controls only; original source code read, no case alignment or case depth accessed.',tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),skipped=len(result.skipped),source_sha256={n:v['sha256'] for n,v in sources.items()},source_pins=sources,original_qualified_source_pins=original_sources,elapsed_seconds=time.monotonic()-start,logs={n:pin(DIRECTORY/n) for n in ['stdout.txt','stderr.txt']})
 (DIRECTORY/'controls.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
 print(json.dumps({k:report[k] for k in ['status','tests_run','errors','failures','skipped','elapsed_seconds']}))
 sys.exit(0 if success else 1)
