#!/usr/bin/env python3
"""Snapshot only the new depth binding guard and run fictional controls."""
import argparse
import datetime
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import platform
import shutil
import sys
import unittest

FILES=['scripts/launch_secondary_depth_v2.py','tests/test_secondary_depth_guard_v2.py',
       'scripts/qualify_secondary_depth_guard_v2.py','scripts/vendor/guard_helpers.py']

def pins(root):
    return {name:dict(bytes=(root/name).stat().st_size,sha256=hashlib.sha256((root/name).read_bytes()).hexdigest()) for name in FILES}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--version',required=True);a=p.parse_args()
    if Path(a.version).name!=a.version or a.version in ('.','..'):raise ValueError('new single version directory required')
    base=Path(__file__).resolve().parents[2];out=base/a.version;out.mkdir();before=pins(base/'workflow')
    for name in FILES:
        dst=out/'workflow'/name;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(base/'workflow'/name,dst)
    if pins(out/'workflow')!=before:raise ValueError('snapshot differs')
    spec=importlib.util.spec_from_file_location('qualified_secondary_depth_guard',out/'workflow/tests/test_secondary_depth_guard_v2.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    log=io.StringIO();result=unittest.TextTestRunner(stream=log,verbosity=2).run(unittest.defaultTestLoader.loadTestsFromModule(module));(out/'test-output.txt').write_text(log.getvalue())
    if pins(base/'workflow')!=before or pins(out/'workflow')!=before:raise ValueError('source changed')
    passed=result.wasSuccessful() and not result.skipped and result.testsRun==21
    receipt=dict(status='complete' if passed else 'failed',tests_run=result.testsRun,errors=len(result.errors),failures=len(result.failures),skipped=[dict(test=str(t),reason=r) for t,r in result.skipped],exit_code=0 if passed else 1,source_files=before,completed_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),scope='Fictional metadata and wrapper refusal controls only; no patient inputs, alignment extraction, depth computation, or re-execution of unchanged producer controls.',runtime=dict(host=platform.node(),machine=platform.machine(),python=sys.version,executable=sys.executable),log_sha256=hashlib.sha256((out/'test-output.txt').read_bytes()).hexdigest())
    (out/'controls.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps({k:receipt[k] for k in ['status','tests_run','errors','failures','skipped']}|dict(receipt_sha256=hashlib.sha256((out/'controls.json').read_bytes()).hexdigest())))
    return receipt['exit_code']

if __name__=='__main__':raise SystemExit(main())
