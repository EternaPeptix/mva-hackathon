#!/usr/bin/env python3
"""Snapshot and run only fictional secondary-adapter controls; never a case configuration."""
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


def pins(root):
    return {str(p.relative_to(root)): dict(bytes=p.stat().st_size, sha256=hashlib.sha256(p.read_bytes()).hexdigest())
            for p in sorted(root.rglob('*')) if p.is_file() and '__pycache__' not in p.parts}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version',required=True)
    parser.add_argument('--require-no-skips',action='store_true')
    args=parser.parse_args()
    if Path(args.version).name != args.version or args.version in ('.','..'):
        raise ValueError('one new version directory required')
    base=Path(__file__).resolve().parents[2];output=base/args.version;output.mkdir()
    before=pins(base/'workflow')
    shutil.copytree(base/'workflow',output/'workflow',ignore=shutil.ignore_patterns('__pycache__'))
    if before != pins(output/'workflow'):raise ValueError('source snapshot differs')
    path=output/'workflow/tests/test_native_secondary.py'
    spec=importlib.util.spec_from_file_location('fictional_native_secondary_controls',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    log=io.StringIO();result=unittest.TextTestRunner(stream=log,verbosity=2).run(unittest.defaultTestLoader.loadTestsFromModule(module))
    (output/'test-output.txt').write_text(log.getvalue())
    if before != pins(base/'workflow') or before != pins(output/'workflow'):raise ValueError('source changed during controls')
    passed=result.wasSuccessful() and (not args.require_no_skips or not result.skipped)
    receipt=dict(status='complete' if passed else 'failed',exit_code=0 if passed else 1,
        completed_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),tests_run=result.testsRun,
        failures=len(result.failures),errors=len(result.errors),skipped=[dict(test=str(test),reason=reason) for test,reason in result.skipped],
        runtime=dict(host=platform.node(),system=platform.system(),machine=platform.machine(),python=sys.version,executable=sys.executable),
        scope='Fictional fixtures only. No real-input configuration, case screening, caller, new database or publication executed.',
        source_sha256={name:value['sha256'] for name,value in before.items()},source_snapshot_files=before,
        log_sha256=hashlib.sha256((output/'test-output.txt').read_bytes()).hexdigest())
    (output/'controls.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:receipt[k] for k in ['status','tests_run','failures','errors']}|dict(skipped=len(result.skipped),receipt=str(output/'controls.json'),receipt_sha256=hashlib.sha256((output/'controls.json').read_bytes()).hexdigest())))
    return receipt['exit_code']


if __name__=='__main__':raise SystemExit(main())
