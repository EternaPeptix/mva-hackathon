#!/usr/bin/env python3
"""Snapshot new public-target wrapper and execute fictional controls only."""
import argparse
import datetime
import hashlib
import io
import importlib.util
import json
from pathlib import Path
import platform
import shutil
import sys
import unittest


def pins(root):
    return {str(p.relative_to(root)):dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
            for p in sorted(root.rglob('*')) if p.is_file() and '__pycache__' not in p.parts}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--version',required=True);parser.add_argument('--require-no-skips',action='store_true');a=parser.parse_args()
    if Path(a.version).name!=a.version or a.version in ('.','..'):raise ValueError('one new version directory required')
    base=Path(__file__).resolve().parents[2];out=base/a.version;out.mkdir()
    before=pins(base/'workflow');shutil.copytree(base/'workflow',out/'workflow',ignore=shutil.ignore_patterns('__pycache__'))
    if before!=pins(out/'workflow'):raise ValueError('source snapshot differs')
    spec=importlib.util.spec_from_file_location('fictional_public_targets',out/'workflow/tests/test_public_targets.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    log=io.StringIO();result=unittest.TextTestRunner(stream=log,verbosity=2).run(unittest.defaultTestLoader.loadTestsFromModule(module));(out/'test-output.txt').write_text(log.getvalue())
    if before!=pins(base/'workflow') or before!=pins(out/'workflow'):raise ValueError('source changed during qualification')
    passed=result.wasSuccessful() and (not a.require_no_skips or not result.skipped)
    receipt=dict(status='complete' if passed else 'failed',exit_code=0 if passed else 1,tests_run=result.testsRun,errors=len(result.errors),failures=len(result.failures),skipped=[dict(test=str(t),reason=r) for t,r in result.skipped],
        completed_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),source_files=before,scope='Fictional public references and intentional corruptions only; no real 84-gene target generation, alignment, index, case read, variant source or download.',runtime=dict(host=platform.node(),machine=platform.machine(),python=sys.version,executable=sys.executable),log_sha256=hashlib.sha256((out/'test-output.txt').read_bytes()).hexdigest())
    (out/'controls.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps({k:receipt[k] for k in ['status','tests_run','errors','failures','skipped']}|dict(receipt_sha256=hashlib.sha256((out/'controls.json').read_bytes()).hexdigest())))
    return receipt['exit_code']


if __name__=='__main__':raise SystemExit(main())
