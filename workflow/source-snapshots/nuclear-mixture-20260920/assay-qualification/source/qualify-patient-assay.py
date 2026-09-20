#!/usr/bin/env python3
"""Freeze candidate sources and run finite fictional controls only."""
import argparse,datetime,hashlib,importlib.util,json,os,platform,shutil,sys,time,unittest
from pathlib import Path
HERE=Path(__file__).resolve().parent
NAMES=['patient_common.py','run-patient-assay.py','audit-patient-output.py','guard-patient-assay.py','guard-fictional-pilot.py','fictional_controls.py','test_patient_assay.py','qualify-patient-assay.py']
EXECUTED=set(NAMES[:5])
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();out=a.output.resolve()
    if platform.node().split('.')[0]!='512S2':raise ValueError('native Mac Studio512S2 qualification required')
    out.mkdir();source=out/'source';source.mkdir()
    for name in NAMES:shutil.copyfile(HERE/name,source/name)
    frozen={name:sha(source/name) for name in NAMES};sys.path.insert(0,str(source));os.environ['MVA_ASSAY_TEST_OUTPUT']=str(out/'fictional-inputs-outputs')
    spec=importlib.util.spec_from_file_location('test_patient_assay',source/'test_patient_assay.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    suite=unittest.defaultTestLoader.loadTestsFromModule(module);start=time.monotonic()
    with (out/'tests.stderr.log').open('x') as log:result=unittest.TextTestRunner(stream=log,verbosity=2).run(suite)
    import patient_common as c
    runtime={k:c.fingerprint(str(c.ROOT/v),c.RUNTIME[k]) for k,v in {'java':'tools/bio/lib/jvm/bin/java','gatk':'gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar'}.items()}
    changed={n:sha(source/n) for n in NAMES}!=frozen
    passed=result.wasSuccessful() and not result.skipped and not changed
    record=dict(status='complete' if passed else 'failed',completed_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),skipped=[str(x) for x in result.skipped],exit_code=0 if passed else 1,elapsed_seconds=time.monotonic()-start,host=platform.node(),machine=platform.machine(),python=c.fingerprint(sys.executable),python_version=sys.version,source_sha256=frozen,executed_source_sha256={n:frozen[n] for n in sorted(EXECUTED)},runtime=runtime,runtime_sha256={k:v['sha256'] for k,v in runtime.items()},all_frozen_sources_preserved=not changed,scope='Fictional wrapper, read-table and saved-output controls only. No patient alignment or callset opened. Two complete native fixture pipelines, not donor-mixture accuracy qualification.')
    with (out/'controls.json').open('x') as f:f.write(json.dumps(record,indent=2,sort_keys=True,allow_nan=False)+'\n')
    print(json.dumps(dict(status=record['status'],tests_run=result.testsRun,controls_sha256=sha(out/'controls.json'))));return record['exit_code']
if __name__=='__main__':raise SystemExit(main())
