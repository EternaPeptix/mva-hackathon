#!/usr/bin/env python3
"""Recheck completed native output after a wrapper failure, preserving the failure."""
import argparse
import json
from pathlib import Path
import time

from build_structural_exclusions import fingerprint, save, stamp
from run_svaba_pilot import EXE_SHA, OPTIONS, inspect_outputs


def review(original, output):
    if output.resolve().is_relative_to(Path(__file__).resolve().parents[2]):raise ValueError('native review stays outside Git')
    output.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    command=json.loads((original/'command.json').read_text());state=json.loads((original/'status.json').read_text())
    inputs=json.loads((original/'inputs.json').read_text())['files']
    if command['exit_code']!=0 or state['status']!='failed':raise ValueError('expected a completed native call and separate failed wrapper')
    if state.get('error')!='start out of range (-1)':raise ValueError('this recovery is restricted to the reviewed boundary-coordinate failure')
    if inputs['native']['sha256']!=EXE_SHA:raise ValueError('unreviewed native build')
    for name in ('alignment','reference','alignment_index','native'):
        if fingerprint(Path(inputs[name]['path']))['sha256']!=inputs[name]['sha256']:raise ValueError('native input changed: '+name)
    executed=original/'executed-run_svaba_pilot.py'
    if fingerprint(executed)['sha256']!=inputs['code']['sha256']:raise ValueError('original executed wrapper source not preserved')
    for p in sorted(original.glob('svaba.*')):
        if p.is_file():(output/p.name).symlink_to(p.resolve())
    checked=inspect_outputs(output,Path(inputs['reference']['path']),Path(inputs['alignment']['path']))
    artifacts={p.name:fingerprint(p) for p in output.iterdir() if p.is_file()}
    result={'status':'complete','mode':state['mode'],'completed_at':stamp(),'review_elapsed_seconds':time.monotonic()-start,
            'native_command':command,'inputs':inputs,'native_outputs':checked,'artifacts':artifacts,
            'original_failed_status':fingerprint(original/'status.json'),'original_executed_wrapper':fingerprint(executed),
            'review_code':fingerprint(Path(__file__)),'current_inspection_code':fingerprint(Path(__file__).with_name('run_svaba_pilot.py')),
            'reference_boundary_code':fingerprint(Path(__file__).with_name('convert_svaba_pilot.py')),
            'scope':'Readback of the exact completed native output; no native caller rerun, coordinate repair or output filtering',
            'boundary_policy':'VCF permits virtual breakends at 0 or length+1. These N-reference boundaries have no base to verify and remain explicitly unverified; this does not establish a biological telomere event.'}
    save(output/'summary.json',result);save(output/'status.json',{'status':'complete','completed_at':result['completed_at'],'scope':result['scope']});return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('original','output'):p.add_argument('--'+name,required=True,type=Path)
    a=p.parse_args();r=review(a.original,a.output);print(json.dumps(r['native_outputs']))
