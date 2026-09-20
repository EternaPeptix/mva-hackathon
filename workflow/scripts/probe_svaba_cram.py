#!/usr/bin/env python3
"""Check the native caller on the previously verified synthetic CRAM encoding."""
import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

import pysam

from build_structural_exclusions import fingerprint, save, stamp
from run_svaba_pilot import EXE_SHA, OPTIONS, SYNTH_REF_SHA, inspect_outputs


def canonical(path, sample):
    rows=[]
    with pysam.VariantFile(str(path)) as reader:
        if list(reader.header.samples)!=[sample]:raise ValueError('caller sample path mismatch')
        for r in reader:
            rows.append({'chrom':r.chrom,'pos':r.pos,'ref':r.ref,'alt':r.alts,'qual':r.qual,'filter':list(r.filter),
                         'info':{k:v for k,v in r.info.items() if k not in ('MATEID','SCTG')},
                         'format':dict(r.samples[sample]),'phased':r.samples[sample].phased})
    return Counter(json.dumps(r,sort_keys=True,allow_nan=False) for r in rows)


def run(root,output):
    output.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    state={'status':'running','pid':os.getpid(),'started_at':stamp(),'scope':'synthetic same-read BAM/CRAM comparison'};save(output/'status.json',state)
    try:
        native=root/'tools/svaba-1.2.0/bin/svaba';fixture=root/'results/structural-synthetic-20260919-v2'
        reference=fixture/'fixture/reference.fa';bam=fixture/'aligned.bam';cram=fixture/'cram-check/aligned.cram'
        checks={native:EXE_SHA,reference:SYNTH_REF_SHA,cram:'eba11eed579cde5164905ba02628a8449cc438e1fe340dbe52e0c1ca0ac0d147',
                Path(str(cram)+'.crai'):'3bae45edf4cee7a900d38182bc4183e6132a4e139b17a367c1cad1d6168955c1'}
        inputs={}
        for p,expected in checks.items():
            identity=fingerprint(p)
            if identity['sha256']!=expected:raise ValueError('synthetic CRAM/runtime input checksum mismatch')
            inputs[str(p)]=identity
        # Each SQ URI must identify this verified local reference. An explicit
        # empty local cache path prevents a default public MD5-reference lookup.
        with pysam.AlignmentFile(str(cram),reference_filename=str(reference)) as reader,pysam.FastaFile(str(reference)) as ref:
            if list(zip(reader.references,reader.lengths))!=list(zip(ref.references,ref.lengths)):raise ValueError('CRAM reference dictionary mismatch')
            for sq in reader.header.to_dict()['SQ']:
                if sq.get('UR')!=str(reference) or sq.get('M5')!=hashlib.md5(ref.fetch(sq['SN']).upper().encode()).hexdigest():
                    raise ValueError('CRAM local reference URI or sequence MD5 mismatch')
        empty_cache=output/'empty-reference-cache';empty_cache.mkdir()
        env=dict(os.environ);env['REF_PATH']=str(empty_cache);env['REF_CACHE']=str(empty_cache/'%2s/%2s/%s');env['OMP_NUM_THREADS']='2'
        args=[str(native),'run','-t',str(cram),'-G',str(reference),'-a','svaba',*OPTIONS]
        tick=time.monotonic()
        with (output/'native.log').open('w') as log:p=subprocess.run(args,cwd=output,env=env,stdout=log,stderr=subprocess.STDOUT)
        command={'command':args,'exit_code':p.returncode,'elapsed_seconds':time.monotonic()-tick,'reference_environment':{k:env[k] for k in ('REF_PATH','REF_CACHE')}}
        save(output/'command.json',command)
        if p.returncode:raise RuntimeError('native synthetic CRAM caller failed; log retained')
        inspection=inspect_outputs(output,reference,cram);comparisons={}
        baseline=root/'results/svaba-synthetic-20260919'
        for kind in ('indel','sv'):
            for filter_prefix in ('','unfiltered.'):
                name='svaba.svaba.'+filter_prefix+kind+'.vcf'
                old=canonical(baseline/name,str(bam));new=canonical(output/name,str(cram))
                comparisons[name]={'BAM_records':sum(old.values()),'CRAM_records':sum(new.values()),'same_canonical_fields':old==new,
                                   'BAM_only':[json.loads(x) for x in (old-new).elements()],
                                   'CRAM_only':[json.loads(x) for x in (new-old).elements()]}
        result={'completed_at':stamp(),'elapsed_seconds':time.monotonic()-start,'scope':state['scope'],'inputs':inputs,
                'command':command,'native_output_checks':inspection,'comparisons':comparisons,
                'all_canonical_call_fields_identical':all(x['same_canonical_fields'] for x in comparisons.values()),
                'comparison_policy':'Ignore record ID, MATEID, SCTG and input sample-column path only; compare all remaining coordinates, alleles, FILTER, QUAL, INFO and genotype/sample fields as multisets',
                'code':fingerprint(Path(__file__)),'caller_wrapper':fingerprint(Path(__file__).with_name('run_svaba_pilot.py')),
                'limitations':['Same fictional reads re-encoded; not new biological evidence','File-format agreement does not repair genotype mistakes or the missing inversion junction']}
        save(output/'summary.json',result);state.update(status='complete',completed_at=result['completed_at']);save(output/'status.json',state);return result
    except BaseException as error:
        state.update(status='failed',failed_at=stamp(),error_type=type(error).__name__,error=str(error));save(output/'status.json',state);raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('root','output'):p.add_argument('--'+name,required=True,type=Path)
    a=p.parse_args();r=run(a.root,a.output)
    print(json.dumps({'elapsed_seconds':r['elapsed_seconds'],'all_canonical_call_fields_identical':r['all_canonical_call_fields_identical'],'comparisons':r['comparisons']}))
