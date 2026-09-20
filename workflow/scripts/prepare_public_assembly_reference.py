#!/usr/bin/env python3
"""Prepare an independent, sequence-identical BWA index for public assembly tests."""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

import pysam

from build_structural_exclusions import fingerprint, save, stamp

SOURCE_SHA = '3c8def6d325c5d1e934b2dd530c4d1709f027677a15859467071fddf3fff2026'
BWA_SHA = '34f338e639956e9d14d582283c44a8b52ad0fff8944ccd1dabaeb8d651bcf04e'


def prepare(source, bwa, output):
    if output.resolve().is_relative_to(Path(__file__).resolve().parents[2]):raise ValueError('reference assets stay outside Git')
    output.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    state={'status':'running','pid':os.getpid(),'started_at':stamp(),'scope':'public HG002 assembly reference only'}
    save(output/'status.json',state)
    try:
        if shutil.disk_usage(output).free<350*2**30:raise ValueError('preserve disk headroom for active patient workflow')
        inputs={'compressed_reference':fingerprint(source),'bwa':fingerprint(bwa),'code':fingerprint(Path(__file__))}
        if inputs['compressed_reference']['sha256']!=SOURCE_SHA or inputs['bwa']['sha256']!=BWA_SHA:raise ValueError('public reference or BWA checksum mismatch')
        state['stage']='decompress unchanged reference';save(output/'status.json',state)
        reference=output/'GRCh38_no_alt.fa';digest=hashlib.sha256()
        with gzip.open(source,'rb') as reader, reference.open('wb') as writer:
            for block in iter(lambda:reader.read(8*1024*1024),b''):writer.write(block);digest.update(block)
        pysam.faidx(str(reference))
        with pysam.FastaFile(str(source)) as old,pysam.FastaFile(str(reference)) as new:
            if list(zip(old.references,old.lengths))!=list(zip(new.references,new.lengths)):raise ValueError('reference dictionary changed')
            for name in old.references:
                if old.fetch(name)!=new.fetch(name):raise ValueError('reference sequence changed')
        state['stage']='build public BWA index';save(output/'status.json',state)
        tick=time.monotonic();args=[str(bwa),'index',str(reference)]
        with (output/'bwa-index.log').open('w') as log:p=subprocess.run(args,stdout=log,stderr=subprocess.STDOUT)
        command={'command':args,'exit_code':p.returncode,'elapsed_seconds':time.monotonic()-tick}
        save(output/'command.json',command)
        if p.returncode:raise RuntimeError('public BWA indexing failed')
        files={suffix:fingerprint(Path(str(reference)+suffix)) for suffix in ('','.fai','.amb','.ann','.bwt','.pac','.sa')}
        if files['']['sha256']!=digest.hexdigest():raise ValueError('reference changed during indexing')
        result={'status':'complete','completed_at':stamp(),'elapsed_seconds':time.monotonic()-start,'inputs':inputs,'outputs':files,'command':command,'all_reference_sequences_identical':True,'scope':state['scope']}
        save(output/'summary.json',result);state.update(status='complete',completed_at=result['completed_at']);save(output/'status.json',state);return result
    except BaseException as error:
        state.update(status='failed',failed_at=stamp(),error_type=type(error).__name__,error=str(error));save(output/'status.json',state);raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('source','bwa','output'):p.add_argument('--'+name,required=True,type=Path)
    a=p.parse_args();r=prepare(a.source,a.bwa,a.output);print(json.dumps({'status':r['status'],'elapsed_seconds':r['elapsed_seconds']}))
