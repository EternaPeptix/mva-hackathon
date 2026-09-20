#!/usr/bin/env python3
"""Repeat the retained synthetic structural experiment through explicit CRAM 3.0."""
import argparse
from collections import Counter
import itertools
import json
import os
from pathlib import Path
import subprocess
import time

import pysam
from build_structural_exclusions import fingerprint, save, stamp
from probe_structural_variants import read_calls


def require(value, message):
    if not value:
        raise ValueError(message)


def identity(read):
    return (read.query_name,read.flag,read.reference_id,read.reference_start,read.mapping_quality,
            read.cigarstring,read.next_reference_id,read.next_reference_start,read.template_length,
            read.query_sequence,read.query_qualities,sorted(read.get_tags(with_value_type=True)))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('experiment','output','samtools','delly'):
        p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--expected-summary-sha256',required=True)
    a=p.parse_args();os.nice(10);start=time.monotonic()
    require(fingerprint(a.experiment/'summary.json')['sha256']==a.expected_summary_sha256,'synthetic summary differs')
    old=json.loads((a.experiment/'summary.json').read_text())
    truth=json.loads((a.experiment/'fixture/truth.json').read_text())
    require(truth['scope']=='Synthetic diploid simple-SV reference; no patient sequence or observed biological effect.','synthetic fixture required')
    bam=a.experiment/'aligned.bam';ref=a.experiment/'fixture/reference.fa'
    for path,expected in [(bam,old['output_sha256']['aligned.bam']),(ref,truth['files']['reference.fa']['sha256']),
                          (a.delly,old['tools']['delly']['sha256']),(a.samtools,old['tools']['samtools']['sha256'])]:
        require(fingerprint(path)['sha256']==expected,'synthetic input/tool differs')
    a.output.mkdir(parents=True,exist_ok=False);commands=[]
    def run(cmd):
        cmd=list(map(str,cmd));log=a.output/('step-%02d.log'%len(commands));t=time.monotonic()
        with log.open('x') as stream:r=subprocess.run(cmd,stdout=stream,stderr=subprocess.STDOUT)
        commands.append({'command':cmd,'exit_code':r.returncode,'elapsed_seconds':time.monotonic()-t,'log':fingerprint(log)})
        save(a.output/'commands.json',commands);require(r.returncode==0,'native synthetic command failed')
    cram=a.output/'aligned.cram'
    run([a.samtools,'view','-C','--output-fmt-option','version=3.0','--output-fmt-option','store_md=1',
         '--output-fmt-option','store_nm=1','-T',ref,'-o',cram,bam])
    run([a.samtools,'index',cram]);run([a.samtools,'quickcheck','-v',cram])
    with cram.open('rb') as f:require(f.read(6)==b'CRAM\x03\x00','output is not CRAM 3.0')
    counts=Counter()
    with pysam.AlignmentFile(str(bam)) as before,pysam.AlignmentFile(str(cram),reference_filename=str(ref)) as after:
        require(before.references==after.references and before.lengths==after.lengths,'dictionary differs')
        require(before.header.to_dict()['RG']==after.header.to_dict()['RG'],'read groups differ')
        for x,y in itertools.zip_longest(before,after):
            require(x is not None and y is not None and identity(x)==identity(y),'decoded read evidence differs')
            counts['records']+=1;counts['duplicate_records']+=x.is_duplicate;counts['supplementary_records']+=x.is_supplementary
    native=[row['command'] for row in old['commands'] if row['stage']=='native structural discovery']
    require(len(native)==1,'prior native command not unique')
    cmd=list(native[0]);require(cmd[0]==str(a.delly) and cmd[-1]==str(bam),'unexpected prior command')
    cmd[-1]=str(cram);cmd[cmd.index('-o')+1]=str(a.output/'calls.bcf');run(cmd)
    calls=json.loads(json.dumps(read_calls(a.output/'calls.bcf','SYNTH_SV')))
    require(calls==old['native_calls'],'structural calls differ from retained synthetic BAM')
    result={'completed_at':stamp(),'elapsed_seconds':time.monotonic()-start,'script':fingerprint(Path(__file__)),
        'prior_summary':fingerprint(a.experiment/'summary.json'),'commands':commands,'read_identity_counts':dict(counts),
        'all_decoded_fields_and_typed_tags_equal':True,'cram_version':[3,0],'native_call_count':len(calls),
        'all_compared_call_fields_equal':True,'native_calls':calls,'delly':fingerprint(a.delly),
        'outputs':{n:fingerprint(a.output/n) for n in ('aligned.cram','aligned.cram.crai','calls.bcf','calls.bcf.csi')},
        'scope':'Bounded synthetic format/caller check, not whole-genome accuracy or completion of patient discovery.'}
    save(a.output/'summary.json',result);print(json.dumps({k:result[k] for k in ('completed_at','elapsed_seconds','read_identity_counts','native_call_count','all_compared_call_fields_equal')}))


if __name__=='__main__':main()
