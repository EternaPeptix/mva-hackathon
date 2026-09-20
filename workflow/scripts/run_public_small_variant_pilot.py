#!/usr/bin/env python3
"""Run the unchanged genome caller on a pinned public chr20-only alignment.

The complete reference dictionary and ordinary contig sharding are retained.
Only chr20 contains mapped reads, so this is not whole-genome accuracy evidence.
This entry point never starts or changes the patient alignment/calling controllers.
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

import pysam

from build_structural_exclusions import fingerprint, save, stamp
from prepare_public_sv_benchmark import ALIGNMENT


PINS = {
    'alignment': '34ac157739e1feeb590f6eb7e11046ccc2aa3277fd55a3ce0942e774d931ed81',
    'alignment_index': '245dcdcc474b1971eeb36fc9a8af37b76ac8b4c7ce06a7bdee64e0c42bd93f34',
    'reference': '9cce8b926416dd96b152deea85188495b75f7ac8d634cc723a017067be8702b7',
    'reference_index': '502a1b8fb73ccd53285c28a0f12df90c818b4fe3de1e862ef47c593ef1a0a4b4',
    'genome_caller': '5af91153e9926b0d8b985d0575b0a6ef51e1496d16ff3614c4cb4d8b2b4ac358',
    'gatk': 'b2824999fafe4dc83036e334b9765d4a7de4f4ff4bec1e929d77bea04a1560c7',
    'python': '101eed1d9e3db137495f9a54720826eb3d7ac25b8c5a0e82d95372d8e7377e8f',
    'java': 'ec64ac719290ad080c7f8a4f5d83402e4152f98f486fe14150cf734a46ebcbb1',
    'samtools': '7518cbb230a35393d7d3a9f02fc8b4519e4a37b5a17fbed06b9752ca1d962a91',
    'bcftools': 'c427681749daabb97be4a010c31c6b24349aaea6abc7bea27d6ce9126874b260',
}


def run(root, output):
    if output.resolve().is_relative_to(Path(__file__).resolve().parents[2]):
        raise ValueError('genome-scale calls must remain outside the repository')
    output.mkdir(parents=True,exist_ok=False);tick=time.monotonic()
    state={'status':'running','stage':'verify pinned public inputs','pid':os.getpid(),'started_at':stamp()}
    save(output/'status.json',state)
    try:
        alignment=root/'ref/giab-sv-20260919/inputs'/ALIGNMENT
        reference=root/'ref/giab-sv-20260919/assembly-reference/GRCh38_no_alt.fa'
        runtime=root/'tools/bio'
        paths={'alignment':alignment,'alignment_index':Path(str(alignment)+'.bai'),
               'reference':reference,'reference_index':Path(str(reference)+'.fai'),
               'genome_caller':root/'repo/workflow/scripts/call_genome.py',
               'gatk':root/'gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar',
               'python':runtime/'bin/python','java':runtime/'lib/jvm/bin/java',
               'samtools':runtime/'bin/samtools','bcftools':runtime/'bin/bcftools'}
        inputs={}
        for name,path in paths.items():
            inputs[name]=fingerprint(path)
            if inputs[name]['sha256']!=PINS[name]:raise ValueError('public pilot input/build changed: '+name)
        with pysam.FastaFile(str(reference)) as ref, pysam.AlignmentFile(str(alignment)) as bam:
            if list(zip(ref.references,ref.lengths))!=list(zip(bam.references,bam.lengths)):
                raise ValueError('reference/alignment dictionaries differ')
            if bam.header.to_dict().get('HD',{}).get('SO')!='coordinate' or not bam.check_index():
                raise ValueError('sorted indexed public BAM required')
            groups=bam.header.to_dict().get('RG',[])
            if len(groups)!=1 or groups[0].get('SM')!='HG002' or groups[0].get('ID')!='HG002.pcr-free':
                raise ValueError('unexpected public sample/read group')
            mapped={s.contig:s.mapped for s in bam.get_index_statistics() if s.mapped}
            if mapped!={'chr20':19430746}:raise ValueError('unexpected public mapped-read scope')
            dictionary={'contigs':len(ref.references),'bases':sum(ref.lengths),'chr20_bases':ref.get_reference_length('chr20')}
        save(output/'inputs.json',{'files':inputs,'dictionary':dictionary,'mapped_records_by_contig':mapped,
             'wrapper':fingerprint(Path(__file__)),'pysam_version':pysam.__version__,
             'scope':'Public prealigned HG002 chr20, same raw diploid discovery code intended for the child; not whole-genome accuracy'})
        shutil.copyfile(paths['genome_caller'],output/'executed-call_genome.py')
        shutil.copyfile(Path(__file__),output/'executed-run_public_small_variant_pilot.py')
        args=[paths['python'],paths['genome_caller'],'--root',root,'--reference',reference,
              '--alignment',alignment,'--sample','HG002','--output',output/'calls','--workers','4']
        args=list(map(str,args));save(output/'command.json',{'command':args})
        state['stage']='unchanged GATK and bcftools discovery';start=time.monotonic()
        with (output/'native.log').open('w') as log:
            child=subprocess.Popen(args,stdout=log,stderr=subprocess.STDOUT)
            state['child_pid']=child.pid;save(output/'status.json',state);exit_code=child.wait()
        command={'command':args,'exit_code':exit_code,'elapsed_seconds':time.monotonic()-start,
                 'log':fingerprint(output/'native.log')};save(output/'command.json',command)
        if exit_code:raise RuntimeError('native public small-variant workflow failed; outputs retained')
        complete=json.loads((output/'calls/status.json').read_text())
        if complete['status']!='complete':raise ValueError('native completion receipt missing')
        outputs={}
        for name,entry in complete['outputs'].items():
            p=Path(entry['path']);f=fingerprint(p)
            if f['sha256']!=entry['sha256']:raise ValueError('completed public callset changed')
            outputs[name]={**f,'counts':entry['counts']}
        result={'status':'complete','completed_at':stamp(),'elapsed_seconds':time.monotonic()-tick,
                'inputs':fingerprint(output/'inputs.json'),'command':command,'outputs':outputs,
                'native_definition':fingerprint(output/'calls/run-definition.json'),
                'native_completion':fingerprint(output/'calls/status.json'),
                'code':fingerprint(Path(__file__)),'accuracy_status':'Not established; independent truth comparison required',
                'limitations':['One prealigned public chromosome, different upstream mapping/duplicate-marking versions',
                               'Raw discovery without BQSR or post-call quality filtering',
                               'Empty non-chr20 shards preserve the ordinary pipeline and do not validate those chromosomes',
                               'No patient genotype, clinical sensitivity or causal probability follows']}
        save(output/'summary.json',result);state.update(status='complete',completed_at=result['completed_at']);save(output/'status.json',state)
        return result
    except BaseException as error:
        state.update(status='failed',failed_at=stamp(),error_type=type(error).__name__,error=str(error));save(output/'status.json',state);raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('root','output'):p.add_argument('--'+name,required=True,type=Path)
    a=p.parse_args();print(json.dumps(run(a.root,a.output)))
