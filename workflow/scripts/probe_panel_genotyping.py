#!/usr/bin/env python3
"""Native constrained genotyping controls, with missingness and third-allele evidence."""
import argparse
from collections import Counter
import datetime
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

import pysam
from panel_marker_evidence import assess_marker


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def create_fixture(root):
    reference=root/'synthetic.fa'
    reference.write_text('>1\n'+'A'*10000+'\n');pysam.faidx(str(reference))
    cases=[('homref',1001,{'A':30},40,0),('heterozygous',2001,{'A':15,'G':15},40,0),
           ('homalt',3001,{'G':30},40,0),('no_coverage',4001,{},40,0),
           ('low_depth',5001,{'A':3},40,0),('third_allele_mixture',6001,{'A':20,'C':10},40,0),
           ('third_allele_only',7001,{'C':30},40,0),('low_base_quality',8001,{'A':30},10,0),
           ('duplicate_only',9001,{'A':30},40,0x400)]
    header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'1','LN':10000}],
            'RG':[{'ID':'SYNTH.1','SM':'SYNTH','PL':'ILLUMINA'}]}
    bam=root/'synthetic.bam'
    with pysam.AlignmentFile(str(bam),'wb',header=header) as writer:
        for label,pos,counts,bq,flag in cases:
            for base,n in counts.items():
                for i in range(n):
                    read=pysam.AlignedSegment(writer.header);read.query_name=f'{label}_{base}_{i}'
                    sequence=list('A'*150);sequence[75]=base
                    read.query_sequence=''.join(sequence);read.query_qualities=[bq]*150
                    read.reference_id=0;read.reference_start=pos-1-75;read.mapping_quality=60
                    read.cigarstring='150M';read.flag=flag;read.set_tag('RG','SYNTH.1');writer.write(read)
    pysam.index(str(bam))
    (root/'sites.bed').write_text(''.join(f'1\t{pos-1}\t{pos}\n' for _,pos,_,_,_ in cases))
    with pysam.BGZFile(str(root/'alleles.tsv.gz'),'wb') as stream:
        stream.write(''.join(f'1\t{pos}\tA,G\n' for _,pos,_,_,_ in cases).encode())
    pysam.tabix_index(str(root/'alleles.tsv.gz'),seq_col=0,start_col=1,end_col=1)
    return cases,reference,bam


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bcftools',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False);os.nice(10)
    started=time.monotonic();cases,reference,bam=create_fixture(a.output)
    commands=[
        [str(a.bcftools),'mpileup','-Ob','-o',str(a.output/'pileup.bcf'),'-f',str(reference),
         '-R',str(a.output/'sites.bed'),'-a','AD,DP,SP','-q','20','-Q','20','-d','10000',
         '--ff','0xF04','--skip-indels',str(bam)],
        [str(a.bcftools),'call','-m','--keep-alts','--insert-missed','--constrain','alleles',
         '--targets-file',str(a.output/'alleles.tsv.gz'),'--ploidy','2','--annotate','GQ',
         '-Ob','-o',str(a.output/'calls.bcf'),str(a.output/'pileup.bcf')]]
    runs=[]
    for i,command in enumerate(commands):
        tick=time.monotonic();r=subprocess.run(command,text=True,capture_output=True)
        log=a.output/f'step-{i}.stderr';log.write_text(r.stderr)
        runs.append({'command':command,'exit_code':r.returncode,'elapsed_seconds':time.monotonic()-tick,'log_sha256':sha(log)})
        if r.returncode:raise RuntimeError('native genotyping control failed')
    observed={}
    for filename,key in [('pileup.bcf','pileup'),('calls.bcf','constrained_call')]:
        with pysam.VariantFile(str(a.output/filename)) as reader:
            if list(reader.header.samples)!=['SYNTH']:raise ValueError('synthetic sample identity lost')
            for row in reader:
                if row.pos in observed and key in observed[row.pos]:raise ValueError('duplicate result marker')
                sample=row.samples['SYNTH']
                observed.setdefault(row.pos,{})[key]={'alleles':row.alleles,
                    'sample':{field:sample.get(field) for field in ('GT','AD','DP','GQ','PL')}}
    expected_positions={pos for _,pos,_,_,_ in cases}
    if {pos for pos,r in observed.items() if 'constrained_call' in r}!=expected_positions:
        raise AssertionError('requested marker accounting failed')
    results=[]
    expected_gt={'homref':(0,0),'heterozygous':(0,1),'homalt':(1,1)}
    for label,pos,counts,bq,flag in cases:
        row=observed[pos];call=row['constrained_call'];sample=call['sample']
        if call['alleles']!=('A','G'):raise AssertionError('requested alleles were not retained')
        if label in expected_gt and sample['GT']!=expected_gt[label]:raise AssertionError('positive genotype control failed')
        if label in ('no_coverage','low_base_quality','duplicate_only') and sample['GT'] is not None and any(x is not None for x in sample['GT']):
            raise AssertionError('unobserved marker was converted into a called genotype')
        pileup=row.get('pileup');nonpanel_depth=None
        if pileup is not None and pileup['sample']['AD'] is not None:
            nonpanel_depth=sum(depth for allele,depth in zip(pileup['alleles'],pileup['sample']['AD'])
                               if allele not in ('A','G','<*>','<NON_REF>') and depth is not None)
        if label.startswith('third_allele') and nonpanel_depth!=counts['C']:
            raise AssertionError('non-panel read evidence was not retained in the unconstrained pileup')
        eligibility=assess_marker(('A','G'),pileup,call)
        if eligibility['usable_for_declared_ROH_analysis'] != (label in expected_gt):
            raise AssertionError('marker qualification did not distinguish positive and excluded controls')
        results.append({'label':label,'position':pos,'simulated_base_counts':counts,
                        'nonpanel_depth_in_original_pileup':nonpanel_depth,'eligibility':eligibility,**row})
    result={'completed_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'elapsed_seconds':time.monotonic()-started,'commands':runs,'controls':results,
        'positive_genotypes_and_missingness_checks_passed':True,
        'bcftools_version':subprocess.check_output([str(a.bcftools),'--version'],text=True),
        'script_sha256':sha(Path(__file__)),'bcftools_sha256':sha(a.bcftools),
        'evidence_filter_sha256':sha(Path(__file__).with_name('panel_marker_evidence.py')),
        'files':{p.name:sha(p) for p in a.output.iterdir() if p.is_file()},
        'scope':'Synthetic mapped-read genotyping controls only; not end-to-end alignment truth, patient genotypes or clinical sensitivity.'}
    (a.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
