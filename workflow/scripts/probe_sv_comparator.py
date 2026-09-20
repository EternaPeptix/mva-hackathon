#!/usr/bin/env python3
"""Execute native SV comparator controls, including an equivalent split allele."""
import argparse
import json
import os
from pathlib import Path
import random
import subprocess
import time

import pysam
import truvari

from build_structural_exclusions import fingerprint, save, stamp


def run_command(args, log, env):
    start = time.monotonic()
    with log.open('w') as handle: p = subprocess.run(list(map(str, args)), stdout=handle, stderr=subprocess.STDOUT, env=env)
    receipt = {'command': list(map(str,args)), 'exit_code': p.returncode,
               'elapsed_seconds': time.monotonic()-start, 'log': fingerprint(log)}
    if p.returncode: raise RuntimeError('comparator failed; log retained at ' + str(log))
    return receipt


def environment(root):
    env = dict(os.environ)
    env['PATH'] = ':'.join(str(root/'tools'/name/'bin') for name in ('truvari-5.5.0','mafft-7.526','bio')) + ':' + env.get('PATH','')
    env['OMP_NUM_THREADS'] = '2'; env['OPENBLAS_NUM_THREADS'] = '2'
    return env


def compare(root, baseline, comp, reference, bed, output, config, refine=True):
    exe = root/'tools/truvari-5.5.0/bin/truvari'; env = environment(root)
    args = [exe,'bench','-b',baseline,'-c',comp,'-f',reference,'--includebed',bed,'-o',output,*config['bench_options']]
    commands = [run_command(args, output.with_suffix('.bench.log'), env)]
    initial = json.loads((output/'summary.json').read_text())
    if refine:
        commands.append(run_command([exe,'refine',output,'--use-original-vcfs','--align','mafft','--threads','2',
                        '--mafft-params','--auto --thread 1'],output.with_suffix('.refine.log'),env))
    refined = None; refine_status = 'not_requested'
    if refine:
        summary_path = output/'refine.variant_summary.json'
        if summary_path.exists():
            refined = json.loads(summary_path.read_text()); refine_status = 'completed'
        else:
            log = output.with_suffix('.refine.log').read_text()
            candidates = (output/'candidate.refine.bed').read_text().strip()
            if candidates or 'No regions to be refined. Exiting' not in log:
                raise ValueError('native refinement produced no summary without verified empty candidates')
            refined = {key:initial[key] for key in ('TP-base','TP-comp','FP','FN','precision','recall','f1','base cnt','comp cnt')}
            refine_status = 'no_candidate_regions; original_counts_preserved'
    return {'bench': initial, 'refined': refined, 'refinement_status': refine_status,
            'commands': commands}


def write_vcf(path, seq, rows):
    header = pysam.VariantHeader(); header.add_line('##contig=<ID=chr1,length=12000>')
    header.add_line('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Synthetic control type">')
    header.add_line('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Synthetic control length">')
    header.add_line('##INFO=<ID=END,Number=1,Type=Integer,Description="End position">')
    header.formats.add('GT',1,'String','Genotype'); header.add_sample('CONTROL')
    with pysam.VariantFile(str(path),'wz',header=header) as out:
        for name,start,end,alt in rows:
            ref = seq[start-1:end] if end>start else seq[start-1]
            if alt=='<DEL>': ref=seq[start-1]
            record=out.new_record(contig='chr1',start=start-1,stop=end if end>start else start,id=name,alleles=(ref,alt))
            record.filter.add('PASS'); record.samples['CONTROL']['GT']=(1,0);record.samples['CONTROL'].phased=True
            if alt=='<DEL>':record.info['SVTYPE']='DEL';record.info['SVLEN']=start-end
            out.write(record)
    pysam.tabix_index(str(path),preset='vcf')


def probe(root, output, config_path):
    output.mkdir(parents=True,exist_ok=False);start=time.monotonic();config=json.loads(config_path.read_text())
    rng=random.Random(20260919);seq=''.join(rng.choices('ACGT',k=12000));ref=output/'control.fa'
    ref.write_text('>chr1\n'+seq+'\n');pysam.faidx(str(ref));bed=output/'regions.bed';bed.write_text('chr1\t0\t12000\n')
    deletion=('DEL',2000,2120,seq[1999]);insertion=('INS',5000,5000,seq[4999]+'A'*120)
    baseline=[deletion,insertion];base=output/'base.vcf.gz';write_vcf(base,seq,baseline)
    cases=[('exact',baseline,{'TP-base':2,'TP-comp':2,'FN':0,'FP':0}),
           ('empty',[],{'TP-base':0,'TP-comp':0,'FN':2,'FP':0}),
           ('wrong_locus',[('WRONG',9000,9120,seq[8999])],{'TP-base':0,'TP-comp':0,'FN':2,'FP':1}),
           ('wrong_insert_sequence',[deletion,('WRONG_INS',5000,5000,seq[4999]+'T'*120)],{'TP-base':1,'TP-comp':1,'FN':1,'FP':1}),
           ('symbolic_deletion',[('SYMBOLIC',2000,2120,'<DEL>')],{'TP-base':1,'TP-comp':1,'FN':1,'FP':0})]
    results=[]
    for name,rows,expect in cases:
        comp=output/(name+'.vcf.gz');write_vcf(comp,seq,rows)
        got=compare(root,base,comp,ref,bed,output/name,config,refine=name=='exact')
        passed=all(got['bench'][k]==v for k,v in expect.items())
        if name=='exact':
            passed=passed and got['refinement_status']=='no_candidate_regions; original_counts_preserved' and got['refined']['FP']==got['refined']['FN']==0
        results.append({'name':name,'expected':expect,'passed':passed,**got})
        save(output/'partial-results.json',results)
        if not passed:raise ValueError('native comparator control failed: '+name)
    # A 125-base span with five retained central bases equals two phased 60-base deletions.
    base2=output/'complex-base.vcf.gz';comp2=output/'split-comp.vcf.gz'
    write_vcf(base2,seq,[('COMPLEX',2000,2125,seq[1999]+seq[2060:2065])])
    write_vcf(comp2,seq,[('PART_A',2000,2060,seq[1999]),('PART_B',2065,2125,seq[2064])])
    got=compare(root,base2,comp2,ref,bed,output/'split_allele',config)
    passed=(got['bench']['FN']==1 and got['bench']['FP']==2 and got['refined']['FN']==0 and got['refined']['FP']==0)
    results.append({'name':'split_allele','passed':passed,'expected':'initial FN=1, FP=2; after full-context MAFFT refinement FN=0, FP=0',**got})
    receipt={'completed_at':stamp(),'elapsed_seconds':time.monotonic()-start,'all_passed':all(x['passed'] for x in results),
             'cases':results,'code':fingerprint(Path(__file__)),'configuration':fingerprint(config_path),
             'scope':'synthetic comparator controls only; no public-read accuracy or patient conclusion'}
    save(output/'summary.json',receipt)
    if not receipt['all_passed']:raise ValueError('complex-allele refinement control failed')
    return receipt


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('root','output','config'):p.add_argument('--'+name,required=True,type=Path)
    a=p.parse_args();r=probe(a.root,a.output,a.config);print(json.dumps({'all_passed':r['all_passed'],'cases':len(r['cases']),'elapsed_seconds':r['elapsed_seconds']}))
