#!/usr/bin/env python3
"""Use the unchanged public comparator and report SvABA/DELLY truth overlap."""
import argparse
from collections import Counter
import importlib.metadata
import json
from pathlib import Path
import time

import pysam

from audit_public_sv_comparison import read_window, records, size_bin
from build_structural_exclusions import fingerprint, save, stamp
from prepare_public_sv_benchmark import ALIGNMENT, REFERENCE
from probe_sv_comparator import compare


def key(row): return (row['chrom'],row['pos'],row['id'],row['type'],row['size'])


def run(root, output, config_path, controls=None, converted=None):
    output.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    state={'status':'running','started_at':stamp()};save(output/'status.json',state)
    try:
        config=json.loads(config_path.read_text())
        if importlib.metadata.version('truvari')!=config['truvari_version']:raise ValueError('comparator version changed')
        control=controls or root/'results/svaba-synthetic-comparison-20260919/summary.json'
        checked=json.loads(control.read_text())
        if checked['configuration']['sha256']!=fingerprint(config_path)['sha256']:raise ValueError('synthetic comparison used different settings')
        if any(checked['native_comparator']['refined'][k]!=v for k,v in {'TP-base':3,'TP-comp':3,'FP':0,'FN':0}.items()):raise ValueError('native conversion control failed')
        if checked['conversion_code']['sha256']!=fingerprint(Path(__file__).with_name('convert_svaba_pilot.py'))['sha256']:
            raise ValueError('conversion code changed after native synthetic controls')
        prior=root/'results/structural-public-comparison-20260919-v2';prior_receipt=json.loads((prior/'summary.json').read_text())
        if prior_receipt['config']['sha256']!=fingerprint(config_path)['sha256']:raise ValueError('public comparator settings changed')
        baseline=prior/'truth-chr20.no-star.full-context.vcf.gz'
        bed=root/'results/structural-public-pilot-20260919/official-chr20-benchmark.bed'
        reference=root/'ref/giab-sv-20260919/inputs'/REFERENCE
        for name,path in [('baseline',baseline),('benchmark_bed',bed),('reference',reference)]:
            if fingerprint(path)['sha256']!=prior_receipt['inputs'][name]['sha256']:raise ValueError('public comparator input changed: '+name)
        conversion=converted or root/'results/svaba-public-conversion-20260919';conversion_receipt=json.loads((conversion/'summary.json').read_text())
        comp=conversion/'comparison.full-context.vcf.gz'
        if conversion_receipt['unaccounted_records'] or fingerprint(comp)['sha256']!=conversion_receipt['output_vcf']['sha256']:
            raise ValueError('public conversion incomplete or changed')
        if conversion_receipt['code']['sha256']!=checked['conversion_code']['sha256']:raise ValueError('public conversion code not covered by controls')
        state['stage']='unchanged native benchmark and full-context refinement';save(output/'status.json',state)
        comparison=compare(root,baseline,comp,reference,bed,output/'benchmark',config)
        bench=output/'benchmark';base=records(bench/'refine.base.vcf.gz');called=records(bench/'refine.comp.vcf.gz')
        native=comparison['refined']
        if Counter(r['decision'] for r in base)!=Counter({'TP':native['TP-base'],'FN':native['FN']}):raise ValueError('baseline count reconciliation failed')
        if Counter(r['decision'] for r in called)!=Counter({'TP':native['TP-comp'],'FP':native['FP']}):raise ValueError('call count reconciliation failed')
        by_class={};by_size={}
        for label,rows in [('baseline',base),('comparison',called)]:
            by_class[label]={':'.join(k):v for k,v in sorted(Counter((r['comparison_class'],r['decision']) for r in rows).items())}
            by_size[label]={':'.join(k):v for k,v in sorted(Counter((r['comparison_class'],size_bin(r['size']),r['decision']) for r in rows).items())}
        old=records(prior/'benchmark/refine.base.vcf.gz');old_keys={key(r) for r in old};new_keys={key(r) for r in base}
        if len(old_keys)!=len(old) or len(new_keys)!=len(base) or old_keys!=new_keys:raise ValueError('public truth universes differ')
        previous_tp={key(r) for r in old if r['decision']=='TP'};current_tp={key(r) for r in base if r['decision']=='TP'}
        sets={'matched_by_both':previous_tp & current_tp,'DELLY_only':previous_tp-current_tp,
              'SvABA_only':current_tp-previous_tp,'missed_by_both':new_keys-(previous_tp|current_tp)}
        overlap={name:{'total':len(s),'classes':dict(Counter(k[3] for k in s))} for name,s in sets.items()}
        selections=[]
        rules=[('first additional SvABA deletion',base,lambda r:key(r) in sets['SvABA_only'] and r['type']=='DEL'),
               ('first additional SvABA insertion',base,lambda r:key(r) in sets['SvABA_only'] and r['type']=='INS'),
               ('first deletion missed by both',base,lambda r:key(r) in sets['missed_by_both'] and r['type']=='DEL'),
               ('first insertion missed by both',base,lambda r:key(r) in sets['missed_by_both'] and r['type']=='INS'),
               ('first unmatched SvABA deletion',called,lambda r:r['decision']=='FP' and r['type']=='DEL'),
               ('first unmatched SvABA insertion',called,lambda r:r['decision']=='FP' and r['type']=='INS')]
        with pysam.AlignmentFile(str(root/'ref/giab-sv-20260919/inputs'/ALIGNMENT)) as reads:
            for label,rows,predicate in rules:
                eligible=sorted((r for r in rows if predicate(r)),key=key)
                if not eligible:selections.append({'selection':label,'status':'no eligible variant'});continue
                r=eligible[0];selections.append({'selection':label,'variant':r,'local_reads':read_window(reads,r)})
        result={'completed_at':stamp(),'elapsed_seconds':time.monotonic()-start,'comparison':comparison,
                'by_class':by_class,'by_size':by_size,'truth_overlap_with_DELLY':overlap,
                'truth_union_matched':len(previous_tp|current_tp),'identical_truth_universe':len(new_keys),
                'selected_examples':selections,'configuration':fingerprint(config_path),'conversion':fingerprint(conversion/'summary.json'),
                'synthetic_controls':fingerprint(control),'code':fingerprint(Path(__file__)),
                'limitations':['Single prealigned public chromosome and draft GIAB v5.0q truth; no patient or whole-genome sensitivity conclusion',
                               'Truth overlap is a set of matched baseline records from two separate comparisons, not a validated combined callset',
                               'Reference-copy and assembly-breakend representation assumptions remain visible in conversion records',
                               'Approximate local CIGAR support neither establishes nor refutes complete allele identity',
                               'Genotype agreement is separate from detection and never a causal probability']}
        save(output/'summary.json',result);state.update(status='complete',completed_at=result['completed_at']);save(output/'status.json',state);return result
    except BaseException as error:
        state.update(status='failed',failed_at=stamp(),error_type=type(error).__name__,error=str(error));save(output/'status.json',state);raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('root','output','config'):p.add_argument('--'+name,required=True,type=Path)
    p.add_argument('--controls',type=Path);p.add_argument('--converted',type=Path)
    a=p.parse_args();r=run(a.root,a.output,a.config,a.controls,a.converted)
    print(json.dumps({'comparison':r['comparison']['refined'],'by_class':r['by_class'],'truth_overlap':r['truth_overlap_with_DELLY']}))
