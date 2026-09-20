#!/usr/bin/env python3
"""Reconcile public SV benchmark counts and inspect fixed illustrative read windows.

Local CIGAR support is approximate and does not validate a complete allele.
Selection is by coordinate within declared categories, not by favorable outcome.
"""
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import time

import pysam
import truvari

from build_structural_exclusions import fingerprint, save, stamp
from prepare_public_sv_benchmark import ALIGNMENT, REFERENCE


def size_bin(size):
    for end,label in ((100,'50-99'),(500,'100-499'),(1000,'500-999'),(10000,'1000-9999')):
        if size<end:return label
    return '10000-50000'


def records(path):
    out=[]
    with truvari.VariantFile(str(path)) as reader:
        sample=list(reader.header.samples)[0]
        for e in reader:
            kind=e.var_type().name;d=dict(e.samples[sample]);info=dict(e.info)
            out.append({'chrom':e.chrom,'pos':e.pos,'start':e.start,'end':e.stop,'id':e.id,
                        'type':kind,'comparison_class':'INS_like' if kind in ('DUP','INS') else kind,
                        'size':e.var_size(),'decision':d.get('BD'),'GT':d.get('GT'),
                        'repeat_annotated':'TRF' in info,'low_complexity':info.get('LCR'),
                        'refined':info.get('IsRefined',False),'raw_info':info})
    return out


def read_window(reader, row):
    start=row['start'];end=row['end'];size=row['size']
    # Keep bounded local read diagnostics even for large symbolic deletions/DUPs.
    windows=[(max(0,start-250),start+250)]
    if end>start+500:windows.append((end-250,end+250))
    seen=set();passing=set();cigar_support=set();clips=set();directions=Counter();observed=Counter();n=Counter()
    for a,b in windows:
        for read in reader.fetch(row['chrom'],a,b):
            key=(read.query_name,read.flag,read.reference_start,read.cigarstring)
            if key in seen:continue
            seen.add(key);n['alignment_records']+=1
            if read.flag & 0xF04 or not 20<=read.mapping_quality<255:continue
            passing.add(read.query_name);n['qualifying_alignment_records']+=1
            if any(op==4 and length>=10 for op,length in read.cigartuples or ()):clips.add(read.query_name)
            rp=read.reference_start;qp=0
            for op,length in read.cigartuples or ():
                if op in (1,2) and abs(rp-start)<=100:
                    observed[('INS' if op==1 else 'DEL')+':'+str(length)]+=1
                    expected_op=2 if row['type']=='DEL' else 1
                    if op==expected_op and abs(length-size)<=max(2,round(size*.1)):
                        quals=read.query_qualities
                        nearby=quals[max(0,qp-5):min(len(quals),qp+5)] if quals is not None else ()
                        if nearby and min(nearby)>=20:
                            cigar_support.add(read.query_name);directions['reverse' if read.is_reverse else 'forward']+=1
                if op in (0,2,3,7,8):rp+=length
                if op in (0,1,4,7,8):qp+=length
    return {**dict(n),'qualifying_fragments':len(passing),'clipped_fragments':len(clips),
            'approximately_size_matched_CIGAR_fragments':len(cigar_support),'supporting_alignment_directions':dict(directions),
            'near_start_CIGAR_indels':dict(observed),
            'policy':'Each mate counted as an alignment; support collapsed by fragment name. Exclude 0xF04 and MQ<20/MQ255, local flanking BQ>=20 for approximate CIGAR support. No read names saved.',
            'limitation':'Within 100 bp and max(2 bp, 10%) length tolerance; repetitive allele identity, haplotype, copy number and clinical significance are not established.'}


def audit(root, comparison, pilot, output):
    output.mkdir(parents=True,exist_ok=False);start=time.monotonic();bench=comparison/'benchmark'
    base=records(bench/'refine.base.vcf.gz');comp=records(bench/'refine.comp.vcf.gz')
    stat=json.loads((bench/'refine.variant_summary.json').read_text())
    if Counter(r['decision'] for r in base)!=Counter({'TP':stat['TP-base'],'FN':stat['FN']}):raise ValueError('base count mismatch')
    if Counter(r['decision'] for r in comp)!=Counter({'TP':stat['TP-comp'],'FP':stat['FP']}):raise ValueError('comparison count mismatch')
    categories={};bins={}
    for group,rows in [('baseline',base),('comparison',comp)]:
        c=Counter((r['comparison_class'],r['decision']) for r in rows)
        categories[group]={':'.join(k):v for k,v in sorted(c.items())}
        c=Counter((r['comparison_class'],size_bin(r['size']),r['decision']) for r in rows)
        bins[group]={':'.join(k):v for k,v in sorted(c.items())}
    query=[(int(f[1]),int(f[2])) for l in (pilot/'query-intervals.bed').read_text().splitlines() if (f:=l.split('\t'))]
    outside=[r for r in base if not any(a<=r['start'] and r['end']<=b for a,b in query)]
    initial_fp={r['id'] for r in records(bench/'fp.vcf.gz')}
    initial_tp={r['id'] for r in records(bench/'tp-comp.vcf.gz')}
    final_tp={r['id'] for r in comp if r['decision']=='TP'}
    rescued=final_tp-initial_tp
    if not initial_tp.issubset(final_tp) or not rescued.issubset(initial_fp):
        raise ValueError('unexpected refinement reclassification')
    select=[]
    rules=[('first missed insertion',base,lambda r:r['decision']=='FN' and r['type']=='INS'),
           ('first missed deletion',base,lambda r:r['decision']=='FN' and r['type']=='DEL'),
           ('first missed deletion without TRF annotation',base,lambda r:r['decision']=='FN' and r['type']=='DEL' and not r['repeat_annotated']),
           ('first matched deletion',comp,lambda r:r['decision']=='TP' and r['type']=='DEL'),
           ('first matched insertion',comp,lambda r:r['decision']=='TP' and r['type']=='INS'),
           ('first unmatched duplication',comp,lambda r:r['decision']=='FP' and r['type']=='DUP'),
           ('first unmatched deletion',comp,lambda r:r['decision']=='FP' and r['type']=='DEL'),
           ('first comparison changed from FP to TP by refinement',comp,lambda r:r['id'] in rescued)]
    with pysam.AlignmentFile(str(root/'ref/giab-sv-20260919/inputs'/ALIGNMENT)) as reads:
        for name,rows,predicate in rules:
            eligible=sorted((r for r in rows if predicate(r)),key=lambda r:(r['pos'],r['size']))
            if not eligible:select.append({'selection':name,'status':'no eligible variant'});continue
            row=eligible[0]
            nearest=sorted(base,key=lambda x:abs(x['pos']-row['pos']))[:3]
            select.append({'selection':name,'variant':row,'local_reads':read_window(reads,row),
                           'nearest_benchmark_variants':[{'pos':x['pos'],'type':x['type'],'size':x['size'],'decision':x['decision']} for x in nearest]})
    receipt={'completed_at':stamp(),'elapsed_seconds':time.monotonic()-start,'scope':'public HG002 chr20 illustrative diagnostics',
             'counts_reconcile_with_native_refined_summary':True,'counts':categories,'by_size':bins,
             'baseline_events_not_fully_in_native_query_windows':len(outside),
             'outside_query_counts':dict(Counter(r['decision']+':'+r['type'] for r in outside)),
             'comparison_IDs_changed_from_FP_to_TP':sorted(rescued),
             'curated_examples':select,'selection_policy':'First by coordinate in each named outcome/type/repeat category; examples may overlap and are not additional independent observations. Refinement rescue requires an actual initial-FP to final-TP change; IsRefined alone is not sufficient.',
             'limitations':['Native FP/FN labels are comparison results, not adjudicated truth about every discordance.',
                            'The local read diagnostics neither validate the complete allele nor distinguish every mapping, repeat, phasing or benchmark error.',
                            'Whole-genome validation, broader variant methods and independent discordance review remain open.'],
             'code':fingerprint(Path(__file__)),'comparison_summary':fingerprint(comparison/'summary.json')}
    save(output/'summary.json',receipt);return receipt


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('root','comparison','pilot','output'):p.add_argument('--'+name,required=True,type=Path)
    a=p.parse_args();r=audit(a.root,a.comparison,a.pilot,a.output)
    print(json.dumps({'counts':r['counts'],'outside_query':r['outside_query_counts'],'elapsed_seconds':r['elapsed_seconds'],
                      'examples':[{'selection':x['selection'],'pos':x.get('variant',{}).get('pos'),'type':x.get('variant',{}).get('type'),
                                   'size':x.get('variant',{}).get('size'),'fragments':x.get('local_reads',{}).get('qualifying_fragments'),
                                   'approximate_support':x.get('local_reads',{}).get('approximately_size_matched_CIGAR_fragments')} for x in r['curated_examples']]}))
