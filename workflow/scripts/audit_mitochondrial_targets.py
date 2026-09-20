#!/usr/bin/env python3
"""Compare selected mitochondrial SNV evidence before and after mitochondrial-only remapping."""
import argparse
from collections import Counter
import json
from pathlib import Path

import pysam

from audit_targeted_snvs import fragment_evidence, validate_targets
from summarize_targeted_snvs import raw_base_profile
from prepare_gatk_compatible_cram import digest, save, stamp


def quality_profile(bam, target):
    counts = Counter()
    with pysam.AlignmentFile(str(bam), 'rb') as reader:
        for read in reader.fetch(target['chrom'], target['pos']-1, target['pos']):
            if read.flag & 0xF04:
                continue
            query = next((q for q, r in read.get_aligned_pairs() if r == target['pos']-1), None)
            if query is None:
                continue
            base = read.query_sequence[query].upper()
            if base != target['alt']:
                continue
            bq = -1 if read.query_qualities is None else read.query_qualities[query]
            cycle = read.query_length-query if read.is_reverse else query+1
            counts[(read.get_tag('RG'), read.mapping_quality, bq, bool(read.is_reverse),
                    cycle, read.query_length, read.get_tag('NM') if read.has_tag('NM') else -1)] += 1
    return [{'read_group':g,'MQ':mq,'raw_BQ':bq,'reverse':rev,'sequencing_cycle_1based':cycle,
             'read_length':length,'NM':nm,'reads':n} for (g,mq,bq,rev,cycle,length,nm),n in sorted(counts.items())]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('root','run','config','output'):
        p.add_argument('--'+name,type=Path,required=True)
    args=p.parse_args()
    if args.output.exists():raise ValueError('existing target audit preserved')
    state=json.loads((args.run/'status.json').read_text())
    if state['status']!='complete':raise ValueError('mitochondrial discovery incomplete')
    config=json.loads(args.config.read_text())
    targets=config['targets']
    sources={}
    for arm,bam,reference in [('original',args.run/'subset.bam',args.root/'ref/realignment/GRCh38_nochr.fa'),
                              ('regular',args.run/'regular.bam',args.root/'ref/mitochondrial-20260920/Homo_sapiens_assembly38.chrM.fasta'),
                              ('shifted',args.run/'shifted.bam',args.root/'ref/mitochondrial-20260920/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta')]:
        if digest(bam)!=state['output_sha256'][bam.name]:raise ValueError('completed alignment changed')
        sources[arm]={'bam':bam,'reference':reference}
    results=[]
    for target in targets:
        row={'target':target,'evidence':{}}
        remapped='regular' if 576<=target['pos']<=16024 else 'shifted'
        for arm in ('original',remapped):
            selected=dict(target)
            selected['chrom']='MT' if arm=='original' else 'chrM'
            if arm=='shifted':selected['pos']=(target['pos']-1-8000)%16569+1
            selected['key']=f"{selected['chrom']}:{selected['pos']}:{selected['ref']}:{selected['alt']}"
            bam,ref=sources[arm]['bam'],sources[arm]['reference']
            validate_targets([selected],ref)
            values={}
            for mq in (20,60):
                filters={'minimum_mapping_quality':mq,'minimum_base_quality':20,'exclude_flags':'0xF04'}
                values[f'MQ{mq}_BAQ20']=fragment_evidence(bam,ref,selected,filters)
            values['raw_profile']=raw_base_profile(bam,selected)
            values['raw_ALT_quality_profile']=quality_profile(bam,selected)
            row['evidence'][arm]=values
        results.append(row)
    source_files=[Path(__file__),args.config,Path(__file__).with_name('audit_targeted_snvs.py'),
                  Path(__file__).with_name('audit_baseline_reads.py'),Path(__file__).with_name('summarize_targeted_snvs.py'),
                  Path(__file__).with_name('prepare_gatk_compatible_cram.py')]
    result={'completed_at':stamp(),'case_status_sha256':digest(args.run/'status.json'),'config':config,
        'source_sha256':{str(x):digest(x) for x in source_files},'targets':results,
        'scope':'Bounded same-read comparisons using existing fragment/BAQ methods and raw ALT quality/cycle bins. Original whole-genome duplicate/MQ decisions and mitochondrial remapping remain distinct. No validated low-level heteroplasmy, contamination or NUMT exclusion.'}
    save(args.output,result)
    print(json.dumps({'status':'complete','targets':len(results),'output_sha256':digest(args.output)}))


if __name__=='__main__':main()
