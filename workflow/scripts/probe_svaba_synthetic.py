#!/usr/bin/env python3
"""Check SvABA against unchanged engineered truth and native comparator controls."""
import argparse
import json
from pathlib import Path

import pysam

from build_structural_exclusions import fingerprint, save, stamp
from probe_structural_variants import compare_truth
from probe_sv_comparator import compare


def run(root, output, config_path, conversion=None):
    output.mkdir(parents=True,exist_ok=False)
    fixture=root/'results/structural-synthetic-20260919-v2/fixture'
    converted=conversion or root/'results/svaba-synthetic-conversion-20260919'
    truth=json.loads((fixture/'truth.json').read_text())
    if fingerprint(fixture/'truth.json')['sha256']!='b770ea22c1b44a4df9d1fd1d402f46d1774d078231e5484a6859fa526a6906ea':
        raise ValueError('synthetic truth changed')
    events=json.loads((converted/'all-native-events.json').read_text())
    represented=json.loads((converted/'conversion-records.json').read_text())
    pass_events=[e for e in events if e['filter']==['PASS']]
    truth_check=compare_truth(truth,pass_events)
    direct=[]; baseline_rows=[]; reference=fixture/'reference.fa'
    with pysam.FastaFile(str(reference)) as ref:
        for event in truth['events_on_one_haplotype']:
            kind=event['type']
            if kind not in ('DEL','DUP','INS'):continue
            chrom=event['chrom'];left=event['left_cut_0based'];right=event['right_cut_0based']
            if kind=='DEL':pos=left; allele_ref=ref.fetch(chrom,left-1,right);alt=allele_ref[0]
            elif kind=='DUP':
                pos=right;allele_ref=ref.fetch(chrom,right-1,right);alt=allele_ref+ref.fetch(chrom,left,right)
            else:pos=left;allele_ref=ref.fetch(chrom,left-1,left);alt=allele_ref+event['inserted_sequence']
            baseline_rows.append((chrom,pos,event['id'],allele_ref,alt))
            nearby=[r for r in represented if r['filters']==['PASS'] and r['chrom']==chrom and
                    abs(r['pos']-pos)<=150 and (len(r['alt'])-len(r['ref']))*(len(alt)-len(allele_ref))>0]
            for r in nearby:
                start=min(pos-1,r['pos']-1)-100;end=max(pos-1+len(allele_ref),r['pos']-1+len(r['ref']))+100
                expected=ref.fetch(chrom,start,pos-1)+alt+ref.fetch(chrom,pos-1+len(allele_ref),end)
                observed=ref.fetch(chrom,start,r['pos']-1)+r['alt']+ref.fetch(chrom,r['pos']-1+len(r['ref']),end)
                direct.append({'event_id':event['id'],'native_id':r['id'],'same_local_haplotype':observed==expected,
                               'native_GT':r['gt'],'truth_GT':event['expected_genotype'],
                               'same_genotype':sorted(r['gt'])==sorted(event['expected_genotype'])})
        if {x['event_id'] for x in direct}!={'deletion','tandem_duplication','insertion'} or not all(x['same_local_haplotype'] for x in direct):
            raise ValueError('direct synthetic haplotype conversion check failed')
        header=pysam.VariantHeader()
        for chrom,length in zip(ref.references,ref.lengths):header.contigs.add(chrom,length=length)
        header.formats.add('GT',1,'String','Known synthetic genotype');header.add_sample('SYNTH_SV')
        base=output/'linear-truth.vcf.gz'
        with pysam.VariantFile(str(base),'wz',header=header) as writer:
            for chrom,pos,name,r,a in baseline_rows:
                v=writer.new_record(contig=chrom,start=pos-1,id=name,alleles=(r,a));v.filter.add('PASS')
                v.samples['SYNTH_SV']['GT']=(0,1);writer.write(v)
        bed=output/'regions.bed';bed.write_text(''.join(f'{c}\t0\t{n}\n' for c,n in zip(ref.references,ref.lengths)))
    pysam.tabix_index(str(base),preset='vcf')
    config=json.loads(config_path.read_text())
    benchmark=compare(root,base,converted/'comparison.full-context.vcf.gz',reference,bed,output/'benchmark',config)
    counts=benchmark['refined']
    if any(counts[k]!=v for k,v in {'TP-base':3,'TP-comp':3,'FN':0,'FP':0}.items()):
        raise ValueError('native converted linear-event comparator control failed')
    result={'completed_at':stamp(),'scope':'Independent engineered truth; not a human-reference sensitivity estimate',
            'direct_linear_haplotype_checks':direct,'PASS_native_truth_comparison':truth_check,
            'all_native_truth_comparison':compare_truth(truth,events),'native_comparator':benchmark,
            'truth':fingerprint(fixture/'truth.json'),'code':fingerprint(Path(__file__)),
            'conversion_code':fingerprint(Path(__file__).with_name('convert_svaba_pilot.py')),'configuration':fingerprint(config_path),
            'interpretation':['Linear deletion, duplication and 100-base insertion match the planted local haplotypes',
                              'Do not treat detected location as proof that both inversion connections or the correct genotype were recovered',
                              'All misses, filtered calls, unmatched calls and the unchanged chromosome are retained in the event ledger']}
    save(output/'summary.json',result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('root','output','config'):p.add_argument('--'+name,required=True,type=Path)
    p.add_argument('--conversion',type=Path)
    a=p.parse_args();r=run(a.root,a.output,a.config,a.conversion)
    print(json.dumps({'completed_at':r['completed_at'],'direct_linear_haplotype_checks':r['direct_linear_haplotype_checks'],
                      'native_comparator':r['native_comparator']['refined'],'PASS_native_truth_comparison':r['PASS_native_truth_comparison']}))
