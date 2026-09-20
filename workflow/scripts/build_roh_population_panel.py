#!/usr/bin/env python3
"""Select independently sourced common markers for subsequent explicit genotyping."""
import argparse
from array import array
from bisect import bisect_right
from collections import Counter
import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import time

import pysam


def hashes(path, include_md5=False):
    before=path.stat();sha=hashlib.sha256();md5=hashlib.md5() if include_md5 else None
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):
            sha.update(block)
            if md5 is not None:md5.update(block)
    after=path.stat()
    if (before.st_size,before.st_mtime_ns)!=(after.st_size,after.st_mtime_ns):
        raise ValueError('input changed while hashing')
    result={'path':str(path.resolve()),'bytes':after.st_size,'mtime_ns':after.st_mtime_ns,'sha256':sha.hexdigest()}
    if md5 is not None:result['md5']=md5.hexdigest()
    return result


def load_mask(path,lengths):
    intervals={chrom:(array('I'),array('I')) for chrom in lengths}
    counts=Counter()
    with path.open() as stream:
        for line in stream:
            if not line.strip() or line.startswith('#') or line.startswith('track '):continue
            fields=line.split()
            if len(fields)!=4 or fields[3]!='strict':raise ValueError('expected passed-sites BED with strict label')
            chrom=fields[0].removeprefix('chr')
            if chrom not in lengths:
                counts['non_autosomal_intervals']+=1;continue
            start,end=map(int,fields[1:3])
            if not 0<=start<end<=lengths[chrom]:raise ValueError('mask interval outside exact reference')
            starts,ends=intervals[chrom]
            if ends and start<ends[-1]:raise ValueError('overlapping or unsorted mask intervals')
            counts['input_autosomal_intervals']+=1
            if ends and start==ends[-1]:ends[-1]=end
            else:starts.append(start);ends.append(end)
    if any(not starts for starts,ends in intervals.values()):raise ValueError('mask missing a requested autosome')
    return intervals,dict(counts)


def included(intervals,chrom,pos0):
    starts,ends=intervals[chrom];i=bisect_right(starts,pos0)-1
    return i>=0 and pos0<ends[i]


def one_af(record,key):
    values=record.info.get(key)
    if values is None or len(values)!=1:return None
    value=values[0]
    if value is None or not math.isfinite(value) or not 0<=value<=1:return None
    return value


def build(source,reference,mask,config,output):
    output.mkdir(parents=True,exist_ok=False)
    started=time.monotonic()
    inputs={'source':hashes(source,True),'reference':hashes(reference),'mask':hashes(mask),
            'reference_index':hashes(Path(str(reference)+'.fai')),
            'runner':hashes(Path(__file__))}
    if inputs['source']['bytes']!=config['source_bytes'] or inputs['source']['md5']!=config['source_publisher_md5']:
        raise ValueError('population resource does not match publisher manifest')
    if inputs['reference']['sha256']!=config['reference_sha256']:
        raise ValueError('exact reference checksum mismatch')
    if inputs['source']['sha256']!=config['source_sha256'] or inputs['mask']['sha256']!=config['mask_sha256']:
        raise ValueError('population resource or mask differs from pinned source')
    if config['scope']!='public_population_marker_preparation' or config['mask_policy']!='include_IGSR_20160622_GRCh38_StrictMask_passed_BED':
        raise ValueError('unreviewed panel scope or mask policy')
    chroms=config['chromosomes']
    if chroms!=list(map(str,range(1,23))):raise ValueError('all 22 autosomes required')
    pops=config['superpopulations']
    if set(pops)!={'AFR','AMR','EAS','EUR','SAS'}:raise ValueError('five population frequency fields required')
    if config['minimum_spacing_bp']<1 or not 0<config['minimum_global_minor_AF']<0.5 or not 0<config['minimum_superpopulation_minor_AF']<0.5:
        raise ValueError('invalid common-marker selection thresholds')
    counts=Counter();selected=Counter();ns_counts=Counter();an_counts=Counter();last_selected={}
    af_keys=['AF']+[p+'_AF' for p in pops]
    with pysam.FastaFile(str(reference)) as fasta,pysam.VariantFile(str(source)) as reader:
        if len(reader.header.samples):raise ValueError('expected public sites-only resource')
        lengths={c:fasta.get_reference_length(c) for c in chroms}
        intervals,mask_counts=load_mask(mask,lengths)
        source_header={'reference':[str(r).strip() for r in reader.header.records if r.key=='reference'],
                       'AF_fields':{k:reader.header.info[k].description for k in af_keys},
                       'contigs':list(reader.header.contigs),'sample_columns':len(reader.header.samples)}
        header=pysam.VariantHeader()
        header.add_meta('source','MVA independently selected IGSR common-marker panel')
        for c in chroms:header.contigs.add(c,length=lengths[c])
        for key in af_keys:
            header.info.add('POP_'+key,1,'Float','Source population allele frequency: '+key)
        header.info.add('POP_AN',1,'Integer','Source global called-allele count')
        header.info.add('POP_NS',1,'Integer','Source number of samples with data')
        panel_path=output/'panel.vcf.gz';bed_path=output/'sites.bed';allele_path=output/'alleles.tsv.gz'
        with pysam.VariantFile(str(panel_path),'wz',header=header) as writer,bed_path.open('w') as bed,pysam.BGZFile(str(allele_path),'wb') as alleles:
            def process(group):
                if not group:return
                counts['source_records']+=len(group)
                r=group[0];c=r.contig
                if c not in lengths:counts['non_autosomal_records']+=len(group);return
                if len(group)!=1:counts['duplicate_position_records']+=len(group);return
                if not 1<=r.pos<=lengths[c]:raise ValueError('source position outside reference bounds')
                if tuple(r.filter)!=('PASS',):counts['non_PASS']+=1;return
                if len(r.alleles)!=2 or r.ref==r.alts[0] or any(len(a)!=1 or a not in 'ACGT' for a in r.alleles):
                    counts['not_biallelic_SNP']+=1;return
                counts['PASS_autosomal_biallelic_SNP']+=1
                values=[one_af(r,key) for key in af_keys]
                an,ns=r.info.get('AN'),r.info.get('NS')
                if any(v is None for v in values) or an is None or ns is None:
                    counts['missing_frequency_or_sample_count']+=1;return
                if an<config['minimum_AN']:
                    counts['insufficient_population_AN']+=1;return
                eps=config['floating_point_tolerance']
                bounds=[config['minimum_global_minor_AF']]+[config['minimum_superpopulation_minor_AF']]*len(pops)
                if any(min(v,1-v)+eps<low for v,low in zip(values,bounds)):
                    counts['outside_common_AF_bounds']+=1;return
                counts['common_in_all_selected_populations']+=1
                if not included(intervals,c,r.pos-1):counts['outside_passed_mask']+=1;return
                counts['inside_passed_mask']+=1
                if c in last_selected and r.pos-last_selected[c]<config['minimum_spacing_bp']:
                    counts['excluded_spacing']+=1;return
                if fasta.fetch(c,r.pos-1,r.pos).upper()!=r.ref:raise ValueError('selected REF differs from exact reference')
                rec=writer.new_record(contig=c,start=r.pos-1,alleles=r.alleles,id=r.id)
                rec.filter.add('PASS')
                for key,value in zip(af_keys,values):rec.info['POP_'+key]=value
                rec.info['POP_AN']=an;rec.info['POP_NS']=ns
                writer.write(rec);bed.write(f'{c}\t{r.pos-1}\t{r.pos}\n')
                alleles.write(f'{c}\t{r.pos}\t{r.ref},{r.alts[0]}\n'.encode())
                selected[c]+=1;last_selected[c]=r.pos;an_counts[str(an)]+=1;ns_counts[str(ns)]+=1
            order={c:i for i,c in enumerate(reader.header.contigs)}
            previous=None;group=[]
            for record in reader:
                key=(order[record.contig],record.pos)
                if previous is not None and key<previous:raise ValueError('source VCF order changed')
                if previous is not None and key!=previous:process(group);group=[]
                group.append(record.copy());previous=key
            process(group)
        if set(selected)!=set(chroms):raise ValueError('selected panel missing an autosome')
        pysam.tabix_index(str(panel_path),preset='vcf',force=False)
        pysam.tabix_index(str(allele_path),seq_col=0,start_col=1,end_col=1,force=False)
        mask_bases={c:sum(end-start for start,end in zip(*intervals[c])) for c in chroms}
    for identity in inputs.values():
        st=Path(identity['path']).stat()
        if (st.st_size,st.st_mtime_ns)!=(identity['bytes'],identity['mtime_ns']):raise ValueError('input changed during panel build')
    result={'completed_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'elapsed_seconds':time.monotonic()-started,'config':config,'inputs':inputs,'source_header':source_header,
        'counts':dict(counts),'selected_markers_by_chromosome':dict(selected),'selected_markers':sum(selected.values()),
        'selected_source_AN_distribution':dict(an_counts),'selected_source_NS_distribution':dict(ns_counts),
        'mask_counts':mask_counts,'passed_mask_bases_by_chromosome':mask_bases,
        'outputs':{name:hashes(output/name) for name in ['panel.vcf.gz','panel.vcf.gz.tbi','sites.bed','alleles.tsv.gz','alleles.tsv.gz.tbi']},
        'scope':'Public genotype-independent panel preparation only; all selected REF bases verified. Not patient genotypes, ancestry assignment, linkage independence, callability or a completed ROH analysis.'}
    (output/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ['source','reference','mask','config','output']:p.add_argument('--'+arg,type=Path,required=True)
    args=p.parse_args();os.nice(10)
    result=build(args.source,args.reference,args.mask,json.loads(args.config.read_text()),args.output)
    print(json.dumps({k:result[k] for k in ['completed_at','elapsed_seconds','selected_markers','selected_markers_by_chromosome','counts','selected_source_AN_distribution','selected_source_NS_distribution','scope']},indent=2))


if __name__=='__main__':main()
