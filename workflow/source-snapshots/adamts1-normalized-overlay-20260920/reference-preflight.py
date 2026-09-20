#!/usr/bin/env python3
"""Read variant/reference inputs and construct a bounded preparation; no reads."""
import argparse, gzip, hashlib, importlib.util, json, os, sys, time
from pathlib import Path

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('normalized_overlay',HERE/'normalized-overlay.py')
overlay=importlib.util.module_from_spec(spec);spec.loader.exec_module(overlay)

def run(args):
    import pysam
    root=args.root.resolve();output=args.output.resolve();method=root/'methods/adamts1-joint-review-20260920'
    overlay.require(output.is_relative_to(method) and not output.exists(),'new owned preparation output required')
    config=json.loads(args.config.read_text());selection=json.loads(args.spec.read_text());pins={}
    def pin(path,expected=None):
        value=overlay.fingerprint(path,expected);pins[value['path']]=value;return value
    config_pin=pin(args.config);pin(args.spec);pin(__file__);pin(HERE/'normalized-overlay.py')
    scripts=root/'methods/joint-haplotype-review-20260920/qualification-v1/workflow/scripts'
    for name,digest in config['source_sha256'].items():
        p=pin(scripts/name);overlay.require(p['sha256']==digest,'qualified source differs')
    sys.path.insert(0,str(scripts))
    import review_joint_haplotypes as native
    from native_coding_records import source_context
    from joint_haplotype_helpers import enumerate_haplotypes
    target=config['targets'][0]
    overlay.require(len(config['targets'])==1 and target['id']==selection['target_id'],'target selection differs')
    overlay.same({k:target[k] for k in ('chrom','core_start0','core_end0')},{k:selection[k] for k in ('chrom','core_start0','core_end0')},'core selection differs')
    overlay.same(config['window_flanks_bp'],[3,10,25],'flanks differ')
    overlay.same(config['mapping_quality_thresholds'],[0,20,60],'MQ views differ')
    overlay.same(config['base_quality_thresholds'],[20,30],'BQ views differ')
    overlay.require(config['raw_base_qualities_no_BAQ'] is True and config['maximum_hypothesis_states']==256,'quality or state policy differs')
    input_pins={}
    for name,expected in config['inputs'].items():
        if name not in ('alignment','alignment_index'): input_pins[name]=pin(root/expected['path'],expected)
    inventory_pin=pin(root/selection['normalized_inventory_input']['path'],selection['normalized_inventory_input'])
    normalized=[]
    with gzip.open(inventory_pin['path'],'rt') as stream:
        for line in stream:
            row=json.loads(line)
            if row['key']==selection['normalized_key']: normalized.append(row)
    overlay.require(len(normalized)==1,'normalized inventory identity missing or duplicated')
    normalized=normalized[0]
    overlay.same(normalized['allele'],selection['normalized_deletion'],'normalized allele differs')
    core_start,core_end=target['core_start0'],target['core_end0'];entries=[];outer=[];census=[]
    for source in config['caller_sources']:
        source_pin=input_pins[source['input']];digest=hashlib.sha256();count=0
        with pysam.VariantFile(source_pin['path']) as reader:
            overlay.same(list(reader.header.samples),[config['sample']],'native sample differs')
            for ordinal,record in enumerate(reader,1):
                count+=1;digest.update(str(record).encode());begin=record.pos-1;end=max(record.stop,begin+len(record.ref))
                if record.contig!=target['chrom'] or begin>=core_end+500 or end<=core_start-500: continue
                context=source_context(record,config['sample'],source['caller'],source_pin['sha256'],ordinal)
                row={'context':context,'reported_end0':record.stop};outer.append(row)
                if begin<core_end+25 and end>core_start-25: entries.append(row)
        overlay.require(count==source['expected_records'] and digest.hexdigest()==source['expected_rendered_record_payload_sha256'],'source count or ordered payload differs')
        census.append({'caller':source['caller'],'input':source_pin,'records':count,'rendered_record_payload_sha256':digest.hexdigest()})
    contexts={row['context']['source_record_id']:row['context'] for row in entries}
    overlay.require(len(contexts)==len(entries)==6,'six distinct full native records required')
    overlay.same(sorted(contexts),sorted(target['required_source_record_ids']),'required six source identities differ')
    deletion_id=selection['roles']['raw_deletion']['source_record_id']
    overlay.same(normalized['source_call']['native_context'],contexts[deletion_id],'normalized/original full context binding differs')
    overlay.same(normalized['source_call']['source_alt_index'],selection['roles']['raw_deletion']['source_alt_index'],'normalized/original ALT index differs')
    windows=[]
    with pysam.FastaFile(input_pins['reference']['path']) as fasta:
        for flank in config['window_flanks_bp']:
            start,end=core_start-flank,core_end+flank;reference=fasta.fetch(target['chrom'],start,end).upper()
            edits,omitted=native.window_edits(entries,start,end)
            family=enumerate_haplotypes(target['chrom'],start,reference,edits,callers=['gatk','bcftools'],max_states=256,max_haplotype_bases=config['maximum_haplotype_bases'])
            overlay.same(family['counters'],{'attempted_combinations':16,'compatible_combinations':11,'conflicting_REF_span_combinations':5,'distinct_sequences':9},'bounded native enumeration differs')
            overlay.require(len(omitted)==1 and omitted[0]['alt']=='*' and omitted[0]['source_alt_index']==1 and omitted[0]['reason']=='nonconcrete_ALT_or_REF','spanning-star omission differs')
            lineages=overlay.overlay_lineages(reference,start,selection,contexts)
            original=[dict(row,id='preparation_F%d_H%d'%(flank,index)) for index,row in enumerate(family['haplotypes'])]
            classes=overlay.union_classes(original,lineages['lineages'],'preparation_F'+str(flank))
            windows.append({'chrom':target['chrom'],'start0':start,'end0':end,'flank_bp':flank,'reference_sequence':reference,
                'native_family':family,'omitted_native_alleles':omitted,'overlay':lineages,'union_classes':list(classes.values()),
                'reference_bases':len(reference),'hypothesis_length_range':[min(map(len,classes)),max(map(len,classes))]})
    for value in list(pins.values()): overlay.same(overlay.fingerprint(value['path']),value,'preparation input changed')
    result={'status':'complete','mode':'reference_only_preparation','completed_at':overlay.stamp(),'native_config_pin':config_pin,
        'native_config_body':config,'overlay_spec':selection,'verified_files':pins,'native_source_census':census,
        'contexts':contexts,'outer_context':outer,'outer_scope':{'chrom':target['chrom'],'start0':core_start-500,'end0':core_end+500},
        'normalized_inventory_row':normalized,'windows':windows,
        'not_opened_inputs':['alignment','alignment_index'],
        'limitations':['Reference and native-variant preparation only; no alignment read extraction, observer run, annotation, caller or patient-support count.',
            'Alignment and alignment-index pins are inherited for later fresh verification, not opened in this preparation.',
            'Constructed sequence identity does not establish genotype, phase, origin, function, pathogenicity or causality.']}
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('x') as stream:json.dump(result,stream,indent=2,sort_keys=True);stream.write('\n')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--root',type=Path,required=True);parser.add_argument('--config',type=Path,required=True)
    parser.add_argument('--spec',type=Path,required=True);parser.add_argument('--output',type=Path,required=True);run(parser.parse_args())
