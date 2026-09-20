#!/usr/bin/env python3
"""Independent exact-allele inventory audit. No producer or vendor imports.

Repeats the pinned bcftools normalization engine: technical reproduction, not an
independent normalization algorithm or biological validation.
"""
import argparse
from collections import Counter, defaultdict
from contextlib import contextmanager
import datetime
import gzip
import hashlib
import itertools
import json
import math
from pathlib import Path
import platform
import re
import subprocess
import sys
import time

import pysam

CHROMS = list(map(str, range(1, 23)))
CALLERS = ['gatk', 'bcftools']
PLP = {'Pathogenic', 'Likely_pathogenic', 'Pathogenic/Likely_pathogenic'}
STARS = {'practice_guideline': 4, 'reviewed_by_expert_panel': 3,
 'criteria_provided,multiple_submitters,no_conflicts': 2,
 'criteria_provided,single_submitter': 1,
 'criteria_provided,conflicting_classifications': 1,
 'criteria_provided,conflicting_interpretations': 1,
 'no_assertion_criteria_provided': 0, 'no_classification_provided': 0,
 'no_assertion_provided': 0, 'no_classification_for_the_individual_variant': 0,
 'no_classification_for_the_single_variant': 0,
 'no_classifications_from_unflagged_records': 0,
 'no_assertion_for_the_individual_variant': 0, '': 0}
FROZEN_PUBLIC = {'clinvar': '47ea3626cd90154d82304afe3b8f717aab2557a1330f4ea2c809892386e08d49',
 'clinvar_index': 'b220738d2c086062f136600cdd31a95bec54586890cb7d92868ad78bdc632c03',
 'panel': '8089748ab8645336eedab56e4f548f7c902e2f7d163fa3a559b49c9fbd93151d'}


def check(condition, message):
    if not condition: raise ValueError(message)


def canon(value): return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def equal(left, right, message):
    check(canon(left) == canon(right), message)


def parse_json(text):
    def pairs(items):
        result = {}
        for key, value in items:
            check(key not in result, 'duplicate JSON key')
            result[key] = value
        return result
    def bad(value): raise ValueError('nonfinite JSON: ' + value)
    return json.loads(text, object_pairs_hook=pairs, parse_constant=bad)


def read_json(path): return parse_json(Path(path).read_text())


def write_json(path, value):
    with Path(path).open('x') as stream: stream.write(json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n')


def now(): return datetime.datetime.now(datetime.timezone.utc).isoformat()


class Files:
    def __init__(self): self.pins = {}
    def hash(self, path, expected=None):
        path = Path(path).resolve(); key = str(path)
        if expected and 'path' in expected:
            equal(key, str(Path(expected['path']).resolve()), 'pin path mismatch')
        if key not in self.pins:
            old = path.stat(); digest = hashlib.sha256()
            with path.open('rb') as stream:
                for block in iter(lambda: stream.read(8*1024*1024), b''): digest.update(block)
            new = path.stat()
            equal([old.st_size, old.st_mtime_ns], [new.st_size, new.st_mtime_ns], 'changed during hash')
            self.pins[key] = dict(path=key, bytes=new.st_size, mtime_ns=new.st_mtime_ns, sha256=digest.hexdigest())
        result = self.pins[key]
        if expected:
            for field in ['sha256', 'bytes']:
                if field in expected: equal(result[field], expected[field], 'pin '+field+' mismatch: '+key)
        return result
    def get(self, pin):
        self.hash(pin['path'], pin)
        return read_json(pin['path'])
    def stable(self):
        for pin in self.pins.values():
            stat = Path(pin['path']).stat()
            equal([stat.st_size, stat.st_mtime_ns], [pin['bytes'], pin['mtime_ns']], 'input changed during audit')


def proof_bind(pin, entries, files, root=None):
    lookup = {}
    for item in entries:
        path = Path(item['path'])
        if not path.is_absolute():
            check(root is not None, 'relative proof without root')
            path = root / path
        key = str(path.resolve())
        check(key not in lookup, 'duplicate proof path')
        lookup[key] = dict(item, path=key)
    key = str(Path(pin['path']).resolve())
    check(key in lookup, 'missing proof member')
    files.hash(key, lookup[key]); files.hash(key, pin)


def plain(value):
    if isinstance(value, (tuple, list)): return [plain(v) for v in value]
    if isinstance(value, dict): return {str(k): plain(v) for k,v in value.items()}
    check(value is None or type(value) in (str,int,bool,float), 'unsupported field type')
    if type(value) is float: check(math.isfinite(value), 'nonfinite field')
    return value


def text(info, name):
    value = info.get(name)
    return '' if value is None else ','.join(map(str,value)) if isinstance(value,(tuple,list)) else str(value)


def symbols(info): return {part.rsplit(':',1)[0] for part in text(info,'GENEINFO').split('|') if ':' in part}


def signal(value):
    tokens = re.split('[|,]', value)
    for token in tokens:
        token = re.sub(r'^\d+:', '', token.strip('_ '))
        token = re.sub(r'_?\(\d+\)$', '', token).strip('_ ')
        if token in PLP: return True
    return False


def classification(info):
    value = text(info,'CLNSIG'); conflict = text(info,'CLNSIGCONF'); included = text(info,'CLNSIGINCL')
    clean = bool(value) and set(re.split('[|,]',value)).issubset(PLP)
    if not (clean or any(map(signal,[value,conflict,included]))): return None
    status = text(info,'CLNREVSTAT').replace(',_',',')
    check(status in STARS, 'unknown assertion review status')
    stars = STARS[status]
    if clean and stars >= 2: name = 'aggregate_plp_2plus_stars'
    elif clean and stars == 1 and 'conflict' not in status: name = 'aggregate_plp_1star'
    elif clean and stars == 0: name = 'aggregate_plp_0star_review'
    elif signal(value) or signal(conflict): name = 'conflicting_or_mixed_classification_review'
    else: name = 'haplotype_or_genotype_only_review'
    return name, stars


def concrete(value): return isinstance(value,str) and re.fullmatch('[ACGT]+',value) is not None


def allele(record):
    check(len(record.alts or ()) == 1, 'not biallelic')
    return (record.contig, record.pos, record.ref, record.alts[0])


def original_index(normalized, context):
    raw = normalized.info.get('SF_ORIG')
    if raw is None:
        equal([normalized.contig,normalized.pos,normalized.ref,list(normalized.alts or ())],
              [context['chrom'],context['pos'],context['ref'],context['alts']], 'untagged source differs')
        return 1
    if isinstance(raw,tuple):
        check(all(type(v) is str for v in raw), 'bad old-record tuple'); raw=','.join(raw)
    check(type(raw) is str, 'bad old-record type')
    parts=raw.split('|')
    if len(parts)==4 and ',' not in parts[3]: parts.append('1')
    check(len(parts)==5, 'old-record field count')
    chrom,pos,ref,alts,idx=parts
    idx=int(idx)
    equal([chrom,int(pos),ref,alts.split(',')], [context['chrom'],context['pos'],context['ref'],context['alts']], 'old-record lineage differs')
    check(1<=idx<=len(context['alts']), 'old-record index outside ALTs')
    return idx


def literal(fasta, normalized, context, index):
    check(normalized.contig==context['chrom'], 'normalized chromosome differs')
    start=min(normalized.pos,context['pos'])-1
    end=max(normalized.pos-1+len(normalized.ref),context['pos']-1+len(context['ref']))
    ref=fasta.fetch(normalized.contig,start,end).upper()
    check(len(ref)==end-start,'short reference window')
    result=[]
    for pos,r,a in [(normalized.pos,normalized.ref,normalized.alts[0]),(context['pos'],context['ref'],context['alts'][index-1])]:
        offset=pos-1-start; equal(ref[offset:offset+len(r)],r,'REF mismatch')
        result.append(ref[:offset]+a+ref[offset+len(r):])
    equal(*result,'different literal edit')


@contextmanager
def norm(executable, reference, source, log):
    command=[str(executable),'norm','-f',str(reference),'-c','e','-m','-any','--old-rec-tag','SF_ORIG','-Ou',str(source)]
    with Path(log).open('xb') as stderr:
        process=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=stderr)
        try:
            with pysam.VariantFile(process.stdout) as reader: yield reader,command
            process.stdout.close()
            check(process.wait()==0,'normalization failed: '+str(log))
        finally:
            if process.poll() is None:
                process.terminate()
                try: process.wait(timeout=10)
                except subprocess.TimeoutExpired: process.kill();process.wait()


def regions_from_model(model, fasta):
    equal([r['chrom'] for r in model['regions']],['X','Y'],'X/Y order')
    for region in model['regions']:
        equal(region['length'],fasta.get_reference_length(region['chrom']),'X/Y reference extent')
        end=0
        for start,stop,ploidy in region['intervals_0based_half_open']:
            check(all(type(v) is int for v in (start,stop,ploidy)),'interval types')
            check(start==end and stop>start and ploidy in (1,2),'ploidy tiling')
            if region['chrom']=='Y': equal(ploidy,1,'Y ploidy')
            end=stop
        equal(end,region['length'],'ploidy extent')
    check({v[2] for v in model['regions'][0]['intervals_0based_half_open']}=={1,2},'X PAR model')
    check(bool(model['required_N_masks']),'missing Y masks')
    for mask in model['required_N_masks']:
        check(mask['chrom']=='Y' and type(mask['start']) is int and type(mask['end']) is int,'mask coordinates')
        check(0<=mask['start']<mask['end']<=fasta.get_reference_length('Y'),'mask interval')
        check(set(fasta.fetch('Y',mask['start'],mask['end']).upper())=={'N'},'Y mask content')
    return model['regions']


def model_ploidy(record,source,regions):
    check(record.contig in source['contigs'],'record outside source')
    if source['kind']=='autosome': return 2
    choices=[(end,p) for region in regions if region['chrom']==record.contig
             for start,end,p in region['intervals_0based_half_open'] if start<=record.start<end]
    check(len(choices)==1 and record.stop<=choices[0][0],'REF crosses model boundary')
    return choices[0][1]


def sources_from_proofs(config,files,fasta):
    inputs=config['inputs']; root=Path(config['root']).resolve()
    handoff=files.get(inputs['handoff']); verified=files.get(inputs['autosome_verification'])
    equal(handoff['status'],'complete','handoff incomplete'); equal(verified['status'],'complete','native audit incomplete')
    for scope in [handoff,verified['scope']]:
        equal(scope['selected_autosomes'],CHROMS,'autosome membership')
        equal(scope['pending_autosomes'],[],'pending autosomes')
        equal(scope['all_autosomes_complete'],True,'completion flag')
    equal(handoff['sample'],config['sample'],'sample binding')
    proof_bind(inputs['handoff'],verified['files'],files)
    for label in ['reference','reference_index','bcftools']: files.hash(inputs[label]['path'],handoff['inputs'][label])
    equal([s['shard_id'] for s in handoff['shards']],list(range(22)),'shard completeness')
    equal([c for s in handoff['shards'] for c,_ in s['contigs']],CHROMS,'chromosome completeness')
    sources=[]
    for shard in handoff['shards']:
        check(len(shard['contigs'])==1,'ambiguous shard')
        chrom,length=shard['contigs'][0]; equal(length,fasta.get_reference_length(chrom),'autosome extent')
        equal(sorted(shard['methods']),sorted(CALLERS),'caller membership')
        for caller in CALLERS:
            method=shard['methods'][caller]
            equal(method['replayed_split_payload_identical'],True,'missing native payload proof')
            proof_bind(method['original_raw_pin'],verified['files'],files)
            sources.append(dict(caller=caller,kind='autosome',contigs=[chrom],pin=files.hash(method['original_raw_pin']['path'],method['original_raw_pin']),expected_census=method['original_raw'],id=caller+'-'+chrom))
    status=files.get(inputs['fixed_xy_status']); review=files.get(inputs['fixed_xy_review']); proof=files.get(inputs['fixed_xy_verification'])
    equal(status['status'],'complete','fixed status incomplete');equal(review['status'],'complete','fixed review incomplete')
    equal(review['fixed_run_status_sha256'],inputs['fixed_xy_status']['sha256'],'fixed status link')
    a=datetime.datetime.fromisoformat(review['verified_at'].replace('Z','+00:00'));b=datetime.datetime.fromisoformat(proof['verified_at'].replace('Z','+00:00'))
    check(a.tzinfo is not None and b.tzinfo is not None and b>=a,'fixed terminal times')
    check(re.fullmatch('[a-f0-9]{64}',review['script_sha256']) is not None,'fixed review source hash')
    equal(proof['script_sha256'],review['script_sha256'],'fixed review source binding')
    equal(review['fixed_outputs_verified'],len(status['output_sha256']),'fixed output count')
    check(type(review['original_files_preserved']) is int and review['original_files_preserved']>0,'fixed original preservation')
    equal(review['bcftools_original_reused_byte_for_byte'],True,'bcftools fixed preservation')
    equal(proof['fixed_outputs'],status['output_sha256'],'fixed output map')
    proof_bind(inputs['fixed_xy_review'],proof['files'],files,root)
    equal(status['config']['sample'],config['sample'],'fixed sample')
    equal(status['config']['reference_sha256'],inputs['reference']['sha256'],'fixed reference')
    regions=regions_from_model(status['config'],fasta)
    folder=Path(inputs['fixed_xy_status']['path']).parent
    for caller in CALLERS:
        path=folder/(caller+'.original.vcf.gz')
        pin=files.hash(path,{'sha256':status['output_sha256'][path.name]})
        files.hash(str(path)+'.tbi',{'sha256':status['output_sha256'][path.name+'.tbi']})
        count=review['raw_record_accounting']['fixed_gatk' if caller=='gatk' else 'bcftools']['records']
        sources.append(dict(caller=caller,kind='fixed_xy',contigs=['X','Y'],pin=pin,expected_census={'counts':{'records':count}},id=caller+'-fixed-xy'))
    equal(len({s['pin']['path'] for s in sources}),46,'distinct native source paths')
    return sources,regions


def public_inventory(config,run,out,files,fasta):
    inputs=config['inputs']; panel=files.get(inputs['panel'])
    genes=[row['symbol'] for row in panel['rows']]
    equal(len(genes),config['expected_genes'],'panel length'); equal(len(set(genes)),len(genes),'duplicate panel gene')
    if 'total' in panel: equal(panel['total'],len(genes),'panel total')
    counts=Counter({k:0 for k in ['original_records','panel_gene_records','selected_source_records','selected_biallelic_source_records','selected_multiallelic_source_records','normalized_lineages','signal_nonprimary_or_nonconcrete_records','selected_nonconcrete_normalized_lineages']})
    contexts={}
    with pysam.VariantFile(inputs['clinvar']['path']) as raw, pysam.VariantFile(run/'public-selected.vcf.gz') as selected:
        check('##reference=GRCh38\n' in str(raw.header),'public assembly')
        equal([s.split('=',1)[1] for s in str(raw.header).splitlines() if s.startswith('##fileDate=')],[config['public_snapshot']['clinvar_fileDate']],'public date')
        equal((run/'original-clinvar-header.txt').read_text(),str(raw.header),'public header preservation')
        for ordinal,record in enumerate(raw,1):
            counts['original_records']+=1
            if not symbols(record.info).intersection(genes): continue
            counts['panel_gene_records']+=1
            if classification(record.info) is None: continue
            if record.contig not in CHROMS+['X','Y'] or not concrete(record.ref) or not any(map(concrete,record.alts or ())):
                counts['signal_nonprimary_or_nonconcrete_records']+=1;continue
            copy=next(selected,None); check(copy is not None,'missing public selected record')
            equal(copy.info['NS_PUBLIC_ORD'],ordinal,'public selection ordinal')
            del copy.info['NS_PUBLIC_ORD'];equal(str(copy),str(record),'public selected original payload')
            contexts[ordinal]=dict(source_record_id='clinvar:'+inputs['clinvar']['sha256']+':'+str(ordinal),source_sha256=inputs['clinvar']['sha256'],source_ordinal=ordinal,chrom=record.contig,pos=record.pos,ref=record.ref,alts=list(record.alts),rendered_record=str(record),INFO=plain(dict(record.info)),id=record.id,QUAL=plain(record.qual),FILTER=list(record.filter))
            counts['selected_source_records']+=1
            counts['selected_biallelic_source_records' if len(record.alts)==1 else 'selected_multiallelic_source_records']+=1
        check(next(selected,None) is None,'extra public selected record')
    lookup=defaultdict(list);seen=set()
    with norm(inputs['bcftools']['path'],inputs['reference']['path'],run/'public-selected.vcf.gz',out/'public-renormalization.stderr') as (reader,command):
        for record in reader:
            ordinal=record.info['NS_PUBLIC_ORD'];check(type(ordinal) is int and ordinal in contexts,'public normalized ordinal')
            if not all(map(concrete,record.alleles)):
                counts['selected_nonconcrete_normalized_lineages']+=1;continue
            context=contexts[ordinal];idx=original_index(record,context);literal(fasta,record,context,idx)
            key=allele(record);unique=(ordinal,idx,key);check(unique not in seen,'duplicate public lineage');seen.add(unique)
            source_assertion=classification(context['INFO']);check(source_assertion is not None,'original assertion absent')
            source_stratum,stars=source_assertion;normalized_assertion=classification(record.info);ambiguous=len(context['alts'])>1
            if not ambiguous: equal(normalized_assertion,source_assertion,'biallelic assertion changed')
            lookup[key].append(dict(source=context,source_alt_index=idx,normalized_record=str(record),normalized_allele=dict(zip(['chrom','pos','ref','alt'],key)),stratum='multiallelic_assertion_mapping_review' if ambiguous else source_stratum,source_record_stratum=source_stratum,stars=stars,normalized_record_stratum=normalized_assertion[0] if normalized_assertion else None,normalized_record_stars=normalized_assertion[1] if normalized_assertion else None,assertion_mapping=dict(status='ambiguous_multiallelic_source' if ambiguous else 'unambiguous_biallelic_source',original_CLNALLE=context['INFO'].get('CLNALLE'),normalized_CLNALLE=plain(record.info['CLNALLE']) if 'CLNALLE' in record.info else None),genes=sorted(symbols(record.info).intersection(genes)),classification=text(record.info,'CLNSIG'),review_status=text(record.info,'CLNREVSTAT'),public_original_multiallelic=ambiguous,interpretation='Exact dated assertion context; condition/variant/inheritance reportability remains pending.'))
            counts['normalized_lineages']+=1
    counts['unique_normalized_identities']=len(lookup)
    return lookup,dict(counts),sorted(genes),command


def native_context(record,source,sample,ordinal):
    calls=record.samples[sample];fields=plain(dict(calls));gt=fields.get('GT') or []
    check(all(v is None or type(v) is int and 0<=v<=len(record.alts or ()) for v in gt),'invalid GT allele')
    return dict(source_record_id=source['caller']+':'+source['pin']['sha256']+':'+str(ordinal),source_namespace=json.dumps([source['caller'],sample,source['pin']['sha256']],separators=(',',':')),caller=source['caller'],sample=sample,source_sha256=source['pin']['sha256'],source_ordinal=ordinal,chrom=record.contig,pos=record.pos,ref=record.ref,alts=list(record.alts or ()),FORMAT=fields,GT_phased=bool(calls.phased),FILTER=list(record.filter),QUAL=plain(record.qual),rendered_record=str(record))


def shape(context,ploidy):
    gt=context['FORMAT'].get('GT') or [];informative=any(v is not None for v in gt)
    if informative: equal(len(gt),ploidy,'informative GT ploidy')
    state='absent' if 'GT' not in context['FORMAT'] else 'fully_missing' if not informative else 'partial' if None in gt else 'fully_called'
    return dict(state=state,expected_model_ploidy=ploidy,encoded_GT_slots=len(gt),observed_GT_ploidy=len(gt) if informative else None,interpretation='Missing GT encoding is retained; encoded slot count does not establish ploidy.')


def native_call(context,idx):
    f=context['FORMAT'];gt=f.get('GT') or [];ad=f.get('AD');separator='|' if context['GT_phased'] else '/'
    call=dict(source_allele=dict(chrom=context['chrom'],pos=context['pos'],ref=context['ref'],alt=context['alts'][idx-1]),source_all_alts=context['alts'],source_alt_index=idx,GT=separator.join('.' if v is None else str(v) for v in gt),GT_values=gt,GT_phased=context['GT_phased'],alt_copies=gt.count(idx),called_ploidy=len(gt),fully_called=bool(gt) and None not in gt,alt_depth=ad[idx] if isinstance(ad,list) and len(ad)>idx else None,FILTER=context['FILTER'],QUAL=context['QUAL'],source_record_id=context['source_record_id'],caller=context['caller'],source_namespace=context['source_namespace'],native_context=context)
    call.update({k:f.get(k) for k in ['DP','GQ','AD','PL','PGT','PID','PS']})
    return call


def axis(context,name,threshold):
    value=context['FORMAT'].get(name);present=name in context['FORMAT']
    state='missing' if value is None else 'invalid_review' if type(value) not in (int,float) or not math.isfinite(value) or value<0 else 'explicit_zero' if value==0 else 'below_historical_threshold' if value<threshold else 'at_or_above_historical_threshold'
    return dict(state=state,value=value,historical_threshold=threshold,FORMAT_presence='absent' if not present else 'present_missing' if value is None else 'present_value',interpretation='Descriptive historical bin only; not an eligibility or validated quality threshold.')


def quality(context,call,ploidy):
    filters=context['FILTER'];gt=call['GT_values'];called=call['alt_copies']>0
    return dict(filter_state='explicit_PASS' if filters==['PASS'] else 'unfiltered' if not filters else 'named_filters',filter_values=filters,DP=axis(context,'DP',10),GQ=axis(context,'GQ',20),allele_participation='observed_in_GT' if called else 'not_observed_in_partial_or_missing_GT' if not call['fully_called'] else 'ALT_absent_from_GT',nonconcrete_in_GT=any(v is not None and v>0 and not concrete(context['alts'][v-1]) for v in gt),genotype=shape(context,ploidy),called_ALT_present=called,fully_called_ALT_present=called and call['fully_called'],historical_bins_are_eligibility_gate=False)


def native_inventory(source,sample,regions,lookup,config,run,out,fasta):
    tagged=run/source['id']/'ordinal-tagged-original.bcf';digest=hashlib.sha256();count=0
    with pysam.VariantFile(source['pin']['path']) as raw,pysam.VariantFile(tagged) as copied:
        equal(list(raw.header.samples),[sample],'native sample header')
        equal(list(copied.header.samples),[sample],'copied sample header')
        for count,record in enumerate(raw,1):
            copy=next(copied,None);check(copy is not None,'native copy omitted record')
            equal(copy.info['NS_SOURCE_ORD'],count,'native copy ordinal');del copy.info['NS_SOURCE_ORD']
            equal(str(copy),str(record),'native copy payload differs')
            ctx=native_context(record,source,sample,count);shape(ctx,model_ploidy(record,source,regions))
            digest.update(str(record).encode())
        check(next(copied,None) is None,'extra copied native record')
    equal(count,source['expected_census']['counts']['records'],'raw record census')
    if 'rendered_record_payload_sha256' in source['expected_census']:
        equal(digest.hexdigest(),source['expected_census']['rendered_record_payload_sha256'],'raw payload census')
    chosen=defaultdict(list);counts=Counter();inputs=config['inputs']
    with norm(inputs['bcftools']['path'],inputs['reference']['path'],tagged,out/(source['id']+'-renormalization.stderr')) as (reader,command):
        for record in reader:
            counts['normalized_records']+=1;model_ploidy(record,source,regions)
            ordinal=record.info['NS_SOURCE_ORD'];check(type(ordinal) is int and 1<=ordinal<=count,'normalized native ordinal')
            if not all(map(concrete,record.alleles)):
                counts['nonconcrete_normalized_records']+=1;continue
            key=allele(record)
            if key in lookup: chosen[ordinal].append(record.copy())
    expected=[];seen=set();recovered=set()
    with pysam.VariantFile(source['pin']['path']) as raw:
        for ordinal,record in enumerate(raw,1):
            if ordinal not in chosen: continue
            recovered.add(ordinal);context=native_context(record,source,sample,ordinal);p=model_ploidy(record,source,regions)
            for normalized in chosen[ordinal]:
                idx=original_index(normalized,context);literal(fasta,normalized,context,idx)
                equal(model_ploidy(normalized,source,regions),p,'normalized model class')
                key=allele(normalized);identity=(ordinal,idx,key);check(identity not in seen,'duplicate native normalization lineage');seen.add(identity)
                call=native_call(context,idx)
                expected.append(dict(caller=source['caller'],source_scope=source['id'],normalized_allele=dict(zip(['chrom','pos','ref','alt'],key)),normalized_record=str(normalized),source_call=call,quality=quality(context,call,p),clinvar=lookup[key],reportability='pending_disease_variant_and_inheritance_review'))
    equal(sorted(recovered),sorted(chosen),'unrecovered source ordinal')
    counts['exact_matching_occurrences']=len(expected)
    return expected,dict(records=count,rendered_record_payload_sha256=digest.hexdigest(),copied_original_payload_preserved=True),dict(counts),command


def grouped(rows):
    buckets=defaultdict(list)
    for row in rows: buckets[canon(dict(caller=row['caller'],allele=row['normalized_allele']))].append(row)
    groups=[];counts={c:Counter() for c in CALLERS}
    excluded={'source_record_id','source_namespace','caller','sample','source_sha256','source_ordinal'}
    for key,rows in sorted(buckets.items()):
        description=json.loads(key);caller=description['caller']
        any_alt=any(row['quality']['called_ALT_present'] for row in rows)
        full=any(row['quality']['fully_called_ALT_present'] for row in rows)
        # Source identity is excluded; every remaining original field is compared.
        evidence={canon({k:v for k,v in row['source_call']['native_context'].items() if k not in excluded}) for row in rows}
        groups.append(description|dict(source_lineages=[dict(source_record_id=r['source_call']['source_record_id'],source_alt_index=r['source_call']['source_alt_index']) for r in rows],source_occurrences=len(rows),duplicate_normalized_identity=len(rows)>1,differing_native_evidence=len(evidence)>1,called_ALT_present=any_alt,fully_called_ALT_present=full,interpretation='No consensus genotype is assigned; every source occurrence is retained.'))
        count=counts[caller];count['unique_exact_matching_identities']+=1;count['exact_matching_occurrences']+=len(rows)
        count['unique_called_ALT_identities']+=int(any_alt);count['unique_fully_called_ALT_identities']+=int(full)
        count['unique_identities_without_called_ALT']+=int(not any_alt);count['duplicate_normalized_identities']+=int(len(rows)>1)
        if any_alt:
            for stratum in {c['stratum'] for row in rows for c in row['clinvar']}:count['called_ALT_unique_identity_stratum/'+stratum]+=1
    fields=['unique_exact_matching_identities','exact_matching_occurrences','unique_called_ALT_identities','unique_fully_called_ALT_identities','unique_identities_without_called_ALT','duplicate_normalized_identities']
    return groups,{c:{k:counts[c][k] for k in fields}|dict(counts[c]) for c in CALLERS}


def public_contract(config):
    equal(config['expected_genes'],84,'real panel size')
    snapshot=config['public_snapshot']
    equal(snapshot['clinvar_fileDate'],'2026-09-13','real literal header date')
    equal(snapshot['clinvar_snapshot_label'],'20260913','real snapshot label')
    equal(snapshot['clinvar_resource_filename'],'clinvar_20260913.vcf.gz','real resource filename')
    equal(Path(config['inputs']['clinvar']['path']).name,snapshot['clinvar_resource_filename'],'real resource basename')
    for k,v in FROZEN_PUBLIC.items():equal(config['inputs'][k]['sha256'],v,'real public pin')


def audit(config_path,config_sha,run,status_sha,qualification,qualification_sha,output):
    config_path=Path(config_path).resolve();run=Path(run).resolve();output=Path(output).resolve();qualification=Path(qualification).resolve()
    check(not output.exists(),'existing audit output preserved')
    check(not output.is_relative_to(run),'audit outputs must be outside producer output tree')
    files=Files();config=files.get(dict(path=str(config_path),sha256=config_sha));status=files.get(dict(path=str(run/'status.json'),sha256=status_sha))
    equal(status['status'],'complete','producer not complete');equal(config['ready'],True,'config not ready');equal(config['schema_version'],1,'config schema')
    check(type(config['fixture_only']) is bool,'fixture-only type')
    equal(config['quality_policy'],'descriptive_native_axes_no_eligibility_gate','quality policy')
    qualification_data=files.get(dict(path=str(qualification),sha256=qualification_sha))
    equal(qualification_data['status'],'complete','qualification incomplete')
    for k in ['errors','failures','exit_code']:equal(qualification_data[k],0,'qualification '+k)
    equal(qualification_data['skipped'],[],'qualification skips')
    check(type(qualification_data['tests_run']) is int and qualification_data['tests_run']>=30,'qualification controls')
    expected_outputs={}
    for pin in status['outputs']:
        p=Path(pin['path']).resolve();check(p.is_relative_to(run),'output pin outside run');check(str(p) not in expected_outputs,'duplicate output pin')
        expected_outputs[str(p)]=files.hash(p,pin)
    equal(sorted(expected_outputs),sorted(str(p.resolve()) for p in run.rglob('*') if p.is_file() and p!=run/'status.json'),'producer output inventory')
    summary=read_json(run/'summary.json');files.hash(run/'summary.json',{'sha256':status['summary_sha256']})
    equal(summary['status'],'complete','summary incomplete');equal(summary['fixture_only'],config['fixture_only'],'fixture declaration');equal(status['fixture_only'],config['fixture_only'],'status fixture declaration')
    equal(summary['sample'],config['sample'],'summary sample');equal(summary['public_snapshot'],config['public_snapshot'],'public snapshot')
    for key,pin in summary['inputs'].items():
        equal(str(Path(key).resolve()),str(Path(pin['path']).resolve()),'summary input key');files.hash(pin['path'],pin)
    proof_bind(dict(path=str(config_path),sha256=config_sha),list(summary['inputs'].values()),files)
    source_pins=qualification_data['source_sha256']
    for rel in ['scripts/review_native_secondary.py','scripts/vendor/native_coding_records.py','scripts/vendor/screen_secondary_findings.py']:
        matches=[p for p in summary['inputs'].values() if p['path'].endswith('/'+rel)]
        check(len(matches)==1,'executed source member ambiguous: '+rel)
        equal(matches[0]['sha256'],source_pins[rel],'executed source not qualified')
        files.hash(qualification.parent/'workflow'/rel,{'sha256':source_pins[rel]})
    for rel,pin in qualification_data['source_snapshot_files'].items():
        files.hash(qualification.parent/'workflow'/rel,pin)
        equal(pin['sha256'],source_pins[rel],'qualification manifest consistency')
    equal(summary['helper_sha256'],{Path(k).name:v for k,v in source_pins.items() if '/vendor/' in k},'helper manifest')
    required={'reference','reference_index','bcftools','clinvar','clinvar_index','panel','handoff','autosome_verification','fixed_xy_status','fixed_xy_review','fixed_xy_verification'}
    equal(sorted(config['inputs']),sorted(required),'configured input membership')
    for pin in config['inputs'].values():
        check(Path(pin['path']).is_absolute() and type(pin['bytes']) is int and pin['bytes']>=0,'input path/size')
        files.hash(pin['path'],pin);proof_bind(pin,list(summary['inputs'].values()),files)
    for name,index in [('reference','reference_index'),('clinvar','clinvar_index')]:
        suffix='.fai' if name=='reference' else '.tbi'
        equal(str(Path(config['inputs'][name]['path']+suffix).resolve()),str(Path(config['inputs'][index]['path']).resolve()),'index binding')
    if not config['fixture_only']:
        public_contract(config)
    output.mkdir(parents=True)
    try:
        with pysam.FastaFile(config['inputs']['reference']['path']) as fasta:
            sources,regions=sources_from_proofs(config,files,fasta)
            for source in sources: proof_bind(source['pin'],list(summary['inputs'].values()),files)
            equal(summary['source_scope'],{c:CHROMS+['X','Y'] for c in CALLERS},'summary complete scope');equal(summary['source_files'],len(sources),'summary source count')
            lookup,public_counts,genes,public_command=public_inventory(config,run,output,files,fasta)
            equal(summary['public_counts'],public_counts,'public census');equal(summary['panel_genes'],genes,'panel genes')
            saved=read_json(run/'source-censuses.json');equal(len(saved),len(sources),'source census count')
            expected=[];commands=[public_command]
            for source,census in zip(sources,saved):
                # Preserve original metadata fields while verifying current identity.
                equal(census['source']['id'],source['id'],'source order')
                for key in ['caller','kind','contigs','expected_census','id']:equal(census['source'][key],source[key],'source declaration '+key)
                files.hash(source['pin']['path'],census['source']['pin'])
                rows,counts,normalization,command=native_inventory(source,config['sample'],regions,lookup,config,run,output,fasta)
                equal(census['census'],counts,'source copy census');equal(census['counts'],normalization,'source normalization census')
                expected.extend(rows);commands.append(command)
        with gzip.open(run/'exact-matching-occurrences.jsonl.gz','rt') as stream:
            actual=[parse_json(line) for line in stream]
        equal(actual,expected,'complete exact occurrence payload/order/multiplicity')
        groups,counts=grouped(expected)
        equal(read_json(run/'exact-identity-groups.json'),groups,'identity grouping')
        equal(summary['counts_by_caller'],counts,'caller counts')
        axes={c:{axis:dict(Counter((row['quality'][axis]['state'] if axis in ['DP','GQ'] else row['quality'][axis]) for row in expected if row['caller']==c)) for axis in ['filter_state','DP','GQ']} for c in CALLERS}
        equal(summary['exact_occurrence_quality_axes'],axes,'quality-axis counts')
        files.stable();files.hash(__file__)
        result=dict(status='complete',verified_at=now(),fixture_only=config['fixture_only'],source_files=len(sources),public_counts=public_counts,counts_by_caller=counts,exact_occurrences=len(expected),independent_source_reconstruction=True,normalization_engine='same pinned bcftools; independent invocation, not an independent normalization algorithm',normalization_commands=commands,verified_files=list(files.pins.values()),limitations=['Same specimen and same normalization engine; no independent biological validation.','Exact dated occurrence/GT membership only; no clinical negative, phase, reportability, callability or quality eligibility inference.'])
        write_json(output/'verification.json',result)
        write_json(output/'status.json',dict(status='complete',completed_at=now(),verification_sha256=hashlib.sha256((output/'verification.json').read_bytes()).hexdigest()))
        return result
    except BaseException as error:
        write_json(output/'status.json',dict(status='failed',failed_at=now(),error_type=type(error).__name__,error=str(error)))
        raise


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['config','run','qualification','output']:p.add_argument('--'+name,type=Path,required=True)
    for name in ['config','status','qualification']:p.add_argument('--expected-'+name+'-sha256',required=True)
    a=p.parse_args();r=audit(a.config,a.expected_config_sha256,a.run,a.expected_status_sha256,a.qualification,a.expected_qualification_sha256,a.output)
    print(json.dumps({k:r[k] for k in ['status','source_files','exact_occurrences','counts_by_caller']}))


if __name__=='__main__':main()
