#!/usr/bin/env python3
"""Independent ADAMTS1 overlay audit; stdlib only, no producer imports/read access."""
import argparse
from collections import Counter, defaultdict
import copy
import datetime
import hashlib
import json
from pathlib import Path
import re


def require(value,message):
    if not value:raise ValueError(message)


def canonical(value):return json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False)


def same(a,b,message):require(canonical(a)==canonical(b),message)


def parse(text):
    def pairs(items):
        d={}
        for k,v in items:require(k not in d,'duplicate JSON key');d[k]=v
        return d
    def bad(v):raise ValueError('nonfinite JSON '+v)
    return json.loads(text,object_pairs_hook=pairs,parse_constant=bad)


def sha(text):return hashlib.sha256(text.encode()).hexdigest()


def stamp():return datetime.datetime.now(datetime.timezone.utc).isoformat()


def write(path,value):
    with path.open('x') as f:f.write(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n')


class Pins:
    def __init__(self):self.items={}
    def add(self,path,expected=None):
        path=Path(path).resolve();key=str(path)
        require(not any(key.endswith(s) for s in ['.bam','.cram','.sam','.bam.bai','.cram.crai','.fastq','.fastq.gz','.fq','.fq.gz']),'read/alignment access forbidden')
        if key not in self.items:
            before=path.stat();h=hashlib.sha256()
            with path.open('rb') as stream:
                for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
            after=path.stat();same([before.st_size,before.st_mtime_ns],[after.st_size,after.st_mtime_ns],'changed input while hashing')
            self.items[key]=dict(path=key,bytes=after.st_size,mtime_ns=after.st_mtime_ns,sha256=h.hexdigest())
        result=self.items[key]
        if expected:
            if 'path' in expected:same(key,str(Path(expected['path']).resolve()),'pin path mismatch')
            for field in ['sha256','bytes']:
                if field in expected:same(result[field],expected[field],'pin '+field+' mismatch')
        return result
    def load(self,path,expected=None):self.add(path,expected);return parse(Path(path).read_text())
    def stable(self):
        for p in self.items.values():
            st=Path(p['path']).stat();same([st.st_size,st.st_mtime_ns],[p['bytes'],p['mtime_ns']],'input changed during audit')


def compose(reference,start0,edits):
    """Validate original spans once, then substitute from right to left."""
    require(type(start0) is int and start0>=0 and re.fullmatch('[ACGT]+',reference) is not None,'bad reference window')
    spans=[]
    for edit in edits:
        require(type(edit['pos']) is int and all(re.fullmatch('[ACGT]+',edit[k] or '') for k in ['ref','alt']),'bad edit')
        left=edit['pos']-1-start0;right=left+len(edit['ref'])
        require(0<=left<right<=len(reference),'edit outside window');same(reference[left:right],edit['ref'],'literal REF mismatch')
        spans.append((left,right,edit['alt']))
    spans.sort()
    require(all(a[1]<=b[0] for a,b in zip(spans,spans[1:])),'overlapping original spans')
    result=reference
    for left,right,alt in reversed(spans):result=result[:left]+alt+result[right:]
    return result


def original_edit(contexts,role):
    context=contexts[role['source_record_id']];idx=role['source_alt_index']
    same(context['source_record_id'],role['source_record_id'],'context identity')
    require(type(idx) is int and 1<=idx<=len(context['alts']),'original ALT index')
    return dict(chrom=context['chrom'],pos=context['pos'],ref=context['ref'],alt=context['alts'][idx-1],caller=context['caller'],source_record_id=context['source_record_id'],source_alt_index=idx,called=idx in (context['FORMAT'].get('GT') or []))


def overlay_proof(reference,start,spec,contexts):
    roles=spec['roles'];raw=original_edit(contexts,roles['raw_deletion']);snv=original_edit(contexts,roles['disjoint_substitution']);overlap=original_edit(contexts,roles['overlapping_substitution'])
    comparisons=[original_edit(contexts,r) for r in roles['comparison_deletions']];star=original_edit(contexts,roles['spanning_call']);star_context=contexts[star['source_record_id']]
    same([star['source_alt_index'],star_context['alts'],star_context['FORMAT']['GT']],[2,['*',star['alt']],[1,2]],'spanning star source representation')
    require(len(star_context['FORMAT']['PL'])==6 and re.fullmatch('[ACGT]+',star['alt']) is not None,'spanning star PL/concrete ALT')
    require(len(comparisons)==2 and len({r['source_record_id'] for r in comparisons})==2 and len({r['caller'] for r in comparisons})==1,'comparison lineage')
    require(len({raw['source_record_id'],snv['source_record_id'],overlap['source_record_id']})==3,'duplicate derived source identity')
    require(raw['caller']==snv['caller']==overlap['caller'],'derived callers')
    require(all(e['chrom']==spec['chrom'] for e in [raw,snv,overlap,star]+comparisons),'source chromosomes')
    normal=spec['normalized_deletion'];same(normal['chrom'],spec['chrom'],'normalized chromosome')
    key=':'.join(str(normal[k]) for k in ['chrom','pos','ref','alt']);same(key,spec['normalized_key'],'normalization key')
    first=compose(reference,start,[normal]);same(first,compose(reference,start,[raw]),'raw-normalized sequence equality')
    joint=compose(reference,start,[normal,snv]);same(joint,compose(reference,start,comparisons),'joint-double-deletion sequence equality')
    rejected=[]
    for label,edits in [('raw_deletion_plus_disjoint_substitution',[raw,snv]),('normalized_deletion_plus_overlapping_substitution',[normal,overlap]),('normalized_deletion_plus_both_substitutions',[normal,overlap,snv])]:
        try:compose(reference,start,edits)
        except ValueError as e:
            same(str(e),'overlapping original spans','overlap rejected for wrong reason')
            rejected.append(dict(label=label,edits=edits,reason='overlapping REF edits'))
        else:raise ValueError('expected overlap rejection missing')
    common=dict(caller=raw['caller'],normalization_key=key,original_deletion_edit=raw,normalized_deletion=normal,raw_normalized_sequence_equal=True,interpretation='Constructed normalized-representation lineage; not native genotype or biological phase.')
    derived=[dict(sequence=first,lineage=common|dict(label='normalized_equivalent_deletion',substitution_edits=[])),dict(sequence=joint,lineage=common|dict(label='normalized_deletion_plus_disjoint_substitution',substitution_edits=[snv],exact_sequence_comparison_lineage=comparisons))]
    return dict(lineages=derived,rejected_compositions=rejected,preserved_spanning_call=star_context,spanning_concrete_ALT_edit=star,raw_normalized_equal=True,joint_sequence_equal=True)


def sequence_classes(old,derived,prefix,reference,start,expected_count):
    require(len(old)==expected_count,'wrong original class count')
    classes={};ids=set()
    for item in old:
        require(item['id'] not in ids and item['sequence'] not in classes,'duplicate native class or ID');ids.add(item['id'])
        for lineage in item['source_combinations']:
            same(compose(reference,start,lineage['edits']),item['sequence'],'original lineage literal sequence')
        classes[item['sequence']]=dict(sequence=item['sequence'],original_hypothesis_ids=[item['id']],original_lineages=item['source_combinations'],derived_lineages=[])
    for item in derived:
        require(item['sequence'] in classes,'new unsupported sequence class')
        classes[item['sequence']]['derived_lineages'].append(item['lineage'])
    for i,seq in enumerate(sorted(classes)):
        classes[seq].update(id=prefix+'_U%04d'%i,sequence_sha256=sha(seq),sequence_length=len(seq))
    return classes


def count_arm(original,classes):
    observations=original['observations'];groups=defaultdict(list);observation_counts=Counter();seen=set();accepted_count=0
    for row in observations:
        key=row['observation_id'];require(type(key) is int and key not in seen,'duplicate/bad observation id');seen.add(key)
        require(row['status'] in ['accepted','rejected'],'unknown observation status')
        identity=(row['RG'],row['QNAME']);valid=all(type(x) is str and bool(x) for x in identity)
        if valid:groups[identity].append(row)
        if row['status']=='accepted':
            require(valid and type(row['sequence']) is str,'accepted observation lacks identity/sequence')
            accepted_count+=1;observation_counts[row['sequence']]+=1
    saved=original['fragments'];same([(v['RG'],v['QNAME']) for v in saved],sorted(groups),'original fragment identity/order')
    counts=dict(raw_observations=len(observations),raw_fragments=len(groups),accepted_observations=accepted_count,conflicting_fragments=0,nonconflicting_spanning_fragments=0,matched_fragments=0,unmatched_fragments=0,fragments_in_derived_lineage_classes=0)
    assignments=[];support=Counter()
    for saved_fragment,(identity,rows) in zip(saved,sorted(groups.items())):
        accepted=[r for r in rows if r['status']=='accepted'];sequences=sorted({r['sequence'] for r in accepted})
        ids=[r['observation_id'] for r in rows];accepted_ids=[r['observation_id'] for r in accepted]
        same(saved_fragment['observation_ids'],ids,'original fragment observations');same(saved_fragment['accepted_observation_ids'],accepted_ids,'original accepted membership');same(saved_fragment['sequences'],sequences,'original fragment sequences')
        item=dict(RG=identity[0],QNAME=identity[1],observation_ids=ids,accepted_observation_ids=accepted_ids,class_id=None,status='no_eligible_span')
        if not sequences:old_status='no_eligible_spanning_observation'
        elif len(sequences)>1:old_status='mate_or_repeated_observation_conflict';counts['conflicting_fragments']+=1;item['status']='conflict'
        else:
            counts['nonconflicting_spanning_fragments']+=1;sequence=sequences[0]
            if sequence in classes:
                old_status='matched';counts['matched_fragments']+=1;support[sequence]+=1;item.update(status='matched',class_id=classes[sequence]['id'])
                if classes[sequence]['derived_lineages']:counts['fragments_in_derived_lineage_classes']+=1
            else:old_status='unmatched';counts['unmatched_fragments']+=1;item['status']='unmatched'
        same(saved_fragment['status'],old_status,'original fragment status');assignments.append(item)
    for k in ['raw_observations','accepted_observations','conflicting_fragments','nonconflicting_spanning_fragments','matched_fragments','unmatched_fragments']:same(counts[k],original['counts'][k],'original counter '+k)
    same(counts['raw_fragments'],original['counts']['fragments_with_identity'],'original group count')
    perclass={v['id']:dict(observations=observation_counts[sequence],fragments=support[sequence]) for sequence,v in classes.items()}
    same(sum(v['fragments'] for v in perclass.values()),counts['matched_fragments'],'lineage multiplication of fragments')
    return dict(counts=counts,class_support=perclass,fragment_assignments=assignments,original_observation_evidence=original)


def audit_content(native,context,selection,overlay,summary,expected_classes=9):
    spec=selection['overlay_spec'];config=selection['native_config_body'];same(native['config'],config,'native config payload')
    original_contexts={r['context']['source_record_id']:r['context'] for r in context['targets'][spec['target_id']]}
    same(len(original_contexts),len(context['targets'][spec['target_id']]),'duplicate native context')
    same(original_contexts,selection['contexts'],'prepared/native context equality');same(overlay['native_contexts'],original_contexts,'overlay full native context preservation')
    preflights={r['flank_bp']:r for r in selection['windows']};same(len(preflights),len(selection['windows']),'duplicate prepared window')
    same([w['flank_bp'] for w in native['windows']],config['window_flanks_bp'],'native flank scope')
    same(len(native['windows']),len(overlay['windows']),'overlay window count');same(len(native['windows']),len(summary['windows']),'summary window count')
    arm_total=0;tables=[]
    for window,observed,bounded in zip(native['windows'],overlay['windows'],summary['windows']):
        pre=preflights[window['flank_bp']]
        for k in ['chrom','start0','end0','reference_sequence']:same(window[k],pre[k],'reference window '+k)
        old=[native['hypotheses'][i] for i in window['hypothesis_ids']]
        same({h['sequence']:h['source_combinations'] for h in old},{h['sequence']:h['source_combinations'] for h in pre['native_family']['haplotypes']},'prepared original lineages')
        same(window['enumeration_metadata']['rejected_combinations'],pre['native_family']['rejected_combinations'],'native rejection preservation')
        same(window['omitted_overlapping_native_alleles'],pre['omitted_native_alleles'],'native omission preservation')
        proof=overlay_proof(window['reference_sequence'],window['start0'],spec,original_contexts)
        classes=sequence_classes(old,proof['lineages'],'ADAMTS1_F'+str(window['flank_bp']),window['reference_sequence'],window['start0'],expected_classes)
        same(observed['overlay'],proof,'derived source lineages and equality proof')
        same(observed['classes'],list(classes.values()),'complete original/derived class payload')
        for k in ['chrom','start0','end0','flank_bp']:same(observed[k],window[k],'overlay geometry '+k)
        for k in ['original_class_count','union_class_count']:same(observed[k],expected_classes,'class count');same(bounded[k],expected_classes,'bounded class count')
        expected_gates=[(mq,bq) for mq in config['mapping_quality_thresholds'] for bq in config['base_quality_thresholds']]
        same([(a['minimum_mapping_quality'],a['minimum_base_quality']) for a in window['arms']],expected_gates,'native arm scope')
        arms=[]
        for arm in window['arms']:
            mq,bq=arm['minimum_mapping_quality'],arm['minimum_base_quality'];source=arm['evidence']
            same(source['window'],{'start0':window['start0'],'end0':window['end0']},'native observation window')
            same(source['thresholds'],dict(min_MQ=mq,min_raw_BQ=bq,excluded_flags=0xF04,reject_MQ_255=True,BAQ_applied=False),'native admission policy')
            arms.append(dict(minimum_mapping_quality=mq,minimum_base_quality=bq,**count_arm(source,classes)))
        same(observed['arms'],arms,'complete observation/fragment preservation and reassignment')
        want=dict(flank_bp=window['flank_bp'],original_class_count=expected_classes,union_class_count=expected_classes,arms=[dict(MQ=a['minimum_mapping_quality'],raw_BQ=a['minimum_base_quality'],counts=a['counts'],class_support=a['class_support']) for a in arms])
        same(bounded,want,'bounded summary counters')
        tables.append(want);arm_total+=len(arms)
    return dict(windows=len(tables),arms=arm_total,tables=tables,original_classes_per_window=expected_classes,derived_lineages_per_window=2,new_sequence_classes=0)


def run(args):
    root=args.root.resolve();native=args.native_run.resolve();overlay_run=args.overlay_run.resolve();out=args.output.resolve()
    require(not out.exists(),'existing independent output preserved')
    require(not out.is_relative_to(native) and not out.is_relative_to(overlay_run),'independent output separate')
    pins=Pins()
    def load(path,digest=None):return pins.load(path,{'sha256':digest} if digest else None)
    status=load(overlay_run/'status.json',args.expected_overlay_status_sha256);same(status['status'],'complete','overlay incomplete')
    summary=load(overlay_run/'summary.json',status['summary_sha256']);same(summary['status'],'complete','overlay summary incomplete')
    overlay=pins.load(overlay_run/'evidence.json',summary['evidence']);same(overlay['status'],'complete','overlay evidence incomplete')
    same(summary['verified_files'],overlay['verified_files'],'overlay manifest copies')
    for path,pin in overlay['verified_files'].items():same(str(Path(path).resolve()),str(Path(pin['path']).resolve()),'overlay input key');pins.add(path,pin)
    def declared(path,digest):
        path=Path(path).resolve();require(str(path) in overlay['verified_files'],'input absent from overlay manifest');pins.add(path,overlay['verified_files'][str(path)]);return load(path,digest)
    selection=declared(args.selection,args.expected_selection_sha256);same(selection['status'],'complete','selection incomplete');same(selection['mode'],'reference_only_preparation','selection scope')
    controls=declared(args.controls,args.expected_controls_sha256)
    same(controls['status'],'complete','overlay controls incomplete');same(controls['tests_run'],21,'overlay control count')
    for k in ['errors','failures','skipped']:same(controls[k],0,'overlay controls '+k)
    same(controls['source_sha256']['normalized-overlay.py'],args.expected_overlay_source_sha256,'qualified overlay source')
    candidates=[p for p in overlay['verified_files'].values() if p['path'].endswith('/normalized-overlay.py')]
    require(bool(candidates),'missing executed overlay source')
    for p in candidates:same(p['sha256'],args.expected_overlay_source_sha256,'executed overlay source')
    verified=declared(args.native_verification,args.expected_native_verification_sha256);same(verified['status'],'complete','original native audit incomplete')
    index={str(Path(v['path']).resolve()):v for v in verified['verified_files'].values()}
    require(len(index)==len(verified['verified_files']),'duplicate native audit path')
    original={}
    for name in ['status.json','summary.json','evidence.json','regional-native-context.json']:
        p=native/name;require(str(p) in index,'native artifact not independently verified')
        original[name]=pins.load(p,index[str(p)]);require(str(p) in overlay['verified_files'],'native artifact absent from overlay inputs');pins.add(p,overlay['verified_files'][str(p)])
    same(original['status.json']['status'],'complete','native status incomplete');same(original['summary.json']['status'],'complete','native summary incomplete')
    same(original['status.json']['summary_sha256'],pins.items[str(native/'summary.json')]['sha256'],'native summary/status link')
    for name in ['evidence.json','regional-native-context.json']:pins.add(native/name,original['summary.json']['outputs'][name])
    same(original['summary.json']['inputs']['config']['sha256'],selection['native_config_pin']['sha256'],'native config pin')
    same(selection['native_config_body']['window_flanks_bp'],[3,10,25],'actual flank scope')
    same(selection['native_config_body']['mapping_quality_thresholds'],[0,20,60],'actual MQ scope')
    same(selection['native_config_body']['base_quality_thresholds'],[20,30],'actual BQ scope')
    require(selection['overlay_spec']['target_id']=='ADAMTS1','actual target')
    out.mkdir(parents=True)
    try:
        result=audit_content(original['evidence.json'],original['regional-native-context.json'],selection,overlay,summary)
        same(result['arms'],18,'complete actual arm count');pins.add(__file__);pins.stable()
        result.update(status='complete',verified_at=stamp(),verified_files=pins.items,scope='Original already-audited observations only; no read extraction, alignment, genotyper or annotation.',limitations=['Admission/quality geometry relies on the separately pinned native independent audit; this checker independently regroups its retained observations.','Exact sequence equivalence and counts do not establish genotype, phase, locus origin, biological function or clinical interpretation.'])
        write(out/'verification.json',result);write(out/'status.json',dict(status='complete',completed_at=stamp(),verification_sha256=hashlib.sha256((out/'verification.json').read_bytes()).hexdigest()))
        return result
    except BaseException as e:
        write(out/'status.json',dict(status='failed',failed_at=stamp(),error_type=type(e).__name__,error=str(e)));raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['root','native-run','overlay-run','selection','controls','native-verification','output']:p.add_argument('--'+name,type=Path,required=True)
    for name in ['overlay-status','selection','controls','native-verification','overlay-source']:p.add_argument('--expected-'+name+'-sha256',required=True)
    result=run(p.parse_args());print(json.dumps({k:result[k] for k in ['status','windows','arms','original_classes_per_window','derived_lineages_per_window','new_sequence_classes']}))
