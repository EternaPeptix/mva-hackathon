#!/usr/bin/env python3
"""Verify the explicit 64-ID food-overlay admission; never confer new reading.

check(root=...) supports compact repository checks. Full food-evidence replay
additionally requires work_root and corpus. No source annotations are rewritten.
"""
import argparse
import hashlib
import importlib.util
import json
import sys
from pathlib import Path
sys.dont_write_bytecode = True

SOURCE_SHA='e5b8b017bb68af38818f1cc8c1bee47adb16afe0367981104a6b17886d2954d9'
PARENT_SHA='a9a12a95285b7456be572811a0a03891748535abe263706184ef517069938fb2'
FOOD_SOURCE_SHA='f4395a70e1e499c8e0971bf8a884c74ed09ce2a6a1efccb1138e2e02974a0fb1'
FOOD_REPO_SHA='9dcffb1be97f6559f42e4b4b3268df60a6cfe4d403c1fd9729ecf2c46ef68183'
COMMON_SHA='1eaa78d49dfd9aa4a17ea1bee1e339668cf07f686c0b6311582bce11cab3675e'
COMMON_VALIDATION_SHA='c273e833b680b34dec9aed0d4e837e9b3f65ba11435417b4e344c7e2a345f999'
BASE_VERIFIER_SHA='c8424d6e574c30a91af6a4bae19f9d9a23eae1d1372e0113c5caeb5bd0677100'
PRIOR={
 'outer-revalidation-review.json':'0bdc47535b70046b0d0a115da26339ec845aacb4d38827a9cec65e88e0dd81bb',
 'outer-revalidation-next128-review.json':'26102e1a6e5bd6cd5204e4d38d92e53b5757709cb4cfda8a839fe7c9769fda7c',
 'outer-revalidation-batch02-review.json':'222ef9edb446dae6c87782eb5a672a904ba45342b5866e49fe3a854f8444ef03',
 'outer-revalidation-batch03-review.json':'28ad16f471d3951ec911498b4285ab023306239c547d4d294d4eec3a237daedb',
 'outer-revalidation-batch04-review.json':'20a7d66a849ea614cf5367b2b02d0d0e0f20f87a61fb46fe19ef103bfef1777c',
 'outer-revalidation-batch05-review.json':'46b5d57e18f6b4207fbbe0e5139b93e3cb0824d2a8e643fcba981d1423ac9bfd',
 'outer-revalidation-batch06-review.json':'843b26bf60409ddbe21ec733c8d6d309b98829d5e0d81253b04b4811be7f3276',
 'outer-revalidation-batch07-review.json':'0132f1a63792eaead21fbb9d9f430d813897abfac1454fd5fffc01eadf1aa840',
 'outer-revalidation-batch08-review.json':'c186a07facdce5cad33e256dfcaa46b65d635885d5e04558c53f9436cbb2f3fc',
}
SAMPLE_FILES={
 'prior-food-admission-samples.json':'f0ec6f432c3cc11527f2bc592bc2eea3df8f413cbcb619a53ab8bd520c319a0e',
 'common-admission-samples.json':'2b9a70038abbad4759dc06355d3cf7588139025c30a8a62830eb425741d9faa2',
}
SAMPLE_LOGS=[
 ('prior-food-sample-a.jsonl','4a628b95363b104a3d0f3386bfeaa63f29400f50e69bccf7019f60a4b6ac7f2e','c296fc',5388),
 ('prior-food-sample-b.jsonl','44208788c61ab7fac6bf7f91a3755467661ea2f295e551c74d79a41d37559709','12998d',4880),
 ('common-sample-a.jsonl','532b03385583074e3b93e65935c29944ccedd08b67e606bfb62abe3c12bc5c7a','3a9433',5337),
 ('common-sample-b.jsonl','a9603f5a1c9c833e8a93d98ab0083a646bb6655019af8a65b63870a44e1e18f9','0f5aa9',3985),
]
FALLBACK=Path('/Users/jeweled/mva-hackathon/artifacts/source-revalidation-20260922/repo/workflow/track2-drug-census')
ADMISSION_NAME='outer-food-admission.json'

def require(condition,message):
    if not condition: raise ValueError(message)

def read_json(path): return json.loads(Path(path).read_text())
def digest(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def object_digest(value): return hashlib.sha256(json.dumps(value,ensure_ascii=False,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def key(row): return row['source']+':'+str(row.get('record_id',row.get('id')))

def default_root():
    for p in [Path(__file__).resolve().parent,*Path(__file__).resolve().parents]:
        if (p/'outer-revalidation-selection.json').is_file(): return p
    return FALLBACK

def module(path,name):
    spec=importlib.util.spec_from_file_location(name,path);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod

def prepare(root,food_review=None,common_review=None):
    root=Path(root)
    require(digest(root/'outer-revalidation-selection.json')==PARENT_SHA,'Original3598 selection bytes changed')
    require(digest(root/'verify_revalidation_batch08.py')==BASE_VERIFIER_SHA,'Baseline verifier changed')
    baseline=module(root/'verify_revalidation_batch08.py','food_admission_baseline').check(root=root)
    require(baseline['result']=='PASS' and baseline['cumulative_accepted_revalidation_records']==1210 and baseline['pending_source_revalidation']==2388,'Baseline gates/counts differ')
    parent=read_json(root/'outer-revalidation-selection.json'); original=parent['records']; original_by_id={key(r):r for r in original}
    require(len(original)==len(original_by_id)==3598,'Original selection ID count differs')
    accepted=[]; prior_manifest=[]
    for name,sha in PRIOR.items():
        require(digest(root/name)==sha,'Prior accepted bytes changed: '+name)
        rows=read_json(root/name)['records'];accepted.extend(rows)
        prior_manifest.append({'file':name,'sha256':sha,'records':len(rows),'record_objects_sha256':object_digest(rows)})
        require(all(key(r) in original_by_id and all(r.get(k)==v for k,v in original_by_id[key(r)].items()) for r in rows),'Prior accepted selection object differs')
    accepted_ids={key(r) for r in accepted};require(len(accepted)==len(accepted_ids)==1210,'Prior accepted IDs duplicate or differ')
    pending=[r for r in original if key(r) not in accepted_ids];require(len(pending)==2388,'Pending queue differs')
    food_path=Path(food_review or root/'food-substances/abstract-review.json')
    common_path=Path(common_review or root/'food-substances/common-abstract-review.json')
    require(digest(food_path)==FOOD_REPO_SHA,'Compact literal-food review changed')
    require(digest(common_path)==COMMON_SHA,'Common-food review changed')
    food=read_json(food_path); common=read_json(common_path)
    require(food['integration']['source_review_sha256']==FOOD_SOURCE_SHA,'Original food review provenance differs')
    require(len(food['records'])==67 and len(common['records'])==154,'Food layer sizes differ')
    layers={'literal_food':food,'common_food':common}; indexes={name:{r['source_id']:(i,r) for i,r in enumerate(d['records'])} for name,d in layers.items()}
    require(all(len(indexes[name])==len(d['records']) for name,d in layers.items()),'Duplicate overlay IDs')
    a=set(indexes['literal_food']);b={r['source_id'] for r in common['records'] if r['abstract_present']}
    old_pending=[r for r in pending if key(r) in a];common_pending=[r for r in pending if key(r) in b]
    common_only=[r for r in common_pending if key(r) not in a];expected=[r for r in pending if key(r) in a|b]
    shared={key(r) for r in old_pending}&{key(r) for r in common_pending}
    require((len(old_pending),len(common_pending),len(common_only),len(shared),len(expected))==(26,40,38,2,64),'26+40-2 admission reconciliation differs')
    require(shared=={'MED:21753893','MED:37096022'},'Food intersection changed')
    require('MED:18082453' not in {key(r) for r in expected},'Missing abstract selected')
    return dict(root=root,parent=parent,original=original,accepted=accepted,accepted_ids=accepted_ids,pending=pending,prior_manifest=prior_manifest,layers=layers,indexes=indexes,old_pending=old_pending,common_pending=common_pending,common_only=common_only,shared=shared,expected=expected)

def evidence_reference(layer,position,row):
    ev=row.get('read_evidence')
    require(isinstance(ev,dict) and row.get('abstract_characters',0)>0,'Missing record reading evidence: '+row['source_id'])
    require(row.get('reading_credit',1)==1 and row.get('abstract_present',True) is True,'No reading credit for source: '+row['source_id'])
    require(ev.get('end_marker_observed') is True and ev.get('untruncated_output_observed') is True and bool(ev.get('tool_chunk_id')),'Incomplete reading declaration: '+row['source_id'])
    require(ev['abstract_sha256']==row['abstract_sha256'] and ev['abstract_characters']==row['abstract_characters'],'Reading evidence identity differs')
    require(row['automatic_outer_corpus_credit']==row['new_appraisal_credit']==0,'Overlay unexpectedly grants automatic credit')
    return {'layer':layer,'review_file':'food-substances/abstract-review.json' if layer=='literal_food' else 'food-substances/common-abstract-review.json','review_file_sha256':FOOD_REPO_SHA if layer=='literal_food' else COMMON_SHA,'record_index':position,'source_id':row['source_id'],'abstract_sha256':row['abstract_sha256'],'record_object_sha256':object_digest(row),'annotation_excluding_source_metadata_sha256':object_digest({k:v for k,v in row.items() if k!='source_metadata'}),'source_metadata_sha256':object_digest(row['source_metadata']),'read_evidence_sha256':object_digest(ev)}

def sample_gate(state):
    groups=[];sample_objects=[]
    for name,population,indices,comparison in [('prior-food-admission-samples.json',state['old_pending'],[0,5,10,15,20,25],'cb0737'),('common-admission-samples.json',state['common_only'],[0,7,15,22,30,37],'2840bd')]:
        chosen=[population[i] for i in indices];sample_objects.append(chosen)
        groups.append({'selection_file':name,'selection_sha256':SAMPLE_FILES[name],'population_count':len(population),'predetermined_zero_based_indices':indices,'source_ids':[key(r) for r in chosen],'selection_objects_sha256':object_digest(chosen),'annotation_comparison_tool_chunk_id':comparison,'annotation_comparison_complete_declared':True})
    logs=[]
    for i,(name,sha,chunk,chars) in enumerate(SAMPLE_LOGS):
        ids=groups[i//2]['source_ids'][(i%2)*3:(i%2+1)*3]
        logs.append({'log_file':name,'log_sha256':sha,'source_ids':ids,'emitted_records':3,'abstract_characters':chars,'tool_chunk_id':chunk,'end_marker_observed_declared':True,'untruncated_literal_read_declared':True})
    return {'status':'COMPLETE_DECLARED_BY_ROOT','reviewer':'/root','sample_records':12,'sample_groups':groups,'literal_logs':logs,'substantive_annotation_corrections':0,'additional_reading_credit':0,'declaration_provenance':'Parent/root completion message in this task on 2026-09-23: all twelve predetermined literal strings read with observed end markers and all twelve annotations compared; no substantive correction needed.','scope_limit':'Twelve bounded spot-checks, not independent re-review of all64 or all212 food sources. Viewing and annotation comparison are root declarations; file checks cannot prove them.'}

def construct_admission(state):
    records=[]
    for selected in state['expected']:
        sid=key(selected);refs=[]
        for layer,index in state['indexes'].items():
            if sid in index:
                i,row=index[sid]
                require(row['abstract_sha256']==selected['abstract_sha256'] and row['title']==selected['title'],'Overlay/source selection mismatch: '+sid)
                refs.append(evidence_reference(layer,i,row))
        records.append({**selected,'source_id':sid,'admission_credit':1,'new_literal_reading_credit':0,'new_appraisal_credit':0,'annotation_policy':'Reference the unchanged scientific annotation; do not duplicate or synthesize it.','evidence_references':refs})
    after=[r for r in state['pending'] if key(r) not in {key(x) for x in records}]
    input_refs={'literal_food_compact_review':{'file':'food-substances/abstract-review.json','sha256':FOOD_REPO_SHA,'original_external_review_sha256':FOOD_SOURCE_SHA},'common_food_review':{'file':'food-substances/common-abstract-review.json','sha256':COMMON_SHA},'common_structural_validation_sha256':COMMON_VALIDATION_SHA}
    return {'schema_version':'explicit-food-outer-admission/1','date':'2026-09-23','reference_commit':'907eeb7219b2bf19dc5ae5affbb81508488e5349','status':'EXPLICIT_RECONCILIATION_ADMITTED','source_sha256':SOURCE_SHA,'parent_selection_sha256':PARENT_SHA,'original_selection_record_objects_sha256':object_digest(state['original']),'prior_accepted_inputs':state['prior_manifest'],'prior_accepted_record_objects_sha256':object_digest(state['accepted']),'evidence_inputs':input_refs,'selection_rule':'Original3598 order; exclude all1210 prior accepted IDs; admit the union of already-read literal-food and available common-food IDs that remain in the2388 queue. No new outcome/title filter.','old_food_pending_count':26,'common_food_pending_count':40,'shared_pending_ids':sorted(state['shared']),'selected_count':64,'reviewed_count':64,'admitted_count':64,'prior_accepted_revalidation_records':1210,'cumulative_accepted_revalidation_records':1274,'previous_outer_ledger_count':3371,'combined_outer_ledger_count':3435,'pending_before':2388,'pending_after':2324,'pending_count':2324,'pending_before_record_objects_sha256':object_digest(state['pending']),'pending_after_record_objects_sha256':object_digest(after),'missing_abstracts_excluded':['MED:18082453'],'new_literal_reading_credit':0,'new_appraisal_credit':0,'food_overlay_automatic_credit_remains':0,'root_sample_gate':sample_gate(state),'verification_limits':['This separate reconciliation grants64 unique outer-ledger admissions once; original overlays retain zero automatic outer credit.','Exact source/annotation/hash integrity and declared reading gates do not establish scientific correctness, independent replication or MVA efficacy.','Prior1210 records and all3598 selection objects are preserved by pinned files and canonical object hashes.','No raw abstracts or duplicated scientific annotations are included.'],'records':records}

def literal_log(path,source_index):
    lines=[json.loads(s) for s in Path(path).read_text().splitlines()]
    require(len(lines)>=3,'Incomplete literal log: '+str(path));rows=lines[1:-1]
    require(1<=len(rows)<=4 and lines[0]=={'corpus_sha256':SOURCE_SHA,'requested_records':len(rows)},'Literal log header/count mismatch')
    require(lines[-1]=={'end_of_batch':True,'emitted_records':len(rows),'total_abstract_characters':sum(r['abstract_characters'] for r in rows)},'Missing literal end marker')
    for row in rows:
        src=source_index[row['source_id']];raw=src.get('abstractText','')
        require(bool(raw) and row['abstractText']==raw and row['title']==src['title'],'Literal source mismatch')
        require(row['abstract_characters']==len(raw) and row['abstract_sha256']==hashlib.sha256(raw.encode()).hexdigest(),'Literal source length/hash mismatch')
    return rows,lines[-1]

def full_replay(state,admission,work_root,corpus):
    work_root=Path(work_root);corpus=Path(corpus);require(digest(corpus)==SOURCE_SHA,'Corpus changed')
    source_index={key(r):r for r in read_json(corpus)['results']}
    b=work_root/'food-common-abstracts-20260923';old=work_root/'food-abstracts-20260923'
    require(digest(old/'review.json')==FOOD_SOURCE_SHA and digest(b/'common-review.json')==COMMON_SHA and digest(b/'common-validation.json')==COMMON_VALIDATION_SHA,'External overlay input changed')
    original_food=read_json(old/'review.json');require(len(original_food['records'])==67,'Original food count differs')
    for full,compact in zip(original_food['records'],state['layers']['literal_food']['records']):
        require({k:v for k,v in full.items() if k!='source_metadata'}=={k:v for k,v in compact.items() if k!='source_metadata'},'Food scientific annotation changed during compaction')
        require(compact['source_metadata']=={k:full['source_metadata'][k] for k in compact['source_metadata']},'Food compact metadata changed')
    checked_logs={};total=0
    for admitted in admission['records']:
        sid=admitted['source_id'];raw=source_index[sid].get('abstractText','')
        require(bool(raw) and hashlib.sha256(raw.encode()).hexdigest()==admitted['abstract_sha256'],'Admitted source missing/changed')
        for ref in admitted['evidence_references']:
            row=state['indexes'][ref['layer']][sid][1];ev=row['read_evidence'];base=old if ref['layer']=='literal_food' else b
            rel=Path(ev['log_file']);require(not rel.is_absolute() and '..' not in rel.parts,'Unsafe log path')
            path=base/rel;log_id=(ref['layer'],str(rel))
            if log_id not in checked_logs:
                require(digest(path)==ev.get('log_sha256',next((x['log_sha256'] for x in state['layers']['literal_food'].get('read_evidence_manifest',[]) if x['log_file']==str(rel)),None)),'Original reading log hash differs')
                rows,end=literal_log(path,source_index);checked_logs[log_id]={x['source_id']:x for x in rows};total+=end['total_abstract_characters']
            require(sid in checked_logs[log_id] and checked_logs[log_id][sid]['abstract_sha256']==admitted['abstract_sha256'],'Record not in referenced literal log')
    gate=admission['root_sample_gate']
    for group in gate['sample_groups']:
        require(digest(b/group['selection_file'])==group['selection_sha256'],'Predetermined sample file changed')
        chosen=read_json(b/group['selection_file']);require(object_digest(chosen)==group['selection_objects_sha256'] and [key(x) for x in chosen]==group['source_ids'],'Sample selection objects differ')
    for log in gate['literal_logs']:
        require(digest(b/log['log_file'])==log['log_sha256'],'Root sample log changed')
        rows,end=literal_log(b/log['log_file'],source_index)
        require([x['source_id'] for x in rows]==log['source_ids'] and end['total_abstract_characters']==log['abstract_characters'],'Root sample literal identities/counts differ')
    return {'referenced_original_literal_logs_checked':len(checked_logs),'characters_in_referenced_logs':total,'root_sample_logs_checked':4,'root_sample_characters':19590,'new_reading_credit':0}

def check(root=None,admission_path=None,food_review=None,common_review=None,work_root=None,corpus=None):
    root=Path(root or default_root());state=prepare(root,food_review,common_review)
    admission=read_json(admission_path or root/ADMISSION_NAME)
    ids=[key(r) for r in admission['records']]
    require(len(ids)==len(set(ids))==64,'Duplicate/missing admission IDs')
    require(not set(ids)&state['accepted_ids'],'Admission overlaps prior accepted1210')
    require(ids==[key(r) for r in state['expected']],'Admission differs from exact pending food intersection')
    require(all(r.get('admission_credit')==1 and r.get('new_literal_reading_credit')==0 for r in admission['records']),'Missing/invalid per-record admission or reading credit')
    expected=construct_admission(state)
    require(admission==expected,'Admission source/hash/reference/gate/count objects differ from exact reconstruction')
    require((work_root is None)==(corpus is None),'Full replay requires work_root and corpus together')
    evidence=full_replay(state,admission,work_root,corpus) if work_root is not None else None
    return {'result':'PASS','mode':'full_food_source_evidence_replay' if evidence else 'compact_integrity','admitted_once':64,'new_literal_reading_credit':0,'prior_accepted_revalidation_records':1210,'cumulative_accepted_revalidation_records':1274,'combined_outer_ledger_records':3435,'pending_source_revalidation':2324,'original3598_selection_preserved':True,'prior1210_preserved':True,'root_sample_checks_declared':12,'root_sample_viewing_machine_proven':False,'scientific_annotations_duplicated':False,'external_evidence':evidence}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=default_root());p.add_argument('--admission',type=Path);p.add_argument('--food-review',type=Path);p.add_argument('--common-review',type=Path);p.add_argument('--work-root',type=Path);p.add_argument('--corpus',type=Path)
    a=p.parse_args();print(json.dumps(check(a.root,a.admission,a.food_review,a.common_review,a.work_root,a.corpus),indent=2))

if __name__=='__main__':main()
