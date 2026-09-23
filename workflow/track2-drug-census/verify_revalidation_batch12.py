#!/usr/bin/env python3
"""Validate batch12 after batch11 and the explicit food64 gate; compact or source replay.

Checks artifact integrity and declared reading, not scientific correctness.
"""
import argparse, hashlib, importlib.util, json, sys
from collections import Counter
from pathlib import Path
sys.dont_write_bytecode=True
SOURCE_SHA='e5b8b017bb68af38818f1cc8c1bee47adb16afe0367981104a6b17886d2954d9'
PARENT_SHA='a9a12a95285b7456be572811a0a03891748535abe263706184ef517069938fb2'
SELECTION_SHA='01f2fd084296906dbab53dfde4a1433567e65362599298471f4e3ec681ec578b'
FOOD_SHA='1fe105032d42ced1d3299e39595a0b3f5e5c7811d0ad48b01205292b2fb5e255'
FOOD_VERIFIER_SHA='7e5a9022833d50fd1676c6f0dda7474506fdabd2ac8e84fed2b7776d5e5222a5'
READER_SHA='354ac74ad402b41a25239a17400aec4f12aff8f4caf18bdec18a82059866dada'
FALLBACK=Path('/Users/jeweled/mva-hackathon/artifacts/source-revalidation-20260922/repo/workflow/track2-drug-census')
SAMPLE_INDICES=[0,15,32,63,64,79,96,127,128,143,160,191,192,207,224,255]
SAMPLE_LOGS=[('root-sample-a.jsonl', '92ba4382c046f6e755e8e7b5919cec2d9a9b5047c767d8b06308e2ac15c3a4d3', '4b4981', 4483), ('root-sample-b.jsonl', 'e5b33d7a244a5742d11a5d310adc1b58c5a0fba0cd9831d9526ee39340e24d0a', '652410', 6724), ('root-sample-c.jsonl', '8590f0c70ef95eef04b348754d321a405e3c7e0b15a6dd2bf8a51a6952d66d72', '77b72a', 5818), ('root-sample-d.jsonl', '9eca46fc0c5b9e5504e4b9f953579987beb56f617a9778c2957b43a1e0b1872e', '3da1c2', 4341)]
PREVIOUS_PINS={'outer-revalidation-batch11-review.json': '976f82f82b7aff76d26a828ee4b1d59025388d51ec369e74d4471133f1b57116', 'outer-revalidation-batch11-selection.json': '1bb2cd624b218b7bcb99e2e1fbe647b93620d53ad18a777b59f9c2d2d47b5d0a', 'verify_revalidation_batch11.py': '6cc0f840743c7d304ff87e0d7085ef9a6b4d82626bd0b589aa0d7033cc54b791'}
FIELD_MAP={'n':'batch_review_index','kind':'saved_content_kind','model':'model','intervention':'intervention_type','outcome':'outcome_direction','uncertainty':'uncertainty','mentioned':'mentioned_nonadministered_or_unresolved_agents','followup':'followup_reason','bub':'bub1b_relevance','source_issue':'source_issue'}

def require(x,msg):
    if not x: raise ValueError(msg)
def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def object_sha(x):return hashlib.sha256(json.dumps(x,ensure_ascii=False,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def key(x):return x['source']+':'+str(x.get('record_id',x.get('id')))
def default_root():
    for p in [Path(__file__).resolve().parent,*Path(__file__).resolve().parents]:
        if (p/'outer-revalidation-selection.json').exists():return p
    return FALLBACK
def import_file(path,name):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def baseline(root,selection):
    root=Path(root)
    require(sha(root/'outer-revalidation-selection.json')==PARENT_SHA,'Original3598 selection bytes changed')
    require(sha(root/'outer-food-admission.json')==FOOD_SHA and sha(root/'verify_food_admission.py')==FOOD_VERIFIER_SHA,'Food admission gate changed')

    for name,digest in PREVIOUS_PINS.items():require(sha(root/name)==digest,'Prior batch11 gate changed: '+name)
    previous=import_file(root/'verify_revalidation_batch11.py','batch12_previous_gate').check(root=root)
    require(previous['result']=='PASS' and previous['cumulative_accepted_revalidation_records']==1786,'Prior1786 gate not satisfied')
    parent=load(root/'outer-revalidation-selection.json');seen=set()
    for item in selection['excluded_reviews']:
        require(sha(root/item['file'])==item['sha256'],'Prior accepted file changed: '+item['file'])
        rows=load(root/item['file'])['records'];ids={key(x) for x in rows}
        require(len(ids)==len(rows)==item['record_count'] and not seen&ids,'Prior exclusion duplicate/count error')
        seen|=ids
    require(len(seen)==1786,'Expected all1786 exclusions, including food64')
    pending=[x for x in parent['records'] if key(x) not in seen]
    require(len(parent['records'])==3598 and len(pending)==1812,'Queue before batch12 differs')
    require(selection['records']==pending[:256],'Batch12 is not first256 exact remaining selection objects')
    require(selection['source_sha256']==SOURCE_SHA and selection['parent_selection_sha256']==PARENT_SHA,'Selection source provenance wrong')
    require(selection['prior_accepted_revalidation_records']==1786 and selection['pending_before']==1812 and selection['pending_after']==1556,'Selection counts wrong')
    return seen

def original_note(row):
    note={}
    for name in row['annotation_source_fields']:
        if name=='agents':note[name]=[[x['name'],x['role']] for x in row['administered_agents']]
        else:note[name]=row[FIELD_MAP[name]]
    return note

def publication_counts(rows):
    types=Counter();links=Counter();ret=[];pre=[];err=[];versions=[];lineage=[];reviews=[]
    for r in rows:
        m=r['publication_metadata_flags'];p=m.get('pubTypeList',{}).get('pubType',[]);c=m.get('commentCorrectionList',{}).get('commentCorrection',[])
        types.update(p);links.update(x.get('type','') for x in c)
        if 'Retracted Publication' in p:ret.append(key(r))
        if r['source']=='PPR' or any(x.lower()=='preprint' for x in p):pre.append(key(r))
        if any('erratum' in x.get('type','').lower() or 'correction' in x.get('type','').lower() for x in c):err.append(key(r))
        if m.get('versionList') or 'versionNumber' in m:versions.append(key(r))
        if any(x.get('type') in ['Preprint in','Preprint of','Update of','Update in'] for x in c):lineage.append(key(r))
        if any('review' in x.lower() for x in p):reviews.append(key(r))
    return {'publication_type_counts':dict(sorted(types.items())),'comment_correction_link_types':dict(sorted(links.items())),'retracted_records':ret,'preprint_records':pre,'records_with_erratum_or_correction_links':err,'records_with_version_metadata':versions,'records_with_preprint_or_update_links':lineage,'review_publication_records':reviews}

def expected_audit(selection):
    ids=[key(selection['records'][i]) for i in SAMPLE_INDICES]
    logs=[]
    for i,(name,h,chunk,chars) in enumerate(SAMPLE_LOGS):logs.append({'log_file':name,'log_sha256':h,'source_ids':ids[i*4:(i+1)*4],'total_abstract_characters':chars,'tool_output_chunk_id':chunk,'end_marker_observed':True,'tool_output_untruncated':True})
    return {'status':'COMPLETE_DECLARED_BY_ROOT','reviewer':'/root','predetermined_zero_based_indices':SAMPLE_INDICES,'source_ids':ids,'checked_count':16,'literal_read_count':16,'literal_evidence':logs,'annotation_comparison_tool_chunk_ids':['c4d113','6bbf4d','b601fa','8dabee'],'annotation_comparison_complete_declared':True,'substantive_correction_count':0,'additional_unique_reading_credit':0,'scope':'Sixteen predetermined literal source/annotation spot-checks only; no independent full 256-record scientific review or error-rate claim. Actual viewing/comparison is a root declaration, not machine-proven by hashes.'}


def check_log(path,expected_ids,source_index=None):
    lines=[json.loads(x) for x in Path(path).read_text().splitlines()]
    require(len(lines)==6,'Expected header, four literal records and end marker: '+str(path))
    rows=lines[1:-1]
    require(lines[0]=={'corpus_sha256':SOURCE_SHA,'requested_records':4},'Invalid log source/header')
    require([x['source_id'] for x in rows]==expected_ids,'Log record order differs')
    total=0
    for r in rows:
        raw=r.get('abstractText');require(isinstance(raw,str) and bool(raw),'Missing literal string')
        require(len(raw)==r['abstract_characters'] and hashlib.sha256(raw.encode()).hexdigest()==r['abstract_sha256'],'Literal length/hash mismatch')
        if source_index is not None:
            src=source_index[r['source_id']];require(raw==src['abstractText'] and r['title']==src['title'],'Literal differs from source')
        total+=len(raw)
    require(lines[-1]=={'end_of_batch':True,'emitted_records':4,'total_abstract_characters':total},'Missing/invalid end marker')
    return rows,total

def check(root=None,corpus=None,evidence_root=None,selection_path=None,review_path=None):
    root=Path(root or default_root());selection_path=Path(selection_path or root/'outer-revalidation-batch12-selection.json');review_path=Path(review_path or root/'outer-revalidation-batch12-review.json')
    require(sha(selection_path)==SELECTION_SHA,'Frozen batch12 selection changed')
    selection=load(selection_path);review=load(review_path);prior=baseline(root,selection)
    require((corpus is None)==(evidence_root is None),'Full replay requires corpus and evidence_root together')
    full=corpus is not None;source_index=None
    if full:
        require(sha(corpus)==SOURCE_SHA,'Corpus bytes changed');source_index={key(x):x for x in load(corpus)['results']};evidence_root=Path(evidence_root)
    rows=review['records'];require(len(rows)==len({key(x) for x in rows})==256,'Missing/duplicate review records')
    require(not {key(x) for x in rows}&prior,'New records overlap prior1786')
    require(review['selection_sha256']==SELECTION_SHA and review['source_sha256']==SOURCE_SHA and review['parent_selection_sha256']==PARENT_SHA,'Review provenance differs')
    require(review['excluded_reviews']==selection['excluded_reviews'],'Prior acceptance manifest differs')
    require(review['reader_script_sha256']==READER_SHA and sha(root/'read_abstracts.py')==READER_SHA,'Reader changed')
    counts={'selected_count':256,'reviewed_count':256,'pending_count':1556,'prior_accepted_revalidation_records':1786,'cumulative_accepted_revalidation_records':2042,'previous_outer_ledger_count':3947,'combined_outer_ledger_count':4203,'batch_count':64,'total_abstract_characters':370003}
    require(all(review.get(k)==v for k,v in counts.items()),'Review coverage counts differ')
    require(review['independent_audit']==expected_audit(selection),'Root sample gate differs')
    require(review['independent_audit']['status']=='COMPLETE_DECLARED_BY_ROOT' and review['independent_audit']['annotation_comparison_complete_declared'] is True,'Root annotation comparison still pending')
    manifest=review['batch_read_evidence'];require(len(manifest)==64,'Expected64 bounded batch declarations')
    total=0;chunk_ids=[];notes_rebuilt=[]
    for i,(chosen,row) in enumerate(zip(selection['records'],rows),1):
        require(all(row.get(k)==v for k,v in chosen.items()),'Original selected object changed')
        require(row['batch_review_index']==i and row['review_depth']=='complete_saved_abstract' and row['annotation_status']=='source_revalidated' and row['reading_credit']==1,'Missing complete-record credit')
        require(row['full_text_read'] is False,'Unexpected full-text credit')
        expected_reader='/root/ms023_primary' if i<=128 else '/root/bromosporine_primary'
        require(row['readers']==[expected_reader],'Per-record reader attribution differs')
        note=original_note(row);require(object_sha(note)==row['original_annotation_sha256'],'Original scientific annotation changed');notes_rebuilt.append(note)
        require(row['note']==row['outcome_direction']+' '+row['uncertainty'],'Derived note differs')
        require(row['requires_primary_followup']==('followup' in row['annotation_source_fields']),'Follow-up flag differs')
        require(object_sha(row['publication_metadata_flags'])==row['publication_metadata_sha256'],'Metadata hash differs')
        require('abstractText' not in row and 'abstractText' not in row['publication_metadata_flags'],'Raw abstract in ledger')
        if full:
            src=source_index[key(row)];require(row['publication_metadata_flags']=={k:v for k,v in src.items() if k not in ['abstractText','id','source','title']},'Missing/changed source metadata')
            require(hashlib.sha256(src['abstractText'].encode()).hexdigest()==row['abstract_sha256'] and row['title']==src['title'],'Source abstract/title mismatch')
    for n,batch in enumerate(manifest,1):
        chosen=rows[(n-1)*4:n*4];ids=[key(x) for x in chosen]
        expected_reader='/root/ms023_primary' if n<=32 else '/root/bromosporine_primary'
        require(batch['readers']==[expected_reader],'Batch reader attribution differs')
        require(batch['batch']==n and batch['source_ids']==ids and batch['log_file']==f'batch-{n:02d}.jsonl','Batch order/IDs differ')
        require(batch['emitted_records']==batch['requested_records']==4 and batch['end_of_batch'] is True and batch['end_marker_observed'] is True and batch['tool_output_untruncated'] is True,'Incomplete batch declaration')
        require(bool(batch['tool_output_chunk_id']),'Missing observed chunk declaration');chunk_ids.append(batch['tool_output_chunk_id'])
        require(batch['total_abstract_characters']==sum(x['read_evidence']['abstract_characters'] for x in chosen),'Batch character total differs')
        for row in chosen:
            ev=row['read_evidence'];require(ev['batch']==n and ev['log_file']==batch['log_file'] and ev['log_sha256']==batch['log_sha256'] and ev['tool_output_chunk_id']==batch['tool_output_chunk_id'] and ev['abstract_sha256']==row['abstract_sha256'] and ev['end_marker_observed'] is True and ev['tool_output_untruncated'] is True,'Per-record reading reference differs')
        if full:
            require(sha(evidence_root/batch['log_file'])==batch['log_sha256'],'Batch log changed')
            literal,t=check_log(evidence_root/batch['log_file'],ids,source_index);require(t==batch['total_abstract_characters'],'Source batch count differs')
            item=review['annotation_inputs'][n-1];require(item['file']==f'notes-batch-{n:02d}.json' and sha(evidence_root/item['file'])==item['sha256'],'Notes input changed')
            rawnotes=load(evidence_root/item['file']);require(rawnotes['notes']==notes_rebuilt[(n-1)*4:n*4] and rawnotes['chunk_id']==batch['tool_output_chunk_id'] and rawnotes['end_marker_observed'] is True and rawnotes['untruncated_output_observed'] is True,'Original notes/declaration differ')
        total+=batch['total_abstract_characters']
    require(total==370003 and len(set(chunk_ids))==64,'Literal total/chunk identities differ')
    require(review['publication_metadata_counts']==publication_counts(rows),'Publication status counts differ')
    require(review['source_content_kinds']==dict(sorted(Counter(x['saved_content_kind'] for x in rows).items())),'Content-kind counts differ')
    require(review['followup_count']==sum(x['requires_primary_followup'] for x in rows),'Bounded follow-up count differs')
    require(review['record_objects_sha256']==object_sha(rows),'Review record-object digest differs')
    if full:
        ep=evidence_root/review['evidence_manifest']['file'];require(sha(ep)==review['evidence_manifest']['sha256'],'Evidence manifest changed')
        require(load(ep)['batches']==manifest,'Evidence manifest batches differ')
        for log in review['independent_audit']['literal_evidence']:
            require(sha(evidence_root/log['log_file'])==log['log_sha256'],'Root sample log changed');_,t=check_log(evidence_root/log['log_file'],log['source_ids'],source_index);require(t==log['total_abstract_characters'],'Root sample characters differ')
    return {'result':'PASS','mode':'full_source_evidence_replay' if full else 'compact_integrity','reviewed_records':256,'prior_accepted_revalidation_records':1786,'cumulative_accepted_revalidation_records':2042,'combined_outer_ledger_records':4203,'pending_source_revalidation':1556,'batch_count':64,'total_abstract_characters':370003,'independent_cross_reader_checks_declared':16,'root_viewing_machine_proven':False,'scientific_interpretations_machine_validated':False,'prior1786_and_original3598_preserved':True,'followup_count':review['followup_count'],'publication_metadata_counts':publication_counts(rows)}

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=default_root());p.add_argument('--corpus',type=Path);p.add_argument('--evidence-root',type=Path);p.add_argument('--selection',type=Path);p.add_argument('--review',type=Path);a=p.parse_args();print(json.dumps(check(a.root,a.corpus,a.evidence_root,a.selection,a.review),indent=2))
