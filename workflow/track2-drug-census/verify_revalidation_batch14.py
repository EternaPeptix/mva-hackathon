#!/usr/bin/env python3
"""Validate batch14 after batch13 and the explicit food64 gate; compact or source replay.

Checks artifact integrity and declared reading, not scientific correctness.
"""
import argparse, hashlib, importlib.util, json, sys
from collections import Counter
from pathlib import Path
sys.dont_write_bytecode=True
SOURCE_SHA='e5b8b017bb68af38818f1cc8c1bee47adb16afe0367981104a6b17886d2954d9'
PARENT_SHA='a9a12a95285b7456be572811a0a03891748535abe263706184ef517069938fb2'
SELECTION_SHA='a0c4c8ea38a0af241123c472bef40f3545a425841481741db691e4e0d373accc'
FOOD_SHA='1fe105032d42ced1d3299e39595a0b3f5e5c7811d0ad48b01205292b2fb5e255'
FOOD_VERIFIER_SHA='7e5a9022833d50fd1676c6f0dda7474506fdabd2ac8e84fed2b7776d5e5222a5'
READER_SHA='354ac74ad402b41a25239a17400aec4f12aff8f4caf18bdec18a82059866dada'
FALLBACK=Path('/Users/jeweled/mva-hackathon/artifacts/source-revalidation-20260922/repo/workflow/track2-drug-census')
SAMPLE_INDICES=[0,31,64,127,128,159,192,255,256,287,320,383,384,415,448,511]
SAMPLE_LOGS=[('root-sample-a.jsonl', '874b759b4aaa237bcae5f4e84d9ee760cd5102cdeeec6c84c195d0acafc6d828', '77ea9d', 5822), ('root-sample-b.jsonl', '7c53d0f896509b7fd6e7a80104ec681dbe61314074c7facabbfb473edd92002d', '06b9b6', 3479), ('root-sample-c.jsonl', '9fd3a4f1534ae9045f490e9da812dde1ce2761fac6c27b48a22d117923f91021', 'e83a04', 4422), ('root-sample-d.jsonl', '14e9aba26738120985b7f86958504817ac1e5a99be7881f0f7d275be654c7ce1', '382245', 6711)]
PREVIOUS_PINS={'outer-revalidation-batch13-review.json': '05fadd536221b75c2d8d31472d27f6aeb806c67bd772c2dfe21f92fcdf615a65', 'outer-revalidation-batch13-selection.json': '8e70719097bd9a38a76b5d81a7ee0000d890c0fa8ca79135b0e0ce9fa592f328', 'verify_revalidation_batch13.py': '387402df81012b8e372513a9592eba6c754a783466c378ece7362c758dc83ca4'}
ROOT_COMPARISONS=[('root-compared-notes-1.json', 'c59bc5ee64a2730a9ac3eae611ab62f868788a51f1cb5dea2dc41b6e76bf63e9', 'f57dd5'), ('root-compared-notes-2.json', '0d83a48ed2ec0709e5e9e797241c5152ec6b24b465cd8f0279009fdd9e0701d6', '496116'), ('root-compared-notes-3.json', '9f164968c121958d885d6a990b9192de5137a2e0cbd94d7b3409cd5b2ac233da', '738551'), ('root-compared-notes-4.json', 'c03083964bf70e0a9cd0f61b42e3dc2e0ac67b4f54ee740ddc9fc4487964cf28', 'b9e693')]
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

    for name,digest in PREVIOUS_PINS.items():require(sha(root/name)==digest,'Prior batch13 gate changed: '+name)
    previous=import_file(root/'verify_revalidation_batch13.py','batch14_previous_gate').check(root=root)
    require(previous['result']=='PASS' and previous['cumulative_accepted_revalidation_records']==2554,'Prior2554 gate not satisfied')
    parent=load(root/'outer-revalidation-selection.json');seen=set()
    for item in selection['excluded_reviews']:
        require(sha(root/item['file'])==item['sha256'],'Prior accepted file changed: '+item['file'])
        rows=load(root/item['file'])['records'];ids={key(x) for x in rows}
        require(len(ids)==len(rows)==item['record_count'] and not seen&ids,'Prior exclusion duplicate/count error')
        seen|=ids
    require(len(seen)==2554,'Expected all2554 exclusions, including food64')
    pending=[x for x in parent['records'] if key(x) not in seen]
    require(len(parent['records'])==3598 and len(pending)==1044,'Queue before batch14 differs')
    require(selection['records']==pending[:512],'Batch14 is not first512 exact remaining selection objects')
    require(selection['source_sha256']==SOURCE_SHA and selection['parent_selection_sha256']==PARENT_SHA,'Selection source provenance wrong')
    require(selection['prior_accepted_revalidation_records']==2554 and selection['pending_before']==1044 and selection['pending_after']==532,'Selection counts wrong')
    return seen

def original_note(row):
    note={}
    for name in row['annotation_source_fields']:
        if name=='agents':note[name]=[[x['name'],x['role']] for x in row['administered_agents']]
        else:note[name]=row[FIELD_MAP[name]]
    return note

def publication_counts(rows):
    types=Counter();links=Counter();ret=[];pre=[];err=[];versions=[];lineage=[];reviews=[];concerns=[]
    for r in rows:
        m=r['publication_metadata_flags'];p=m.get('pubTypeList',{}).get('pubType',[]);c=m.get('commentCorrectionList',{}).get('commentCorrection',[])
        types.update(p);links.update(x.get('type','') for x in c)
        if 'Retracted Publication' in p:ret.append(key(r))
        if r['source']=='PPR' or any(x.lower()=='preprint' for x in p):pre.append(key(r))
        if any('erratum' in x.get('type','').lower() or 'correction' in x.get('type','').lower() for x in c):err.append(key(r))
        if m.get('versionList') or 'versionNumber' in m:versions.append(key(r))
        if any('expression of concern' in x.get('type','').lower() for x in c):concerns.append(key(r))
        if any(x.get('type') in ['Preprint in','Preprint of','Update of','Update in'] for x in c):lineage.append(key(r))
        if any('review' in x.lower() for x in p):reviews.append(key(r))
    return {'publication_type_counts':dict(sorted(types.items())),'comment_correction_link_types':dict(sorted(links.items())),'retracted_records':ret,'preprint_records':pre,'records_with_erratum_or_correction_links':err,'records_with_version_metadata':versions,'records_with_preprint_or_update_links':lineage,'review_publication_records':reviews,'records_with_expression_of_concern_links':concerns}

def expected_audit(selection):
    ids=[key(selection['records'][i]) for i in SAMPLE_INDICES]
    logs=[]
    for i,(name,h,chunk,chars) in enumerate(SAMPLE_LOGS):logs.append({'log_file':name,'log_sha256':h,'source_ids':ids[i*4:(i+1)*4],'total_abstract_characters':chars,'tool_output_chunk_id':chunk,'end_marker_observed':True,'tool_output_untruncated':True})
    return {'status':'COMPLETE_DECLARED_BY_ROOT','reviewer':'/root','predetermined_zero_based_indices':SAMPLE_INDICES,'source_ids':ids,'checked_count':16,'literal_read_count':16,'literal_evidence':logs,'annotation_comparison_tool_chunk_ids':[x[2] for x in ROOT_COMPARISONS],'annotation_comparison_evidence':[{'file':x[0],'sha256':x[1],'tool_output_chunk_id':x[2]} for x in ROOT_COMPARISONS],'annotation_comparison_complete_declared':True,'substantive_correction_count':0,'additional_unique_reading_credit':0,'scope':'Sixteen predetermined literal source/annotation spot-checks only; no independent full 512-record scientific review or error-rate claim. Actual viewing/comparison is a root declaration, not machine-proven by hashes.'}


def check_log(path,expected_ids,source_index=None):
    text=Path(path).read_text();require(text.endswith('\n'),'JSONL missing final newline: '+str(path))
    # JSON strings may legally contain Unicode line separators; only LF delimits JSONL.
    lines=[json.loads(x) for x in text[:-1].split('\n')]
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
    root=Path(root or default_root());selection_path=Path(selection_path or root/'outer-revalidation-batch14-selection.json');review_path=Path(review_path or root/'outer-revalidation-batch14-review.json')
    require(sha(selection_path)==SELECTION_SHA,'Frozen batch14 selection changed')
    selection=load(selection_path);review=load(review_path);prior=baseline(root,selection)
    require((corpus is None)==(evidence_root is None),'Full replay requires corpus and evidence_root together')
    full=corpus is not None;source_index=None
    if full:
        require(sha(corpus)==SOURCE_SHA,'Corpus bytes changed');source_index={key(x):x for x in load(corpus)['results']};evidence_root=Path(evidence_root)
    rows=review['records'];require(len(rows)==len({key(x) for x in rows})==512,'Missing/duplicate review records')
    require(not {key(x) for x in rows}&prior,'New records overlap prior2554')
    require(review['selection_sha256']==SELECTION_SHA and review['source_sha256']==SOURCE_SHA and review['parent_selection_sha256']==PARENT_SHA,'Review provenance differs')
    require(review['excluded_reviews']==selection['excluded_reviews'],'Prior acceptance manifest differs')
    require(review['reader_script_sha256']==READER_SHA and sha(root/'read_abstracts.py')==READER_SHA,'Reader changed')
    counts={'selected_count':512,'reviewed_count':512,'pending_count':532,'prior_accepted_revalidation_records':2554,'cumulative_accepted_revalidation_records':3066,'previous_outer_ledger_count':4715,'combined_outer_ledger_count':5227,'batch_count':128,'total_abstract_characters':764721}
    require(all(review.get(k)==v for k,v in counts.items()),'Review coverage counts differ')
    require(review['independent_audit']==expected_audit(selection),'Root sample gate differs')
    require(review['independent_audit']['status']=='COMPLETE_DECLARED_BY_ROOT' and review['independent_audit']['annotation_comparison_complete_declared'] is True,'Root annotation comparison still pending')
    require(len(ROOT_COMPARISONS)==4,'Root comparison evidence incomplete')
    manifest=review['batch_read_evidence'];require(len(manifest)==128,'Expected 128 bounded batch declarations')
    require(len(review['annotation_inputs'])==128,'Expected 128 annotation inputs')
    total=0;chunk_ids=[];notes_rebuilt=[]
    for i,(chosen,row) in enumerate(zip(selection['records'],rows),1):
        require(all(row.get(k)==v for k,v in chosen.items()),'Original selected object changed')
        require(row['batch_review_index']==i and row['review_depth']=='complete_saved_abstract' and row['annotation_status']=='source_revalidated' and row['reading_credit']==1,'Missing complete-record credit')
        require(row['full_text_read'] is False,'Unexpected full-text credit')
        expected_reader='/root/ms023_primary' if i<=256 else '/root/bromosporine_primary'
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
        expected_reader='/root/ms023_primary' if n<=64 else '/root/bromosporine_primary'
        require(batch['readers']==[expected_reader],'Batch reader attribution differs')
        require(batch['batch']==n and batch['source_ids']==ids and batch['log_file']==f'batch-{n:02d}.jsonl','Batch order/IDs differ')
        require(batch['emitted_records']==batch['requested_records']==4 and batch['end_of_batch'] is True and batch['end_marker_observed'] is True and batch['tool_output_untruncated'] is True,'Incomplete batch declaration')
        require(bool(batch['tool_output_chunk_id']),'Missing observed chunk declaration');chunk_ids.append(batch['tool_output_chunk_id'])
        require(batch['total_abstract_characters']==sum(x['read_evidence']['abstract_characters'] for x in chosen),'Batch character total differs')
        require(batch['abstract_sha256s']==[x['abstract_sha256'] for x in chosen],'Batch abstract digests differ')
        item=review['annotation_inputs'][n-1];require(item['record_count']==4 and item['annotation_objects_sha256']==object_sha(notes_rebuilt[(n-1)*4:n*4]),'Annotation input object digest differs')
        for row in chosen:
            ev=row['read_evidence'];require(ev['batch']==n and ev['log_file']==batch['log_file'] and ev['log_sha256']==batch['log_sha256'] and ev['tool_output_chunk_id']==batch['tool_output_chunk_id'] and ev['abstract_sha256']==row['abstract_sha256'] and ev['end_marker_observed'] is True and ev['tool_output_untruncated'] is True,'Per-record reading reference differs')
        if full:
            require(sha(evidence_root/batch['log_file'])==batch['log_sha256'],'Batch log changed')
            literal,t=check_log(evidence_root/batch['log_file'],ids,source_index);require(t==batch['total_abstract_characters'],'Source batch count differs')
            item=review['annotation_inputs'][n-1];require(item['file']==f'notes-batch-{n:02d}.json' and sha(evidence_root/item['file'])==item['sha256'],'Notes input changed')
            rawnotes=load(evidence_root/item['file']);require(rawnotes['notes']==notes_rebuilt[(n-1)*4:n*4] and rawnotes['chunk_id']==batch['tool_output_chunk_id'] and rawnotes['end_marker_observed'] is True and rawnotes['untruncated_output_observed'] is True,'Original notes/declaration differ')
        total+=batch['total_abstract_characters']
    require(total==764721 and len(set(chunk_ids))==128,'Literal total/chunk identities differ')
    require(review['publication_metadata_counts']==publication_counts(rows),'Publication status counts differ')
    require(review['source_content_kinds']==dict(sorted(Counter(x['saved_content_kind'] for x in rows).items())),'Content-kind counts differ')
    require(review['followup_count']==sum(x['requires_primary_followup'] for x in rows),'Bounded follow-up count differs')
    require(review['record_objects_sha256']==object_sha(rows),'Review record-object digest differs')
    if full:
        ep=evidence_root/review['evidence_manifest']['file'];require(sha(ep)==review['evidence_manifest']['sha256'],'Evidence manifest changed')
        evidence=load(ep);require(evidence['batches']==manifest,'Evidence manifest batches differ')
        require(evidence['source_sha256']==SOURCE_SHA and evidence['selection_sha256']==SELECTION_SHA and evidence['reader_script_sha256']==READER_SHA,'Evidence source provenance differs')
        require(evidence['record_count']==512 and evidence['batch_count']==128 and evidence['total_abstract_characters']==764721,'Evidence coverage differs')
        require(evidence['reader_scope']=={'/root/ms023_primary':{'first_batch':1,'last_batch':64,'record_count':256},'/root/bromosporine_primary':{'first_batch':65,'last_batch':128,'record_count':256}},'Evidence reader scope differs')
        for group,item in enumerate(review['independent_audit']['annotation_comparison_evidence'],1):
            ap=evidence_root/item['file'];require(sha(ap)==item['sha256'],'Root annotation snapshot changed');comp=load(ap)
            require(comp['selection_sha256']==SELECTION_SHA and comp['sample_group']==group and comp['end_of_annotation_group'] is True,'Root annotation snapshot header/end differs')
            chosen_indices=SAMPLE_INDICES[(group-1)*4:group*4]
            require([x['zero_based_index'] for x in comp['annotations']]==chosen_indices,'Root annotation sample order differs')
            for a in comp['annotations']:
                row=rows[a['zero_based_index']]
                require(a['source_id']==key(row) and a['title']==row['title'] and a['abstract_sha256']==row['abstract_sha256'] and a['original_note']==original_note(row),'Root compared annotation no longer matches')
        for log in review['independent_audit']['literal_evidence']:
            require(sha(evidence_root/log['log_file'])==log['log_sha256'],'Root sample log changed');_,t=check_log(evidence_root/log['log_file'],log['source_ids'],source_index);require(t==log['total_abstract_characters'],'Root sample characters differ')
    return {'result':'PASS','mode':'full_source_evidence_replay' if full else 'compact_integrity','reviewed_records':512,'prior_accepted_revalidation_records':2554,'cumulative_accepted_revalidation_records':3066,'combined_outer_ledger_records':5227,'pending_source_revalidation':532,'batch_count':128,'total_abstract_characters':764721,'independent_cross_reader_checks_declared':16,'root_viewing_machine_proven':False,'scientific_interpretations_machine_validated':False,'prior2554_and_original3598_preserved':True,'followup_count':review['followup_count'],'publication_metadata_counts':publication_counts(rows)}

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=default_root());p.add_argument('--corpus',type=Path);p.add_argument('--evidence-root',type=Path);p.add_argument('--selection',type=Path);p.add_argument('--review',type=Path);a=p.parse_args();print(json.dumps(check(a.root,a.corpus,a.evidence_root,a.selection,a.review),indent=2))
