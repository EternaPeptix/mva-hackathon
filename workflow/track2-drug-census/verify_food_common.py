#!/usr/bin/env python3
"""Single-scope common-food consolidation and replayable structural verifier.

This module verifies saved artifacts. It cannot prove that a person or model viewed
tool output, establish scientific correctness, or confer outer-corpus admission.
"""
from pathlib import Path
from collections import Counter
import argparse
import csv
import hashlib
import html
import json
import re

CORPUS_SHA = 'e5b8b017bb68af38818f1cc8c1bee47adb16afe0367981104a6b17886d2954d9'
READER_SHA = '354ac74ad402b41a25239a17400aec4f12aff8f4caf18bdec18a82059866dada'
MENTIONS_SHA = '9408dd04ea11f09299b7928c5d88e1b874c0041d847ea4c8ed0208159daf33a6'
PINNED = {
 'selection.json':'3ea6d4bda3faad87024cce42f827d8ece31ee5763ce868f6b568d9ffc0815e73',
 'part-a-selection.json':'39cef55c6d0e86308d0de84fa14bdc563bf30589282ac33e70f5e4ad7669c961',
 'part-b-selection.json':'9b04adbc4d28ca77ddd65735ee8148d34d13e7be4e7fbf81cf02478cd6089065',
 'part-a-review.json':'264e7d747af02418c10662419399906d0e58f2721f85f24918e4f977753d5aa2',
 'part-b-review.json':'369af34cb426912f5588f6e3577b3fa01118010e2264dcaba9a7231e20215d33',
 'part-a-batch-evidence.json':'0e05ff0af01c8ffeed08b4d65f8014f0de6345b0e301810f8f3bf0ff665e396b',
 'part-b-batch-evidence.json':'12f82e9e418a5e0a81c00d28783275b2f480afe20897844faa41772a0b69daef',
}
LITERAL_SELECTION_SHA = 'b676239a155b56d781b4f1f3165067ab646d67b8b9a36cacb34a918247270807'
LITERAL_REVIEW_SHA = 'f4395a70e1e499c8e0971bf8a884c74ed09ce2a6a1efccb1138e2e02974a0fb1'
META_KEYS = ['id','source','pmid','pmcid','doi','title','pubTypeList','publicationStatus','commentCorrectionList','versionList','versionNumber','firstPublicationDate']
EXPECTED_PARTS = {'a':{'records':77,'queries':81,'reads':76,'batches':19,'characters':128599}, 'b':{'records':77,'queries':80,'reads':77,'batches':20,'characters':131919}}
REPO = Path(__file__).resolve().parents[2]
TITLE_FORMAT_TAG = re.compile(
    r"</?(?:abstract|title|p|sec|label|italic|bold|i|b|em|strong|u|sup|sub|span|div|br|xref|ext-link|styled-content|named-content|inline-formula|disp-formula|article-title|list|list-item|h[1-6]|table|tr|td|th|tbody|thead|footnote|fn|sc|monospace)"
    r"(?:\s+[A-Za-z_:][\w:.-]*(?:\s*=\s*(?:\"[^\"]*\"|'[^']*'|[^\s<>]+))?)*\s*/?>", re.IGNORECASE)

def inventory_title_only(title):
    # Exact inventory display-title convention from broad_hub_screen.clean.
    # This is NEVER applied to abstract strings or used for reading credit.
    for _ in range(3):
        decoded=html.unescape(title)
        if decoded==title: break
        title=decoded
    return TITLE_FORMAT_TAG.sub(' ',title)

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def load(path):
    return json.loads(Path(path).read_text())

def canonical_sha(obj):
    return hashlib.sha256(json.dumps(obj,ensure_ascii=False,sort_keys=True,separators=(',',':')).encode()).hexdigest()

def string_leaves(obj):
    if isinstance(obj,str):
        yield obj
    elif isinstance(obj,dict):
        require('abstractText' not in obj,'Raw-abstract field in output annotations')
        for value in obj.values(): yield from string_leaves(value)
    elif isinstance(obj,list):
        for value in obj: yield from string_leaves(value)

def require(condition, message):
    if not condition:
        raise ValueError(message)

def sid(row):
    return row['source'] + ':' + str(row.get('record_id',row.get('id')))

def query_keys(records):
    return [(x['source_row_id'],x['query'],r['source'],r['record_id']) for r in records for x in r['food_matches']]

def status_summary(source_records):
    retracted, preprints, reviews, cases, errata, preprint_links, comments, versions = [],[],[],[],[],[],[],[]
    for row in source_records:
        key=sid(row); pubs=row.get('pubTypeList',{}).get('pubType',[])
        links=row.get('commentCorrectionList',{}).get('commentCorrection',[])
        if row['source']=='PPR' or any('preprint' in p.lower() for p in pubs): preprints.append(key)
        if any('review' in p.lower() for p in pubs): reviews.append(key)
        if any('case report' in p.lower() for p in pubs): cases.append(key)
        if any('retract' in p.lower() for p in pubs) or any('retract' in p.get('type','').lower() for p in links): retracted.append(key)
        for link in links:
            item={'source_id':key,'saved_link':link}
            if link.get('type')=='Erratum in': errata.append(item)
            if link.get('type')=='Preprint in': preprint_links.append(item)
            if link.get('type') in ['Comment in','Comment on']: comments.append(item)
        if 'versionList' in row or 'versionNumber' in row:
            versions.append({'source_id':key,**{k:row[k] for k in ['versionNumber','versionList'] if k in row}})
    return {'represented_source_ids':len(source_records),'retracted_record_count':len(retracted),'retracted_source_ids':retracted,'preprint_record_count':len(preprints),'preprint_source_ids':preprints,'review_publication_count':len(reviews),'review_source_ids':reviews,'case_report_publication_count':len(cases),'case_report_source_ids':cases,'records_with_erratum_links':len({x['source_id'] for x in errata}),'erratum_link_count':len(errata),'erratum_links':errata,'records_with_preprint_links':len({x['source_id'] for x in preprint_links}),'preprint_link_count':len(preprint_links),'preprint_links':preprint_links,'comment_link_count':len(comments),'comment_links':comments,'records_with_explicit_version_metadata':len(versions),'explicit_version_metadata':versions,'status_scope':'Frozen corpus metadata only; linked notices and versions are not newly read or independently validated by this consolidation.'}

def validate_common_overlay(workdir, corpus_path, mentions_path, reader_path, literal_selection_path, literal_review_path):
    """Return consolidated review and structural report; do not mutate any inputs."""
    b=Path(workdir)
    for name, expected in PINNED.items(): require(sha(b/name)==expected,f'Changed pinned input: {name}')
    for path, expected in [(corpus_path,CORPUS_SHA),(mentions_path,MENTIONS_SHA),(reader_path,READER_SHA),(literal_selection_path,LITERAL_SELECTION_SHA),(literal_review_path,LITERAL_REVIEW_SHA)]:
        require(sha(path)==expected,f'Changed source or prior-overlay input: {path}')
    selection=load(b/'selection.json'); corpus=load(corpus_path)
    source_index={(r['source'],str(r['id'])):r for r in corpus['results']}
    require(selection['source_sha256']==CORPUS_SHA and selection['mentions_sha256']==MENTIONS_SHA,'Selection source fingerprints differ')
    require((selection['record_count'],selection['mention_rows'],selection['available_abstracts'])==(154,161,153),'Wrong frozen selection size')
    with Path(mentions_path).open() as f: all_mentions=list(csv.DictReader(f,delimiter='\t'))
    common_mentions=[r for r in all_mentions if r['retrieval_policy']=='COMMON_OR_SHORT_NAME_CONTEXT_AMBIGUOUS']
    common_keys=[(r['source_row_id'],r['query'],r['source'],r['record_id']) for r in common_mentions]
    require(len(common_keys)==len(set(common_keys))==161,'Expected 161 unique FDA common-name mention rows')
    require(Counter(query_keys(selection['records']))==Counter(common_keys),'Selection does not exactly cover common-name rows')
    mention_index={(r['source_row_id'],r['query'],r['source'],r['record_id']):r for r in common_mentions}
    literal_selection=load(literal_selection_path); literal_review=load(literal_review_path)
    old_ids={sid(r) for r in literal_selection['records']}
    require(len(old_ids)==len(literal_selection['records'])==len(literal_review['records'])==67,'Wrong prior literal overlay size')
    require({r['source_id'] for r in literal_review['records']}==old_ids,'Prior literal selection/review IDs differ')
    require(sum(len(r['food_matches']) for r in literal_selection['records'])==73,'Wrong prior literal query count')
    require(literal_review['selection_sha256']==LITERAL_SELECTION_SHA and literal_review['corpus_sha256']==CORPUS_SHA,'Prior literal provenance differs')
    for old in literal_selection['records']:
        raw=source_index[(old['source'],old['record_id'])].get('abstractText','')
        require(bool(raw) and hashlib.sha256(raw.encode()).hexdigest()==old['abstract_sha256'],f'Prior literal source unavailable or changed: {sid(old)}')

    records=[]; part_provenance=[]; batch_checks=[]; source_checks=[]; followups=[]; selected_parts=[]; all_chunks=[]
    for part in ['a','b']:
        expect=EXPECTED_PARTS[part]
        sel=load(b/f'part-{part}-selection.json'); review=load(b/f'part-{part}-review.json'); manifest=load(b/f'part-{part}-batch-evidence.json')
        selected_parts.extend(sel['records'])
        require(sel['parent_selection_sha256']==PINNED['selection.json'],'Part selection has wrong parent')
        require(sel['source_sha256']==CORPUS_SHA and sel['mentions_sha256']==MENTIONS_SHA,'Part source fingerprint differs')
        require(review['selection_sha256']==PINNED[f'part-{part}-selection.json']==manifest['selection_sha256'],'Part selection hash inconsistent')
        require(review['selection_metadata']=={k:v for k,v in sel.items() if k!='records'},'Part selection metadata changed')
        require(review['corpus_sha256']==CORPUS_SHA,'Part corpus hash wrong')
        require(manifest.get('corpus_sha256',manifest.get('source_sha256'))==CORPUS_SHA,'Manifest corpus hash wrong')
        require(review['reader_script_sha256']==manifest['reader_script_sha256']==READER_SHA,'Reader hash wrong')
        require(review['selection_modified_by_reviewer'] is False,'Selection edit declared')
        require(len(review['records'])==len(sel['records'])==expect['records'],'Part record count wrong')
        require(sum(len(r['food_matches']) for r in sel['records'])==expect['queries'],'Part query count wrong')
        require(len(manifest['batches'])==expect['batches'],'Part batch count wrong')
        require([x['batch'] for x in manifest['batches']]==list(range(1,expect['batches']+1)),'Part batch numbering/order wrong')
        manifest_ref=review['read_evidence_manifest']
        if isinstance(manifest_ref,dict):
            require(manifest_ref['path']==f'part-{part}-batch-evidence.json' and manifest_ref['sha256']==sha(b/manifest_ref['path']),'Referenced manifest hash wrong')
        else: require(manifest_ref==f'part-{part}-batch-evidence.json','Referenced manifest path wrong')
        seen_logs={}; ordered_log_ids=[]; chars=0
        for batch in manifest['batches']:
            log_name=batch['log_file']; require(Path(log_name).name==log_name,'Unsafe/nonlocal log path')
            require(sha(b/log_name)==batch['log_sha256'],f'Log hash wrong: {log_name}')
            log_text=(b/log_name).read_text()
            if 'log_bytes' in batch: require(batch['log_bytes']==(b/log_name).stat().st_size,'Manifest log byte count wrong')
            if 'output_characters' in batch: require(batch['output_characters']==len(log_text),'Manifest output character count wrong')
            if 'corpus_sha256' in batch: require(batch['corpus_sha256']==CORPUS_SHA,'Batch corpus fingerprint wrong')
            if 'end_of_batch' in batch: require(batch['end_of_batch'] is True,'Manifest end marker declaration wrong')
            lines=[json.loads(s) for s in log_text.splitlines()]
            require(len(lines)>=3,f'Log has no complete record: {log_name}')
            literals=lines[1:-1]; count=len(literals)
            require(1<=count<=4,f'Batch not bounded to four: {log_name}')
            if 'requested_records' in batch: require(batch['requested_records']==count,'Manifest requested count wrong')
            require(lines[0]=={'corpus_sha256':CORPUS_SHA,'requested_records':count},f'Log header wrong: {log_name}')
            last={'end_of_batch':True,'emitted_records':count,'total_abstract_characters':sum(x['abstract_characters'] for x in literals)}
            require(lines[-1]==last,f'Missing/incorrect literal end marker: {log_name}')
            require(batch['emitted_records']==count and batch['total_abstract_characters']==last['total_abstract_characters'],'Manifest batch totals wrong')
            require(batch['source_ids']==[x['source_id'] for x in literals],'Manifest source order wrong')
            require(batch.get('end_marker_observed') is True and batch.get('untruncated_output_observed') is True,'Missing reviewer output declarations')
            require(isinstance(batch.get('tool_chunk_id'),str) and bool(batch['tool_chunk_id']),'Missing declared tool chunk ID')
            all_chunks.append((manifest['reader'],batch['tool_chunk_id']))
            for item in literals:
                require(set(item)=={'source_id','title','abstract_characters','abstract_sha256','abstractText'},'Unexpected literal reader schema')
                src_key=item['source_id'].split(':',1); src=source_index[tuple(src_key)]; raw=src.get('abstractText','')
                require(bool(raw) and item['abstractText']==raw,f'Literal mismatch: {item["source_id"]}')
                require(item['abstract_characters']==len(raw) and item['abstract_sha256']==hashlib.sha256(raw.encode()).hexdigest(),'Literal length/hash mismatch')
                require(item['title']==src['title'],'Literal title mismatch')
                require(item['source_id'] not in seen_logs,'Duplicate literal ID in part')
                seen_logs[item['source_id']]=(batch,item); ordered_log_ids.append(item['source_id'])
            chars+=last['total_abstract_characters']
            batch_checks.append({'part':part,'batch':batch['batch'],'log_file':log_name,'log_sha256':batch['log_sha256'],'source_ids':batch['source_ids'],'emitted_records':count,'abstract_characters':last['total_abstract_characters'],'log_exact_to_source':True,'literal_end_marker_valid':True,'tool_chunk_id_as_declared':batch['tool_chunk_id'],'reviewer':manifest['reader'],'viewing_and_untruncated_status':'Declared by original reader; not machine-proven by this check.'})
        require(len(seen_logs)==expect['reads'] and chars==expect['characters'],'Part read/character total wrong')
        require(ordered_log_ids==[sid(x) for x in sel['records'] if x['abstract_present']],'Literal order/available coverage wrong')
        for selected,out in zip(sel['records'],review['records']):
            key=sid(selected); src=source_index[(selected['source'],selected['record_id'])]; raw=src.get('abstractText','')
            require(all(out[k]==v for k,v in selected.items()),f'Selected field modified: {key}')
            require(out['source_id']==key and selected['title']==src['title'],'Review identity mismatch')
            metadata={k:src[k] for k in META_KEYS if k in src}
            require(out['source_metadata']==metadata,f'Missing/extra/changed compact metadata: {key}')
            require(all(k in metadata for k in ['id','source','title','pubTypeList']),f'Required bibliographic fields absent: {key}')
            require(selected['abstract_present'] is bool(raw),'Availability flag differs')
            require(selected['already_in_literal_food_overlay'] is (key in old_ids),'Prior-literal overlap flag differs')
            require(out['automatic_outer_corpus_credit']==out['new_appraisal_credit']==0,'Unexpected credit')
            require(out['full_text_read_in_this_layer'] is False,'Unexpected full-text credit')
            require('abstractText' not in out and 'abstractText' not in out['source_metadata'],'Raw abstract copied into review')
            require(len(out['food_query_reviews'])==len(selected['food_matches']),'Missing query annotation')
            for match,q in zip(selected['food_matches'],out['food_query_reviews']):
                require((q['source_row_id'],q['food_query'])==(match['source_row_id'],match['query']),'Query role assigned to wrong row')
                require(bool(q['query_role']),'Empty query role')
                mk=(match['source_row_id'],match['query'],selected['source'],selected['record_id'])
                require(mention_index[mk]['title']==inventory_title_only(selected['title']),'Mention display title mismatch')
                require(mention_index[mk]['abstract_sha256']==hashlib.sha256((raw or '').encode()).hexdigest(),'Mention abstract hash mismatch')
            pubs=src.get('pubTypeList',{}).get('pubType',[]); links=src.get('commentCorrectionList',{}).get('commentCorrection',[])
            require(out['publication_types_as_saved']==pubs and out['correction_and_version_links_as_saved']==links,'Status annotations differ from saved metadata')
            require(out['preprint_as_saved'] is (src['source']=='PPR' or any('preprint' in p.lower() for p in pubs)),'Preprint flag differs')
            require(out['retraction_flag_as_saved'] is (any('retract' in p.lower() for p in pubs) or any('retract' in p.get('type','').lower() for p in links)),'Retraction flag differs')
            if raw:
                require(selected['abstract_sha256']==hashlib.sha256(raw.encode()).hexdigest(),'Selected abstract hash differs')
                require(out['abstract_characters']==len(raw),'Review character count differs')
                require(key in seen_logs,'Available source not read in logs')
                batch,item=seen_logs[key]; ev=out['read_evidence']
                require(ev is not None and ev['batch']==batch['batch'] and ev['log_file']==batch['log_file'] and ev['log_sha256']==batch['log_sha256'],'Record log reference wrong')
                require(ev['abstract_sha256']==item['abstract_sha256'] and ev['abstract_characters']==len(raw),'Record reading identity wrong')
                require(ev['end_marker_observed'] is True and ev['untruncated_output_observed'] is True and ev['tool_chunk_id']==batch['tool_chunk_id'],'Record declarations do not match manifest')
                require(out['readers']==[manifest['reader']] and out['read_date']==manifest.get('read_date',manifest.get('date')),'Reader/date mismatch')
                require(out.get('reading_credit',1)==1,'Available record reading credit contradicts declaration')
            else:
                require(key=='MED:18082453','Unexpected missing abstract')
                require(all(json.loads(mention_index[(m['source_row_id'],m['query'],selected['source'],selected['record_id'])]['matched_fields'])==['title'] for m in selected['food_matches']),'Missing source has non-title match claim')
                require(selected['abstract_sha256'] is None and out['abstract_characters']==0,'Missing record has nonempty abstract identity')
                require(out['read_evidence'] is None and out['readers']==[] and out['read_date'] is None and out['reading_credit']==0,'Missing record has reading claim')
                require(out['annotation_status']=='no_saved_abstract' and out['review_depth']=='not_read_no_saved_abstract','Missing record is not explicitly marked')
                require(all(q['query_role']=='not_assessed_no_saved_abstract' and q['efficacy_for_bub1b_mva']=='NOT_ASSESSED_NO_SAVED_ABSTRACT' for q in out['food_query_reviews']),'Missing record receives substance/result inference')
                require(key not in seen_logs,'Missing record appears in literal log')
            records.append(out)
            source_checks.append({'source_id':key,'part':part,'abstract_present':bool(raw),'abstract_sha256':selected['abstract_sha256'],'part_annotation_canonical_sha256':canonical_sha(out),'selection_fields_preserved':True,'compact_metadata_exact':True,'query_rows':len(selected['food_matches']),'literal_log_verified':bool(raw),'missing_source_has_zero_reading':not bool(raw)})
        followups.extend(review['primary_followup_flags'])
        part_provenance.append({'part':part,'selection_file':f'part-{part}-selection.json','selection_sha256':PINNED[f'part-{part}-selection.json'],'review_file':f'part-{part}-review.json','review_sha256':PINNED[f'part-{part}-review.json'],'manifest_file':f'part-{part}-batch-evidence.json','manifest_sha256':PINNED[f'part-{part}-batch-evidence.json'],'reader_as_declared':manifest['reader'],'counts':expect,'records_preserved_without_annotation_changes':True})
    require(selected_parts==selection['records'],'Parts do not reproduce parent selection/order')
    require(len(records)==len({r['source_id'] for r in records})==154,'Duplicate/missing consolidated IDs')
    output_strings=set(string_leaves(records))
    require(all(source_index[(r['source'],r['record_id'])].get('abstractText','') not in output_strings for r in selection['records'] if r['abstract_present']),'Raw abstract copied into another annotation field')
    require(len(all_chunks)==len(set(all_chunks))==39,'Duplicate/missing declared reader chunks')
    available={r['source_id'] for r in records if r['abstract_present']}; all_ids={r['source_id'] for r in records}
    overlap=all_ids & old_ids; new_ids=all_ids-old_ids
    require((len(available),len(overlap),len(new_ids),len(available-old_ids))==(153,8,146,145),'Common/previous overlap arithmetic differs')
    global_ids=all_ids|old_ids; global_available=available|old_ids
    require((len(global_ids),len(global_available))==(213,212),'Global source arithmetic differs')
    combined_query_keys=query_keys(selection['records'])+query_keys(literal_selection['records'])
    require(len(combined_query_keys)==len(set(combined_query_keys))==234,'Global query coverage duplicates/misses a row')
    require(Counter(combined_query_keys)==Counter((r['source_row_id'],r['query'],r['source'],r['record_id']) for r in all_mentions),'Global FDA mention coverage not exact')
    roles=Counter(q['query_role'] for r in records for q in r['food_query_reviews']); queries=Counter(q['food_query'] for r in records for q in r['food_query_reviews'])
    common_status=status_summary([source_index[(r['source'],r['record_id'])] for r in selection['records']])
    global_status=status_summary([source_index[tuple(k.split(':',1))] for k in sorted(global_ids)])
    summary={'represented_source_ids':154,'query_rows_annotated':161,'saved_abstract_reads_declared':153,'missing_abstract_source_ids':['MED:18082453'],'title_metadata_only_query_rows':1,'literal_batches':39,'literal_abstract_characters':260518,'prior_literal_food_overlap_ids':sorted(overlap),'prior_literal_food_overlap_count':8,'source_ids_not_in_prior_literal_overlay':146,'available_ids_not_in_prior_literal_overlay':145,'exact_query_role_counts':dict(sorted(roles.items())),'distinct_original_query_role_labels':len(roles),'query_name_counts':dict(sorted(queries.items())),'bounded_followup_count':len(followups),'bounded_followup_source_ids':[x['source_id'] for x in followups],'automatic_outer_corpus_credit':0,'new_appraisal_credit':0}
    summary['review_abstract_classifications_by_part_readers']=sum(r['review_by_abstract_reading'] for r in records)
    summary['reader_review_classification_differs_from_publication_type_ids']=[r['source_id'] for r in records if r['review_by_abstract_reading'] != any('review' in p.lower() for p in r['publication_types_as_saved'])]
    coverage={'prior_literal_source_ids':67,'prior_literal_query_rows':73,'common_source_ids':154,'common_available_abstracts':153,'overlap_source_ids':8,'unique_saved_abstract_ids':212,'unique_all_hit_source_ids':213,'all_query_rows_annotated':234,'missing_abstract_source_ids':['MED:18082453'],'arithmetic':'67 + 153 - 8 = 212 available/read-declared source IDs; 67 + 154 - 8 = 213 represented source IDs; 73 + 161 = 234 distinct query rows.','source_ids_are_not_independent_study_counts':True,'outer_admission_credit':0}
    provenance={'corpus_sha256':CORPUS_SHA,'mentions_sha256':MENTIONS_SHA,'reader_script_sha256':READER_SHA,'selection_sha256':PINNED['selection.json'],'pinned_input_hashes':PINNED,'previous_literal_selection_sha256':LITERAL_SELECTION_SHA,'previous_literal_review_sha256':LITERAL_REVIEW_SHA,'part_provenance':part_provenance}
    common={'schema_version':'food-common-name-consolidated-review/1','date':'2026-09-23','scope':'Exact consolidation of both common-name parts; original annotations preserved unchanged.','selection_metadata':{k:v for k,v in selection.items() if k!='records'},'provenance':provenance,'summary':summary,'global_food_overlay_coverage':coverage,'publication_status_as_saved':common_status,'global_food_publication_status_as_saved':global_status,'reading_boundary':'Original readers declare complete literal viewing. This consolidation verifies saved strings, hashes, manifests and declared chunk IDs; those machine checks cannot prove viewing and add no new reading credit.','root_scientific_sample_audit_status':'NOT_PROVIDED_TO_THIS_CONSOLIDATION','outer_admission_status':'NOT_PERFORMED: any proposed 64-ID admission belongs to separate root-owned eligibility and audit gates.','primary_followup_flags':followups,'records':records}
    validation={'schema_version':'food-common-name-structural-validation/1','date':'2026-09-23','result':'PASS','provenance':provenance,'summary':summary,'global_food_overlay_coverage':coverage,'publication_status_as_saved':common_status,'global_food_publication_status_as_saved':global_status,'checks':{'frozen_input_fingerprints':True,'all_161_common_FDA_mention_rows_exact':True,'all_234_total_FDA_mention_rows_exact':True,'all_154_source_identities_exact':True,'all_153_literal_strings_exact':True,'all_39_log_hashes_and_end_markers_valid':True,'all_declared_chunk_ids_present_unique_by_reader':True,'all_record_reading_references_match_manifests':True,'selection_fields_and_original_annotations_preserved':True,'compact_bibliographic_status_metadata_exact':True,'missing_abstract_explicit_with_zero_reading':True,'previous_overlap_flags_exact':True,'raw_abstracts_absent_from_consolidated_review':True,'zero_automatic_outer_or_appraisal_credit':True},'machine_check_limit':'Artifact integrity and consistency only. Actual viewing, untruncated UI delivery and scientific accuracy are original-reader declarations, not machine-proven facts. No new abstract reading or root sample audit is claimed.','source_checks':source_checks,'batch_checks':batch_checks}
    validation['representation_conventions']={'missing_abstract':{'source_id':'MED:18082453','corpus_abstract_field':'absent','selection_abstract_sha256':None,'mention_abstract_sha256':hashlib.sha256(b'').hexdigest(),'reading_credit':0},'inventory_display_title':'FDA mentions use the inventory display-title formatter; corpus, selection, review and literal-log titles remain unchanged. No abstract normalization is used.','formatted_title_source_ids':[sid(r) for r in selection['records'] if inventory_title_only(r['title'])!=r['title']],'formatted_title_query_rows':sum(len(r['food_matches']) for r in selection['records'] if inventory_title_only(r['title'])!=r['title'])}
    return common, validation

def make_findings(common):
    s=common['summary']; st=common['publication_status_as_saved']; gst=common['global_food_publication_status_as_saved']
    text=['Both common-name parts are consolidated without changing their record annotations. **154 source IDs and all 161 FDA common-name query rows are represented: 153 saved abstracts are declared read, and MED:18082453 remains explicitly unread because the saved corpus contains no abstract.** Its title and Comment-on metadata receive no experimental inference or reading credit.',
    'Structural checks passed against the frozen corpus, selection, FDA mention rows, both reviews, both manifests and all 39 literal logs. The logs contain 260,518 abstract characters. Exact strings, hashes, source order, end markers, per-record references and declared chunk IDs match. These checks verify artifacts; they do not prove actual viewing or scientific correctness. The original reading declarations belong to the two part reviewers. Root scientific sample-audit results have not been supplied to this consolidation.',
    '**Coverage remains separate from admission.** Eight IDs overlap the previous literal-food overlay. The common layer therefore represents 146 IDs absent from that overlay, of which 145 have saved abstracts. Across both food layers, 67 + 153 − 8 = **212 unique saved-abstract IDs**, while 67 + 154 − 8 = **213 total hit IDs**. All **234 distinct query rows** are annotated, including one title/metadata-only row. Source-ID counts are not independent-study counts. This package awards **zero automatic outer-corpus credit and zero new appraisal credit**; the proposed 64-ID outer admission must pass separate root-owned gates.',
    f'Saved common-layer status: **{st["retracted_record_count"]} retraction flags, {st["preprint_record_count"]} preprint records, {st["review_publication_count"]} review publications, {st["case_report_publication_count"]} case-report publications, {st["erratum_link_count"]} erratum links on {st["records_with_erratum_links"]} records, and {st["preprint_link_count"]} journal-to-preprint links**. The errata are attached to MED:22457777 (saved reference only, no linked ID) and MED:40744332 (MED:40883118 and MED:41807160). The preprints are PPR178023 and PPR624548; only PPR178023 has explicit saved version metadata (version 1 and its version list). Linked notices and versions were not newly read, so their contents/reasons are not inferred.',
    f'The original readers classified {s["review_abstract_classifications_by_part_readers"]} abstracts as reviews/syntheses, while publication metadata labels {st["review_publication_count"]} as reviews. MED:25844087 accounts for that distinction: its reader classification is preserved alongside its saved research-article/Journal Article types. The combined 213-ID food union has {gst["retracted_record_count"]} saved retraction flag, {gst["preprint_record_count"]} selected preprints, {gst["erratum_link_count"]} erratum links on {gst["records_with_erratum_links"]} records, and {gst["records_with_explicit_version_metadata"]} records with explicit version metadata. The full exact metadata inventory is retained in both JSON outputs; retraction reasons and linked-version equivalence are not invented.',
    f'The original role vocabulary is retained: **{s["distinct_original_query_role_labels"]} distinct query-role labels across 161 rows**. Exact counts for every role are in `summary.exact_query_role_counts`; counts are not collapsed into a new efficacy taxonomy. Common-name query totals are:',
    '| Query | Rows |\n|---|---:|\n'+'\n'.join(f'| {q} | {n} |' for q,n in s['query_name_counts'].items()),
    f'The **{s["bounded_followup_count"]} bounded follow-up IDs** below are questions, not new drug appraisals or nominations. The missing comment is included as a source-recovery question. Original model, compound/derivative/mixture/reagent distinctions, functional positives, nulls and BUB1B limitations remain in the unchanged source records.',
    ', '.join(s['bounded_followup_source_ids']) + '.',
    'The most consequential preserved distinctions include actual unsaturated-fatty-acid hepatocyte rescue with unnamed species at the abstract level; Sod1-oocyte segregation-error reduction by ascorbic acid or NAC without explicit drug-restored BUBR1 localization; folic-acid cell protection versus cow expression markers; collagen biomaterial/islet function versus dietary collagen; and Wy-14643 mouse liver recovery where fatty acids are downstream measurements. Several source-specific BUB1B effects instead involve genetic perturbation or cancer progression. These observations do not establish MVA rescue or justify promotion over provisional niacin/apitegromab; everolimus remains conditional.',
    'Replay the structural checks with `python3 build-common-review.py --check-only`. Optional path arguments support relocation. The verifier rejects missing or changed identity/metadata/log/end-marker evidence; its sole allowed absent abstract is the explicit zero-reading MED:18082453 record. The saved package includes no raw abstracts; the original literal strings remain only in the part logs.']
    return '\n\n'.join(text)+'\n'

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workdir',type=Path,default=Path(__file__).resolve().parent)
    parser.add_argument('--corpus',type=Path)
    parser.add_argument('--mentions',type=Path,default=REPO/'workflow/track2-drug-census/food-substances/mentions.tsv')
    parser.add_argument('--reader',type=Path,default=REPO/'workflow/track2-drug-census/read_abstracts.py')
    parser.add_argument('--literal-selection',type=Path)
    parser.add_argument('--literal-review',type=Path)
    parser.add_argument('--check-only',action='store_true')
    args=parser.parse_args(); b=args.workdir
    common,validation=validate_common_overlay(b,args.corpus or b.parent/'literature-gene-unrestricted.json',args.mentions,args.reader,args.literal_selection or b.parent/'food-abstracts-20260923/selection.json',args.literal_review or b.parent/'food-abstracts-20260923/review.json')
    findings=make_findings(common)
    if args.check_only:
        require(load(b/'common-review.json')==common,'Saved common review differs from exact recomputation')
        require(load(b/'common-validation.json')==validation,'Saved validation differs from exact recomputation')
        require((b/'common-findings.md').read_text()==findings,'Saved findings differ from exact recomputation')
    else:
        for name,obj in [('common-review.json',common),('common-validation.json',validation)]:
            (b/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n')
        (b/'common-findings.md').write_text(findings)
    print(json.dumps({'result':'PASS','mode':'check-only' if args.check_only else 'build','represented_ids':154,'query_rows':161,'literal_read_declarations':153,'missing':['MED:18082453'],'batches':39,'characters':260518,'unique_saved_abstract_ids_across_food_layers':212,'all_hit_ids_across_food_layers':213,'all_food_query_rows':234,'outer_credit':0,'artifacts':{name:sha(b/name) for name in ['common-review.json','common-validation.json','common-findings.md']}},indent=2))

if __name__=='__main__': main()
