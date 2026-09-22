#!/usr/bin/env python3
"""Reconcile the saved census; fail on missing rows/IDs or unsafe alias-policy regressions."""
import csv,json,hashlib
from collections import Counter
from pathlib import Path
P=Path(__file__).resolve().parent

def read(name):return list(csv.DictReader((P/name).open(),delimiter='\t'))
def main():
    c=json.loads((P/'coverage.json').read_text());l=json.loads((P/'literature-coverage.json').read_text())
    catalog=read('catalog.tsv');supplement=read('drugcentral-status-unverified.tsv');clinical=read('chembl-clinical-records.tsv')
    ids={r['catalog_id'] for r in catalog};all_ids=ids|{'DC'+r['drugcentral_id'] for r in supplement}
    assert len(catalog)==len(ids)==c['catalog_counts']['parent_or_unmatched_dc_entities']
    assert len(all_ids)==len(ids)+len(supplement)==l['catalog_entities_screened']
    assert len(clinical)==len({r['molecule_chembl_id'] for r in clinical})==c['source_counts']['chembl_unique_clinical_or_withdrawn_records']
    assert sum(len(json.loads(r['chembl_ids'])) for r in catalog)==len(clinical)
    assert sum(r['withdrawn_flag']=='True' for r in clinical)==c['source_counts']['chembl_withdrawn_flag']
    assert len(read('drugcentral-interactions.tsv'))==c['source_counts']['drugcentral_interactions_all']==c['source_counts']['drugcentral_interactions_screened']
    coverage=read('literature-drug-coverage.tsv');assert {r['catalog_id'] for r in coverage}==all_ids and len(coverage)==len(all_ids)
    assert all(r['literature_evidence_state']=='UNKNOWN_NO_SPECIFIC_ALIAS_TITLE_ABSTRACT_MATCH' for r in coverage if int(r['specific_alias_paper_count'])==0)
    for filename in ['literature-mentions.tsv','literature-ambiguous-mentions.tsv','drugcentral-interactions.tsv','catalog-appraisals.tsv']:
        assert all(r['catalog_id'] in all_ids for r in read(filename)),filename
    policy={r['normalized_alias']:r for r in read('literature-alias-policy.tsv')}
    for key in ['search','equal','success','murine','serine','atp']:
        if key in policy:assert policy[key]['policy']=='ambiguous_manual_review',key
    for key in ['rapamycin','decitabine','azd1152','everolimus']:
        assert policy[key]['policy']=='specific_alias_comention_only',key
    rows=read('bub1b-clinical-drug-activities.tsv')
    das=[r for r in rows if r['molecule_pref_name']=='DASATINIB'];assert len(das)==1
    assert das[0]['standard_relation']=='>' and float(das[0]['standard_value'])==75 and das[0]['standard_units']=='%'
    assert all(r['target_chembl_id'] in c['coverage']['bub1b_annotated_targets'] for r in rows)
    if (P/'target-assay-source-checks.json').exists():
        target=read('bub1b-all-target-activities.tsv')
        checks=json.loads((P/'target-assay-source-checks.json').read_text())
        assert hashlib.sha256((P/'bub1b-all-target-activities.tsv').read_bytes()).hexdigest()==checks['source_api_rows_sha256']
        activity={r['activity_id']:r for r in target}
        assert len(target)==len(activity)==checks['row_count']
        assert len({r['molecule_chembl_id'] for r in target})==checks['molecule_count']
        grouped=[aid for group in checks['document_groups'] for aid in group['activity_ids']]
        assert len(grouped)==len(set(grouped))==len(target) and set(grouped)==set(activity)
        for group in checks['document_groups']:
            assert all(activity[aid]['document_chembl_id']==group['document_chembl_id'] for aid in group['activity_ids'])
        for fix in checks['corrections']:
            source=activity[fix['activity_id']]
            assert source['standard_type']==fix['recorded_standard_type']=='Inhibition'
            assert source['molecule_chembl_id']==fix['recorded_molecule_chembl_id']
            assert source['standard_value']==fix['recorded_value']
            assert source['assay_description']==fix['recorded_assay_description']
            assert fix['source_endpoint']=='residual activity percent of control'
            assert fix['compound_and_concentration_error']==(fix['source_molecule_chembl_id']!=source['molecule_chembl_id'])
        # Four retained source cells and two assignment corrections, not rewritten API data.
        assert {x['source_cell'] for x in checks['corrections']}=={'F53','G53','H53','I53'}
        assert sum(x['compound_and_concentration_error'] for x in checks['corrections'])==2
    if (P/'literature-triage.tsv').exists():
        triage=read('literature-triage.tsv'); mentions=read('literature-mentions.tsv')
        key=lambda r:(r['catalog_id'],r['source'],r['record_id'])
        assert len(triage)==len(mentions) and {key(r) for r in triage}=={key(r) for r in mentions}
        assert sum(r['gene_literal_in_title_abstract']=='True' for r in triage)==l['specific_alias_direct_gene_comention_rows']
        for record_id,drug in [('41700350','SILVER'),('33244299','CAFFEIC ACID'),('17464181','DOXECITINE')]:
            hits=[r for r in triage if r['record_id']==record_id and r['drug_name']==drug]
            assert hits and all(r['triage_category']=='lexical_false_match' for r in hits),(record_id,drug)
        summary=json.loads((P/'literature-triage-summary.json').read_text())
        if 'final_unresolved_pass' in summary:
            final_ids=set(summary['final_unresolved_pass']['record_ids'])
            final_rows=[r for r in triage if r['record_id'] in final_ids and r['gene_literal_in_title_abstract']=='True']
            assert final_rows and all(r['abstract_available']=='True' for r in final_rows)
            assert all(r['triage_method']=='manual_abstract_review' for r in final_rows)
            assert all(not r['triage_category'].startswith('unresolved') for r in final_rows)
        assert all(r['actual_BUB1B_variant_muscle_rescue']=='NOT_ESTABLISHED_BY_THIS_ROW' for r in triage)
    if (P/'pubmed-supplement-summary.json').exists():
        ps=json.loads((P/'pubmed-supplement-summary.json').read_text())
        pr=read('pubmed-supplement-records.tsv'); pm={r['pmid']:r for r in pr}
        assert len(pr)==len(pm)==ps['fetched_supplement_records']
        assert dict(Counter(r['record_relation'] for r in pr))==ps['record_relations']
        extra={r['pmid'] for r in pr if r['same_pmid_in_base']=='False'}
        assert extra==set(ps['additional_record_ids'])
        assert ps['union_record_count']==ps['base_records']+len(extra)
        assert sum(r['abstract_available']=='False' for r in pr)==ps['missing_abstracts_in_supplement']
        assert sum(r['adds_previously_missing_abstract']=='True' for r in pr)==ps['newly_available_abstracts']
        for name,count in [('pubmed-supplement-mentions.tsv','specific_alias_rows'),('pubmed-supplement-ambiguous.tsv','ambiguous_alias_rows')]:
            hits=read(name);assert len(hits)==ps[count]
            for hit in hits:
                assert hit['catalog_id'] in all_ids and hit['record_id'] in pm
                assert hit['record_relation']==pm[hit['record_id']]['record_relation']
                assert hit['journal_update_in_base']==pm[hit['record_id']]['journal_update_in_base']
    if (P/'gene-record-screen-summary.json').exists():
        gs=json.loads((P/'gene-record-screen-summary.json').read_text());gr=read('gene-record-screen.tsv')
        assert len(gr)==len({(r['source'],r['record_id']) for r in gr})==gs['records']
        assert sum(int(r['abstract_available']) for r in gr)==gs['abstracts_available']
        assert sum(int(r['selected_title_screen']) for r in gr)==gs['selected_title_records']
        assert sum(int(r['manual_title_seen']) for r in gr)==gs['manual_title_seen_including_all_missing_abstract_titles']
        assert sum(int(r['manual_full_abstract_read']) for r in gr)==gs['manual_full_abstracts_read']
        for row in gr:
            if row['manual_full_abstract_read']=='1':assert row['abstract_available']=='1' and row['manual_appraisal']
            else:assert row['verification_level']=='not_primary_adjudicated' and not row['manual_appraisal']
        assert all('abstractText' not in row for row in gr)
    if (P/'gene-abstract-review-summary.json').exists():
        ar=read('gene-abstract-review.tsv')
        summary=json.loads((P/'gene-abstract-review-summary.json').read_text())
        notes=json.loads((P/'abstract-review-notes.json').read_text())
        assert summary['corpus_sha256']==notes['corpus_sha256']
        assert len(ar)==len({(r['source'],r['record_id']) for r in ar})==summary['source_records']
        assert dict(Counter(r['review_depth'] for r in ar))==summary['review_depth_counts']
        assert sum(r['abstract_available']=='True' for r in ar)==summary['available_saved_abstracts']
        actual={(r['source'],r['record_id']):r for r in ar}
        for note in notes['records']:
            row=actual[(note['source'],note['record_id'])]
            assert row['review_depth']==note['review_depth'] and row['note']==note['note']
            if note['review_depth']=='complete_saved_abstract':
                assert row['abstract_available']=='True' and row['abstract_sha256']==note['abstract_sha256']
                assert json.loads(row['administered_agents'])==note['administered_agents']
        assert summary['available_saved_abstracts_read']+summary['available_saved_abstracts_unread']==summary['available_saved_abstracts']
        assert summary['all_available_saved_abstracts_read']==(summary['available_saved_abstracts_unread']==0)
        assert all('abstractText' not in row for row in notes['records'])
    if (P/'missing-records-followup.json').exists():
        follow=json.loads((P/'missing-records-followup.json').read_text())
        historical=P/'abstract-review-notes.json'
        assert hashlib.sha256(historical.read_bytes()).hexdigest()==follow['historical_source_notes_sha256']
        notes=json.loads(historical.read_text())
        assert follow['corpus_sha256']==notes['corpus_sha256']
        old={(r['source'],r['record_id']):r for r in notes['records'] if r['review_depth']=='title_metadata_only'}
        recovered={(r['source'],r['record_id']):r for r in follow['records']}
        assert len(recovered)==len(follow['records'])==follow['reviewed_records'] and set(recovered)==set(old)
        depths=Counter()
        for key,row in recovered.items():
            assert row['title']==old[key]['title'] and row['historical_review_depth']==old[key]['review_depth']
            assert row['historical_abstract_sha256']==old[key].get('abstract_sha256','')
            assert row['status']=='reviewed' and row['source_url'] and row['interpretation'] and row['access_limitation']
            assert row['supplemental_depth'] in {'complete_body','complete_new_abstract','partial','title_only'}
            depths[row['supplemental_depth']]+=1
        assert dict(depths)==follow['supplemental_depth_counts']
        combined=Counter(complete_accessible_body=sum(r['review_depth']=='complete_accessible_body' for r in notes['records']))
        for depth,count in depths.items():combined['complete_accessible_body' if depth=='complete_body' else depth]+=count
        assert dict(combined)==follow['combined_depth_among_28_originally_missing_abstracts']
    if (P/'outer-corpus-review.json').exists():
        outer=json.loads((P/'outer-corpus-review.json').read_text())
        selection=json.loads((P/'outer-corpus-selection.json').read_text())
        ranking=json.loads((P/'outer-corpus-candidate-ranking.json').read_text())
        overlap=json.loads((P/'outer-corpus-overlap-audit.json').read_text())
        selected={(r['source'],r['record_id']):r for r in selection['records']}
        reviewed={(r['source'],r['record_id']):r for r in outer['records']}
        assert len(selected)==len(selection['records'])==selection['selected_count']
        assert len(reviewed)==len(outer['records'])==outer['reviewed_count']==outer['selected_count']
        assert set(selected)==set(reviewed) and len(reviewed)==selection['reviewed_count']
        assert selection['corpus_sha256']==l['corpus_sha256']
        assert selection['unrestricted_count']==l['corpus_records_scanned']
        assert selection['inner_title_abstract_count']==l['gene_title_abstract_query_count']
        assert selection['outer_difference_count']==selection['unrestricted_count']-selection['inner_title_abstract_count']
        assert len(reviewed)<=selection['outer_with_saved_abstract']<=selection['outer_difference_count']
        assert ranking['count']==len(ranking['candidates'])
        gene_ids={(r['source'],r['record_id']) for r in json.loads((P/'abstract-review-notes.json').read_text())['records']}
        assert not set(reviewed)&gene_ids
        assert overlap['selected_count']==len(reviewed)
        assert not overlap['inner_corpus_overlap'] and not overlap['authored_gene_review_overlap']
        for key,row in reviewed.items():
            assert row['title']==selected[key]['title'] and row['abstract_sha256']==selected[key]['abstract_sha256']
            assert len(row['abstract_sha256'])==64 and row['review_depth']=='complete_saved_abstract'
            assert row['note'] and row['model'] and isinstance(row['requires_primary_followup'],bool)
            assert 'abstractText' not in row and 'abstract' not in row
    if (P/'outer-next120-review.json').exists():
        batch=json.loads((P/'outer-next120-review.json').read_text())
        selection_path=P/'outer-next120-selection.json'
        selection=json.loads(selection_path.read_text())
        ranking_path=P/'outer-corpus-candidate-ranking.json'
        ranking=json.loads(ranking_path.read_text())
        assert hashlib.sha256(selection_path.read_bytes()).hexdigest()==batch['selection_sha256']
        assert hashlib.sha256(ranking_path.read_bytes()).hexdigest()==selection['ranking_sha256']==batch['ranking_sha256']
        assert selection['source_sha256']==batch['source_sha256']==l['corpus_sha256']
        excluded=set()
        for name,meta in selection['excluded_sets'].items():
            path=P/name
            assert hashlib.sha256(path.read_bytes()).hexdigest()==meta['sha256'],name
            rows=json.loads(path.read_text())['records'] if path.suffix=='.json' else meta['source_ids']
            keys={(r['source'],r['record_id']) for r in rows}
            assert len(keys)==meta['count'],name
            excluded.update(keys)
        eligible=[(i,r) for i,r in enumerate(ranking['candidates'],1) if (r['source'],r['record_id']) not in excluded]
        chosen=selection['records'];reviewed=batch['records']
        assert len(eligible)>=len(chosen)
        assert len(chosen)==len(reviewed)==selection['selected_count']==selection['reviewed_count']==batch['reviewed_count']==batch['selected_count']
        assert len({(r['source'],r['record_id']) for r in chosen})==len(chosen)
        for selected,note,(position,ranked) in zip(chosen,reviewed,eligible[:len(chosen)]):
            assert (selected['source'],selected['record_id'])==(ranked['source'],ranked['record_id'])
            assert selected['ranking_position']==position
            for field in ['source','record_id','title','abstract_sha256']:
                assert selected[field]==note[field],(note['record_id'],field)
            assert len(note['abstract_sha256'])==64 and note['review_depth']=='complete_saved_abstract'
            assert note['note'] and note['model'] and note['uncertainty']
            assert isinstance(note['requires_primary_followup'],bool)
            assert not {'abstractText','abstract'}&set(note)
    if (P/'outer-ranked-completion-review.json').exists():
        batch=json.loads((P/'outer-ranked-completion-review.json').read_text())
        selection_path=P/'outer-ranked-completion-selection.json'
        selection=json.loads(selection_path.read_text())
        assert hashlib.sha256(selection_path.read_bytes()).hexdigest()==batch['selection_sha256']
        ranking_path=P/'outer-corpus-candidate-ranking.json'
        ranking=json.loads(ranking_path.read_text())
        assert hashlib.sha256(ranking_path.read_bytes()).hexdigest()==batch['ranking_sha256']==selection['ranking_sha256']
        assert batch['source_sha256']==selection['source_sha256']==l['corpus_sha256']
        excluded=set()
        for name,digest in selection['excluded_ledgers'].items():
            path=P/name
            assert hashlib.sha256(path.read_bytes()).hexdigest()==digest,name
            excluded.update((r['source'],r['record_id']) for r in json.loads(path.read_text())['records'])
        supplement_path=P/'outer-next120-selection.json'
        assert hashlib.sha256(supplement_path.read_bytes()).hexdigest()==selection['supplement_exclusion_selection_sha256']
        supplement_exclusion=json.loads(supplement_path.read_text())['excluded_sets']['pubmed-supplement-review.md']
        excluded.update((r['source'],r['record_id']) for r in supplement_exclusion['source_ids'])
        eligible=[(i,r) for i,r in enumerate(ranking['candidates'],1) if (r['source'],r['record_id']) not in excluded]
        chosen=selection['records'];reviewed=batch['records']
        assert len(eligible)==len(chosen)==len(reviewed)==selection['selected_count']==selection['reviewed_count']==batch['selected_count']==batch['reviewed_count']
        assert len({(r['source'],r['record_id']) for r in chosen})==len(chosen)
        for selected,note,(position,ranked) in zip(chosen,reviewed,eligible):
            assert (selected['source'],selected['record_id'])==(ranked['source'],ranked['record_id'])
            assert selected['ranking_position']==position
            for field in ['source','record_id','title','abstract_sha256','ranking_position']:
                assert selected[field]==note[field],(note['record_id'],field)
            assert note['review_depth']=='complete_saved_abstract' and len(note['abstract_sha256'])==64
            assert note['note'] and note['model'] and note['uncertainty']
            assert isinstance(note['requires_primary_followup'],bool)
            assert not {'abstractText','abstract'}&set(note)
        assert sum(r['requires_primary_followup'] for r in reviewed)==batch['followup_count']
        historical=[json.loads((P/name).read_text())['records'] for name in ['outer-corpus-review.json','outer-next120-review.json']]
        combined=[r for rows in historical+[reviewed] for r in rows]
        combined_ids={(r['source'],r['record_id']) for r in combined}
        ranked_ids={(r['source'],r['record_id']) for r in ranking['candidates']}
        assert len(combined_ids)==len(combined) and ranked_ids<=combined_ids
        assert len(combined_ids)==433 and len(ranked_ids)==419
    if (P/'outer-unfiltered-review.json').exists():
        from verify_unfiltered import check
        check()
    if (P/'outer-unfiltered-next576-review.json').exists():
        from verify_next576 import check as check_next576
        check_next576()
    if (P/'broad-hub-summary.json').exists():
        broad=json.loads((P/'broad-hub-summary.json').read_text())
        broad_rows=read('broad-hub-screen.tsv');broad_mentions=read('broad-hub-mentions.tsv')
        by_id={r['hub_row_id']:r for r in broad_rows}
        assert len(broad_rows)==len(by_id)==broad['source_names_screened']==7540
        assert len(broad_mentions)==broad['lexical_mention_rows']==2148
        assert len({r['mention_id'] for r in broad_mentions})==len(broad_mentions)
        assert len({r['normalized_preferred_name'] for r in broad_rows})==broad['normalized_names_screened']
        mentions_by_id=Counter(r['hub_row_id'] for r in broad_mentions)
        assert set(mentions_by_id)<=set(by_id)
        assert len(mentions_by_id)==broad['matched_source_names']==544
        assert len(by_id)-len(mentions_by_id)==broad['unmatched_source_names']
        for r in broad_rows:
            assert int(r['literature_record_count'])==mentions_by_id[r['hub_row_id']]
            candidates=json.loads(r['name_candidate_ids'])+json.loads(r['structure_candidate_ids'])
            assert set(candidates)<=all_ids
            if r['crosswalk_status']=='NO_MATCH':assert not candidates
            if not mentions_by_id[r['hub_row_id']]:assert r['evidence_state']=='UNKNOWN_NO_PREFERRED_NAME_MATCH'
        for r in broad_mentions:
            assert r['evidence_state']=='PREFERRED_NAME_LEXICAL_MENTION_ONLY_UNREVIEWED'
            assert len(r['source_record_sha256'])==64
            assert json.loads(r['matched_fields_token_spans'])
            assert not {'abstractText','abstract','body_text'}&r.keys()
        assert sum(r['crosswalk_status']=='NO_MATCH' for r in broad_rows)==broad['no_catalog_match_source_names']==2713
        for name in ['screen_script','catalog','aliases','supplement','literature_match_script']:
            meta=broad['inputs'][name]
            assert hashlib.sha256((P/meta['file']).read_bytes()).hexdigest()==meta['sha256']
    decisions=json.loads((P/'decisions.json').read_text())['decisions'];joins=read('decision-alias-joins.tsv')
    assert len(decisions)==len(joins) and {r['decision_id'] for r in decisions}=={r['decision_id'] for r in joins}
    print(json.dumps({'result':'PASS','clinical_source_records':len(clinical),'main_catalog_entities':len(ids),'status_unverified_entities':len(supplement),'literature_screened_entities':len(all_ids),'curated_appraisals':len(decisions),'checks':'source ID/count reconciliation; complete interaction and literature coverage; UNKNOWN semantics; ambiguity regressions; original assay relation/units; triage coverage and compound identity guards; supplemental version/alias reconciliation; catalog-independent review depth; full saved-abstract and selected outer-corpus assessment reconciliation, including fixed-order exclusions and all419 ranked records across433 outer records; unfiltered384 plus next576 selection/review reconciliation, totaling1393 outer records; original missing143 depth reconciliation; Broad7540 preferred-name and2148 mention reconciliation with input hashes; all-target source preservation and correction mappings; complete decision joins'},indent=2))
if __name__=='__main__':main()
