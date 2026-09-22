#!/usr/bin/env python3
"""Reconcile the saved census; fail on missing rows/IDs or unsafe alias-policy regressions."""
import csv,json
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
    print(json.dumps({'result':'PASS','clinical_source_records':len(clinical),'main_catalog_entities':len(ids),'status_unverified_entities':len(supplement),'literature_screened_entities':len(all_ids),'checks':'source ID/count reconciliation; complete interaction and literature coverage; UNKNOWN semantics; ambiguity regressions; original assay relation/units; triage coverage and compound identity guards'},indent=2))
if __name__=='__main__':main()
