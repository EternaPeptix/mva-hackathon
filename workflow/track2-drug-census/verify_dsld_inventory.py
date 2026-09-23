#!/usr/bin/env python3
"""Verify the committed NIH discovery snapshot, not interpretation or efficacy."""
import csv,hashlib,json
from pathlib import Path

def check(root=None):
    root=Path(root or Path(__file__).resolve().parent);p=root/'dsld-ingredients'
    def js(name):return json.loads((p/name).read_text())
    def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
    def rows(name):return list(csv.DictReader((p/name).open(),delimiter='\t'))
    s=js('summary.json');counts=s['counts']
    for name,h in s['output_sha256'].items():assert Path(name).name==name and sha(p/name)==h,name
    inputs=s['input_sha256']
    for name in ['build_dsld_inventory.py','broad_hub_screen.py','catalog.tsv','aliases.tsv','drugcentral-status-unverified.tsv','decisions.json']:
        hits=[h for k,h in inputs.items() if Path(k).name==name];assert len(hits)==1 and sha(root/name)==hits[0],name
    assert sha(root/'food-substances/substances.tsv')==next(h for k,h in inputs.items() if k.endswith('/food-substances/substances.tsv'))
    assert sha(p/'source-retrieval-manifest.json')==next(h for k,h in inputs.items() if Path(k).name=='retrieval-manifest.json')
    groups=rows('groups.tsv');cells=rows('alias-cells.tsv');mentions=rows('mentions.tsv');records=rows('source-records.tsv')
    assert len(groups)==len({x['group_id'] for x in groups})==counts['source_groups']==6467
    assert len(cells)==counts['all_source_alias_cells']==85255
    assert len(mentions)==counts['group_alias_record_mentions']==41496
    assert len(records)==counts['corpus_source_records']==7362
    assert all(x['evidence_state']=='UNKNOWN_NOT_ASSESSED_BY_INVENTORY' for x in groups)
    assert all(x['interpretation']=='LEXICAL_DISCOVERY_CONTEXT_UNREVIEWED' for x in mentions)
    ids={x['group_id'] for x in groups};assert all(x['group_id'] in ids for x in cells+mentions)
    assert sum(v is None for x in groups for v in json.loads(x['category_as_reported']))==592
    assert all(json.loads(x['nutrient_info_as_reported'])==[] for x in groups)
    assert len({x['group_name_as_reported'] for x in groups})==6080
    assert counts['automatic_manual_abstract_reads']==counts['automatic_primary_reads']==counts['automatic_new_appraisals']==0
    assert len(json.loads((root/'decisions.json').read_text())['decisions'])==counts['target_inventory']['appraisals']==170
    assert len(rows('name-crossrefs.tsv'))==counts['name_crossref_rows']==6803
    parent=js('parent-reconciliation.json');assert parent['current_summary_sha256']==sha(p/'summary.json') and parent['current_builder_sha256']==sha(root/'build_dsld_inventory.py')
    assert parent['historical_ledger_count']==168 and parent['current_ledger_count']==170 and parent['new_manual_reading_or_efficacy_credit']==0
    overlay=js('first-role-review.json');assert overlay['source_snapshot']['mentions_sha256']==sha(p/'mentions.tsv')
    assert overlay['reading_manifest_sha256']==sha(p/'first-role-reading-manifest.json')
    chosen=[x for record in overlay['records'] for x in record['nih_unflagged_rows']]
    keys={tuple(sorted(x.items())) for x in mentions}
    assert len(chosen)==len({tuple(sorted(x.items())) for x in chosen})==6
    assert all(tuple(sorted(x.items())) in keys for x in chosen)
    assert len(overlay['records'])==4 and overlay['credit']['new_query_role_annotations']==6
    assert all(overlay['credit'][k]==0 for k in ['new_unique_abstract_reads','new_outer_corpus_reads','new_primary_reads','new_drug_appraisals'])
    return {'result':'PASS','scope':'Committed hashes, source-table relationships and unknown states; full replay needs external source cache. No scientific-role validation.','source_groups':6467,'source_cells':85255,'lexical_mentions':41496,'added_reading_or_efficacy_credit':0}

if __name__=='__main__':print(json.dumps(check(),indent=2))
