#!/usr/bin/env python3
"""Verify the CPDat v4.0/v4.1 Formulation and Article overlay snapshots."""
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OVERLAYS = (
    ROOT / 'cpdat-v4.0',
    ROOT / 'cpdat-v4.0' / 'articles',
    ROOT / 'cpdat-v4.1',
    ROOT / 'cpdat-v4.1' / 'articles',
)
REFINEMENTS = ROOT / 'cpdat-v4.1' / 'ledger-refinements-20261008.json'
REFINEMENTS_SHA256 = '3e55c137b70fd396574e7e0d45092fd7831afe23ee60abf0df76d92805a1650f'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def tsv(path):
    with Path(path).open(newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def decision_snapshot(ledger_bytes, recorded_count, recorded_sha256, refinements=None):
    """Authenticate a dated ledger, allowing appends and exact audited refinements."""
    ledger = json.loads(ledger_bytes)
    decisions = ledger['decisions']
    assert 0 <= recorded_count <= len(decisions), 'Recorded decision snapshot is truncated'
    if recorded_count == len(decisions):
        snapshot_bytes = ledger_bytes
    else:
        # Preserve metadata and key order. This reconstruction is accepted only
        # when it reproduces the exact bytes authenticated by the recorded hash.
        ledger['decisions'] = decisions[:recorded_count]
        snapshot_bytes = (json.dumps(ledger, indent=2, ensure_ascii=False) + '\n').encode('utf-8')
    if hashlib.sha256(snapshot_bytes).hexdigest() == recorded_sha256:
        return decisions[:recorded_count]
    assert refinements is not None, 'Recorded decision snapshot changed'
    assert refinements['recorded_decision_count'] == recorded_count, 'Refinement snapshot count differs'
    assert refinements['recorded_ledger_sha256'] == recorded_sha256, 'Refinement snapshot hash differs'
    ledger['decisions'] = decisions[:recorded_count]
    for key, values in refinements['metadata_refinements'].items():
        assert key != 'decisions', 'Decision objects require separately audited refinements'
        assert key in ledger and ledger[key] in (values['original'], values['accepted_current']), 'Unrecorded ledger metadata refinement'
        ledger[key] = values['original']
    seen = set()
    for item in refinements['decision_refinements']:
        index = item['index']
        assert 0 <= index < recorded_count and index not in seen, 'Invalid refinement index'
        seen.add(index)
        current = ledger['decisions'][index]
        original = item['original_object']
        assert current['decision_id'] == original['decision_id'] == item['decision_id'], 'Refined decision identity changed'
        if current != original:
            current_hash = hashlib.sha256(json.dumps(current, ensure_ascii=False, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
            assert current_hash == item['accepted_current_object_sha256'], 'Unrecorded decision refinement'
        ledger['decisions'][index] = original
    snapshot_bytes = (json.dumps(ledger, indent=2, ensure_ascii=False) + '\n').encode('utf-8')
    assert hashlib.sha256(snapshot_bytes).hexdigest() == recorded_sha256, 'Recorded decision snapshot changed'
    return ledger['decisions']


def check_overlay(overlay, ledger_bytes, refinements=None):
    manifest = json.loads((overlay / 'source-manifest.json').read_text())
    summary = json.loads((overlay / 'screen-summary.json').read_text())
    for item in manifest['outputs']:
        path = overlay / item['path']
        assert path.stat().st_size == item['bytes'], (overlay.name, item['path'])
        assert digest(path) == item['sha256'], (overlay.name, item['path'])
    assert manifest['screening_script_sha256'] == digest(ROOT / 'build_cpdat_bub1b_screen.py')
    assert summary['decision_ledger_sha256'] == manifest['decision_ledger_sha256']
    assert summary['decision_ledger_count'] == manifest['decision_ledger_count']
    decisions = decision_snapshot(ledger_bytes, summary['decision_ledger_count'], manifest['decision_ledger_sha256'], refinements)
    decision_ids = {d['decision_id'] for d in decisions}
    assert summary['cpdat_archive_sha256'] == manifest['archive_sha256']
    assert summary['gene_corpus_sha256'] == manifest['target_corpus_sha256']
    assert summary['puc_kind'] == manifest['puc_kind']

    chemicals = tsv(overlay / manifest['output_table'])
    categories = tsv(overlay / 'product-category-summary.tsv')
    mentions_path = overlay / 'bub1b-exact-name-co-mentions.tsv'
    mentions = tsv(mentions_path)
    assert all(not line.endswith(b'\t') for line in mentions_path.read_bytes().splitlines())
    assert len(chemicals) == summary['unique_selected_kind_dtxsids']
    assert len({r['dtxsid'] for r in chemicals}) == len(chemicals)
    assert summary['selected_kind_composition_rows_with_dtxsid'] == sum(int(r['composition_rows']) for r in chemicals)
    assert sum(int(r['composition_rows_with_dtxsid']) for r in categories) == summary['selected_kind_composition_rows_with_dtxsid']
    assert len(categories) == summary['selected_kind_categories']
    assert all(r['puc_kind'] == manifest['puc_kind'] for r in categories)
    assert len(mentions) == summary['exact_unique_name_matches_written']
    source_ids = {(r['source'], r['id']) for r in mentions}
    mention_dtxsids = {sid for row in mentions for sid in row['dtxsids'].split(';')}
    assert len(source_ids) == summary['source_id_records_with_review_class_exact_match']
    assert len(mention_dtxsids) == summary['review_class_exact_match_dtxsids']
    assert all(row['candidate_status'] == 'REVIEW_EXACT_NAME' for row in mentions)
    for row in mentions:
        if row['current_decision_ids'] != 'NO_EXACT_LEDGER_ALIAS':
            assert set(row['current_decision_ids'].split(';')) <= decision_ids

    version = overlay.parents[0].name if overlay.name == 'articles' else overlay.name
    if manifest['puc_kind'] == 'Formulation' and version == 'cpdat-v4.0':
        bpa = [r for r in chemicals if r['dtxsid'] == 'DTXSID7020182']
        assert len(bpa) == 1
        assert int(bpa[0]['composition_rows']) == int(bpa[0]['distinct_product_keys']) == 42
        assert set(bpa[0]['product_categories'].split(';')) == {'Food and drug', 'Home maintenance'}
    if manifest['puc_kind'] == 'Formulation' and version == 'cpdat-v4.1':
        ibuprofen = [r for r in chemicals if r['dtxsid'] == 'DTXSID5020732']
        naproxen_sodium = [r for r in chemicals if r['dtxsid'] == 'DTXSID7045576']
        assert len(ibuprofen) == len(naproxen_sodium) == 1
        assert int(ibuprofen[0]['composition_rows']) == int(ibuprofen[0]['distinct_product_keys']) == 2
        assert int(naproxen_sodium[0]['composition_rows']) == 6
        assert int(naproxen_sodium[0]['distinct_product_keys']) == 5
        assert summary['product_identity_basis'].startswith('Data Document ID + Product Name')
        assert any(r['current_decision_ids'] == 'new-ibuprofen-bubr1-downregulation' for r in mentions if 'DTXSID5020732' in r['dtxsids'])
    if manifest['puc_kind'] == 'Article':
        ozone_rows = [r for r in mentions if r['id'] == '24336422' and 'DTXSID0021098' in r['dtxsids']]
        assert ozone_rows
        assert all('new-ozone-bubr1b-downregulation' in r['current_decision_ids'].split(';') for r in ozone_rows)

    return {
        'overlay': str(overlay.relative_to(ROOT)),
        'result': 'PASS',
        'puc_kind': manifest['puc_kind'],
        'cpdat_release': manifest['source'],
        'screened_decisions': len(decisions),
        'appended_decisions_since_screen': len(json.loads(ledger_bytes)['decisions']) - len(decisions),
        'composition_rows_with_dtxsid': summary['selected_kind_composition_rows_with_dtxsid'],
        'dtxsids': len(chemicals),
        'review_class_exact_name_rows': len(mentions),
        'review_class_source_id_records': len(source_ids),
        'review_class_dtxsids': len(mention_dtxsids),
    }


def check():
    assert digest(REFINEMENTS) == REFINEMENTS_SHA256, 'Dated refinement audit changed'
    refinements = json.loads(REFINEMENTS.read_text())
    ledger_bytes = (ROOT / 'decisions.json').read_bytes()
    decisions = json.loads(ledger_bytes)['decisions']
    decision_ids = {d['decision_id'] for d in decisions}
    assert {
        'new-bisphenol-a-bubr1-exposure',
        'new-naproxen-bubr1-downregulation',
        'new-ibuprofen-bubr1-downregulation',
        'new-mnbp-bubr1b-downregulation',
        'new-ozone-bubr1b-downregulation',
        'new-quercetin-aglycone-bubr1b-screen',
    } <= decision_ids
    return {'result': 'PASS', 'current_decisions': len(decisions), 'dated_refinement_audit': str(REFINEMENTS.relative_to(ROOT)), 'overlays': [check_overlay(o, ledger_bytes, refinements) for o in OVERLAYS]}


if __name__ == '__main__':
    print(json.dumps(check(), indent=2))
