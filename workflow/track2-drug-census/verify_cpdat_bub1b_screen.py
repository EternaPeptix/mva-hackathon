#!/usr/bin/env python3
"""Verify the checked-in, dated CPDat v4.0 discovery overlay and its hashes."""
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OVERLAY = ROOT / 'cpdat-v4.0'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def tsv(path):
    with Path(path).open(newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def check():
    manifest = json.loads((OVERLAY / 'source-manifest.json').read_text())
    summary = json.loads((OVERLAY / 'screen-summary.json').read_text())
    for item in manifest['outputs']:
        path = OVERLAY / item['path']
        assert path.stat().st_size == item['bytes'], item['path']
        assert digest(path) == item['sha256'], item['path']
    assert manifest['screening_script_sha256'] == digest(ROOT / 'build_cpdat_bub1b_screen.py')
    assert manifest['decision_ledger_sha256'] == digest(ROOT / 'decisions.json')
    assert summary['decision_ledger_sha256'] == manifest['decision_ledger_sha256']
    assert summary['cpdat_archive_sha256'] == manifest['archive_sha256']
    assert summary['gene_corpus_sha256'] == manifest['target_corpus_sha256']

    chemicals = tsv(OVERLAY / 'consumer-formulations.tsv')
    categories = tsv(OVERLAY / 'product-category-summary.tsv')
    mentions = tsv(OVERLAY / 'bub1b-exact-name-co-mentions.tsv')
    assert all(not line.endswith(b'\t') for line in (OVERLAY / 'bub1b-exact-name-co-mentions.tsv').read_bytes().splitlines())
    assert len(chemicals) == summary['unique_formulation_dtxsids'] == 7009
    assert len({r['dtxsid'] for r in chemicals}) == len(chemicals)
    assert sum(int(r['formulation_rows']) for r in categories) == summary['formulation_composition_rows_with_dtxsid']
    assert len(mentions) == summary['exact_unique_name_matches_written'] == 1150
    assert len({(r['source'], r['id']) for r in mentions}) == summary['source_id_records_with_review_class_exact_match'] == 725
    assert len({sid for row in mentions for sid in row['dtxsids'].split(';')}) == summary['review_class_exact_match_dtxsids'] == 233
    assert all(row['candidate_status'] == 'REVIEW_EXACT_NAME' for row in mentions)

    bpa = [r for r in chemicals if r['dtxsid'] == 'DTXSID7020182']
    assert len(bpa) == 1
    assert int(bpa[0]['composition_rows']) == int(bpa[0]['unique_product_ids']) == 42
    assert set(bpa[0]['product_categories'].split(';')) == {'Food and drug', 'Home maintenance'}
    decisions = json.loads((ROOT / 'decisions.json').read_text())['decisions']
    assert len(decisions) == 284
    assert sum(d['decision_id'] == 'new-bisphenol-a-bubr1-exposure' for d in decisions) == 1
    return {
        'result': 'PASS',
        'cpdat_release': manifest['source'],
        'formulation_dtxsids': len(chemicals),
        'review_class_exact_name_rows': len(mentions),
        'review_class_source_id_records': len({(r['source'], r['id']) for r in mentions}),
        'current_decisions': len(decisions),
        'scope': 'Hash and count consistency only; does not validate chemical identity, co-mention context, exposure, or efficacy.',
    }


if __name__ == '__main__':
    print(json.dumps(check(), indent=2))
