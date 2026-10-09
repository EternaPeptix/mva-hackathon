#!/usr/bin/env python3
"""Check review coverage; optionally reproduce selection from original retrievals.

No check here establishes scientific validity or proves that a human read a source.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def sha(data):
    return hashlib.sha256(data).hexdigest()


def key(row):
    return row['source'], row.get('record_id', row.get('id'))


def check(corpus=None, gene_corpus=None):
    def read(name):
        return json.loads((ROOT / name).read_text())

    selection = read('outer-unfiltered-selection.json')
    review = read('outer-unfiltered-review.json')
    missing = read('outer-missing-review.json')
    assert sha((ROOT / 'outer-unfiltered-selection.json').read_bytes()) == review['selection_sha256']
    assert review['source_sha256'] == selection['source_sha256']
    selected, reviewed = selection['records'], review['records']
    selected_keys = {key(r) for r in selected}
    assert len(selected) == len(reviewed) == len(selected_keys) == 384
    assert selection['selected_count'] == review['selected_count'] == review['reviewed_count'] == 384
    previous = []
    for name in ('outer-corpus-review.json', 'outer-next120-review.json',
                 'outer-ranked-completion-review.json'):
        previous.extend(read(name)['records'])
    previous_keys = {key(r) for r in previous}
    gene_keys = {key(r) for r in read('abstract-review-notes.json')['records']}
    assert len(previous) == len(previous_keys) == 433
    assert not selected_keys & (previous_keys | gene_keys)
    for index, (s, r) in enumerate(zip(selected, reviewed), 1):
        assert all(s[field] == r[field] for field in s), key(s)
        assert r['unfiltered_index'] == index
        expected_order = sha((selection['seed'] + '\0' + r['source'] + '\0' + r['record_id']).encode())
        assert r['order_key'] == expected_order
        assert r['review_depth'] == 'complete_saved_abstract'
        assert len(r['abstract_sha256']) == 64
        assert r['note'] and r['uncertainty'] and r['model']
        assert isinstance(r['requires_primary_followup'], bool)
        assert not {'abstractText', 'abstract', 'body_text'} & r.keys()
    assert [r['order_key'] for r in selected] == sorted(r['order_key'] for r in selected)
    assert sum(r['requires_primary_followup'] for r in reviewed) == review['followup_count'] == 65
    assert selection['remaining_outer_records'] == 5469
    assert selection['remaining_available_abstracts'] == 5326
    assert selection['remaining_missing_abstracts'] == 143

    missing_keys = {key(r) for r in missing['records']}
    assert len(missing['records']) == len(missing_keys) == missing['record_count'] == 143
    assert not missing_keys & (selected_keys | previous_keys | gene_keys)
    counts = Counter(r['review_depth'] for r in missing['records'])
    assert dict(counts) == missing['reading_counts']
    assert counts == Counter(complete_recovered_abstract=18, complete_available_body=53,
                             retrieval_only_unread_long_body=11, unavailable=56, no_body=5)
    assert counts['complete_recovered_abstract'] + counts['complete_available_body'] == missing['completed_readings'] == 71
    stubs = [r for r in missing['records'] if r['review_depth'] == 'complete_available_body'
             and r['intervention_type'] == 'availability_stub_no_scientific_content']
    assert len(stubs) == 14
    assert sum(r['intervention_type'] == 'keyword_index_no_efficacy_content'
               for r in missing['records']) == 1
    assert sum(r['available_scientific_content'] for r in missing['records']) == 53
    for r in missing['records']:
        assert len(r['source_record_sha256']) == 64
        assert r['note'] and r['uncertainty'] and r['identity_status']
        assert not {'abstractText', 'abstract', 'body_text'} & r.keys()
        if r['review_depth'] == 'unavailable':
            assert r['source_sha256'] is None
        else:
            assert len(r['source_sha256']) == 64

    if (corpus is None) != (gene_corpus is None):
        raise ValueError('Supply both --corpus and --gene-corpus for source replay')
    if corpus is not None:
        assert sha(corpus.read_bytes()) == selection['source_sha256']
        raw = json.loads(corpus.read_text())['results']
        inner = json.loads(gene_corpus.read_text())['results']
        assert len(raw) == 7362 and len(inner) == 1460
        assert {key(r) for r in inner} == gene_keys
        remaining = [r for r in raw if key(r) not in gene_keys | previous_keys]
        abstracts = [r for r in remaining if r.get('abstractText')]
        assert len(remaining) == 5469 and len(abstracts) == 5326
        assert {key(r) for r in remaining if not r.get('abstractText')} == missing_keys
        ordered = sorted(abstracts, key=lambda r: sha((selection['seed'] + '\0' +
                          r['source'] + '\0' + r['id']).encode()))[:384]
        for raw_row, reviewed_row in zip(ordered, reviewed):
            assert key(raw_row) == key(reviewed_row)
            assert raw_row['title'] == reviewed_row['title']
            assert sha(raw_row['abstractText'].encode()) == reviewed_row['abstract_sha256']
    return {'result': 'PASS', 'new_complete_abstracts': 384,
            'combined_outer_selected_records': 817, 'original_saved_abstracts_unread': 4942,
            'missing_records_accounted': 143, 'recovered_readings': 71,
            'availability_only_bodies': 14, 'source_selection_replayed': corpus is not None}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--corpus', type=Path)
    parser.add_argument('--gene-corpus', type=Path)
    args = parser.parse_args()
    print(json.dumps(check(args.corpus, args.gene_corpus), indent=2))
