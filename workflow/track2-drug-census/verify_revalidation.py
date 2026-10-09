#!/usr/bin/env python3
"""Reconcile partial, source-bound re-readings without validating their science.

Hash equality identifies the text associated with a declaration of reading. It
does not prove reading fidelity or efficacy. Unvalidated rows stay pending.
"""
import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def digest(data):
    return hashlib.sha256(data).hexdigest()


def key(row):
    return row['source'], row.get('record_id', row.get('id'))


def check(root=ROOT, corpus=None, gene_corpus=None):
    def read(name):
        return json.loads((root / name).read_text())
    selection = read('outer-revalidation-selection.json')
    review = read('outer-revalidation-review.json')
    assert review['status'] == 'partial_source_revalidation'
    assert review['selection_sha256'] == digest((root / 'outer-revalidation-selection.json').read_bytes())
    assert review['source_sha256'] == selection['source_sha256']
    prior = []
    for name, expected in selection['prior_review_sha256'].items():
        assert digest((root / name).read_bytes()) == expected, name
        prior.extend(read(name)['records'])
    previous = {key(r) for r in prior}
    gene = {key(r) for r in read('abstract-review-notes.json')['records']}
    assert len(prior) == len(previous) == review['previous_outer_ledger_count'] == 2161
    assert len(gene) == 1460
    selected = selection['records']
    lookup = {key(r): r for r in selected}
    assert len(lookup) == len(selected) == selection['selected_count'] == review['selected_count'] == 3598
    assert not set(lookup) & (previous | gene)
    order = lambda r: digest((selection['seed'] + '\0' + r['source'] + '\0' + key(r)[1]).encode())
    assert [r['order_key'] for r in selected] == sorted(order(r) for r in selected)
    for index, record in enumerate(selected, 1729):
        assert record['order_key'] == order(record)
        assert record['unfiltered_index'] == index
    rows = review['records']
    assert len(rows) == len({key(r) for r in rows}) == review['reviewed_count']
    assert review['pending_count'] == len(selected) - len(rows)
    assert review['combined_outer_ledger_count'] == len(previous) + len(rows)
    assert review['original_outer_records_without_abstract'] == selection['remaining_missing_abstracts'] == 143
    assert review['followup_count'] == sum(r['requires_primary_followup'] for r in rows)
    assert review['source_content_kinds'] == dict(Counter(r['saved_content_kind'] for r in rows))
    for row in rows:
        assert all(row[k] == v for k, v in lookup[key(row)].items()), key(row)
        assert row['annotation_status'] == 'source_revalidated'
        assert row['quality_validation'] == 'complete_literal_saved_string_re_read'
        assert row['review_depth'] == 'complete_saved_abstract'
        assert row['readers'] and row['read_date'] and row['saved_content_kind']
        assert row['note'] and row['uncertainty'] and row['model'] and row['source_abstract_quality']
        assert isinstance(row['requires_primary_followup'], bool)
        assert not {'abstractText', 'abstract', 'body_text'} & row.keys()
        assert all(a.get('name') and a.get('role') for a in row['administered_agents'])
    replayed = False
    if (corpus is None) != (gene_corpus is None):
        raise ValueError('Supply both original retrieval files')
    if corpus is not None:
        assert digest(corpus.read_bytes()) == selection['source_sha256']
        assert digest(gene_corpus.read_bytes()) == selection['gene_source_sha256']
        raw = json.loads(corpus.read_text())['results']
        assert {key(r) for r in json.loads(gene_corpus.read_text())['results']} == gene
        remaining = [r for r in raw if key(r) not in previous | gene]
        assert len(remaining) == selection['remaining_outer_records'] == 3741
        available = sorted([r for r in remaining if r.get('abstractText')], key=order)
        assert len(available) == len(selected)
        for actual, saved in zip(available, selected):
            assert key(actual) == key(saved) and actual['title'] == saved['title']
            assert digest(actual['abstractText'].encode()) == saved['abstract_sha256']
        replayed = True
    return dict(result='PASS', checked='identity and declared-depth reconciliation only',
                additional_literal_source_re_readings=len(rows),
                pending_source_revalidation=review['pending_count'],
                combined_outer_ledger_records=review['combined_outer_ledger_count'],
                source_content_kinds=review['source_content_kinds'],
                original_empty_abstract_records=143, source_selection_replayed=replayed)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    parser.add_argument('--corpus', type=Path)
    parser.add_argument('--gene-corpus', type=Path)
    args = parser.parse_args()
    print(json.dumps(check(args.root, args.corpus, args.gene_corpus), indent=2))
