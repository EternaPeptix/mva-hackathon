#!/usr/bin/env python3
"""Reconcile corrected screen; optionally replay original sources independently of saved outputs.

Checks identities and transformations, not scientific reading or clinical validity.
"""
import argparse
import json
import tempfile
from collections import Counter
from pathlib import Path
from literature_rescreen import HERE, digest, read_tsv, row_key, record_key, run


def check(catalog_dir=HERE, out=None, corpus=None, title_abstract_corpus=None):
    out = out or catalog_dir/'literature-v2'
    summary = json.loads((out/'summary.json').read_text())
    for name, expected in summary['input_sha256'].items():
        assert digest((catalog_dir/name).read_bytes()) == expected, name
    for name, expected in summary['script_sha256'].items():
        assert digest((HERE/name).read_bytes()) == expected, name
    for name, expected in summary['output_sha256'].items():
        assert digest((out/name).read_bytes()) == expected, name
    mentions, ambiguous, coverage, triage, records = [read_tsv(out/name) for name in (
        'mentions.tsv', 'ambiguous-mentions.tsv', 'entity-coverage.tsv', 'triage.tsv', 'source-records.tsv')]
    ids = {r['catalog_id'] for r in read_tsv(catalog_dir/'catalog.tsv')}
    ids |= {'DC'+r['drugcentral_id'] for r in read_tsv(catalog_dir/'drugcentral-status-unverified.tsv')}
    assert len(coverage) == len(ids) == summary['catalog_entities']
    assert {r['catalog_id'] for r in coverage} == ids
    assert len(mentions) == len({row_key(r) for r in mentions}) == summary['specific_rows']
    assert len(ambiguous) == len({row_key(r) for r in ambiguous}) == summary['ambiguous_rows']
    counts, amb_counts = Counter(r['catalog_id'] for r in mentions), Counter(r['catalog_id'] for r in ambiguous)
    for r in coverage:
        assert int(r['specific_alias_record_count']) == counts[r['catalog_id']]
        assert int(r['ambiguous_alias_record_count']) == amb_counts[r['catalog_id']]
        expected = 'CO_MENTION_REQUIRES_MANUAL_REVIEW' if counts[r['catalog_id']] else 'UNKNOWN_NO_SPECIFIC_ALIAS_TITLE_ABSTRACT_MATCH'
        assert r['literature_evidence_state'] == expected
    assert len(counts) == summary['matched_entities']
    assert len(ids)-len(counts) == summary['unknown_no_specific_match_entities']
    record_map = {record_key(r): r for r in records}
    assert len(records) == len(record_map) == summary['corpus_records']
    assert sum(r['abstract_available'] == 'True' for r in records) == summary['records_with_abstracts']
    assert sum(r['gene_literal_in_title_abstract'] == 'True' for r in records) == summary['records_with_gene_literal']
    policy = {r['normalized_alias']: r for r in read_tsv(catalog_dir/'literature-alias-policy.tsv')}
    for rows, policy_type in [(mentions, 'specific_alias_comention_only'), (ambiguous, 'ambiguous_manual_review')]:
        for r in rows:
            assert r['catalog_id'] in ids and record_key(r) in record_map
            for alias in json.loads(r['matched_normalized_aliases']):
                assert policy[alias]['policy'] == policy_type
                assert r['catalog_id'] in json.loads(policy[alias]['catalog_ids'])
            source = record_map[record_key(r)]
            assert r['raw_abstract_sha256'] == source['raw_abstract_sha256']
            assert r['direct_gene_name_in_title_or_abstract'] == source['gene_literal_in_title_abstract']
            assert r['evidence_state'] == 'CO_MENTION_ONLY_NOT_EFFICACY_OR_RESCUE'
    assert len(triage) == len(mentions) and {row_key(r) for r in triage} == {row_key(r) for r in mentions}
    assert dict(Counter(r['triage_method'] for r in triage)) == summary['triage_methods']
    assert dict(Counter(r['triage_category'] for r in triage)) == summary['triage_categories']
    assert all(r['actual_BUB1B_variant_muscle_rescue'] == 'NOT_ESTABLISHED_BY_THIS_ROW' for r in triage)
    reviewed = {row_key(r): r for r in triage if r['triage_method'] == 'source_revalidated_recovered_alias_role'}
    review = json.loads((catalog_dir/'decoder-source-role-review.json').read_text())
    expected = {}
    for note in review['records']:
        for role in note['recovered_alias_roles']:
            for gid in role['catalog_ids']:
                key = (gid, *record_key(note))
                expected.setdefault(key, set()).add(role['normalized_alias'])
                assert reviewed[key]['raw_abstract_sha256'] == note['abstract_sha256']
                assert reviewed[key]['source_revalidated_requires_primary_followup'] == str(note['requires_primary_followup'])
    assert set(reviewed) == set(expected) and len(reviewed) == summary['recovered_role_rows']
    for key, aliases in expected.items():
        assert set(json.loads(reviewed[key]['source_revalidated_aliases'])) == aliases
    assert summary['manual_reading_coverage_increment'] == 0
    restorative = read_tsv(out/'restorative-routing.tsv')
    assert len(restorative) == summary['restorative_routing_rows']
    assert {row_key(r) for r in restorative} <= {row_key(r) for r in triage}
    followups = read_tsv(out/'source-role-followups.tsv')
    assert len(followups) == summary['source_role_primary_followup_rows']
    assert {row_key(r) for r in followups} == {key for key, r in reviewed.items()
        if r['source_revalidated_requires_primary_followup'] == 'True'}
    if bool(corpus) != bool(title_abstract_corpus):
        raise ValueError('Supply both corpora to replay')
    if corpus:
        with tempfile.TemporaryDirectory(prefix='mva-rescreen-') as directory:
            replay = run(corpus, title_abstract_corpus, catalog_dir, Path(directory))
            assert replay == summary, 'Source replay differs from saved outputs'
    result = dict(result='PASS', source_replay=bool(corpus), specific_rows=len(mentions),
                  ambiguous_rows=len(ambiguous), entities=len(ids), recovered_role_rows=len(reviewed),
                  limits='Identity, policy and transformation reconciliation; not proof of scientific conclusions.')
    print(json.dumps(result, indent=2))
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--catalog-dir', type=Path, default=HERE)
    parser.add_argument('--out', type=Path)
    parser.add_argument('--corpus', type=Path)
    parser.add_argument('--title-abstract-corpus', type=Path)
    args = parser.parse_args()
    check(args.catalog_dir, args.out, args.corpus, args.title_abstract_corpus)
