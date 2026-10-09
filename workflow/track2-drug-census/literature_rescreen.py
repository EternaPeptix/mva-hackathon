#!/usr/bin/env python3
"""Versioned, inequality-preserving replay of the saved catalog alias screen.

Requires the original retained corpora and immutable alias policy. Writes a new
layer; does not overwrite historical matching, triage or manual-reading ledgers.
No retrieval match or automatic routing label is an efficacy conclusion.
"""
import argparse
import csv
import hashlib
import io
import json
import types
from collections import Counter, defaultdict
from pathlib import Path

from broad_hub_screen import clean, normalize, GENE
from catalog import write_json, write_tsv
import literature_triage as legacy_triage

HERE = Path(__file__).resolve().parent
INPUT_NAMES = (
    'catalog.tsv', 'drugcentral-status-unverified.tsv',
    'literature-alias-policy.tsv', 'literature-mentions.tsv',
    'literature-ambiguous-mentions.tsv', 'literature-coverage.json',
    'decoder-source-role-review.json',
)
SCRIPT_NAMES = ('literature_rescreen.py', 'broad_hub_screen.py',
                'literature_triage.py', 'literature_match.py', 'catalog.py')
# Copy the routing function's globals so even its intermediate normalization uses
# the corrected decoder, without changing historical modules or their bytes.
route_context = types.FunctionType(
    legacy_triage.context.__code__,
    dict(legacy_triage.__dict__, normalize=normalize), 'route_context')


def digest(data):
    return hashlib.sha256(data).hexdigest()


def read_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def assert_unchanged(inputs):
    for path, payload in inputs.items():
        if digest(path.read_bytes()) != digest(payload):
            raise ValueError(f'Input changed during replay: {path.name}')


def record_key(record):
    return record['source'], str(record.get('record_id', record.get('id')))


def row_key(row):
    return row['catalog_id'], row['source'], row['record_id']


def make_trie(policies):
    trie = {}
    for alias in policies:
        if sum(c.isalpha() for c in alias) < 2 or len(alias.replace(' ', '')) < 3:
            continue
        node = trie
        for token in alias.split():
            node = node.setdefault(token, {})
        node[''] = alias
    return trie


def matches(text, trie):
    tokens, found = text.split(), set()
    for start in range(len(tokens)):
        node = trie
        for token in tokens[start:]:
            if token not in node:
                break
            node = node[token]
            if '' in node:
                found.add(node[''])
    return found


def role_index(review, records, corpus_hash):
    """Reject stale annotations instead of silently reusing them on new text."""
    if review.get('status') != 'source_role_review_complete':
        raise ValueError('Source-role review has not been accepted')
    if review['corpus_sha256'] != corpus_hash:
        raise ValueError('Source-role review refers to another corpus')
    result = defaultdict(list)
    seen = set()
    for note in review['records']:
        key = record_key(note)
        if key in seen:
            raise ValueError(f'Duplicate source-role note: {key}')
        seen.add(key)
        if note.get('annotation_status') != 'source_revalidated' or note.get('review_depth') != 'complete_saved_abstract':
            raise ValueError(f'Source-role note is not a complete revalidated reading: {key}')
        record = records[key]
        raw = record.get('abstractText', '') or ''
        if digest(raw.encode()) != note['abstract_sha256']:
            raise ValueError(f"Source-role abstract changed: {record_key(note)}")
        for role in note['recovered_alias_roles']:
            for gid in role['catalog_ids']:
                prior = result[(gid, *key)]
                if any(r['normalized_alias'] == role['normalized_alias'] for _, r in prior):
                    raise ValueError(f'Duplicate source-role alias: {key}')
                prior.append((note, role))
    return result


def triage(row, record, roles):
    abstract = ' '.join(clean(record.get('abstractText', '')).split())
    role_rows = roles.get(row_key(row), [])
    if role_rows:
        aliases = set(row['matched_normalized_aliases'])
        if not all(role['normalized_alias'] in aliases for _, role in role_rows):
            raise ValueError(f'Reviewed alias missing in new match: {row_key(row)}')
        notes = [f"{role['normalized_alias']}: {role['role_evidence']}"
                 for _, role in role_rows]
        note = role_rows[0][0]
        category = '+'.join(sorted({r['role'] for _, r in role_rows}))
        supported_note = ' '.join(notes) + ' ' + note['uncertainty']
        followup = note['followup_reason']
        method = 'source_revalidated_recovered_alias_role'
        reviewed_aliases = sorted({r['normalized_alias'] for _, r in role_rows})
        requires_primary = note['requires_primary_followup']
    else:
        old_row = dict(row, matched_normalized_aliases=json.dumps(
            row['matched_normalized_aliases']),
            direct_gene_name_in_title_or_abstract=str(
                row['direct_gene_name_in_title_or_abstract']))
        category, supported_note, followup, method = route_context(old_row, abstract)
        if method.startswith('manual_'):
            method = 'inherited_' + method + '_not_revalidated_here'
        else:
            method = 'recomputed_' + method
        reviewed_aliases = []
        requires_primary = ''
    result = {k: row[k] for k in (
        'catalog_id', 'drug_name', 'source', 'record_id', 'pmid', 'doi', 'title',
        'matched_normalized_aliases', 'source_url', 'raw_abstract_sha256')}
    result.update(
        gene_literal_in_title_abstract=row['direct_gene_name_in_title_or_abstract'],
        abstract_available=bool(abstract), triage_category=category,
        triage_method=method, abstract_supported_note=supported_note,
        source_revalidated_aliases=reviewed_aliases, follow_up=followup,
        source_revalidated_requires_primary_followup=requires_primary,
        full_primary_verification='NOT_COMPLETED_BY_THIS_RESCREEN',
        actual_BUB1B_variant_muscle_rescue='NOT_ESTABLISHED_BY_THIS_ROW')
    routed = bool((legacy_triage.RESTORATIVE.search(abstract)
                   and legacy_triage.BENEFIT.search(abstract)
                   and not legacy_triage.CANCER.search(clean(record.get('title', ''))))
                  or 'unreviewed-primary lead' in followup)
    return result, routed


def compare(old_rows, new_rows):
    old = {row_key(r): set(json.loads(r['matched_normalized_aliases'])) for r in old_rows}
    new = {row_key(r): set(r['matched_normalized_aliases']) for r in new_rows}
    deltas = []
    for key in sorted(old.keys() | new.keys()):
        added, removed = new.get(key, set()) - old.get(key, set()), old.get(key, set()) - new.get(key, set())
        if added or removed:
            deltas.append(dict(catalog_id=key[0], source=key[1], record_id=key[2],
                               row_added=key not in old, row_removed=key not in new,
                               aliases_added=sorted(added), aliases_removed=sorted(removed)))
    return dict(old_rows=len(old), new_rows=len(new),
                added_rows=len(new.keys() - old.keys()), removed_rows=len(old.keys() - new.keys()),
                changed_rows=len(deltas), rows=deltas)


def run(corpus_path, ta_path, catalog_dir, out):
    out, catalog_dir = out.resolve(), catalog_dir.resolve()
    if out == catalog_dir:
        raise ValueError('Output must be a separate layer; historical tables are protected')
    if out.exists() and any(out.iterdir()):
        raise ValueError('Output must be empty; use a fresh directory for a replay')
    source_bytes = {catalog_dir/name: (catalog_dir/name).read_bytes() for name in INPUT_NAMES}
    script_bytes = {HERE/name: (HERE/name).read_bytes() for name in SCRIPT_NAMES}
    raw, ta_raw = corpus_path.read_bytes(), ta_path.read_bytes()
    protected_bytes = {**source_bytes, **script_bytes, corpus_path: raw, ta_path: ta_raw}
    def saved_json(name):
        return json.loads(source_bytes[catalog_dir/name])
    def saved_tsv(name):
        return list(csv.DictReader(io.StringIO(source_bytes[catalog_dir/name].decode()), delimiter='\t'))
    corpus, ta = json.loads(raw), json.loads(ta_raw)
    corpus_hash = digest(raw)
    historical = saved_json('literature-coverage.json')
    if corpus_hash != historical['corpus_sha256'] or digest(ta_raw) != historical['title_abstract_corpus_sha256']:
        raise ValueError('This replay requires the saved original corpus versions')
    records = {record_key(r): r for r in corpus['results']}
    if len(records) != len(corpus['results']):
        raise ValueError('Duplicate source/record IDs')
    ta_ids = {record_key(r) for r in ta['results']}
    catalog = saved_tsv('catalog.tsv')
    main_count = len(catalog)
    for r in saved_tsv('drugcentral-status-unverified.tsv'):
        catalog.append(dict(catalog_id='DC'+r['drugcentral_id'], preferred_name=r['name'],
                            status_summary='STATUS_UNVERIFIED_SUPPLEMENT'))
    names = {r['catalog_id']: r['preferred_name'] for r in catalog}
    if len(names) != len(catalog):
        raise ValueError('Duplicate catalog entity')
    policy_rows = saved_tsv('literature-alias-policy.tsv')
    policies = {r['normalized_alias']: dict(r, catalog_ids=json.loads(r['catalog_ids'])) for r in policy_rows}
    if len(policies) != len(policy_rows):
        raise ValueError('Duplicate saved alias policy')
    if any(p['policy'] not in {'specific_alias_comention_only', 'ambiguous_manual_review'} for p in policies.values()):
        raise ValueError('Unknown alias-policy class')
    if any(not set(p['catalog_ids']) <= names.keys() for p in policies.values()):
        raise ValueError('Alias refers to absent catalog ID')
    trie = make_trie(policies)
    review = saved_json('decoder-source-role-review.json')
    roles = role_index(review, records, corpus_hash)
    mentions, ambiguous, ledger = [], [], []
    for key, record in records.items():
        title, abstract = normalize(record.get('title', '')), normalize(record.get('abstractText', ''))
        direct = bool(GENE.search(title+' '+abstract))
        abstract_hash = digest((record.get('abstractText', '') or '').encode())
        source_hash = digest(json.dumps(record, ensure_ascii=False, sort_keys=True, separators=(',', ':')).encode())
        ledger.append(dict(source=key[0], record_id=key[1], raw_abstract_sha256=abstract_hash,
                           canonical_record_sha256=source_hash,
                           abstract_available=bool(record.get('abstractText')),
                           gene_literal_in_title_abstract=direct))
        hits = defaultdict(set)
        # Match within fields separately: a drug name cannot span title/abstract.
        for alias in matches(title, trie) | matches(abstract, trie):
            policy = policies[alias]
            for gid in policy['catalog_ids']:
                hits[policy['policy'], gid].add(alias)
        for (policy, gid), aliases in sorted(hits.items()):
            row = dict(catalog_id=gid, drug_name=names[gid], source=key[0], record_id=key[1],
                       pmid=record.get('pmid', ''), pmcid=record.get('pmcid', ''), doi=record.get('doi', ''),
                       title=' '.join(clean(record.get('title', '')).split()),
                       matched_normalized_aliases=sorted(aliases),
                       preferred_name_matched=normalize(names[gid]) in aliases,
                       direct_gene_name_in_title_or_abstract=direct,
                       member_of_gene_title_abstract_query=key in ta_ids,
                       abstract_available=bool(abstract), raw_abstract_sha256=abstract_hash,
                       evidence_state='CO_MENTION_ONLY_NOT_EFFICACY_OR_RESCUE',
                       source_url='https://europepmc.org/article/'+key[0]+'/'+key[1])
            (ambiguous if policy == 'ambiguous_manual_review' else mentions).append(row)
    current_keys = {row_key(r) for r in mentions}
    if not set(roles) <= current_keys:
        raise ValueError('Revalidated source-role rows lost in replay')
    counts = Counter(r['catalog_id'] for r in mentions)
    amb_counts = Counter(r['catalog_id'] for r in ambiguous)
    coverage = [dict(catalog_id=r['catalog_id'], preferred_name=r['preferred_name'],
                     status_summary=r['status_summary'], specific_alias_record_count=counts[r['catalog_id']],
                     ambiguous_alias_record_count=amb_counts[r['catalog_id']],
                     literature_evidence_state='CO_MENTION_REQUIRES_MANUAL_REVIEW' if counts[r['catalog_id']]
                     else 'UNKNOWN_NO_SPECIFIC_ALIAS_TITLE_ABSTRACT_MATCH') for r in catalog]
    triaged, restorative = [], []
    for row in mentions:
        result, route = triage(row, records[record_key(row)], roles)
        triaged.append(result)
        if route:
            restorative.append(result)
    delta = {kind: compare(saved_tsv(old), rows) for kind, old, rows in [
        ('specific', 'literature-mentions.tsv', mentions),
        ('ambiguous', 'literature-ambiguous-mentions.tsv', ambiguous)]}
    tables = {'mentions.tsv': mentions, 'ambiguous-mentions.tsv': ambiguous,
              'entity-coverage.tsv': coverage, 'source-records.tsv': ledger,
              'triage.tsv': triaged, 'restorative-routing.tsv': restorative,
              'source-role-followups.tsv': [r for r in triaged if r['source_revalidated_requires_primary_followup'] is True]}
    summary = dict(schema_version=1, status='CURRENT_CORRECTED_LEXICAL_LAYER',
                   corpus_sha256=corpus_hash, title_abstract_corpus_sha256=digest(ta_raw),
                   corpus_records=len(records), records_with_abstracts=sum(r['abstract_available'] for r in ledger),
                   records_with_gene_literal=sum(r['gene_literal_in_title_abstract'] for r in ledger),
                   catalog_entities=len(catalog), main_entities=main_count,
                   supplement_entities=len(catalog)-main_count, saved_normalized_aliases=len(policies),
                   alias_policy_counts=dict(Counter(p['policy'] for p in policies.values())),
                   specific_rows=len(mentions), ambiguous_rows=len(ambiguous),
                   matched_entities=len(counts), unknown_no_specific_match_entities=len(catalog)-len(counts),
                   direct_gene_specific_rows=sum(r['direct_gene_name_in_title_or_abstract'] for r in mentions),
                   triage_methods=dict(Counter(r['triage_method'] for r in triaged)),
                   triage_categories=dict(Counter(r['triage_category'] for r in triaged)),
                   recovered_role_rows=len(roles), recovered_role_source_records=len(review['records']),
                   restorative_routing_rows=len(restorative),
                   source_role_primary_followup_rows=len(tables['source-role-followups.tsv']),
                   manual_reading_coverage_increment=0,
                   limits='Lexical retrieval and routing only. Recovered roles use the separately source-bound abstract audit; other manual annotations are inherited, not revalidated. Restorative routing is a keyword shortlist, not the complete primary-followup queue; source-role followups are retained separately. Ambiguous hits are unadjudicated. Source record counts are not independent studies. No match is UNKNOWN, not inefficacy. This does not repair rejected manual readings or establish worldwide inventory completeness.',
                   input_sha256={name: digest(source_bytes[catalog_dir/name]) for name in INPUT_NAMES},
                   script_sha256={name: digest(script_bytes[HERE/name]) for name in SCRIPT_NAMES})
    assert_unchanged(protected_bytes)
    out.mkdir(parents=True, exist_ok=True)
    for name, rows in tables.items():
        write_tsv(out/name, rows, list(rows[0]) if rows else ['catalog_id'])
    write_json(out/'delta.json', delta)
    summary['output_sha256'] = {name: digest((out/name).read_bytes()) for name in [*tables, 'delta.json']}
    write_json(out/'summary.json', summary)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--corpus', type=Path, required=True)
    parser.add_argument('--title-abstract-corpus', type=Path, required=True)
    parser.add_argument('--catalog-dir', type=Path, default=HERE)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    result = run(args.corpus, args.title_abstract_corpus, args.catalog_dir, args.out)
    print(json.dumps({key: result[key] for key in (
        'status', 'corpus_records', 'catalog_entities', 'specific_rows', 'ambiguous_rows',
        'matched_entities', 'unknown_no_specific_match_entities', 'records_with_gene_literal',
        'recovered_role_rows', 'manual_reading_coverage_increment')}, indent=2))


if __name__ == '__main__':
    main()
