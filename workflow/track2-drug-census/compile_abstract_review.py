#!/usr/bin/env python3
"""Validate authored abstract assessments against their exact source; never infer review."""
import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

from catalog import write_json, write_tsv

ROOT = Path(__file__).resolve().parent
FIELDS = [
    'source', 'record_id', 'pmid', 'doi', 'title', 'source_url',
    'abstract_available', 'abstract_sha256', 'review_depth', 'reviewer',
    'intervention_type', 'administered_agents', 'model', 'outcome_direction',
    'note', 'requires_primary_followup', 'publication_types',
]
TYPES = {'administered_drug', 'genetic', 'biomarker_association',
         'damaging_exposure', 'basic_mechanism', 'clinical_case', 'review', 'other'}


def digest(data):
    return hashlib.sha256(data).hexdigest()


def compile_review(corpus, notes):
    raw = corpus.read_bytes()
    if digest(raw) != notes['corpus_sha256']:
        raise ValueError('Source snapshot changed; historical reading cannot be inherited')
    records = json.loads(raw)['results']
    source = {(r['source'], r['id']): r for r in records}
    if len(source) != len(records):
        raise ValueError('Duplicate source identities')
    reviewed = {}
    for row in notes['records']:
        key = (row['source'], row['record_id'])
        if key not in source or key in reviewed:
            raise ValueError(f'Unknown or duplicate assessment: {key}')
        original = source[key]
        if row['title'] != original['title']:
            raise ValueError(f'Source title differs: {key}')
        abstract = original.get('abstractText') or ''
        if row['review_depth'] == 'complete_saved_abstract':
            if not abstract or row['abstract_sha256'] != digest(abstract.encode()):
                raise ValueError(f'Abstract identity differs: {key}')
            if row['intervention_type'] not in TYPES:
                raise ValueError(f'Unknown assessment type: {key}')
            if not isinstance(row['administered_agents'], list):
                raise ValueError(f'Agent list must be explicit: {key}')
            if not row.get('model') or not row.get('outcome_direction'):
                raise ValueError(f'Model and outcome scope must be explicit: {key}')
        elif row['review_depth'] in {'complete_accessible_body', 'title_metadata_only'}:
            if abstract:
                raise ValueError(f'Missing-abstract assessment attached to abstract: {key}')
        else:
            raise ValueError(f'Unsupported manual review depth: {key}')
        if not row.get('note') or not isinstance(row['requires_primary_followup'], bool):
            raise ValueError(f'Missing authored assessment: {key}')
        reviewed[key] = row
    out = []
    for original in records:
        key = (original['source'], original['id'])
        abstract = original.get('abstractText') or ''
        row = dict(reviewed.get(key, {}))
        row.update(source=key[0], record_id=key[1], pmid=original.get('pmid', ''),
                   doi=original.get('doi', ''), title=original['title'],
                   abstract_available=bool(abstract),
                   publication_types=original.get('pubTypeList', {}).get('pubType', []))
        row.setdefault('source_url', f'https://europepmc.org/article/{key[0]}/{key[1]}')
        row.setdefault('abstract_sha256', digest(abstract.encode()) if abstract else '')
        row.setdefault('review_depth', 'not_reviewed')
        row.setdefault('intervention_type', 'UNKNOWN')
        row.setdefault('administered_agents', [])
        row.setdefault('note', 'No manual assessment; efficacy and full-text intervention role remain UNKNOWN.')
        row.setdefault('requires_primary_followup', True)
        out.append(row)
    depth = Counter(r['review_depth'] for r in out)
    available = sum(r['abstract_available'] for r in out)
    summary = {
        'corpus_sha256': digest(raw), 'source_records': len(out),
        'available_saved_abstracts': available,
        'missing_saved_abstracts': len(out) - available,
        'review_depth_counts': dict(sorted(depth.items())),
        'available_saved_abstracts_read': depth['complete_saved_abstract'],
        'available_saved_abstracts_unread': available - depth['complete_saved_abstract'],
        'all_available_saved_abstracts_read': available == depth['complete_saved_abstract'],
        'intervention_types_among_read_abstracts': dict(sorted(Counter(
            r['intervention_type'] for r in out if r['review_depth'] == 'complete_saved_abstract').items())),
        'followup_flags': sum(r['requires_primary_followup'] for r in out),
        'scope': 'Complete saved strings, including source truncations; not all original abstracts, full texts, independent studies or drugs. Notes are authored assessments, not generated classifications. Missing evidence is UNKNOWN.',
        'source_order_preserved': True,
    }
    return out, summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--corpus', type=Path, required=True)
    p.add_argument('--notes', type=Path, default=ROOT / 'abstract-review-notes.json')
    p.add_argument('--out', type=Path, default=ROOT)
    args = p.parse_args()
    rows, summary = compile_review(args.corpus, json.loads(args.notes.read_text()))
    args.out.mkdir(parents=True, exist_ok=True)
    write_tsv(args.out / 'gene-abstract-review.tsv', rows, FIELDS)
    write_json(args.out / 'gene-abstract-review-summary.json', summary)
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
