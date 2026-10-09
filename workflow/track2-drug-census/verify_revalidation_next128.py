#!/usr/bin/env python3
"""Reconcile the additive next-128 reading ledger and its external evidence.

This validates identities, hashes, exact selection, preserved input annotations,
and recorded full-string batch evidence. It does not certify that a reader's
scientific interpretation is correct or that any intervention is effective.
Raw abstracts remain in the external source cache, never in the review ledger.
"""
import argparse
import copy
import importlib.util
import json
import hashlib
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SOURCE_SHA = 'e5b8b017bb68af38818f1cc8c1bee47adb16afe0367981104a6b17886d2954d9'
GENE_SHA = '1df717bce779851e8f8f02654cc476de320cf368ad3cbc10dad1499b5a809e90'
PARENT_SHA = 'a9a12a95285b7456be572811a0a03891748535abe263706184ef517069938fb2'
ACCEPTED186_SHA = '0bdc47535b70046b0d0a115da26339ec845aacb4d38827a9cec65e88e0dd81bb'
SELECTION_NAME = 'outer-revalidation-next128-selection.json'
REVIEW_NAME = 'outer-revalidation-next128-review.json'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def key(row):
    return row['source'], row.get('record_id', row.get('id'))


def source_id(row):
    return ':'.join(key(row))


def normalize_record(row):
    row = dict(row)
    if 'originalnote' in row:
        require('note' not in row or row['note'] == row['originalnote'],
                'Conflicting note and originalnote')
        row['note'] = row.pop('originalnote')
    return row


def publication_counts(rows):
    links = Counter()
    with_links = errata = versions = preprints = linked_preprints = 0
    for row in rows:
        flags = row['publication_metadata_flags']
        comments = flags.get('commentCorrectionList', {}).get('commentCorrection', [])
        types = [entry.get('type', '') for entry in comments]
        links.update(types)
        with_links += bool(comments)
        errata += any('erratum' in item.lower() or 'correction' in item.lower()
                      for item in types)
        linked_preprints += any(item.lower() == 'preprint in' for item in types)
        versions += bool(flags.get('versionList', {}).get('version', []))
        preprints += (row['source'] == 'PPR' or any(
            item.lower() == 'preprint' for item in flags.get('pubTypeList', {}).get('pubType', [])))
    return dict(records_with_comment_correction_links=with_links,
                comment_correction_link_types=dict(sorted(links.items())),
                records_with_erratum_or_correction_links=errata,
                records_with_version_lists=versions,
                preprint_records=preprints,
                records_with_linked_preprints=linked_preprints)


def evidence_path(base, name):
    # Evidence names are local basenames, not arbitrary filesystem locations.
    require(isinstance(name, str) and Path(name).name == name and name not in ('', '.', '..'),
            'Invalid evidence filename: ' + str(name))
    path = base / name
    require(path.is_file(), 'Missing external evidence: ' + str(path))
    return path


def normalized_batches(part, payload, evidence_root):
    """Check each literal log, and return a canonical manifest plus its strings."""
    entries = payload.get('batch_read_evidence')
    require(isinstance(entries, list) and len(entries) == 16, f'{part}: expected 16 batches')
    canonical, emitted = [], []
    for batch_number, item in enumerate(entries, 1):
        require(item.get('batch') == batch_number, f'{part}: batch order')
        basename = Path(item['log_file']).name
        expected_name = f'part-{part}-batch-{batch_number:02d}.' + ('jsonl' if part == 'a' else 'log')
        require(basename == expected_name, f'{part}: unexpected batch file')
        path = evidence_path(evidence_root, basename)
        raw = path.read_bytes()
        require(digest(raw) == item['log_sha256'], f'{basename}: log hash mismatch')
        lines = [json.loads(line) for line in raw.decode().splitlines()]
        require(len(lines) == 6, f'{basename}: expected header, four records and end marker')
        header, rows, end = lines[0], lines[1:-1], lines[-1]
        require(header == {'corpus_sha256': SOURCE_SHA, 'requested_records': 4},
                f'{basename}: invalid source header')
        require(end.get('end_of_batch') is True and end.get('emitted_records') == 4,
                f'{basename}: missing/invalid end marker')
        require(item.get('end_of_batch') is True, f'{basename}: missing end declaration')
        observed = item.get('end_of_batch_observed', item.get('reader_observed_end_marker'))
        untruncated = item.get('tool_output_untruncated')
        if untruncated is None and 'tool_output_truncated' in item:
            untruncated = item['tool_output_truncated'] is False
        require(observed is True and untruncated is True,
                f'{basename}: absent complete-output reading declaration')
        ids = [row['source_id'] for row in rows]
        require(ids == item['source_ids'] and len(set(ids)) == 4,
                f'{basename}: record identities differ')
        total_chars = 0
        for row in rows:
            value = row.get('abstractText')
            require(isinstance(value, str) and bool(value), f'{basename}: absent literal string')
            require(len(value) == row['abstract_characters'], f'{basename}: character mismatch')
            require(digest(value.encode()) == row['abstract_sha256'], f'{basename}: abstract hash mismatch')
            total_chars += len(value)
        require(end.get('total_abstract_characters') == total_chars == item['total_abstract_characters'],
                f'{basename}: character total mismatch')
        require(item['requested_records'] == item['emitted_records'] == 4,
                f'{basename}: declaration count mismatch')
        if 'abstract_sha256s' in item:
            require(item['abstract_sha256s'] == [r['abstract_sha256'] for r in rows],
                    f'{basename}: declared abstract hashes differ')
        entry = dict(part=part, batch=batch_number, log_file=basename,
                     log_sha256=digest(raw), source_ids=ids,
                     abstract_sha256s=[r['abstract_sha256'] for r in rows],
                     requested_records=4, emitted_records=4,
                     total_abstract_characters=total_chars, end_of_batch=True,
                     end_of_batch_observed=True, tool_output_untruncated=True,
                     readers=payload['readers'])
        if item.get('tool_output_chunk_id'):
            entry['tool_output_chunk_id'] = item['tool_output_chunk_id']
        canonical.append(entry)
        emitted.extend(rows)
    return canonical, emitted


def check(root=ROOT, corpus=None, gene_corpus=None, selection_path=None,
          review_path=None, evidence_root=None):
    root = Path(root)
    selection_path = Path(selection_path or root / SELECTION_NAME)
    review_path = Path(review_path or root / REVIEW_NAME)
    evidence_root = Path(evidence_root or selection_path.parent)
    require((corpus is None) == (gene_corpus is None), 'Supply both original retrieval files')
    corpus = Path(corpus) if corpus is not None else None
    gene_corpus = Path(gene_corpus) if gene_corpus is not None else None

    parent_path = root / 'outer-revalidation-selection.json'
    accepted_path = root / 'outer-revalidation-review.json'
    require(digest(parent_path.read_bytes()) == PARENT_SHA, 'Original 3598 selection changed')
    require(digest(accepted_path.read_bytes()) == ACCEPTED186_SHA, 'Accepted 186 ledger changed')
    parent = read_json(parent_path)
    accepted = read_json(accepted_path)
    preserved = {name: digest((root / name).read_bytes())
                 for name in parent['prior_review_sha256']}
    preserved.update({'outer-revalidation-selection.json': PARENT_SHA,
                      'outer-revalidation-review.json': ACCEPTED186_SHA})
    verifier = root / 'verify_revalidation.py'
    require(verifier.is_file(), 'Missing historical verifier')
    spec = importlib.util.spec_from_file_location('historical_revalidation', verifier)
    historical = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(historical)
    previous_result = historical.check(root, corpus, gene_corpus)
    require(accepted['reviewed_count'] == len(accepted['records']) == 186,
            'Exactly 186 previously accepted records required')
    excluded = {key(r) for r in accepted['records']}
    require(len(excluded) == 186, 'Duplicate accepted record')
    pending = [r for r in parent['records'] if key(r) not in excluded]
    require(len(parent['records']) == 3598 and len(pending) == 3412, 'Original queue totals differ')

    selection = read_json(selection_path)
    review = read_json(review_path)
    require(selection['parent_selection_sha256'] == PARENT_SHA, 'Parent selection hash mismatch')
    require(selection['excluded_review_sha256'] == ACCEPTED186_SHA, 'Accepted exclusion hash mismatch')
    require(selection['source_sha256'] == SOURCE_SHA and selection['gene_source_sha256'] == GENE_SHA,
            'Selection corpus hashes differ')
    require(selection['records'] == pending[:128], 'Selection is not the exact first 128 pending rows')
    require(selection['selected_count'] == 128 and selection['pending_before'] == 3412
            and selection['prior_outer_ledger_records'] == 2347, 'Selection totals differ')
    selected = selection['records']
    require(review['status'] == 'partial_source_revalidation', 'Unexpected review status')
    require(review['selection_sha256'] == digest(selection_path.read_bytes()), 'Selection file hash mismatch')
    require(review['source_sha256'] == SOURCE_SHA and review['gene_source_sha256'] == GENE_SHA,
            'Review corpus hashes differ')
    require(review['parent_selection_sha256'] == PARENT_SHA and review['excluded_review_sha256'] == ACCEPTED186_SHA,
            'Review prior-input hashes differ')
    rows = review['records']
    require(len(rows) == len({key(r) for r in rows}) == review['reviewed_count'] == 128,
            'Expected 128 distinct reviewed source/ID records')
    require([key(r) for r in rows] == [key(r) for r in selected], 'Review selection order differs')
    require(not ({key(r) for r in rows} & excluded), 'Previously accepted record counted again')
    for row, chosen in zip(rows, selected):
        label = source_id(row)
        require(all(row.get(k) == v for k, v in chosen.items()), f'{label}: selection field changed')
        require(row['annotation_status'] == 'source_revalidated'
                and row['quality_validation'] == 'complete_literal_saved_string_re_read'
                and row['review_depth'] == 'complete_saved_abstract', f'{label}: reading depth absent')
        require(isinstance(row['readers'], list) and row['readers']
                and all(isinstance(r, str) and r.strip() for r in row['readers']), f'{label}: reader absent')
        for field in ('read_date', 'saved_content_kind', 'note', 'uncertainty', 'model',
                      'source_abstract_quality', 'intervention_type', 'outcome_direction'):
            require(isinstance(row.get(field), str) and row[field].strip(), f'{label}: empty {field}')
        require(type(row['requires_primary_followup']) is bool, f'{label}: followup not boolean')
        require(isinstance(row['followup_reason'], str) and
                (not row['requires_primary_followup'] or row['followup_reason'].strip()),
                f'{label}: followup reason absent')
        require(not {'abstractText', 'abstract', 'body_text', 'originalnote'} & row.keys(),
                f'{label}: raw text or unnormalized note key in ledger')
        require(isinstance(row['administered_agents'], list)
                and all(isinstance(a, dict) and a.get('name') and a.get('role')
                        for a in row['administered_agents']), f'{label}: invalid agent roles')
        require(isinstance(row.get('publication_metadata_flags'), dict), f'{label}: metadata absent')

    inputs = review['review_inputs']
    require(isinstance(inputs, list) and len(inputs) == 2, 'Expected two input readers')
    all_input_rows, all_batches, all_strings = [], [], []
    for index, part in enumerate(('a', 'b')):
        item = inputs[index]
        require(item['part'] == part and item['file'] == f'part-{part}-review.json', 'Input part order differs')
        path = evidence_path(evidence_root, item['file'])
        require(digest(path.read_bytes()) == item['sha256'], f'{part}: input review hash mismatch')
        payload = read_json(path)
        part_selection_name = f'part-{part}-selection.json'
        require(item['selection_file'] == part_selection_name, f'{part}: assignment file differs')
        part_selection_path = evidence_path(evidence_root, part_selection_name)
        require(digest(part_selection_path.read_bytes()) == item['selection_sha256'] == payload['selection_sha256'],
                f'{part}: assignment hash mismatch')
        part_selection = read_json(part_selection_path)
        part_rows = [normalize_record(r) for r in payload['records']]
        require(part_selection['records'] == selected[index * 64:(index + 1) * 64], f'{part}: assignment differs')
        require(len(part_rows) == item['record_count'] == 64, f'{part}: expected 64 rows')
        require(payload['source_sha256'] == SOURCE_SHA and payload['parent_selection_sha256'] == PARENT_SHA,
                f'{part}: source hashes differ')
        require(payload['readers'] == item['readers'] and all(r['readers'] == payload['readers'] for r in part_rows),
                f'{part}: reader attribution changed')
        require(payload.get('read_method', payload.get('reading_method', '')).strip(), f'{part}: reading method absent')
        batches, strings = normalized_batches(part, payload, evidence_root)
        require([r['source_id'] for r in strings] == [source_id(r) for r in part_rows],
                f'{part}: actual batch order differs')
        for row, literal in zip(part_rows, strings):
            require(row['abstract_sha256'] == literal['abstract_sha256'] and row['title'] == literal['title'],
                    f'{source_id(row)}: batch source differs')
            if 'read_evidence' in row:
                evidence = row['read_evidence']
                require(evidence['end_of_batch_observed'] is True, 'Per-record end marker absent')
                batch = batches[evidence['batch'] - 1]
                require(evidence['batch_log_file'] == batch['log_file']
                        and evidence['batch_log_sha256'] == batch['log_sha256']
                        and source_id(row) in batch['source_ids']
                        and evidence['abstract_characters'] == literal['abstract_characters'],
                        f'{source_id(row)}: per-record reading evidence differs')
        all_input_rows.extend(part_rows)
        all_batches.extend(batches)
        all_strings.extend(strings)
    revised_rows = copy.deepcopy(all_input_rows)
    by_key = {key(row): row for row in revised_rows}
    revisions = review.get('annotation_revisions', [])
    require(isinstance(revisions, list), 'Invalid annotation revision ledger')
    revised_fields = set()
    revised_ids = set()
    audit_strings = {}
    for evidence in review.get('annotation_revision_evidence', []):
        path = evidence_path(evidence_root, evidence['log_file'])
        raw = path.read_bytes()
        require(digest(raw) == evidence['log_sha256'], 'Annotation audit log hash mismatch')
        lines = [json.loads(line) for line in raw.decode().splitlines()]
        require(len(lines) >= 3, 'Incomplete annotation audit log')
        header, literal_rows, end = lines[0], lines[1:-1], lines[-1]
        require(1 <= len(literal_rows) <= 8 and header['corpus_sha256'] == SOURCE_SHA
                and header['requested_records'] == len(literal_rows), 'Invalid annotation audit header')
        require(end.get('end_of_batch') is True and end['emitted_records'] == len(literal_rows)
                and evidence['end_of_batch_observed'] is True
                and evidence['tool_output_untruncated'] is True, 'Incomplete annotation audit reading declaration')
        require([r['source_id'] for r in literal_rows] == evidence['source_ids'], 'Annotation audit identities differ')
        require(end['total_abstract_characters'] == sum(len(r['abstractText']) for r in literal_rows),
                'Annotation audit character total differs')
        for literal in literal_rows:
            require(literal['abstract_characters'] == len(literal['abstractText'])
                    and digest(literal['abstractText'].encode()) == literal['abstract_sha256'],
                    'Annotation audit literal hash differs')
            require(literal['source_id'] not in audit_strings, 'Duplicate annotation audit record')
            audit_strings[literal['source_id']] = (literal, evidence)
    allowed_revision_fields = {'model', 'note', 'uncertainty', 'followup_reason',
                               'outcome_direction', 'intervention_type', 'source_abstract_quality'}
    for revision in revisions:
        identity = key(revision)
        field = revision['field']
        require(identity in by_key and field in allowed_revision_fields,
                'Revision changes identity/provenance or unknown annotation field')
        require((identity, field) not in revised_fields, 'Duplicate revision of same field')
        revised_fields.add((identity, field))
        revised_ids.add(source_id(revision))
        row = by_key[identity]
        require(row[field] == revision['original_value'], 'Revision does not preserve original annotation value')
        require(isinstance(revision['revised_value'], str) and revision['revised_value'].strip()
                and revision['revised_value'] != revision['original_value'], 'Empty or no-op annotation revision')
        require(revision['reason'].strip() and revision['reviewers'] and revision['read_date'],
                'Annotation revision lacks rationale/attribution')
        literal, evidence = audit_strings[source_id(revision)]
        require(literal['abstract_sha256'] == row['abstract_sha256'] == revision['abstract_sha256']
                and literal['title'] == row['title'], 'Revision audit source differs')
        require(revision['evidence_log_file'] == evidence['log_file']
                and revision['evidence_log_sha256'] == evidence['log_sha256'], 'Revision evidence link differs')
        row[field] = revision['revised_value']
    require(set(audit_strings) == revised_ids, 'Revision audit record set differs')
    require(review['annotation_revision_count'] == len(revisions)
            and review['annotation_revision_record_count'] == len(revised_ids), 'Annotation revision count differs')
    require(rows == revised_rows, 'Consolidation changed annotations outside the explicit revision ledger')
    require(review['batch_read_evidence'] == all_batches, 'Consolidated batch manifest differs')
    require(review['batch_count'] == len(all_batches) == 32, 'Batch total differs')
    require(review['total_abstract_characters'] == sum(r['abstract_characters'] for r in all_strings),
            'Read character total differs')

    require(review['previous_outer_ledger_count'] == accepted['combined_outer_ledger_count'] == 2347,
            'Previous ledger count differs')
    require(review['combined_outer_ledger_count'] == 2347 + len(rows) == 2475, 'Combined count differs')
    require(review['pending_count'] == len(pending) - len(rows) == 3284, 'Pending count differs')
    require(review['selected_count'] == 128 and review['original_outer_records_without_abstract'] == 143,
            'Selection/empty-source totals differ')
    require(review['followup_count'] == sum(r['requires_primary_followup'] for r in rows), 'Followup count differs')
    require(review['source_content_kinds'] == dict(Counter(r['saved_content_kind'] for r in rows)), 'Content kinds differ')
    metadata = publication_counts(rows)
    require(review['publication_metadata_counts'] == metadata, 'Publication counts differ')
    require(review['correction_count'] == metadata['records_with_erratum_or_correction_links']
            and review['version_count'] == metadata['records_with_version_lists'], 'Correction/version count differs')

    if corpus is not None:
        # historical.check has already replayed full selection and both file hashes.
        raw = read_json(corpus)['results']
        lookup = {key(r): r for r in raw}
        require(len(lookup) == len(raw), 'Duplicate source corpus ID')
        for row, literal in zip(rows, all_strings):
            actual = lookup[key(row)]
            require(actual.get('abstractText') == literal['abstractText'], f'{source_id(row)}: literal corpus mismatch')
            require(actual['title'] == row['title'], f'{source_id(row)}: title mismatch')
            flags = row['publication_metadata_flags']
            require(all(actual.get(k) == v for k, v in flags.items()), f'{source_id(row)}: publication metadata mismatch')
            for flag in ('pubTypeList', 'publicationStatus', 'commentCorrectionList', 'versionList', 'versionNumber'):
                require(flag not in actual or flags.get(flag) == actual[flag], f'{source_id(row)}: omitted publication flag {flag}')
            require(all(actual.get(k) == v for k, v in row.get('source_links', {}).items()),
                    f'{source_id(row)}: source link mismatch')
    require(all(digest((root / name).read_bytes()) == sha for name, sha in preserved.items()),
            'Historical input bytes changed during verification')
    return dict(result='PASS', checked='Provenance, exact selection, preserved annotations and declared literal-reading evidence only',
                scientific_reading_fidelity_validated=False, primary_paper_qualification_validated=False,
                additional_literal_source_re_readings=128, prior_accepted_revalidation_records=186,
                prior186_bytes_preserved=True, pending_source_revalidation=3284,
                combined_outer_ledger_records=2475, followup_count=review['followup_count'],
                source_content_kinds=review['source_content_kinds'], publication_metadata_counts=metadata,
                batch_count=len(all_batches), total_abstract_characters=review['total_abstract_characters'],
                external_batch_evidence_verified=True, annotation_revision_count=len(revisions),
                annotation_revision_record_count=len(revised_ids), source_selection_replayed=corpus is not None,
                previous_verifier_result=previous_result['result'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    parser.add_argument('--selection', type=Path)
    parser.add_argument('--review', type=Path)
    parser.add_argument('--evidence-root', type=Path)
    parser.add_argument('--corpus', type=Path)
    parser.add_argument('--gene-corpus', type=Path)
    args = parser.parse_args()
    print(json.dumps(check(args.root, args.corpus, args.gene_corpus, args.selection,
                           args.review, args.evidence_root), indent=2))
