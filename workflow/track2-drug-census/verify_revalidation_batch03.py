#!/usr/bin/env python3
"""Reconcile additive batch03, with structural or explicit full-evidence checks.

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
import sys
sys.dont_write_bytecode = True
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SOURCE_SHA = 'e5b8b017bb68af38818f1cc8c1bee47adb16afe0367981104a6b17886d2954d9'
GENE_SHA = '1df717bce779851e8f8f02654cc476de320cf368ad3cbc10dad1499b5a809e90'
PARENT_SHA = 'a9a12a95285b7456be572811a0a03891748535abe263706184ef517069938fb2'
ACCEPTED186_SHA = '0bdc47535b70046b0d0a115da26339ec845aacb4d38827a9cec65e88e0dd81bb'
NEXT128_SHA = '26102e1a6e5bd6cd5204e4d38d92e53b5757709cb4cfda8a839fe7c9769fda7c'
NEXT128_SELECTION_SHA = '6e110ed66bcf85fe29cd1bcd9f72b9983bd61c012544327722e48d6fdcc43abf'
BATCH02_SHA = '222ef9edb446dae6c87782eb5a672a904ba45342b5866e49fe3a854f8444ef03'
BATCH02_SELECTION_SHA = '028a557db4d83781bf12597cd23f03310d8e2be45b96fe5e7c6cc61155bf9702'
SELECTION_NAME = 'outer-revalidation-batch03-selection.json'
REVIEW_NAME = 'outer-revalidation-batch03-review.json'


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


def publication_counts(rows):
    links = Counter()
    with_links = errata = versions = preprints = linked_preprints = retracted = retraction_links = 0
    for row in rows:
        flags = row['publication_metadata_flags']
        comments = flags.get('commentCorrectionList', {}).get('commentCorrection', [])
        types = [entry.get('type', '') for entry in comments]
        links.update(types)
        with_links += bool(comments)
        errata += any('erratum' in item.lower() or 'correction' in item.lower()
                      for item in types)
        linked_preprints += any(item.lower() == 'preprint in' for item in types)
        retracted += 'Retracted Publication' in flags.get('pubTypeList', {}).get('pubType', [])
        retraction_links += any('retraction' in t.lower() for t in types)
        versions += bool(flags.get('versionList', {}).get('version', []))
        preprints += (row['source'] == 'PPR' or any(
            item.lower() == 'preprint' for item in flags.get('pubTypeList', {}).get('pubType', [])))
    return dict(records_with_comment_correction_links=with_links,
                comment_correction_link_types=dict(sorted(links.items())),
                records_with_erratum_or_correction_links=errata,
                records_with_version_lists=versions,
                preprint_records=preprints,
                records_with_linked_preprints=linked_preprints,
                retracted_publication_records=retracted, records_with_retraction_links=retraction_links)


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
        expected_name = f'part-{part}-batch-{batch_number:02d}.jsonl'
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


def object_digest(value):
    return digest(json.dumps(value, ensure_ascii=False, sort_keys=True,
                             separators=(',', ':')).encode())


def valid_hash(value):
    return isinstance(value, str) and len(value) == 64 and all(c in '0123456789abcdef' for c in value)


def original_rows_from_revisions(review, full=False, evidence_root=None):
    """Undo only explicit, source-linked annotation corrections for input comparison."""
    originals = copy.deepcopy(review['records'])
    lookup = {key(r): r for r in originals}
    revisions = review.get('annotation_revisions', [])
    evidence = review.get('annotation_revision_evidence', [])
    require(isinstance(revisions, list) and isinstance(evidence, list), 'Invalid revision declarations')
    evidence_by_file = {}
    literal_lookup = {}
    for item in evidence:
        require(item['log_file'] not in evidence_by_file and valid_hash(item['log_sha256']), 'Duplicate/invalid audit evidence')
        require(item['end_of_batch_observed'] is True and item['tool_output_untruncated'] is True
                and item['readers'] and item['additional_coverage_credit'] == 0, 'Invalid audit reading declaration')
        require(item['source_ids'] and len(item['source_ids']) == len(set(item['source_ids']))
                and all(s in {source_id(r) for r in originals} for s in item['source_ids']), 'Audit IDs outside assigned batch')
        evidence_by_file[item['log_file']] = item
        if full:
            path = evidence_path(evidence_root, item['log_file'])
            raw = path.read_bytes()
            require(digest(raw) == item['log_sha256'], 'Audit log hash differs')
            lines = [json.loads(line) for line in raw.decode().splitlines()]
            require(len(lines) >= 3, 'Incomplete audit log')
            header, bodies, end = lines[0], lines[1:-1], lines[-1]
            require(1 <= len(bodies) <= 8 and header['corpus_sha256'] == SOURCE_SHA
                    and header['requested_records'] == len(bodies), 'Audit source header differs')
            require(end.get('end_of_batch') is True and end['emitted_records'] == len(bodies)
                    and end['total_abstract_characters'] == sum(len(r['abstractText']) for r in bodies),
                    'Audit end marker/count differs')
            require([r['source_id'] for r in bodies] == item['source_ids'], 'Audit log identities differ')
            for r in bodies:
                require(len(r['abstractText']) == r['abstract_characters']
                        and digest(r['abstractText'].encode()) == r['abstract_sha256'], 'Audit source string hash differs')
                literal_lookup[(item['log_file'], r['source_id'])] = r
    fields_seen, records_seen, used_evidence = set(), set(), set()
    allowed = {'note', 'model', 'outcome_direction', 'uncertainty', 'followup_reason', 'source_abstract_quality'}
    for change in revisions:
        identity, field = key(change), change['field']
        require(identity in lookup and field in allowed and (identity, field) not in fields_seen,
                'Invalid/duplicate annotation revision')
        fields_seen.add((identity, field)); records_seen.add(identity)
        row = lookup[identity]
        require(change['abstract_sha256'] == row['abstract_sha256'], 'Revision source hash differs')
        require(row[field] == change['revised_value'] and change['revised_value'] != change['original_value']
                and isinstance(change['original_value'], str) and change['original_value'].strip(), 'Revision old/new values differ')
        require(change['reason'].strip() and change['reviewers'] and change['read_date'], 'Revision attribution absent')
        item = evidence_by_file[change['evidence_log_file']]
        used_evidence.add(change['evidence_log_file'])
        require(item['log_sha256'] == change['evidence_log_sha256'] and source_id(row) in item['source_ids'],
                'Revision evidence link differs')
        if full:
            literal = literal_lookup[(item['log_file'], source_id(row))]
            require(literal['abstract_sha256'] == row['abstract_sha256'] and literal['title'] == row['title'],
                    'Revision literal source differs')
        row[field] = change['original_value']
    require(used_evidence == set(evidence_by_file), 'Unreferenced revision evidence')
    require(review['annotation_revision_count'] == len(revisions)
            and review['annotation_revision_record_count'] == len(records_seen), 'Revision count differs')
    return originals


def check(root=ROOT, corpus=None, gene_corpus=None, selection_path=None,
          review_path=None, evidence_root=None):
    """Repo-only checks by default; full replay needs three explicit source paths."""
    root = Path(root)
    selection_path = Path(selection_path or root / SELECTION_NAME)
    review_path = Path(review_path or root / REVIEW_NAME)
    provided = [corpus is not None, gene_corpus is not None, evidence_root is not None]
    require(all(provided) or not any(provided),
            'Full replay requires --evidence-root, --corpus and --gene-corpus together')
    full = all(provided)
    if full:
        corpus, gene_corpus, evidence_root = map(Path, (corpus, gene_corpus, evidence_root))
        require(digest(corpus.read_bytes()) == SOURCE_SHA, 'Original source corpus hash differs')
        require(digest(gene_corpus.read_bytes()) == GENE_SHA, 'Original gene corpus hash differs')
    prior_files = {
        'outer-revalidation-selection.json': PARENT_SHA,
        'outer-revalidation-review.json': ACCEPTED186_SHA,
        'outer-revalidation-next128-selection.json': NEXT128_SELECTION_SHA,
        'outer-revalidation-next128-review.json': NEXT128_SHA,
        'outer-revalidation-batch02-selection.json': BATCH02_SELECTION_SHA,
        'outer-revalidation-batch02-review.json': BATCH02_SHA,
    }
    for name, sha in prior_files.items():
        require(digest((root / name).read_bytes()) == sha, 'Accepted historical file changed: ' + name)
    parent = read_json(root / 'outer-revalidation-selection.json')
    accepted = read_json(root / 'outer-revalidation-review.json')
    preserved = {**parent['prior_review_sha256'], **prior_files}
    spec = importlib.util.spec_from_file_location('historical_revalidation', root / 'verify_revalidation.py')
    historical = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(historical)
    previous = historical.check(root, corpus if full else None, gene_corpus if full else None)
    excluded = {key(r) for r in accepted['records']}
    require(len(excluded) == len(accepted['records']) == accepted['reviewed_count'] == 186,
            'Original accepted186 count differs')
    for stem, selsha, outer_count in (
            ('next128', NEXT128_SELECTION_SHA, 2475),
            ('batch02', BATCH02_SELECTION_SHA, 2603)):
        old = read_json(root / f'outer-revalidation-{stem}-review.json')
        old_selection = read_json(root / f'outer-revalidation-{stem}-selection.json')
        old_keys = {key(r) for r in old['records']}
        prior_pending = [r for r in parent['records'] if key(r) not in excluded]
        require(len(old_keys) == len(old['records']) == old['reviewed_count'] == 128,
                stem + ': expected128 distinct records')
        require(not excluded & old_keys and old_selection['records'] == prior_pending[:128],
                stem + ': changed fixed-order queue or overlap')
        for row, chosen in zip(old['records'], old_selection['records']):
            require(all(row.get(k) == v for k, v in chosen.items()), stem + ': source selection differs')
        require(old['selection_sha256'] == selsha and old['source_sha256'] == SOURCE_SHA
                and old['combined_outer_ledger_count'] == outer_count, stem + ': provenance/count differs')
        excluded.update(old_keys)
    pending = [r for r in parent['records'] if key(r) not in excluded]
    require(len(parent['records']) == 3598 and len(excluded) == 442 and len(pending) == 3156,
            'Queue after accepted442 differs')
    selection, review = read_json(selection_path), read_json(review_path)
    exclusion_manifest = [
        {'file': 'outer-revalidation-review.json', 'sha256': ACCEPTED186_SHA, 'record_count': 186},
        {'file': 'outer-revalidation-next128-review.json', 'sha256': NEXT128_SHA, 'record_count': 128},
        {'file': 'outer-revalidation-batch02-review.json', 'sha256': BATCH02_SHA, 'record_count': 128},
    ]
    require(selection['excluded_reviews'] == exclusion_manifest, 'Exclusion manifest differs')
    require(selection['source_sha256'] == SOURCE_SHA and selection['gene_source_sha256'] == GENE_SHA
            and selection['parent_selection_sha256'] == PARENT_SHA, 'Selection source hashes differ')
    require(selection['records'] == pending[:128], 'Batch03 is not exact first128 after excluding442')
    require(selection['selected_count'] == 128 and selection['pending_before'] == 3156
            and selection['pending_after'] == 3028, 'Selection totals differ')
    require(review['status'] == 'partial_source_revalidation', 'Unexpected review status')
    require(review['selection_sha256'] == digest(selection_path.read_bytes()), 'Selection hash differs')
    require(review['excluded_reviews'] == exclusion_manifest
            and review['parent_selection_sha256'] == PARENT_SHA
            and review['source_sha256'] == SOURCE_SHA and review['gene_source_sha256'] == GENE_SHA,
            'Review provenance differs')
    rows, selected = review['records'], selection['records']
    require(len(rows) == len({key(r) for r in rows}) == review['reviewed_count'] == 128,
            'Expected128 distinct source/ID records')
    require([key(r) for r in rows] == [key(r) for r in selected], 'Review assignment order differs')
    require(not {key(r) for r in rows} & excluded, 'New review overlaps accepted442')
    for row, chosen in zip(rows, selected):
        label = source_id(row)
        require(all(row.get(k) == v for k, v in chosen.items()), label + ': changed selection field')
        require(row['annotation_status'] == 'source_revalidated'
                and row['quality_validation'] == 'complete_literal_saved_string_re_read'
                and row['review_depth'] == 'complete_saved_abstract', label + ': incomplete declaration')
        require(isinstance(row['readers'], list) and row['readers']
                and all(isinstance(v, str) and v.strip() for v in row['readers']), label + ': absent readers')
        for field in ('note', 'model', 'uncertainty', 'source_abstract_quality', 'saved_content_kind',
                      'read_date', 'intervention_type', 'outcome_direction'):
            require(isinstance(row.get(field), str) and row[field].strip(), label + ': absent ' + field)
        require(type(row['requires_primary_followup']) is bool, label + ': nonboolean followup')
        require(isinstance(row['followup_reason'], str)
                and (not row['requires_primary_followup'] or row['followup_reason'].strip()), label + ': absent followup reason')
        require(isinstance(row['administered_agents'], list)
                and all(isinstance(a, dict) and a.get('name') and a.get('role') for a in row['administered_agents']),
                label + ': invalid agent roles')
        require(not {'abstractText', 'abstract', 'body_text', 'originalnote'} & row.keys(), label + ': raw/unexpected text key')
        require(isinstance(row.get('publication_metadata_flags'), dict), label + ': metadata missing')
        require(isinstance(row.get('read_evidence'), dict), label + ': reading evidence missing')

    require(review['selected_count'] == 128 and review['previous_outer_ledger_count'] == 2603,
            'Previous/count totals differ')
    require(review['combined_outer_ledger_count'] == 2603 + len(rows) == 2731,
            'Combined outer count differs')
    require(review['pending_count'] == len(pending) - len(rows) == 3028, 'Pending count differs')
    require(review['prior_accepted_revalidation_records'] == 442
            and review['cumulative_accepted_revalidation_records'] == 570, 'Revalidation totals differ')
    require(review['original_outer_records_without_abstract'] == 143, 'Originally empty count differs')
    require(review['followup_count'] == sum(r['requires_primary_followup'] for r in rows), 'Followup total differs')
    require(review['source_content_kinds'] == dict(Counter(r['saved_content_kind'] for r in rows)), 'Content kinds differ')
    require(review['source_quality_counts'] == dict(Counter(r['source_abstract_quality'] for r in rows)), 'Quality declarations differ')
    metadata = publication_counts(rows)
    require(review['publication_metadata_counts'] == metadata, 'Publication metadata totals differ')
    require(review['retracted_publication_count'] == metadata['retracted_publication_records']
            and review['retraction_link_count'] == metadata['records_with_retraction_links'], 'Retraction counts differ')
    require(review['correction_count'] == metadata['records_with_erratum_or_correction_links']
            and review['version_count'] == metadata['records_with_version_lists'], 'Correction/version totals differ')

    inputs, manifests = review['review_inputs'], review['batch_read_evidence']
    original_rows = original_rows_from_revisions(review, full, evidence_root)
    require(len(inputs) == 2 and len(manifests) == review['batch_count'] == 32, 'Input/batch total differs')
    literal_rows = []
    for part_index, part in enumerate(('a', 'b')):
        item = inputs[part_index]
        part_rows = original_rows[part_index * 64:(part_index + 1) * 64]
        require(item['part'] == part and item['file'] == f'part-{part}-review.json'
                and item['selection_file'] == f'part-{part}-selection.json', 'Input names/order differ')
        require(item['record_count'] == len(part_rows) == 64, 'Input part count differs')
        require(all(valid_hash(item[k]) for k in ('sha256', 'selection_sha256', 'record_objects_sha256')),
                'Input hash declaration invalid')
        require(item['record_objects_sha256'] == object_digest(part_rows), 'Source-part object digest differs')
        require(item['readers'] and all(r['readers'] == item['readers'] for r in part_rows), 'Reader attribution differs')
        part_batches = manifests[part_index * 16:(part_index + 1) * 16]
        for index, batch in enumerate(part_batches, 1):
            group = part_rows[(index - 1) * 4:index * 4]
            require(batch['part'] == part and batch['batch'] == index
                    and batch['log_file'] == f'part-{part}-batch-{index:02d}.jsonl', 'Batch filename/order differs')
            require(valid_hash(batch['log_sha256']) and batch['source_ids'] == [source_id(r) for r in group]
                    and batch['abstract_sha256s'] == [r['abstract_sha256'] for r in group], 'Batch identities/hashes differ')
            require(batch['requested_records'] == batch['emitted_records'] == 4
                    and batch['end_of_batch'] is True and batch['end_of_batch_observed'] is True
                    and batch['tool_output_untruncated'] is True, 'Missing complete-batch declaration')
            require(batch['readers'] == item['readers'] and batch.get('tool_output_chunk_id'), 'Batch reader/chunk absent')
            chars = 0
            for row in group:
                evidence = row['read_evidence']
                require(evidence['batch'] == index and evidence.get('batch_log_file', evidence.get('log_file')) == batch['log_file']
                        and evidence.get('batch_log_sha256', evidence.get('log_sha256')) == batch['log_sha256']
                        and evidence.get('end_of_batch_observed', evidence.get('end_marker_observed')) is True
                        and evidence['tool_output_chunk_id'] == batch['tool_output_chunk_id'],
                        source_id(row) + ': per-record evidence differs')
                require(type(evidence['abstract_characters']) is int and evidence['abstract_characters'] > 0,
                        source_id(row) + ': invalid character count')
                chars += evidence['abstract_characters']
            require(chars == batch['total_abstract_characters'], 'Declared batch character count differs')
        if full:
            path = evidence_path(evidence_root, item['file'])
            assignment_path = evidence_path(evidence_root, item['selection_file'])
            require(digest(path.read_bytes()) == item['sha256'], 'Source part-file hash differs')
            require(digest(assignment_path.read_bytes()) == item['selection_sha256'], 'Source assignment hash differs')
            payload, assignment = read_json(path), read_json(assignment_path)
            require(payload['records'] == part_rows, 'Source part objects were altered during consolidation')
            require(payload['selection_sha256'] == item['selection_sha256'], 'Source part assignment declaration differs')
            require(payload['source_sha256'] == SOURCE_SHA and payload['readers'] == item['readers'], 'Source part attribution differs')
            require(payload.get('read_method', payload.get('reading_method', '')).strip(), 'Source reading method absent')
            require(assignment['records'] == selected[part_index * 64:(part_index + 1) * 64]
                    and assignment['source_sha256'] == SOURCE_SHA
                    and assignment['parent_selection_sha256'] == PARENT_SHA
                    and assignment['excluded_reviews'] == exclusion_manifest, 'Source assignment differs')
            canonical, actual_rows = normalized_batches(part, payload, evidence_root)
            require(canonical == part_batches, 'Actual batch manifest differs')
            require([r['source_id'] for r in actual_rows] == [source_id(r) for r in part_rows], 'Actual batch order differs')
            for row, actual in zip(part_rows, actual_rows):
                require(row['abstract_sha256'] == actual['abstract_sha256']
                        and row['title'] == actual['title']
                        and row['read_evidence']['abstract_characters'] == actual['abstract_characters'], 'Actual source text differs')
            literal_rows.extend(actual_rows)
    require(review['total_abstract_characters'] == sum(m['total_abstract_characters'] for m in manifests),
            'Total characters differ')
    if full:
        raw = read_json(corpus)['results']
        lookup = {key(r): r for r in raw}
        require(len(lookup) == len(raw), 'Duplicate source corpus identity')
        for row, literal in zip(rows, literal_rows):
            actual = lookup[key(row)]
            require(actual.get('abstractText') == literal['abstractText']
                    and actual['title'] == row['title'], source_id(row) + ': corpus literal differs')
            flags = row['publication_metadata_flags']
            require(all(actual.get(k) == v for k, v in flags.items()), source_id(row) + ': metadata differs from corpus')
            for field in ('pubTypeList', 'publicationStatus', 'commentCorrectionList', 'versionList', 'versionNumber'):
                require(field not in actual or flags.get(field) == actual[field], source_id(row) + ': omitted metadata ' + field)
            require(all(actual.get(k) == v for k, v in row.get('source_links', {}).items()), 'Source links differ')
    audit = review['independent_audit']
    audit_ids = [source_id(rows[i]) for i in (0,15,32,63,64,79,96,127)]
    require(audit['source_ids'] == audit_ids and audit['record_count'] == 8
            and audit['cross_reader_checks'] == 8 and audit['self_rereads'] == 0
            and audit['additional_coverage_credit'] == 0, 'Independent audit totals/scope differ')
    require(audit['reviewers'] == ['/root'] and audit['sampling_method'].strip()
            and audit['annotation_comparison_tool_chunk_id'] and audit['interpretation_limit'].strip(),
            'Independent audit attribution/scope absent')
    require(len(audit['batch_evidence']) == 2, 'Independent audit batches differ')
    for i, evidence in enumerate(audit['batch_evidence']):
        group_ids = audit_ids[i*4:(i+1)*4]
        require(evidence['source_ids'] == group_ids and valid_hash(evidence['log_sha256'])
                and evidence['end_of_batch_observed'] is True and evidence['tool_output_untruncated'] is True
                and evidence['tool_output_chunk_id'], 'Incomplete independent audit declaration')
        if full:
            data = evidence_path(evidence_root, evidence['log_file']).read_bytes()
            require(digest(data) == evidence['log_sha256'], 'Independent audit log hash differs')
            lines = [json.loads(x) for x in data.decode().splitlines()]
            require(len(lines) == 6 and lines[0] == {'corpus_sha256': SOURCE_SHA, 'requested_records': 4},
                    'Independent audit header/count differs')
            emitted, end = lines[1:-1], lines[-1]
            require([r['source_id'] for r in emitted] == group_ids and end.get('end_of_batch') is True
                    and end['emitted_records'] == 4, 'Independent audit incomplete/IDs differ')
            by_id = {source_id(r): r for r in rows}
            for r in emitted:
                require(r['abstract_sha256'] == digest(r['abstractText'].encode()) == by_id[r['source_id']]['abstract_sha256']
                        and r['title'] == by_id[r['source_id']]['title']
                        and r['abstract_characters'] == len(r['abstractText']), 'Independent audit literal differs')
            require(end['total_abstract_characters'] == sum(r['abstract_characters'] for r in emitted)
                    == evidence['total_abstract_characters'], 'Independent audit character count differs')
    require(all(digest((root / name).read_bytes()) == sha for name, sha in preserved.items()), 'Historical input bytes changed')
    return dict(result='PASS', mode='full_source_evidence_replay' if full else 'structural_only',
                checked='Identity, exact additive queue, annotation objects and reading declarations; literal logs/corpora only in full mode',
                scientific_reading_fidelity_validated=False, primary_paper_qualification_validated=False,
                additional_literal_source_re_readings_declared=128, prior_accepted_revalidation_records=442,
                cumulative_accepted_revalidation_records=570, prior442_bytes_preserved=True,
                pending_source_revalidation=3028, combined_outer_ledger_records=2731,
                followup_count=review['followup_count'], source_content_kinds=review['source_content_kinds'],
                publication_metadata_counts=metadata, batch_count=32,
                total_abstract_characters=review['total_abstract_characters'],
                independent_cross_reader_checks=8, retracted_publication_count=review['retracted_publication_count'],
                annotation_revision_count=review['annotation_revision_count'],
                annotation_revision_record_count=review['annotation_revision_record_count'],
                source_part_objects_verified=full, external_batch_evidence_verified=full,
                source_selection_replayed=full, previous_verifier_result=previous['result'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    parser.add_argument('--selection', type=Path)
    parser.add_argument('--review', type=Path)
    parser.add_argument('--evidence-root', type=Path)
    parser.add_argument('--corpus', type=Path)
    parser.add_argument('--gene-corpus', type=Path)
    args = parser.parse_args()
    print(json.dumps(check(args.root, args.corpus, args.gene_corpus,
                           args.selection, args.review, args.evidence_root), indent=2))
