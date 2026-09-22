#!/usr/bin/env python3
"""EMA source-only product/name coverage layer, never a unique-drug census.

Source: European Medicines Agency (EMA), with acknowledgement required under
https://www.ema.europa.eu/en/about-us/about-website/legal-notice .
No literature retrieval, substance identity merge, salt removal or fuzzy match.
"""
import argparse
import collections
import csv
import datetime
import hashlib
import json
from pathlib import Path
import unicodedata

SOURCE_URL = 'https://www.ema.europa.eu/en/documents/report/medicines-output-medicines_json-report_en.json'
GUIDE_URL = 'https://www.ema.europa.eu/en/about-us/about-website/download-website-data-json-data-format'
TERMS_URL = 'https://www.ema.europa.eu/en/about-us/about-website/legal-notice'
FIELDS = ['category', 'ema_product_number', 'name_of_medicine', 'active_substance',
          'international_non_proprietary_name_common_name', 'medicine_status',
          'opinion_status', 'marketing_authorisation_date',
          'european_commission_decision_date', 'opinion_adopted_date',
          'first_published_date', 'last_updated_date', 'medicine_url']
NAME_FIELDS = ['name_of_medicine', 'active_substance', 'international_non_proprietary_name_common_name']
INPUT_FILES = ['catalog.tsv', 'aliases.tsv', 'drugcentral-status-unverified.tsv',
               'approval-overrides.json', 'broad-hub-screen.tsv']


def norm(s):
    return ' '.join(unicodedata.normalize('NFKC', s or '').casefold().split())


def digest(path):
    b = path.read_bytes()
    return {'bytes': len(b), 'sha256': hashlib.sha256(b).hexdigest()}


def read_tsv(path):
    with path.open(encoding='utf-8', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def build_indexes(catalog_dir):
    indexes = {k: collections.defaultdict(set) for k in ['main', 'status_unverified_supplement', 'fda2026', 'broad_preferred_names']}
    for r in read_tsv(catalog_dir / 'catalog.tsv'):
        indexes['main'][norm(r['preferred_name'])].add(r['catalog_id'])
    for r in read_tsv(catalog_dir / 'aliases.tsv'):
        indexes['main'][norm(r['alias'])].add(r['catalog_id'])
    for r in read_tsv(catalog_dir / 'drugcentral-status-unverified.tsv'):
        for v in [r['name']] + json.loads(r['aliases']):
            indexes['status_unverified_supplement'][norm(v)].add('DC:' + r['drugcentral_id'])
    for r in json.loads((catalog_dir / 'approval-overrides.json').read_text())['records']:
        for v in [r['brand'], r['active_ingredient']] + r['search_aliases']:
            indexes['fda2026'][norm(v)].add(r['brand'])
    for r in read_tsv(catalog_dir / 'broad-hub-screen.tsv'):
        indexes['broad_preferred_names'][norm(r['preferred_source_name'])].add(r['hub_row_id'])
    return indexes


def cross_refs(value, indexes):
    key = norm(value)
    return {layer: sorted(idx[key]) for layer, idx in indexes.items() if key and key in idx}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--ema-json', required=True, type=Path)
    p.add_argument('--catalog-dir', required=True, type=Path)
    p.add_argument('--output-dir', required=True, type=Path)
    p.add_argument('--retrieval-metadata', required=True, type=Path,
                   help='JSON list of original HTTP retrieval records; must contain this feed SHA and URL')
    a = p.parse_args()
    raw = a.ema_json.read_bytes()
    feed_sha = hashlib.sha256(raw).hexdigest()
    feed = json.loads(raw)
    data = feed['data']
    assert isinstance(data, list) and len(data) == feed['meta']['total_records'], 'Record count mismatch'
    assert all(all(k in r for k in FIELDS) for r in data), 'Missing required source fields'
    ids = [r['ema_product_number'] for r in data]
    assert len(ids) == len(set(ids)) and all(ids), 'Duplicate or missing EMA product ID; review before running'
    meta = json.loads(a.retrieval_metadata.read_text())
    retrievals = [m for m in meta if m.get('sha256') == feed_sha and m.get('url') == SOURCE_URL]
    assert len(retrievals) == 1, 'Feed must have one matching original retrieval record'
    retrieval_time = datetime.datetime.fromisoformat(retrievals[0]['retrieved_utc'].replace('Z', '+00:00'))
    attribution = 'European Medicines Agency (EMA), accessed ' + retrieval_time.strftime('%B %Y')
    indexes = build_indexes(a.catalog_dir)
    authorised = []
    for r in data:
        if r['category'] != 'Human' or r['medicine_status'] != 'Authorised':
            continue
        rec = {k: r[k] for k in FIELDS}
        matches = {field: cross_refs(r[field], indexes) for field in NAME_FIELDS}
        rec['exact_name_cross_references_by_field'] = matches
        rec['name_match_state'] = ('EXACT_NAME_CROSS_REFERENCE_ONLY' if any(matches.values())
                                   else 'NO_EXACT_WHOLE_FIELD_MATCH_IDENTITY_UNKNOWN')
        rec['multiple_candidates_in_one_field_layer'] = any(len(ids) > 1 for m in matches.values() for ids in m.values())
        rec['compound_field_has_semicolon'] = any(';' in r[k] for k in NAME_FIELDS[1:])
        rec['substance_or_product_equivalence'] = 'UNKNOWN_NOT_EVALUATED'
        rec['bub1b_efficacy_state'] = 'UNKNOWN_NOT_ASSESSED_IN_EMA_LAYER'
        rec['source_attribution'] = attribution
        rec['source_feed_sha256'] = feed_sha
        authorised.append(rec)
    authorised.sort(key=lambda r: r['ema_product_number'])
    no_match = sum(r['name_match_state'].startswith('NO_EXACT') for r in authorised)
    source_time = datetime.datetime.fromisoformat(feed['meta']['timestamp'].replace('Z', '+00:00'))
    retrieval_time = datetime.datetime.fromisoformat(retrievals[0]['retrieved_utc'].replace('Z', '+00:00'))
    out = a.output_dir
    out.mkdir(parents=True, exist_ok=True)
    summary = {
        'artifact_schema_version': 1,
        'source_attribution': attribution,
        'source_url': SOURCE_URL, 'source_guide_url': GUIDE_URL,
        'source_reuse_terms_url': TERMS_URL,
        'source_reuse_terms_summary': 'EMA-owned public information may be reproduced and distributed for commercial and non-commercial purposes with EMA attribution in each copy. Third-party content and non-EMA images are excluded. No named Creative Commons licence is asserted.',
        'source_feed_sha256': feed_sha, 'source_retrieval': retrievals[0],
        'source_meta_unmodified': feed['meta'],
        'source_explicit_schema_or_version': None,
        'source_schema_or_version_note': 'Source exposes meta.total_records and meta.timestamp; no explicit schema/version field is present.',
        'source_root_keys': list(feed),
        'source_record_fields': sorted(set().union(*(r.keys() for r in data))),
        'source_meta_timestamp_after_retrieval': source_time > retrieval_time,
        'timestamp_note': 'Preserve embedded meta.timestamp, actual retrieval time and HTTP Last-Modified separately. No timestamp correction is inferred.',
        'scope': 'EMA centralised-procedure product inventory, currently Authorised Human rows only in the derived table. Not all national European authorisations or all drugs worldwide.',
        'product_row_unit': 'One EMA product number. Not a unique substance, active moiety, formulation, marketed-country or clinically equivalent-product count.',
        'source_total_product_rows': len(data),
        'source_unique_nonblank_product_ids': len(set(ids)),
        'source_category_status_product_counts': {' | '.join(k): v for k, v in sorted(collections.Counter((r['category'], r['medicine_status']) for r in data).items())},
        'authorised_human_product_rows': len(authorised),
        'authorised_human_product_rows_no_exact_whole_field_match_any_layer': no_match,
        'authorised_human_product_rows_with_exact_name_cross_reference': len(authorised) - no_match,
        'authorised_human_product_rows_missing_initial_authorisation_date': sum(not r['marketing_authorisation_date'] for r in authorised),
        'name_normalization': 'Unicode NFKC, casefold, whitespace collapse only; whole source fields. No fuzzy match, salt removal, combination splitting, brand pruning or identity merge.',
        'identity_limit': 'Every exact match is only a name cross-reference. Multiple IDs are retained. Missing matches are UNKNOWN. Combinations and cell/vaccine descriptions remain product-level fields.',
        'date_limit': 'marketing_authorisation_date is the initial date. european_commission_decision_date can concern a later change. Positive opinion does not equal authorisation. Blank dates remain blank.',
        'catalog_inputs': {f: digest(a.catalog_dir / f) for f in INPUT_FILES},
        'main_census_entity_denominator_changed': False,
        'literature_searched': False,
    }
    (out / 'ema-authorised-human-products.json').write_text(json.dumps(authorised, ensure_ascii=False, indent=2) + '\n')
    with (out / 'ema-authorised-human-products.tsv').open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(authorised[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for r in authorised:
            writer.writerow({k: json.dumps(v, ensure_ascii=False, separators=(',', ':')) if isinstance(v, (dict, list, bool)) else v for k, v in r.items()})
    (out / 'ema-coverage-summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\n')
    assert len(read_tsv(out / 'ema-authorised-human-products.tsv')) == len(authorised)
    print(json.dumps({'authorised_human_product_rows': len(authorised), 'no_exact_whole_field_matches': no_match, 'source_sha256': feed_sha}))


if __name__ == '__main__':
    main()
