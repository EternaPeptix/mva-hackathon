#!/usr/bin/env python3
"""Screen EPA CPDat consumer-formulation chemical names against the fixed BUB1B corpus.

This is lexical discovery only. Product inclusion and a BUB1B co-mention do not
establish exposure, target engagement, efficacy, or a competitive treatment.
"""
import argparse
import csv
import hashlib
import io
import json
import sys
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from literature_match import LOW_SPECIFICITY, normalize  # noqa: E402
from build_dsld_inventory import COMMON  # noqa: E402

CPDAT_GENERIC = {
    'alpha', 'beta', 'gamma', 'delta', 'paper', 'inhibitor', 'inhibitors',
    'modifier', 'modifiers', 'carrier', 'carriers', 'magnetic', 'methyl',
    'ethyl', 'propyl', 'butyl', 'ester', 'esters', 'yeast',
    'saccharomyces cerevisiae', 'not available', 'not specified', 'unknown',
    'none', 'n a', 'na', 'n d', 'nd',
}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def read_cpdat(path):
    chemicals = {}
    row_count = 0
    category_stats = defaultdict(lambda: {'rows': 0, 'dtxsids': set(), 'product_ids': set()})
    with zipfile.ZipFile(path) as archive:
        table = 'cpdat_v4.0/cpdat_v4.0_product_composition_data.csv'
        with archive.open(table) as raw:
            rows = csv.DictReader(io.TextIOWrapper(raw, encoding='utf-8-sig', newline=''))
            for row in rows:
                if (row.get('puc_kind') or '').strip() != 'Formulation':
                    continue
                dtxsid = (row.get('dtxsid') or '').strip()
                if not dtxsid:
                    continue
                row_count += 1
                item = chemicals.setdefault(dtxsid, {
                    'names': set(), 'categories': set(), 'product_ids': set(), 'rows': 0,
                })
                item['rows'] += 1
                category = (row.get('puc_general_category') or '').strip() or 'UNKNOWN'
                category_stats[category]['rows'] += 1
                category_stats[category]['dtxsids'].add(dtxsid)
                for field in ('curated_chemical_name', 'raw_chemical_name'):
                    value = (row.get(field) or '').strip()
                    if value:
                        item['names'].add((field, value))
                if category:
                    item['categories'].add(category)
                product_id = (row.get('cpdat_product_id') or '').strip()
                if product_id:
                    item['product_ids'].add(product_id)
                    category_stats[category]['product_ids'].add(product_id)
    return chemicals, row_count, category_stats


def read_corpus(path):
    data = json.loads(Path(path).read_text())
    records = data.get('results')
    if not isinstance(records, list):
        raise ValueError('Expected a Europe PMC corpus object with a results list')
    return data, records


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cpdat-zip', type=Path, required=True)
    ap.add_argument('--corpus', type=Path, required=True)
    ap.add_argument('--decision-ledger', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    chemicals, formulation_rows, category_stats = read_cpdat(args.cpdat_zip)
    corpus, records = read_corpus(args.corpus)
    ledger = json.loads(args.decision_ledger.read_text())
    ledger_aliases = defaultdict(set)
    for decision in ledger['decisions']:
        for name in decision.get('match_aliases', []):
            alias = normalize(name)
            if alias:
                ledger_aliases[alias].add(decision['decision_id'])
    alias_to_ids = defaultdict(set)
    alias_to_names = defaultdict(set)
    alias_to_categories = defaultdict(set)
    for dtxsid, item in chemicals.items():
        for field, name in item['names']:
            alias = normalize(name)
            if not alias:
                continue
            alias_to_ids[alias].add(dtxsid)
            alias_to_names[alias].add((field, name))
            alias_to_categories[alias].update(item['categories'])

    by_first = defaultdict(set)
    for alias in alias_to_ids:
        words = alias.split()
        by_first[words[0]].add((len(words), alias))

    alias_rows = []
    for dtxsid, item in sorted(chemicals.items()):
        for field, name in sorted(item['names']):
            alias_rows.append({
                'dtxsid': dtxsid,
                'name_field': field,
                'name': name,
                'normalized_name': normalize(name),
                'categories': ';'.join(sorted(item['categories'])),
                'composition_rows': item['rows'],
                'product_ids': len(item['product_ids']),
            })
    with (args.out / 'consumer-formulation-aliases.tsv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(alias_rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(alias_rows)

    chemical_rows = []
    for dtxsid, item in sorted(chemicals.items()):
        curated = sorted(name for kind, name in item['names'] if kind == 'curated_chemical_name')
        raw = sorted(name for kind, name in item['names'] if kind == 'raw_chemical_name')
        chemical_rows.append({
            'dtxsid': dtxsid,
            'curated_names': json.dumps(curated, ensure_ascii=False),
            'raw_names': json.dumps(raw, ensure_ascii=False),
            'product_categories': ';'.join(sorted(item['categories'])),
            'composition_rows': item['rows'],
            'unique_product_ids': len(item['product_ids']),
        })
    with (args.out / 'consumer-formulations.tsv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(chemical_rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(chemical_rows)
    with (args.out / 'product-category-summary.tsv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['product_category', 'formulation_rows', 'unique_dtxsids', 'unique_product_ids'], delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for category, counts in sorted(category_stats.items()):
            writer.writerow({'product_category': category, 'formulation_rows': counts['rows'], 'unique_dtxsids': len(counts['dtxsids']), 'unique_product_ids': len(counts['product_ids'])})

    match_counts = Counter()
    matches = []
    all_source_ids = set()
    review_source_ids = set()
    review_dtxsids = set()
    matched_decisions = set()
    for record in records:
        metadata = {k: record.get(k, '') for k in ('source', 'id', 'pmid', 'pmcid', 'pubYear', 'doi', 'title')}
        for field in ('title', 'abstractText'):
            text = record.get(field) or ''
            words = normalize(text).split()
            found = set()
            for i, word in enumerate(words):
                for size, alias in by_first.get(word, ()):
                    if i + size <= len(words) and ' '.join(words[i:i + size]) == alias:
                        found.add(alias)
            for alias in sorted(found):
                ids = sorted(alias_to_ids[alias])
                low = alias in LOW_SPECIFICITY or alias in COMMON or alias in CPDAT_GENERIC or (len(alias.split()) == 1 and len(alias) <= 4)
                status = 'LOW_SPECIFICITY' if low else ('AMBIGUOUS_DTXSID' if len(ids) > 1 else 'REVIEW_EXACT_NAME')
                match_counts[status] += 1
                source_key = (record.get('source', ''), record.get('id', ''))
                all_source_ids.add(source_key)
                if status != 'REVIEW_EXACT_NAME':
                    continue
                review_source_ids.add(source_key)
                review_dtxsids.update(ids)
                matched_decisions.update(ledger_aliases.get(alias, set()))
                matches.append({
                    **metadata,
                    'matched_field': field,
                    'matched_alias': alias,
                    'matched_names': ' | '.join(f'{kind}:{name}' for kind, name in sorted(alias_to_names[alias])),
                    'dtxsids': ';'.join(ids),
                    'product_categories': ';'.join(sorted(alias_to_categories[alias])),
                    'candidate_status': status,
                    'current_decision_ids': ';'.join(sorted(ledger_aliases.get(alias, set()))) or 'NO_EXACT_LEDGER_ALIAS',
                })
    fields = ['source', 'id', 'pmid', 'pmcid', 'pubYear', 'doi', 'title', 'matched_field', 'matched_alias', 'matched_names', 'dtxsids', 'product_categories', 'candidate_status', 'current_decision_ids']
    with (args.out / 'bub1b-exact-name-co-mentions.tsv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(matches)

    payload = {
        'schema_version': 1,
        'status': 'LEXICAL_DISCOVERY_ONLY_NOT_EFFICACY_REVIEW',
        'cpdat_release': 'CPDat v4.0, Figshare version 5, posted 2025-04-09',
        'cpdat_doi': '10.23645/epacomptox.5352997',
        'cpdat_scope': 'Product composition table, puc_kind=Formulation; includes source product categories such as personal care, household cleaning, food/drug, pesticides, vehicles and others. Product records do not prove current availability, use, exposure, or ingredient concentration in a purchased item.',
        'cpdat_archive_sha256': sha(args.cpdat_zip),
        'gene_corpus_sha256': sha(args.corpus),
        'decision_ledger_sha256': sha(args.decision_ledger),
        'decision_ledger_count': len(ledger['decisions']),
        'gene_corpus_query': corpus.get('query'),
        'gene_corpus_records': len(records),
        'formulation_composition_rows_with_dtxsid': formulation_rows,
        'unique_formulation_dtxsids': len(chemicals),
        'unique_normalized_formulation_names': len(alias_to_ids),
        'cpdat_alias_rows': len(alias_rows),
        'formulation_product_categories': len(category_stats),
        'exact_name_screen_counts': dict(match_counts),
        'exact_unique_name_matches_written': len(matches),
        'source_id_records_with_any_exact_or_low_specificity_match': len(all_source_ids),
        'source_id_records_with_review_class_exact_match': len(review_source_ids),
        'review_class_exact_match_dtxsids': len(review_dtxsids),
        'review_class_matches_to_existing_decision_aliases': len(matched_decisions),
        'limitations': [
            'CPDat v4.0 is a dated source release; the EPA ChemExpo page reports CPDat 4.1 as the current release, which was not obtained in this run.',
            'This is a lexical exact-name screen of curated/raw chemical names against titles and abstracts in the fixed BUB1B/BUBR1 literature corpus.',
            'LOW_SPECIFICITY names use the existing BUB1B corpus common-term list plus explicit CPDat generic labels; they are counted but excluded from the review list. Multi-DTXSID names are also excluded.',
            'A co-mention is not evidence that the chemical was administered, changed BUB1B/BUBR1, or improved MVA phenotype.',
            'No CPDat record establishes current market presence, actual household exposure, dose, route, or product concentration.',
        ],
    }
    (args.out / 'screen-summary.json').write_text(json.dumps(payload, indent=2, ensure_ascii=False) + '\n')
    print(json.dumps({k: payload[k] for k in ('cpdat_release','formulation_composition_rows_with_dtxsid','unique_formulation_dtxsids','unique_normalized_formulation_names','cpdat_alias_rows','exact_name_screen_counts','exact_unique_name_matches_written','source_id_records_with_any_exact_or_low_specificity_match','source_id_records_with_review_class_exact_match','review_class_exact_match_dtxsids','review_class_matches_to_existing_decision_aliases','limitations')},indent=2,ensure_ascii=False))


if __name__ == '__main__':
    main()
