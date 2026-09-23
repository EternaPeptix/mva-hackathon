#!/usr/bin/env python3
"""Preserve FDA food-substance rows and screen literal names in a saved corpus.

Inventory membership is neither current-use authorization nor BUB1B efficacy.
No structure, CAS identity, availability, food safety or drug use is inferred.
"""
import argparse
import json
import os
import re
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import parse_qs, urljoin, urlparse

from broad_hub_screen import clean
from sgc_probe_inventory import Tables, combined_pattern, digest, identity, read_tsv, write_tsv
from literature_match import LOW_SPECIFICITY  # Vocabulary only; never use its old text decoder.

HERE = Path(__file__).resolve().parent
URL = ('https://www.hfpappexternal.fda.gov/scripts/fdcc/index.cfm?order=ASC&search='
       '&set=FoodSubstances&showAll=true&sort=Sortterm_ID&type=basic')
HEADERS = ['CAS Reg. No.* (or other ID)', 'Substance* (sorted A-Z)',
           'Used for*† (Technical Effect)', '21 CFR*']
SUFFIX = re.compile(r'\s*--\s*(PROHIBITED(?: WITH EXCEPTIONS)?|DELISTED|NLFG)\s*$', re.I)
GENE = re.compile(r'\b(?:bub1b|bubr1|bub r1|mosaic variegated aneuploidy)\b')
MENTION_FIELDS = ['source_row_id', 'query', 'retrieval_policy', 'source', 'record_id',
                  'matched_fields', 'gene_literal_in_title_abstract', 'title',
                  'source_url', 'abstract_sha256', 'evidence_state']
COMMON_FOOD_CLASSES = {'casein', 'collagen', 'fatty acids', 'pulp', 'yeasts'}


def retrieval_policy(name):
    # Normalize separators for ambiguity classification only. Identity and
    # retrieval retain the original hyphenation and stereochemical labels.
    policy_name = ' '.join(re.sub(r'[-_]+', ' ', identity(name)).split())
    return ('COMMON_OR_SHORT_NAME_CONTEXT_AMBIGUOUS'
            if policy_name in LOW_SPECIFICITY or policy_name in COMMON_FOOD_CLASSES
            or len(identity(name)) < 5 else 'LITERAL_NAME_CONTEXT_UNASSESSED')


def parse_source(raw):
    text = raw.decode('utf-8')
    declared = {int(x) for x in re.findall(r'Records Found:\s*(\d+)', text)}
    if len(declared) != 1:
        raise ValueError('Missing or inconsistent FDA declared record count')
    parser = Tables()
    parser.feed(text)
    tables = [t for t in parser.tables if len(t) > 1 and [c['text'] for c in t[0]] == HEADERS]
    if len(tables) != 1 or len(tables[0]) - 1 != next(iter(declared)):
        raise ValueError('FDA table header/count mismatch; partial pages cannot pass')
    rows = []
    for number, cells in enumerate(tables[0][1:], 1):
        if len(cells) != 4 or not cells[1]['text'] or len(cells[1]['links']) != 1:
            raise ValueError('Malformed FDA row ' + str(number))
        detail = urljoin(URL, cells[1]['links'][0])
        query = parse_qs(urlparse(detail).query)
        if query.get('set') != ['FoodSubstances'] or len(query.get('id', [])) != 1:
            raise ValueError('Unexpected FDA detail URL')
        label = cells[1]['text']
        match = SUFFIX.search(label)
        retrieval = label[:match.start()].strip() if match else label
        if not retrieval:
            raise ValueError('Empty substance retrieval name')
        rows.append({'source_row_id': f'fda-food:{number:04d}',
                     'source_record_id': query['id'][0], 'detail_url': detail,
                     'cas_or_other_id_as_reported': cells[0]['text'],
                     'substance_as_reported': label,
                     'technical_effect_as_reported': cells[2]['text'],
                     'cfr_as_reported': cells[3]['text'], 'cfr_urls': cells[3]['links'],
                     'explicit_name_status': match.group(1).upper() if match else 'NOT_STATED_IN_NAME',
                     'retrieval_name': retrieval,
                     'retrieval_derivation': 'removed_explicit_status_suffix' if match else 'literal_source_name',
                     'current_permission': 'NOT_DETERMINED_BY_INVENTORY',
                     'BUB1B_efficacy': 'UNKNOWN_NOT_ASSESSED_BY_INVENTORY'})
    # A source can legitimately repeat names/CAS numbers; preserve rows rather
    # than silently merging substances or assigning unique molecular identity.
    updated = re.findall(r'Page Last Updated:\s*(\d{2}/\d{2}/\d{4})', text)
    return rows, sorted(set(updated))


def build(html_path, corpus_path, catalog_dir, output):
    paths = [html_path, corpus_path] + [catalog_dir / n for n in
            ['catalog.tsv', 'aliases.tsv', 'drugcentral-status-unverified.tsv', 'decisions.json']]
    before = {str(p): digest(p.read_bytes()) for p in paths}
    if len({p.name for p in paths}) != len(paths):
        raise ValueError('Input basenames must be distinct for hash metadata')
    rows, updated = parse_source(html_path.read_bytes())
    catalog = read_tsv(catalog_dir / 'catalog.tsv')
    ids = {r['catalog_id'] for r in catalog}
    if len(ids) != len(catalog):
        raise ValueError('Duplicate catalog IDs')
    catalog_index = defaultdict(set)
    for row in read_tsv(catalog_dir / 'aliases.tsv'):
        if row['catalog_id'] not in ids:
            raise ValueError('Dangling catalog alias')
        catalog_index[identity(row['alias'])].add(row['catalog_id'])
    for row in read_tsv(catalog_dir / 'drugcentral-status-unverified.tsv'):
        for alias in set(json.loads(row['aliases']) + [row['name']]):
            catalog_index[identity(alias)].add('DC' + row['drugcentral_id'])
    appraisal_index = defaultdict(set)
    for row in json.loads((catalog_dir / 'decisions.json').read_text())['decisions']:
        for alias in row['match_aliases']:
            appraisal_index[identity(alias)].add(row['decision_id'])
    terms = defaultdict(list)
    for row in rows:
        name = identity(row['retrieval_name'])
        row['catalog_name_candidates'] = sorted(catalog_index[name])
        row['curated_appraisal_name_candidates'] = sorted(appraisal_index[name])
        row['identity_state'] = 'NAME_CROSS_REFERENCE_ONLY_NOT_STRUCTURE_OR_CAS_VERIFIED'
        row['retrieval_policy'] = retrieval_policy(name)
        terms[name].append(row)
    corpus = json.loads(corpus_path.read_bytes())['results']
    if len({(r['source'], str(r['id'])) for r in corpus}) != len(corpus):
        raise ValueError('Duplicate corpus source/ID records')
    pattern = combined_pattern(terms)
    mentions = []
    for record in corpus:
        fields = {f: identity(clean(record.get(f, '') or '')) for f in ['title', 'abstractText']}
        found = defaultdict(set)
        for field, text in fields.items():
            for hit in pattern.finditer(text):
                found[hit.group()].add(field)
        for name, fields_hit in sorted(found.items()):
            for row in terms[name]:
                mentions.append({'source_row_id': row['source_row_id'], 'query': row['retrieval_name'],
                                 'retrieval_policy': row['retrieval_policy'], 'source': record['source'],
                                 'record_id': str(record['id']), 'matched_fields': sorted(fields_hit),
                                 'gene_literal_in_title_abstract': bool(GENE.search(' '.join(fields.values()))),
                                 'title': clean(record.get('title', '') or ''),
                                 'source_url': 'https://europepmc.org/article/' + record['source'] + '/' + str(record['id']),
                                 'abstract_sha256': digest((record.get('abstractText', '') or '').encode()),
                                 'evidence_state': 'LEXICAL_MATCH_ONLY_REVIEW_STATUS_NOT_ASSIGNED'})
    hits = Counter(r['source_row_id'] for r in mentions)
    for row in rows:
        row['lexical_mention_rows'] = hits[row['source_row_id']]
    summary = {'schema_version': 1, 'source_url': URL, 'source_page_update_dates': updated,
               'source_input_sha256': {p.name: before[str(p)] for p in paths},
               'code_sha256': {n: digest((HERE / n).read_bytes()) for n in
                               ['food_substance_inventory.py', 'sgc_probe_inventory.py', 'broad_hub_screen.py', 'literature_match.py']},
               'source_rows': len(rows), 'unique_source_record_ids': len({r['source_record_id'] for r in rows}),
               'unique_literal_labels': len({r['substance_as_reported'] for r in rows}),
               'unique_cas_or_other_id_strings': len({r['cas_or_other_id_as_reported'] for r in rows}),
               'explicit_name_status_counts': dict(Counter(r['explicit_name_status'] for r in rows)),
               'rows_with_catalog_name_candidate': sum(bool(r['catalog_name_candidates']) for r in rows),
               'rows_with_curated_appraisal_name_candidate': sum(bool(r['curated_appraisal_name_candidates']) for r in rows),
               'corpus_source_id_records': len(corpus), 'lexical_mention_rows': len(mentions),
               'mention_rows_by_policy': dict(Counter(r['retrieval_policy'] for r in mentions)),
               'source_rows_with_lexical_hit': len(hits),
               'matched_source_id_records': len({(r['source'], r['record_id']) for r in mentions}),
               'manual_reading_credit': 0, 'new_curated_appraisal_credit': 0,
               'scope': 'Complete retrieved FDA source table, not all foods or unique molecules. Literal preferred names only; no synonyms, CAS matching, individual regulatory review, current stock or efficacy assessment. Explicit suffix flags retained; unflagged does not mean permitted or safe. Co-mentions can be endogenous molecules, media, controls, harmful exposures or background. A lexical miss is unknown. Source rows overlap the existing catalogs; do not add their count to the main entity denominator.'}
    # Do not damage a previously valid package when an input changes while
    # screening. Prepare all products, validate inputs, then replace files.
    # Publishing a directory is not a multi-file transaction; the summary is
    # replaced last, so an interrupted publish is detectable through hashes.
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.fda-food-', dir=output.parent) as temp:
        stage = Path(temp)
        write_tsv(stage / 'substances.tsv', rows, list(rows[0]))
        write_tsv(stage / 'mentions.tsv', mentions, MENTION_FIELDS)
        summary['output_sha256'] = {n: digest((stage / n).read_bytes())
                                    for n in ['substances.tsv', 'mentions.tsv']}
        (stage / 'summary.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False) + '\n')
        for p in paths:
            if digest(p.read_bytes()) != before[str(p)]:
                raise ValueError('Input changed during screen: ' + str(p))
        output.mkdir(parents=True, exist_ok=True)
        for name in ['substances.tsv', 'mentions.tsv', 'summary.json']:
            os.replace(stage / name, output / name)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--html', required=True, type=Path)
    parser.add_argument('--corpus', required=True, type=Path)
    parser.add_argument('--catalog-dir', default=HERE, type=Path)
    parser.add_argument('--out', required=True, type=Path)
    a = parser.parse_args()
    print(json.dumps(build(a.html, a.corpus, a.catalog_dir, a.out), indent=2))


if __name__ == '__main__':
    main()
