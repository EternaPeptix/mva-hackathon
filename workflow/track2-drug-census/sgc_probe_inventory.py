#!/usr/bin/env python3
"""Enumerate SGC probe/control labels and screen a saved literature corpus.

Offline inputs only. Counts are labels/rows, not structure-deduplicated drugs.
No name match, target annotation or literature hit confers efficacy evidence.
"""
import argparse
import csv
import hashlib
import json
import re
import unicodedata
from collections import Counter, defaultdict
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urljoin

from broad_hub_screen import clean

HERE = Path(__file__).resolve().parent
SOURCES = {
    'all': ('https://www.thesgc.org/chemical-probes',
            ['Name', 'Structures', 'Protein Family', 'Target Protein',
             'Negative Control', 'Cellular usage', 'Added']),
    'donated': ('https://www.thesgc.org/donated-chemical-probes',
                ['Probe Name', 'Probe Structure', 'Protein Family',
                 'Target protein', 'Control name', 'Control structure',
                 'Cellular usage']),
}
GENE = re.compile(r'\b(?:bub1b|bubr1|bub r1|mosaic variegated aneuploidy)\b')
EMPTY_SCHEMAS = {
    'controls.tsv': ['control_record_id', 'parent_probe_id', 'parent_probe_label',
                     'source_control_cell', 'identity_scope', 'search_names',
                     'catalog_name_candidates', 'therapeutic_evidence',
                     'lexical_mention_rows'],
    'mentions.tsv': ['entity_record_id', 'role', 'query', 'match_policy', 'source',
                     'record_id', 'matched_fields', 'title', 'source_url',
                     'gene_literal_in_title_abstract', 'abstract_sha256',
                     'evidence_state'],
}


def identity(text):
    """Keep stereochemistry, punctuation and code separators distinct."""
    return ' '.join(unicodedata.normalize('NFKC', text).casefold().split())


def digest(data):
    return hashlib.sha256(data).hexdigest()


def key(prefix, text):
    return prefix + ':' + digest(text.encode())[:16]


def read_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write_tsv(path, rows, fields):
    with path.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t',
                                lineterminator='\n')
        writer.writeheader()
        for row in rows:
            writer.writerow({field: json.dumps(row[field], ensure_ascii=False)
                             if isinstance(row[field], (dict, list)) else row[field]
                             for field in fields})


class Tables(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.tables = []
        self.table = self.row = self.cell = None

    def handle_starttag(self, tag, attributes):
        if tag == 'table':
            if self.table is not None:
                raise ValueError('Nested source table requires explicit handling')
            self.table = []
        if self.table is None:
            return
        if tag == 'tr':
            self.row = []
        if tag in ('td', 'th'):
            self.cell = {'text': [], 'links': []}
        if tag == 'a' and self.cell is not None:
            href = dict(attributes).get('href')
            if href:
                self.cell['links'].append(href)
        if tag == 'br' and self.cell is not None:
            self.cell['text'].append(' ')

    def handle_data(self, data):
        if self.cell is not None:
            self.cell['text'].append(data)

    def handle_endtag(self, tag):
        if tag in ('td', 'th') and self.cell is not None:
            self.cell['text'] = ' '.join(''.join(self.cell['text']).split())
            self.row.append(self.cell)
            self.cell = None
        if tag == 'tr' and self.row is not None:
            self.table.append(self.row)
            self.row = None
        if tag == 'table' and self.table is not None:
            self.tables.append(self.table)
            self.table = None


def parse_source(data, source):
    url, header = SOURCES[source]
    parser = Tables()
    parser.feed(data.decode('utf-8'))
    tables = [table for table in parser.tables if table
              and [cell['text'] for cell in table[0]] == header]
    if len(tables) != 1 or len(tables[0]) < 2:
        raise ValueError('Expected exactly one nonempty SGC probe table')
    result = []
    for number, row in enumerate(tables[0][1:], 1):
        if len(row) != 7 or not row[0]['text'] or len(row[0]['links']) != 1:
            raise ValueError(f'Unexpected source row: {source}:{number}')
        result.append({
            'row_id': f'sgc-{source}:{number:04d}',
            'collection': source, 'source_url': url,
            'source_sha256': digest(data), 'probe_label': row[0]['text'],
            'probe_url': urljoin(url, row[0]['links'][0]),
            'protein_family': row[2]['text'], 'target_text': row[3]['text'],
            'negative_control_cell': row[4]['text'],
            'cellular_usage_as_reported': row[5 if source == 'all' else 6]['text'],
            'added_as_reported': row[6]['text'] if source == 'all' else '',
        })
    return result


def combined_pattern(queries):
    # Bound both ends: BAY-179 must not match BAY-1797. Do not strip + or -.
    terms = sorted({identity(q) for q in queries}, key=lambda q: (-len(q), q))
    if not terms or any(not q for q in terms):
        raise ValueError('Search terms must be nonempty')
    return re.compile(r'(?<![\w])(?:' + '|'.join(re.escape(q) for q in terms) + r')(?![\w])')


def aliases(label):
    result = [(label, 'literal_source_label')]
    if label == 'LLY-507 (multiple off-targets)':
        result = [('LLY-507', 'source_label_without_warning')]
    elif label == 'Bromosporine (BSP)':
        result = [('Bromosporine', 'source_explicit_name'),
                  ('BSP', 'ambiguous_short_source_alias')]
    elif label == '(+)-JQ1':
        result.append(('JQ1', 'ambiguous_unspecified_stereochemistry'))
    return result


def build(main_html, donated_html, corpus_path, catalog_dir, output):
    input_paths = [main_html, donated_html, corpus_path] + [catalog_dir / name for name in
                   ['catalog.tsv', 'aliases.tsv', 'drugcentral-status-unverified.tsv',
                    'decisions.json']]
    hashes = {str(path): digest(path.read_bytes()) for path in input_paths}
    rows = parse_source(main_html.read_bytes(), 'all') + parse_source(donated_html.read_bytes(), 'donated')
    grouped = defaultdict(list)
    for row in rows:
        grouped[identity(row['probe_label'])].append(row)
    catalog_index, appraisal_index = defaultdict(set), defaultdict(set)
    catalog_rows = read_tsv(catalog_dir / 'catalog.tsv')
    catalog_ids = {r['catalog_id'] for r in catalog_rows}
    if len(catalog_ids) != len(catalog_rows):
        raise ValueError('Duplicate main catalog IDs')
    for row in read_tsv(catalog_dir / 'aliases.tsv'):
        if row['catalog_id'] not in catalog_ids:
            raise ValueError('Alias references missing catalog ID: ' + row['catalog_id'])
        catalog_index[identity(row['alias'])].add(row['catalog_id'])
    supplemental_rows = read_tsv(catalog_dir / 'drugcentral-status-unverified.tsv')
    supplemental_ids = {'DC' + r['drugcentral_id'] for r in supplemental_rows}
    if len(supplemental_ids) != len(supplemental_rows) or supplemental_ids & catalog_ids:
        raise ValueError('Duplicate or overlapping supplemental IDs')
    for row in supplemental_rows:
        for alias in set(json.loads(row['aliases']) + [row['name']]):
            catalog_index[identity(alias)].add('DC' + row['drugcentral_id'])
    decisions = json.loads((catalog_dir / 'decisions.json').read_text())['decisions']
    if len({r['decision_id'] for r in decisions}) != len(decisions):
        raise ValueError('Duplicate appraisal IDs')
    for row in decisions:
        for alias in row['match_aliases']:
            appraisal_index[identity(alias)].add(row['decision_id'])
    probes, controls, search_terms = [], [], []
    for normalized, members in sorted(grouped.items()):
        label = members[0]['probe_label']
        pid = key('sgc-probe', normalized)
        queries = aliases(label)
        exact_queries = [q for q, policy in queries if not policy.startswith('ambiguous')]
        ids = sorted({gid for q in exact_queries for gid in catalog_index[identity(q)]})
        appraisals = sorted({did for q in exact_queries for did in appraisal_index[identity(q)]})
        targets = sorted({r['target_text'] for r in members})
        probe = {'probe_id': pid, 'probe_label': label,
                 'source_row_ids': [r['row_id'] for r in members],
                 'source_urls': sorted({r['probe_url'] for r in members}),
                 'protein_families': sorted({r['protein_family'] for r in members}),
                 'target_strings': targets,
                 'literal_BUB1B_target': any(re.search(r'\bBUB1B\b', t, re.I) for t in targets),
                 'catalog_name_candidates': ids, 'appraisal_name_candidates': appraisals,
                 'identity_state': 'NAME_CROSS_REFERENCE_ONLY_NOT_STRUCTURE_VERIFIED',
                 'bub1b_efficacy': 'UNKNOWN_NOT_ASSESSED_BY_INVENTORY',
                 'competitive_state': 'CURATED_APPRAISAL_LINK_ONLY' if appraisals else 'UNASSESSED'}
        probes.append(probe)
        for query, policy in queries:
            search_terms.append((pid, 'probe', query, policy))
        # Keep control fields scoped to their parent probe. A descriptive
        # enantiomer label is not a globally identified molecule.
        for cell in sorted({r['negative_control_cell'] for r in members if r['negative_control_cell']}):
            cid = key('sgc-control', normalized + '\t' + identity(cell))
            scope = 'contextual_stereoisomer_description' if cell == '2S,3R-enantiomer' else 'named_control_cell'
            qs = ['nc_WDR5', 'nc_VHL'] if cell == 'nc_WDR5, nc_VHL' else [cell]
            if scope == 'contextual_stereoisomer_description':
                qs = []
            controls.append({'control_record_id': cid, 'parent_probe_id': pid,
                             'parent_probe_label': label, 'source_control_cell': cell,
                             'identity_scope': scope, 'search_names': qs,
                             'catalog_name_candidates': sorted({gid for q in qs for gid in catalog_index[identity(q)]}),
                             'therapeutic_evidence': 'NOT_INFERRED_FROM_NEGATIVE_CONTROL_ROLE'})
            for q in qs:
                search_terms.append((cid, 'negative_control', q, 'literal_source_control_name'))
    corpus = json.loads(corpus_path.read_bytes())['results']
    keys = [(r['source'], str(r['id'])) for r in corpus]
    if len(keys) != len(set(keys)):
        raise ValueError('Duplicate corpus source/ID records')
    # One combined pattern finds candidate spans; an explicit index preserves
    # ambiguous queries and control roles, including multiple parent probes.
    term_index = defaultdict(list)
    for entry in search_terms:
        term_index[identity(entry[2])].append(entry)
    combined = combined_pattern(term_index)
    mentions = []
    for record in corpus:
        fields = {field: identity(clean(record.get(field, '') or '')) for field in ['title', 'abstractText']}
        found = defaultdict(set)
        for field, text in fields.items():
            for hit in combined.finditer(text):
                found[hit.group()].add(field)
        for matched, fields_hit in sorted(found.items()):
            for entity_id, role, query, policy in term_index[matched]:
                mentions.append({'entity_record_id': entity_id, 'role': role,
                                 'query': query, 'match_policy': policy,
                                 'source': record['source'], 'record_id': str(record['id']),
                                 'matched_fields': sorted(fields_hit),
                                 'title': clean(record.get('title', '')),
                                 'source_url': 'https://europepmc.org/article/' + record['source'] + '/' + str(record['id']),
                                 'gene_literal_in_title_abstract': bool(GENE.search(' '.join(fields.values()))),
                                 'abstract_sha256': digest((record.get('abstractText', '') or '').encode()),
                                 'evidence_state': 'LEXICAL_MATCH_ONLY_CONSULT_REVIEW_OVERLAY'})
    hits = Counter(r['entity_record_id'] for r in mentions)
    for probe in probes:
        probe['lexical_mention_rows'] = hits[probe['probe_id']]
    for control in controls:
        control['lexical_mention_rows'] = hits[control['control_record_id']]
    summary = {
        'schema_version': 1, 'source_pages': {k: v[0] for k, v in SOURCES.items()},
        'source_input_sha256': {path.name: hashes[str(path)] for path in input_paths},
        'script_sha256': digest(Path(__file__).read_bytes()),
        'decoder_sha256': digest((HERE / 'broad_hub_screen.py').read_bytes()),
        'source_rows': dict(Counter(r['collection'] for r in rows)),
        'total_source_rows': len(rows), 'distinct_probe_labels': len(probes),
        'probe_labels_by_page': {s: len({identity(r['probe_label']) for r in rows if r['collection'] == s}) for s in SOURCES},
        'duplicate_labels_within_pages': {s: {n: c for n, c in Counter(r['probe_label'] for r in rows if r['collection'] == s).items() if c > 1} for s in SOURCES},
        'multiple_family_labels': [r['probe_label'] for r in probes if len(r['protein_families']) > 1],
        'donated_labels_absent_from_main': sorted({r['probe_label'] for r in rows if r['collection'] == 'donated'} - {r['probe_label'] for r in rows if r['collection'] == 'all'}),
        'control_records_scoped_to_parent': len(controls),
        'distinct_nonempty_control_cells': len({r['source_control_cell'] for r in controls}),
        'contextual_control_descriptions_not_screened': sum(not r['search_names'] for r in controls),
        'probe_labels_with_catalog_name_candidate': sum(bool(r['catalog_name_candidates']) for r in probes),
        'probe_labels_with_curated_appraisal_name_candidate': sum(bool(r['appraisal_name_candidates']) for r in probes),
        'source_rows_with_literal_BUB1B_target': sum(bool(re.search(r'\bBUB1B\b', r['target_text'], re.I)) for r in rows),
        'corpus_source_id_records': len(corpus), 'mention_rows': len(mentions),
        'mention_rows_by_role': dict(Counter(r['role'] for r in mentions)),
        'distinct_matched_source_ids': len({(r['source'], r['record_id']) for r in mentions}),
        'probe_labels_with_lexical_hit': sum(bool(r['lexical_mention_rows']) for r in probes),
        'manual_abstract_reading_credit': 0, 'new_curated_appraisal_credit': 0,
        'scope': 'Source-table enumeration and conservative name retrieval only. Probe labels and parent-scoped control cells are not structure-deduplicated molecules. No individual detail page, target selectivity, current stock, clinical safety or efficacy is validated by this layer. Blank fields and no matches remain unknown. Not a worldwide compound inventory.'}
    output.mkdir(parents=True, exist_ok=True)
    for name, data in [('source-rows.tsv', rows), ('probes.tsv', probes), ('controls.tsv', controls), ('mentions.tsv', mentions)]:
        write_tsv(output / name, data, list(data[0]) if data else EMPTY_SCHEMAS[name])
    for path in input_paths:
        if digest(path.read_bytes()) != hashes[str(path)]:
            raise ValueError(f'Input changed during run: {path}')
    summary['output_sha256'] = {name: digest((output / name).read_bytes()) for name in ['source-rows.tsv', 'probes.tsv', 'controls.tsv', 'mentions.tsv']}
    (output / 'summary.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False) + '\n')
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--main-html', required=True, type=Path)
    parser.add_argument('--donated-html', required=True, type=Path)
    parser.add_argument('--corpus', required=True, type=Path)
    parser.add_argument('--catalog-dir', default=HERE, type=Path)
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(build(args.main_html, args.donated_html, args.corpus,
                           args.catalog_dir, args.out), indent=2))


if __name__ == '__main__':
    main()
