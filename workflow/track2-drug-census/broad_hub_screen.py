#!/usr/bin/env python3
"""Offline Broad parser/crosswalk audit and preferred-name discovery screen.

Writes only to --output-dir. Never mutates raw inputs or the repository.
No alias expansion, fuzzy match, efficacy inference, or molecule deduplication.
"""
import argparse
import ast
import csv
import hashlib
import html
import io
import json
import re
import unicodedata
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUTPUT = HERE
EXPECTED_REPAIRS = {
    'BRD-K00005246-001-01-9': ('CK-101', 'RX-', '518'),
    'BRD-K00005246-001-02-9': ('CK-101', 'RX-', '518'),
    'BRD-K62685538-001-02-9': ('triptorelin', 'Luteinizing hormone-', 'releasing factor'),
}
ADDITIONAL_LOW_SPECIFICITY = {'propidium iodide','hematoxylin','alizarin','arachidonic acid','guanosine','hypoxanthine','nad','nadh','nadp','nadph','n6 methyladenosine','acrylate','antagonist g'}
GENE = re.compile(r'\b(?:bub1b|bubr1|bub r1|mosaic variegated aneuploidy)\b')
# Restrict removal to syntactically plausible, named HTML/JATS formatting tags.
# In particular, never treat biological inequalities as an arbitrary <...> tag.
FORMAT_TAG = re.compile(
    r"</?(?:abstract|title|p|sec|label|italic|bold|i|b|em|strong|u|sup|sub|span|div|br|xref|ext-link|styled-content|named-content|inline-formula|disp-formula|article-title|list|list-item|h[1-6]|table|tr|td|th|tbody|thead|footnote|fn|sc|monospace)"
    r"(?:\s+[A-Za-z_:][\w:.-]*(?:\s*=\s*(?:\"[^\"]*\"|'[^']*'|[^\s<>]+))?)*\s*/?>",
    flags=re.IGNORECASE,
)



def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(name, value):
    (OUTPUT / name).write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n')


def identity_normalize(text):
    # Exact independent reproduction of the saved crosswalk's imported function.
    text = unicodedata.normalize('NFKC', html.unescape(re.sub('<[^>]+>', ' ', text or ''))).casefold()
    return ' '.join(re.findall(r'[^\W_]+', text, flags=re.UNICODE))


def clean(text):
    # Decode before tag removal so encoded italics do not become literal i tokens.
    text = text or ''
    for _ in range(3):
        decoded = html.unescape(text)
        if decoded == text:
            break
        text = decoded
    return FORMAT_TAG.sub(' ', text)


def normalize(text):
    text = unicodedata.normalize('NFKC', clean(text)).casefold()
    return ' '.join(re.findall(r'[^\W_]+', text, flags=re.UNICODE))


def source_lines(path):
    return [(i, line) for i, line in enumerate(path.read_text().splitlines(), 1)
            if line and not line.startswith('!')]


def read_tsv(path, source=False):
    text = '\n'.join(line for _, line in source_lines(path)) if source else path.read_text()
    result = list(csv.DictReader(io.StringIO(text), delimiter='\t', strict=True))
    assert all(None not in row and all(v is not None for v in row.values()) for row in result)
    return result


def read_samples(path):
    lines = source_lines(path)
    header = lines[0][1].split('\t')
    assert len(header) == 12
    rebuilt = [lines[0][1]]
    manual = []
    repairs = []
    physical_widths = Counter(len(line.split('\t')) for _, line in lines[1:])
    quoted_cells = sum(cell.startswith('"') for _, line in lines[1:] for cell in line.split('\t'))
    i = 1
    while i < len(lines):
        line_number, line = lines[i]
        cells = line.split('\t')
        if len(cells) == 7:
            continuation_number, continuation = lines[i + 1]
            tail = continuation.split('\t')
            assert len(tail) == 6
            assert (cells[1], cells[6], tail[0]) == EXPECTED_REPAIRS[cells[0]]
            repairs.append({'broad_id': cells[0], 'pert_iname': cells[1],
                            'physical_line_numbers': [line_number, continuation_number],
                            'field': 'vendor_name', 'value': cells[6] + '\n' + tail[0]})
            cells = cells[:6] + [cells[6] + '\n' + tail[0]] + tail[1:]
            # Quote only the known unquoted newline field for the independent CSV parser.
            encoded = cells.copy()
            encoded[6] = '"' + encoded[6].replace('"', '""') + '"'
            rebuilt.append('\t'.join(encoded))
            i += 1
        else:
            assert len(cells) == 12, (line_number, len(cells))
            rebuilt.append(line)
        assert not cells[0] or cells[0].startswith('BRD-')
        # Reproduce the original cell-unquoting path only for independent comparison.
        manual.append(dict(zip(header, [v[1:-1].replace('""', '"')
                                        if v.startswith('"') and v.endswith('"') else v
                                        for v in cells])))
        i += 1
    result = list(csv.DictReader(io.StringIO('\n'.join(rebuilt)), delimiter='\t', strict=True))
    assert result == manual, 'Strict CSV parse differs from original split/unquote semantics'
    assert len(result) == 22612 and len(repairs) == 3
    assert {r['broad_id'] for r in repairs} == set(EXPECTED_REPAIRS)
    assert all(None not in row and all(v is not None for v in row.values()) for row in result)
    id_counts = Counter(r['broad_id'] for r in result if r['broad_id'])
    id_names = defaultdict(set)
    for row in result:
        if row['broad_id']:
            id_names[row['broad_id']].add(row['pert_iname'])
    return result, {'physical_data_width_counts': dict(physical_widths),
                    'quoted_source_cells': quoted_cells, 'distinct_nonblank_broad_ids': len(id_counts), 'repeated_nonblank_broad_ids': sum(n > 1 for n in id_counts.values()), 'broad_ids_shared_across_names': {k: sorted(v) for k, v in id_names.items() if len(v) > 1}, 'rows_with_blank_broad_id': sum(not r['broad_id'] for r in result), 'strict_csv_equals_original_cell_parser': True,
                    'repairs': repairs}


def policy_terms(path):
    # Read only the existing explicit constant. Do not execute/import repository code.
    tree = ast.parse(path.read_text())
    node = next(n.value for n in tree.body if isinstance(n, ast.Assign)
                and any(isinstance(t, ast.Name) and t.id == 'LOW_SPECIFICITY' for t in n.targets))
    assert isinstance(node, ast.BinOp) and isinstance(node.op, ast.BitOr)
    words = node.left.args[0].func.value.value.split()
    return set(words) | ast.literal_eval(node.right)


def classify(structure_ids, name_ids, hub_keys, catkeys):
    if len(structure_ids) == 1 and name_ids <= structure_ids:
        return 'SHARED_FULL_INCHIKEY_AND_NAME' if name_ids else 'SHARED_FULL_INCHIKEY_ONLY'
    if structure_ids and name_ids - structure_ids:
        return 'NAME_STRUCTURE_CONFLICT'
    if len(structure_ids) > 1:
        return 'MULTIPLE_STRUCTURE_CANDIDATES'
    if len(name_ids) > 1:
        return 'MULTIPLE_NAME_CANDIDATES'
    if name_ids:
        cid = next(iter(name_ids))
        return 'NAME_ONLY_DIFFERENT_FULL_KEYS' if hub_keys and catkeys[cid] else 'NAME_ONLY_STRUCTURE_UNAVAILABLE'
    return 'NO_MATCH'


def scan_field(tokens, trie):
    found = defaultdict(list)
    for start in range(len(tokens)):
        node, end = trie, start
        while end < len(tokens) and tokens[end] in node:
            node = node[tokens[end]]
            end += 1
            if '' in node:
                found[node['']].append([start, end])
    return found


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source-dir', type=Path, required=True)
    ap.add_argument('--catalog-dir', type=Path, required=True)
    ap.add_argument('--output-dir', type=Path, required=True)
    ap.add_argument('--corpus', type=Path, required=True)
    ap.add_argument('--title-abstract-corpus', type=Path)
    args = ap.parse_args()
    global OUTPUT
    OUTPUT = args.output_dir.resolve()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    paths = {'drug_source': args.source_dir/'repo-drug-annotation-20250818.txt',
             'sample_source': args.source_dir/'repo-sample-annotation-20250818.txt',
             'screen_script': Path(__file__),
             'catalog': args.catalog_dir/'catalog.tsv', 'aliases': args.catalog_dir/'aliases.tsv',
             'supplement': args.catalog_dir/'drugcentral-status-unverified.tsv',
             'literature_match_script': args.catalog_dir/'literature_match.py',
             'corpus': args.corpus}
    if args.title_abstract_corpus:
        paths['title_abstract_corpus'] = args.title_abstract_corpus
    for key, name in [('existing_crosswalk','crosswalk.json'),('existing_summary','summary.json'),('existing_repairs','sample-parse-repairs.json'),('existing_crosswalk_script','crosswalk.py')]:
        if (args.source_dir/name).exists():
            paths[key] = args.source_dir/name
    inputs = {k: {'file': p.name, 'sha256': sha(p)} for k, p in paths.items()}
    drugs = read_tsv(paths['drug_source'], True)
    samples, parser = read_samples(paths['sample_source'])
    catalog = read_tsv(paths['catalog'])
    aliases = read_tsv(paths['aliases'])
    supplement = read_tsv(paths['supplement'])
    previous = json.loads(paths['existing_crosswalk'].read_text())['rows'] if 'existing_crosswalk' in paths else None
    prior_summary = json.loads(paths['existing_summary'].read_text()) if 'existing_summary' in paths else None
    prior_repairs = json.loads(paths['existing_repairs'].read_text()) if 'existing_repairs' in paths else None
    assert len(drugs) == 7540
    if previous is not None:
        assert len(previous) == len(drugs)
    assert len({r['pert_iname'] for r in drugs}) == 7540
    if prior_repairs is not None:
        assert {r['broad_id'] for r in prior_repairs} == set(EXPECTED_REPAIRS)
        for old, new in zip(prior_repairs, parser['repairs']):
            assert old['pert_iname'] == new['pert_iname']
            assert old['source_fragment_before_break'] + '\n' + old['source_fragment_after_break'] == new['value']
    keys, names, sample_by_name = defaultdict(set), defaultdict(set), defaultdict(list)
    catkeys = {}
    for row in catalog:
        cid = row['catalog_id']
        catkeys[cid] = set(json.loads(row['inchikeys']))
        for key in catkeys[cid]:
            keys[key].add(cid)
    for row in aliases:
        assert row['catalog_id'] in catkeys
        names[identity_normalize(row['alias'])].add(row['catalog_id'])
    for row in supplement:
        cid = 'DC' + row['drugcentral_id']
        assert cid not in catkeys
        catkeys[cid] = {row['inchikey']} if row['inchikey'] else set()
        for key in catkeys[cid]:
            keys[key].add(cid)
        for name in json.loads(row['aliases']) + [row['name']]:
            names[identity_normalize(name)].add(cid)
    for row in samples:
        sample_by_name[row['pert_iname']].append(row)
    assert set(sample_by_name) <= {row['pert_iname'] for row in drugs}
    missing_sample_names = sorted({row['pert_iname'] for row in drugs} - set(sample_by_name))
    crosscheck = []
    for index, row in enumerate(drugs, 1):
        cid = f'HUB{index:05d}'
        group = sample_by_name.get(row['pert_iname'], [])
        hub_keys = {s['InChIKey'] for s in group if s['InChIKey']}
        structure_ids = set().union(*(keys[k] for k in hub_keys)) if hub_keys else set()
        name_ids = names[identity_normalize(row['pert_iname'])]
        status = classify(structure_ids, name_ids, hub_keys, catkeys)
        built = {'hub_row_id':cid, **row, 'sample_rows':len(group),
                 'sample_full_inchikeys':sorted(hub_keys), 'structure_candidate_ids':sorted(structure_ids),
                 'name_candidate_ids':sorted(name_ids), 'crosswalk_status':status}
        if previous is not None:
            old = previous[index-1]
            assert all(old[k] == v for k,v in built.items())
        crosscheck.append(built)
    statuses = dict(Counter(r['crosswalk_status'] for r in crosscheck))
    if prior_summary is not None:
        assert statuses == prior_summary['crosswalk_status']
        assert dict(Counter(r['clinical_phase'] for r in drugs)) == prior_summary['source_phases']
        for p in [paths['drug_source'], paths['sample_source']]:
            assert sha(p) == prior_summary['source_hashes'][p.name]
    normalized_to_ids, key_to_ids = defaultdict(list), defaultdict(set)
    for row in crosscheck:
        normalized_to_ids[normalize(row['pert_iname'])].append(row['hub_row_id'])
        for key in row['sample_full_inchikeys']:
            key_to_ids[key].add(row['hub_row_id'])
    shared_keys = {key: sorted(ids) for key, ids in key_to_ids.items() if len(ids) > 1}
    collisions = {key: ids for key, ids in normalized_to_ids.items() if len(ids) > 1}
    audit = {'inputs': inputs, 'parser_audit': parser,
             'drug_rows': len(drugs), 'distinct_source_names': len(drugs), 'sample_rows': len(samples),
             'actual_distinct_sample_names': len(sample_by_name), 'missing_sample_names': missing_sample_names,
             'summary_discrepancies': [{'field':'distinct_sample_names', 'reported':prior_summary['distinct_sample_names'] if prior_summary else None, 'audited':len(sample_by_name), 'cause':'defaultdict access inserted empty groups for 18 source names without sample rows'}] if prior_summary and prior_summary['distinct_sample_names'] != len(sample_by_name) else [],
             'catalog_rows': len(catalog), 'alias_rows': len(aliases), 'supplement_rows': len(supplement),
             'crosswalk_rows_independently_built': len(crosscheck), 'existing_crosswalk_rows_compared':len(crosscheck) if previous is not None else 0, 'crosswalk_status_counts': statuses,
             'source_rows_with_multiple_sample_full_keys': sum(len(r['sample_full_inchikeys']) > 1 for r in crosscheck),
             'shared_full_keys_across_source_names': shared_keys,
             'normalized_preferred_name_collisions': collisions,
             'limits': 'Cross-references are candidates, not deduplicated medicines or efficacy evidence.'}
    low_specificity = {normalize(x) for x in policy_terms(paths['literature_match_script']) | ADDITIONAL_LOW_SPECIFICITY}
    trie = {}
    first_token_candidates = defaultdict(list)
    for key in normalized_to_ids:
        first_token_candidates[key.split()[0]].append((key, key.split()))
        assert key
        node = trie
        for token in key.split():
            node = node.setdefault(token, {})
        node[''] = key
    # Meaningful matching invariants: whole-token boundary, punctuation equivalence,
    # HTML decoding before tag removal, and title/abstract field separation.
    toy = {'paclitaxel': {'': 'paclitaxel'}, 'azd': {'1152': {'': 'azd 1152'}}}
    assert not scan_field(normalize('depaclitaxel').split(), toy)
    assert scan_field(normalize('AZD-1152').split(), toy)['azd 1152'] == [[0, 2]]
    assert normalize('&lt;i&gt;BUB1B&lt;/i&gt;') == 'bub1b'
    assert normalize('p < 0.05 with paclitaxel and x > 2') == 'p 0 05 with paclitaxel and x 2'
    assert scan_field(normalize('p < 0.05 with paclitaxel and x > 2').split(),toy)['paclitaxel'] == [[4,5]]
    assert normalize('&lt;italic&gt;paclitaxel&lt;/italic&gt;') == 'paclitaxel'
    assert normalize('<span class="drug">paclitaxel</span>') == 'paclitaxel'
    assert not scan_field(['azd'], toy) and not scan_field(['1152'], toy)
    corpus = json.loads(paths['corpus'].read_text())
    ta = json.loads(paths['title_abstract_corpus'].read_text()) if 'title_abstract_corpus' in paths else None
    records = corpus['results']
    assert len(records) == corpus['retrievedCount'] == corpus['hitCount'] == 7362
    assert len({(r['source'], r['id']) for r in records}) == len(records)
    ta_ids = {(r['source'], r['id']) for r in ta['results']} if ta else set()
    mentions = []
    by_hub = defaultdict(list)
    abstract_count = 0
    direct_gene_count = 0
    normalized_difference_records = 0
    for record in records:
        fields = {f: normalize(record.get(k, '')) for f, k in [('title', 'title'), ('abstract', 'abstractText')]}
        abstract_count += bool(fields['abstract'])
        direct_gene_fields = [f for f, text in fields.items() if GENE.search(text)]
        direct_gene_count += bool(direct_gene_fields)
        normalized_difference_records += any(fields[f] != identity_normalize(record.get(k, ''))
                                            for f, k in [('title', 'title'), ('abstract', 'abstractText')])
        matches = defaultdict(dict)
        for field, text in fields.items():
            tokens = text.split()
            found = scan_field(tokens, trie)
            independent = defaultdict(list)
            for start, token in enumerate(tokens):
                for key, pattern in first_token_candidates.get(token, []):
                    if tokens[start:start+len(pattern)] == pattern:
                        independent[key].append([start,start+len(pattern)])
            assert found == independent, 'Independent whole-corpus token scan mismatch'
            for name, spans in found.items():
                matches[name][field] = spans
        record_hash = hashlib.sha256(json.dumps(record, sort_keys=True, ensure_ascii=False,
                                               separators=(',', ':')).encode()).hexdigest()
        for key, field_spans in sorted(matches.items()):
            for hid in normalized_to_ids[key]:
                mid = f'M{len(mentions)+1:06d}'
                mention = {'mention_id': mid, 'hub_row_id': hid, 'normalized_preferred_name': key,
                           'source': record['source'], 'record_id': record['id'],
                           'pmid': record.get('pmid'), 'pmcid': record.get('pmcid'),
                           'title': ' '.join(clean(record.get('title', '')).split()),
                           'source_url': 'https://europepmc.org/article/'+record['source']+'/'+record['id'],
                           'source_record_sha256': record_hash, 'matched_fields_token_spans': field_spans,
                           'direct_gene_fields': direct_gene_fields,
                           'member_of_gene_title_abstract_query': ((record['source'],record['id']) in ta_ids) if ta else None,
                           'abstract_available': bool(fields['abstract']),
                           'evidence_state': 'PREFERRED_NAME_LEXICAL_MENTION_ONLY_UNREVIEWED'}
                mentions.append(mention)
                by_hub[hid].append(mention)
    rows = []
    for old in crosscheck:
        hid, key = old['hub_row_id'], normalize(old['pert_iname'])
        lexical_flags = []
        letters = sum(c.isalpha() for c in key)
        if len(key.replace(' ', '')) < 5 or letters < 3:
            lexical_flags.append('short_or_code_ambiguity')
        if key in low_specificity or (len(key.split()) > 1 and key.split()[0] in ['d','l','dl']
                                     and ' '.join(key.split()[1:]) in low_specificity):
            lexical_flags.append('common_biological_chemical_or_assay_term')
        if len(normalized_to_ids[key]) > 1:
            lexical_flags.append('normalized_name_shared_by_multiple_source_rows')
        identity_flags = []
        if len(old['sample_full_inchikeys']) > 1:
            identity_flags.append('multiple_source_sample_full_keys')
        if any(not r['broad_id'] for r in sample_by_name.get(old['pert_iname'], [])):
            identity_flags.append('source_includes_structure_row_without_broad_sample_id')
        duplicate_ids = sorted(set().union(*(set(shared_keys.get(k, [])) for k in old['sample_full_inchikeys'])) - {hid})
        if duplicate_ids:
            identity_flags.append('full_key_also_occurs_under_other_source_names')
        if old['sample_rows'] == 0:
            identity_flags.append('no_source_sample_annotation_rows')
        if old['crosswalk_status'] == 'NO_MATCH':
            identity_flags.append('no_current_catalog_name_or_full_key_candidate_not_proof_of_novelty')
        elif not old['crosswalk_status'].startswith('SHARED_FULL_INCHIKEY'):
            identity_flags.append('catalog_identity_unresolved')
        source_ids = [r['broad_id'] for r in sample_by_name.get(old['pert_iname'], []) if r['broad_id']]
        if len(source_ids) > len(set(source_ids)):
            identity_flags.append('repeated_broad_sample_id_annotation_rows')
        matches = by_hub[hid]
        rows.append({'hub_row_id':hid, 'preferred_source_name':old['pert_iname'], 'normalized_preferred_name':key,
                     'clinical_phase_source_annotation':old['clinical_phase'], 'crosswalk_status':old['crosswalk_status'],
                     'name_candidate_ids':old['name_candidate_ids'], 'structure_candidate_ids':old['structure_candidate_ids'],
                     'sample_rows':old['sample_rows'], 'sample_full_inchikeys':old['sample_full_inchikeys'],
                     'normalized_name_peer_ids':[x for x in normalized_to_ids[key] if x != hid],
                     'shared_full_key_peer_ids':duplicate_ids, 'lexical_ambiguity_flags':lexical_flags,
                     'identity_uncertainty_flags':identity_flags,
                     'mention_policy':'AMBIGUOUS_PREFERRED_NAME_MANUAL_REVIEW' if lexical_flags else 'SPECIFIC_PREFERRED_NAME_DISCOVERY_ONLY',
                     'literature_record_count':len(matches),
                     'direct_gene_record_count':sum(bool(m['direct_gene_fields']) for m in matches),
                     'mention_ids':[m['mention_id'] for m in matches],
                     'evidence_state':'LEXICAL_DISCOVERY_ONLY_UNREVIEWED' if matches else 'UNKNOWN_NO_PREFERRED_NAME_MATCH'})
    unmatched = [r for r in rows if r['crosswalk_status']=='NO_MATCH' and r['literature_record_count']]
    summary = {'source_names_screened':len(rows), 'independent_match_verification':'All title/abstract fields checked by a separate first-token candidate matcher; exact match and span equality passed', 'normalized_names_screened':len(normalized_to_ids),
               'corpus_records_scanned':len(records), 'records_with_abstracts':abstract_count,
               'records_with_direct_gene_literal':direct_gene_count,
               'normalization_differs_from_old_html_order_records':normalized_difference_records,
               'matched_source_names':sum(bool(r['literature_record_count']) for r in rows),
               'unmatched_source_names':sum(not r['literature_record_count'] for r in rows),
               'lexical_mention_rows':len(mentions),
               'ambiguity_policy_counts':dict(Counter(r['mention_policy'] for r in rows)),
               'no_catalog_match_source_names':sum(r['crosswalk_status']=='NO_MATCH' for r in rows),
               'no_catalog_match_names_with_mentions':len(unmatched),
               'no_catalog_match_names_with_direct_gene_mentions':sum(bool(r['direct_gene_record_count']) for r in unmatched),
               'no_catalog_match_specific_names_with_direct_gene_mentions':sum(bool(r['direct_gene_record_count']) and not r['lexical_ambiguity_flags'] for r in unmatched),
               'identity_flag_counts':dict(Counter(flag for r in rows for flag in r['identity_uncertainty_flags'])),
               'audited_distinct_sample_names':len(sample_by_name),
               'source_annotation_rows':len(samples),
               'distinct_nonblank_broad_ids':parser['distinct_nonblank_broad_ids'],
               'sample_rows_without_broad_id':parser['rows_with_blank_broad_id'],
               'scope':'All 7540 preferred source names; no alias expansion, fuzzy matching, full texts, primary followups, or network.',
               'limits':'Mentions are lexical discovery only. NO_MATCH is an unresolved catalog cross-reference, not a new unique medicine. No text-match absence establishes inefficacy. Phase labels are dated source annotations, not verified present approval.'}
    result={'inputs':inputs, 'summary':summary,
            'method':{'normalization':'HTML entity decode up to three passes; strip only whitelisted syntactically plausible HTML/JATS formatting tags, preserving inequality text; Unicode NFKC/casefold; Unicode alphanumeric tokens excluding underscores; punctuation/hyphens become boundaries.',
                      'matching':'Exact contiguous normalized preferred-name tokens, independently within title and abstract. All names included, including short/generic ambiguous names; overlap retained.',
                      'span_convention':'Zero-based [start,end) token spans in the separately normalized field.',
                      'source_record_hash_serialization':'UTF-8 JSON, sorted keys, ensure_ascii=False, compact separators',
                      'ambiguity_terms_source':inputs['literature_match_script'],
                      'additional_explicit_low_specificity_terms':sorted(ADDITIONAL_LOW_SPECIFICITY),
                      'identity':'Retain all source rows; annotate key sharing and collisions without chemical merging.'},
            'rows':rows, 'mentions':mentions}
    assert len(rows)==7540 and len({r['hub_row_id'] for r in rows})==7540
    assert sum(r['literature_record_count'] for r in rows)==len(mentions)
    assert all(sha(paths[k])==v['sha256'] for k,v in inputs.items()), 'An input changed during the audit'
    # Compact row table for portable repository integration.
    fields = ['hub_row_id','preferred_source_name','normalized_preferred_name','clinical_phase_source_annotation','crosswalk_status','name_candidate_ids','structure_candidate_ids','sample_rows','sample_full_inchikeys','mention_policy','literature_record_count','direct_gene_record_count','lexical_ambiguity_flags','identity_uncertainty_flags','evidence_state']
    with (OUTPUT/'preferred-name-screen.tsv').open('w', newline='') as fh:
        writer=csv.DictWriter(fh, fieldnames=fields, delimiter='\t', lineterminator='\n'); writer.writeheader()
        for row in rows:
            writer.writerow({k:json.dumps(row[k],ensure_ascii=False) if isinstance(row[k],list) else row[k] for k in fields})
    mention_fields=['mention_id','hub_row_id','source','record_id','pmid','pmcid','title','source_url','source_record_sha256','matched_fields_token_spans','direct_gene_fields','member_of_gene_title_abstract_query','abstract_available','evidence_state']
    with (OUTPUT/'preferred-name-mentions.tsv').open('w',newline='') as fh:
        writer=csv.DictWriter(fh,fieldnames=mention_fields,delimiter='\t', lineterminator='\n');writer.writeheader()
        for row in mentions:
            writer.writerow({k:json.dumps(row[k],ensure_ascii=False) if isinstance(row[k],(dict,list)) else row[k] for k in mention_fields})
    summary['inputs']=inputs
    summary['source_metadata']={'official_hub_url':'https://repo-hub.broadinstitute.org/repurposing', 'source_banner_url':'http://www.broadinstitute.org/repurposing', 'drug_annotation_url':'https://repo-hub.broadinstitute.org/public/data/repo-drug-annotation-20250818.txt', 'sample_annotation_url':'https://repo-hub.broadinstitute.org/public/data/repo-sample-annotation-20250818.txt', 'source_file_date':'8/18/2025', 'raw_txt_license_notice':'The Drug Repurposing Hub data is provided for non-commercial use only.', 'portal_metadata_license':'CC-BY 4.0 (current portal, root verified 2026-09-21)', 'portal_use_restrictions':'Research use; not for clinical treatment decisions or commercial marketing.', 'license_notice_note':'Retain the dated raw-TXT and current-portal notices separately; do not silently substitute one for the other. This is a research discovery artifact.', 'citation_doi':'10.1038/nm.4306'}
    write_json('audit-results.json', audit)
    write_json('preferred-name-screen.json', result)
    write_json('preferred-name-screen-summary.json', summary)
    write_json('no-match-mentioned-names.json', {'inputs':inputs, 'limits':summary['limits'],
               'rows':sorted(unmatched,key=lambda r:(-r['direct_gene_record_count'],-r['literature_record_count'],r['preferred_source_name'])),
               'mentions':[m for m in mentions if any(r['hub_row_id']==m['hub_row_id'] for r in unmatched)]})
    print(json.dumps(summary,indent=2))

if __name__ == '__main__':
    main()
