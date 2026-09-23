#!/usr/bin/env python3
"""Read-only comparison of saved normalizers under the same saved alias policy."""
import ast
import argparse
import csv
import hashlib
import html
import json
import re
import unicodedata
from collections import Counter
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--corpus', type=Path, required=True)
parser.add_argument('--catalog-dir', type=Path, default=Path(__file__).resolve().parent)
parser.add_argument('--out', type=Path, required=True)
args = parser.parse_args()
REPO, RAW, OUT = args.catalog_dir, args.corpus, args.out
OUT.mkdir(parents=True, exist_ok=True)
INPUTS = [RAW, REPO/'literature_match.py', REPO/'broad_hub_screen.py',
          REPO/'literature-alias-policy.tsv', REPO/'literature-mentions.tsv']

def sha(data):
    return hashlib.sha256(data).hexdigest()

def read_tsv(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter='\t'))

def function_scope(path, names, constants=()):
    """Compile only named function definitions/constants; no module imports/main."""
    tree = ast.parse(path.read_text())
    nodes = [n for n in tree.body
             if (isinstance(n, ast.FunctionDef) and n.name in names)
             or (isinstance(n, ast.Assign) and any(isinstance(t, ast.Name)
                 and t.id in constants for t in n.targets))]
    scope = dict(re=re, html=html, unicodedata=unicodedata)
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(path), 'exec'), scope)
    return scope

def make_trie(policies):
    trie = {}
    for policy in policies:
        key = policy['normalized_alias']
        if not key or sum(c.isalpha() for c in key) < 2 or len(key.replace(' ', '')) < 3:
            continue
        node = trie
        for token in key.split():
            node = node.setdefault(token, {})
        node[''] = key
    return trie

def matches(text, trie):
    tokens, found = text.split(), set()
    for i in range(len(tokens)):
        node, j = trie, i
        while j < len(tokens) and tokens[j] in node:
            node = node[tokens[j]]
            j += 1
            if '' in node:
                found.add(node[''])
    return found

def identity(r):
    return dict(source=r['source'], record_id=r['id'], pmid=r.get('pmid'),
                pmcid=r.get('pmcid'), doi=r.get('doi'))

hash_before = {p.name: sha(p.read_bytes()) for p in INPUTS}
old = function_scope(REPO/'literature_match.py', {'normalize'})['normalize']
newscope = function_scope(REPO/'broad_hub_screen.py', {'normalize', 'clean'}, {'FORMAT_TAG'})
new, clean, format_tag = newscope['normalize'], newscope['clean'], newscope['FORMAT_TAG']
def restricted_tags_only(text):
    # Ablation: retain legacy decode order and one html.unescape pass;
    # change only arbitrary <...> removal to the named-tag expression.
    text = unicodedata.normalize('NFKC', html.unescape(format_tag.sub(' ', text or ''))).casefold()
    return ' '.join(re.findall(r'[^\W_]+', text, flags=re.UNICODE))
policies = read_tsv(REPO/'literature-alias-policy.tsv')
specific = [p for p in policies if p['policy'] == 'specific_alias_comention_only']
policy_by_alias = {p['normalized_alias']:p for p in specific}
assert len(policy_by_alias) == len(specific)
trie = make_trie(specific)
records = json.loads(RAW.read_bytes())['results']
counts = Counter()
changed, deltas, spans = [], [], []
old_rows, new_rows, old_alias_rows, new_alias_rows = set(), set(), set(), set()
title_changed_ids = []
gene_deltas = []
gene = re.compile(r'\b(bub1b|bubr1|bub r1|mosaic variegated aneuploidy)\b')

for r in records:
    a, title = r.get('abstractText', ''), r.get('title', '')
    oa, na, ot, nt = old(a), new(a), old(title), new(title)
    counts['source_records'] += 1
    counts['nonempty_raw_abstracts'] += bool(a)
    counts['missing_raw_abstracts'] += not bool(a)
    if ot != nt:
        title_changed_ids.append(identity(r))
    if oa != na:
        changed.append({**identity(r), 'raw_abstract_sha256':sha(a.encode()),
                        'legacy_normalized_sha256':sha(oa.encode()),
                        'new_normalized_sha256':sha(na.encode()),
                        'legacy_token_count':len(oa.split()), 'new_token_count':len(na.split())})
    removed = []
    for m in re.finditer('<[^>]+>', a or ''):
        if not format_tag.fullmatch(m.group()):
            removed.append({'start_character':m.start(), 'end_character':m.end(),
                            'character_count':len(m.group()), 'span_sha256':sha(m.group().encode()),
                            'new_decoder_retained_token_count':len(new(m.group()).split()),
                            'starts_as_numeric_inequality':bool(re.match(r'<\s*[-+\d.=]', m.group()))})
    if removed:
        spans.append({**identity(r), 'legacy_removed_nonformat_spans':removed})
    af, bf = matches(oa, trie), matches(na, trie)
    full_old = matches(ot+' '+oa, trie)
    abstract_change_only = matches(ot+' '+na, trie)
    full_new = matches(nt+' '+na, trie)
    restricted_only = matches(ot+' '+restricted_tags_only(a), trie)
    if af != bf or full_old != full_new or full_old != abstract_change_only:
        def details(keys):
            return [dict(normalized_alias=k, catalog_ids=json.loads(policy_by_alias[k]['catalog_ids'])) for k in sorted(keys)]
        deltas.append({**identity(r),
                       'abstract_alias_gains':details(bf-af), 'abstract_alias_losses':details(af-bf),
                       'combined_alias_gains_holding_title_legacy':details(abstract_change_only-full_old),
                       'combined_alias_losses_holding_title_legacy':details(full_old-abstract_change_only),
                       'combined_alias_gains_both_fields_new':details(full_new-full_old),
                       'combined_alias_losses_both_fields_new':details(full_old-full_new),
                       'generic_regex_only_repair_reproduces_combined_alias_result':restricted_only==full_new})
    counts['legacy_direct_gene_records'] += bool(gene.search(ot+' '+oa))
    counts['new_direct_gene_records'] += bool(gene.search(nt+' '+na))
    if bool(gene.search(ot+' '+oa)) != bool(gene.search(nt+' '+na)):
        gene_deltas.append({**identity(r), 'legacy_match':bool(gene.search(ot+' '+oa)),
                            'new_match':bool(gene.search(nt+' '+na))})
    counts['all_record_matches_equal_under_regex_only_repair_and_full_new'] += restricted_only == full_new
    for found, rowset, aliasrows in [(full_old, old_rows, old_alias_rows), (full_new, new_rows, new_alias_rows)]:
        for alias in found:
            for gid in json.loads(policy_by_alias[alias]['catalog_ids']):
                rowset.add((r['source'], r['id'], gid))
                aliasrows.add((r['source'], r['id'], gid, alias))

saved = read_tsv(REPO/'literature-mentions.tsv')
saved_rows = {(r['source'],r['record_id'],r['catalog_id']) for r in saved}
saved_alias_rows = {(r['source'],r['record_id'],r['catalog_id'],a) for r in saved
                    for a in json.loads(r['matched_normalized_aliases'])}
hash_after = {p.name: sha(p.read_bytes()) for p in INPUTS}
assert hash_before == hash_after, 'An input changed during the audit; rerun required.'
counts.update(normalized_abstracts_changed=len(changed), normalized_titles_changed=len(title_changed_ids),
              abstracts_with_nonformat_spans_removed_by_legacy=len(spans),
              nonformat_spans_removed_by_legacy=sum(len(r['legacy_removed_nonformat_spans']) for r in spans),
              legacy_nonformat_spans_containing_retained_tokens=sum(s['new_decoder_retained_token_count']>0 for r in spans for s in r['legacy_removed_nonformat_spans']),
              records_with_retained_tokens_in_legacy_nonformat_spans=sum(any(s['new_decoder_retained_token_count']>0 for s in r['legacy_removed_nonformat_spans']) for r in spans),
              new_retained_tokens_inside_legacy_removed_nonformat_spans=sum(s['new_decoder_retained_token_count'] for r in spans for s in r['legacy_removed_nonformat_spans']),
              records_with_specific_alias_delta=len(deltas),
              legacy_specific_catalog_record_rows=len(old_rows), new_specific_catalog_record_rows=len(new_rows),
              specific_catalog_record_rows_gained=len(new_rows-old_rows), specific_catalog_record_rows_lost=len(old_rows-new_rows),
              specific_catalog_record_alias_rows_gained=len(new_alias_rows-old_alias_rows), specific_catalog_record_alias_rows_lost=len(old_alias_rows-new_alias_rows))
result = dict(scope='All saved source-record abstract slots; same saved specific-alias policy. Identity-level outputs only; no raw abstracts.',
              input_sha256=hash_before, input_hashes_unchanged_during_audit=True,
              script_sha256=sha(Path(__file__).read_bytes()),
              policy=dict(rows=len(policies), specific_aliases=len(specific), ambiguous_aliases_excluded=len(policies)-len(specific)),
              counts=dict(counts),
              baseline_reproduction=dict(catalog_record_rows_equal_saved=old_rows==saved_rows,
                                         catalog_record_alias_rows_equal_saved=old_alias_rows==saved_alias_rows,
                                         computed_only_rows=sorted(old_rows-saved_rows), saved_only_rows=sorted(saved_rows-old_rows)),
              specific_alias_deltas=deltas,
              direct_gene_flag_deltas=gene_deltas,
              gained_catalog_record_rows=sorted(new_rows-old_rows), lost_catalog_record_rows=sorted(old_rows-new_rows),
              changed_abstract_identities=changed, changed_title_identities=title_changed_ids,
              legacy_nonformat_removal_identities=spans,
              limitations=['This is a reproducible decoder and lexical-match comparison, not efficacy or drug-intervention qualification.',
                           'The cause of the separate 3,598-record manual annotation errors remains UNKNOWN; this comparison does not prove that those readers used this decoder.',
                           'No historical counts, repository files, catalogs, aliases or saved appraisals were rebuilt or changed.'])
(OUT/'decoder-alias-details.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
compact = {k:v for k,v in result.items() if k not in ['changed_abstract_identities','changed_title_identities','legacy_nonformat_removal_identities','specific_alias_deltas']}
compact['specific_alias_deltas'] = [{**{k:r[k] for k in ['source','record_id','pmid','doi']},
    'gained':r['combined_alias_gains_both_fields_new'], 'lost':r['combined_alias_losses_both_fields_new'],
    'same_with_legacy_titles': r['combined_alias_gains_holding_title_legacy']==r['combined_alias_gains_both_fields_new'] and r['combined_alias_losses_holding_title_legacy']==r['combined_alias_losses_both_fields_new'],
    'generic_regex_only_repair_reproduces_result':r['generic_regex_only_repair_reproduces_combined_alias_result']} for r in deltas]
compact['detail_file'] = 'decoder-alias-details.json'
(OUT/'decoder-alias-delta.json').write_text(json.dumps(compact,indent=2,ensure_ascii=False)+'\n')
print(json.dumps({k:result[k] for k in ['counts','baseline_reproduction','specific_alias_deltas']},indent=2))
