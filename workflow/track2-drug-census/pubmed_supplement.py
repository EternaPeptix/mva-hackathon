#!/usr/bin/env python3
"""Reconcile and screen supplemental PubMed records with the saved alias policy.

Uses external abstract caches; exports metadata and retrieval hits only. A preprint
and its journal update are version-linked, never counted as independent evidence.
"""
import argparse
import csv
import hashlib
import json
import re
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path
from catalog import stamp, write_json, write_tsv
from literature_match import normalize

ROOT = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def parse_articles(path):
    records, metadata = [], []
    for article in ET.parse(path).getroot().findall('PubmedArticle'):
        citation = article.find('MedlineCitation')
        pmid = citation.findtext('PMID')
        title_node = citation.find('Article/ArticleTitle')
        title = ''.join(title_node.itertext()) if title_node is not None else ''
        abstract = ' '.join(''.join(node.itertext()) for node in
                            citation.findall('Article/Abstract/AbstractText'))
        ids = {node.get('IdType'): node.text for node in
               article.findall('PubmedData/ArticleIdList/ArticleId')}
        types = [node.text for node in citation.findall('Article/PublicationTypeList/PublicationType')]
        updates = [node.findtext('PMID') for node in citation.findall('CommentsCorrectionsList/CommentsCorrections')
                   if node.get('RefType') == 'UpdateIn' and node.findtext('PMID')]
        records.append(dict(source='MED', id=pmid, pmid=pmid, title=title,
                            abstractText=abstract, doi=ids.get('doi', ''), pmcid=ids.get('pmc', '')))
        metadata.append(dict(pmid=pmid, title=title, doi=ids.get('doi', ''),
                             publication_types=types, update_in_pmids=updates,
                             abstract_available=bool(abstract),
                             abstract_sha256=hashlib.sha256(abstract.encode()).hexdigest(),
                             source_url='https://pubmed.ncbi.nlm.nih.gov/'+pmid+'/'))
    if not records or len({r['id'] for r in records}) != len(records):
        raise ValueError('Empty or duplicate PubMed article records')
    return records, metadata


def run(corpus, pubmed_xml, pubmed_search, title_abstract, out):
    out.mkdir(parents=True, exist_ok=True)
    broad = json.loads(corpus.read_text())
    by_pmid = {str(r.get('pmid') or r['id']): r for r in broad['results'] if r['source'] == 'MED'}
    records, metadata = parse_articles(pubmed_xml)
    search = json.loads(pubmed_search.read_text())['esearchresult']
    searched_ids = set(search['idlist'])
    assert len(searched_ids) == len(search['idlist']) == int(search['count'])
    fetched_ids = {r['id'] for r in records}
    if not fetched_ids <= searched_ids:
        raise ValueError('Supplement contains records outside the saved PubMed query')
    if not searched_ids <= set(by_pmid) | fetched_ids:
        raise ValueError('Saved PubMed search IDs are not fully represented by the two sources')
    for row, record in zip(metadata, records):
        same_record = row['pmid'] in by_pmid
        known_updates = [p for p in row['update_in_pmids'] if p in by_pmid]
        row['same_pmid_in_base'] = same_record
        row['journal_update_in_base'] = known_updates
        row['record_relation'] = ('SAME_RECORD_IN_BASE' if same_record else
                                  'PREPRINT_WITH_JOURNAL_UPDATE_IN_BASE' if 'Preprint' in row['publication_types'] and known_updates else
                                  'ADDITIONAL_PUBMED_RECORD')
        row['adds_previously_missing_abstract'] = bool(same_record and record['abstractText'] and not by_pmid[row['pmid']].get('abstractText'))
        row['gene_literal_in_title_abstract'] = bool(re.search(r'\b(bub1b|bubr1|bub r1|mosaic variegated aneuploidy)\b',normalize(record['title']+' '+record['abstractText'])))
        row['evidence_status'] = 'SOURCE_RECORD_ONLY_NOT_INDEPENDENT_RESCUE_EVIDENCE'
    # Execute the existing matcher unchanged: the same catalog, aliases and ambiguity
    # policy apply. Temporary full abstract material stays outside the repository.
    with tempfile.TemporaryDirectory(prefix='mva-pubmed-screen-') as folder:
        work = Path(folder)
        for name in ('catalog.tsv', 'aliases.tsv', 'drugcentral-status-unverified.tsv'):
            (work/name).symlink_to(ROOT/name)
        input_path = work/'corpus.json'
        write_json(input_path,dict(query='Supplemental PubMed records from saved independent gene query',
                                   hitCount=len(records),results=records))
        result = subprocess.run([sys.executable, str(ROOT/'literature_match.py'), '--corpus', str(input_path),
                                 '--title-abstract-corpus', str(title_abstract), '--out', str(work)],
                                capture_output=True, text=True, check=True)
        specific = read_tsv(work/'literature-mentions.tsv')
        ambiguous = read_tsv(work/'literature-ambiguous-mentions.tsv')
        relation = {r['pmid']:r for r in metadata}
        for rows, name in ((specific,'pubmed-supplement-mentions.tsv'),
                           (ambiguous,'pubmed-supplement-ambiguous.tsv')):
            for row in rows:
                meta = relation[row['record_id']]
                row['record_relation'] = meta['record_relation']
                row['journal_update_in_base'] = meta['journal_update_in_base']
            fields = list(rows[0]) if rows else ['catalog_id','record_id','record_relation']
            write_tsv(out/name,rows,fields)
        match_summary = json.loads((work/'literature-coverage.json').read_text())
    counts = Counter(r['record_relation'] for r in metadata)
    unique_new = fetched_ids-set(by_pmid)
    summary = dict(built_utc=stamp(), input_sha256={str(p.name):digest(p) for p in
                   (corpus,pubmed_xml,pubmed_search,title_abstract)},
                   matcher_sha256=digest(ROOT/'literature_match.py'),
                   catalog_input_sha256={name:digest(ROOT/name) for name in
                     ('catalog.tsv','aliases.tsv','drugcentral-status-unverified.tsv')},
                   base_records=len(broad['results']), pubmed_query_records=len(searched_ids),
                   fetched_supplement_records=len(records), record_relations=dict(counts),
                   additional_record_ids=sorted(unique_new), union_record_count=len(broad['results'])+len(unique_new),
                   all_saved_pubmed_search_ids_accounted_for=True,
                   missing_abstracts_in_supplement=sum(not r['abstract_available'] for r in metadata),
                   newly_available_abstracts=sum(r['adds_previously_missing_abstract'] for r in metadata),
                   specific_alias_rows=len(specific), ambiguous_alias_rows=len(ambiguous),
                   specific_alias_entities=len({r['catalog_id'] for r in specific}),
                   additional_record_specific_rows=sum(r['record_id'] in unique_new for r in specific),
                   catalog_entities_screened=match_summary['catalog_entities_screened'],
                   interpretation='Additional publication records are not additional independent studies. XML UpdateIn links preserve preprint/journal relationships. All supplemental abstracts are screened with the existing conservative alias policy. Co-mentions and missing evidence establish neither efficacy nor inefficacy; no corpus-wide full-text review is claimed.')
    write_tsv(out/'pubmed-supplement-records.tsv',metadata,list(metadata[0]))
    write_json(out/'pubmed-supplement-summary.json',summary)
    print(json.dumps(summary,indent=2))


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--corpus',type=Path,required=True)
    parser.add_argument('--pubmed-xml',type=Path,required=True)
    parser.add_argument('--pubmed-search',type=Path,required=True)
    parser.add_argument('--title-abstract-corpus',type=Path,required=True)
    parser.add_argument('--out',type=Path,default=ROOT)
    args=parser.parse_args()
    run(args.corpus,args.pubmed_xml,args.pubmed_search,args.title_abstract_corpus,args.out)
