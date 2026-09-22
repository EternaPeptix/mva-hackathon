#!/usr/bin/env python3
"""Reproduce the gene-only literature retrieval used by the 2026-09-21 census.

Reconstructed from the original API calls; retained snapshots remain the
authoritative dated record. Live indexes can change even for a fixed date.
Requires Python 3 and requests. Writes abstracts outside the repository.
"""
import argparse
import datetime
import json
import time
from pathlib import Path

import requests

GENE_QUERY = 'BUB1B OR BUBR1 OR "mosaic variegated aneuploidy"'
TITLE_QUERY = '(TITLE_ABS:BUBR1 OR TITLE_ABS:BUB1B OR TITLE_ABS:"BUB-R1" OR TITLE_ABS:"mosaic variegated aneuploidy")'
EPMC = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search'
NCBI = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'


def get(session, url, params, json_mode=True):
    for attempt in range(6):
        try:
            response = session.get(url, params=params, timeout=60)
            response.raise_for_status()
            if not json_mode:
                if b'<PubmedArticleSet' not in response.content:
                    raise ValueError('Unexpected PubMed XML response')
                return response.content
            result = response.json()
            if url == EPMC and 'resultList' not in result:
                raise ValueError('Europe PMC response has no resultList')
            return result
        except (requests.RequestException, ValueError):
            if attempt == 5:
                raise
            time.sleep(min(2 ** attempt, 16))


def europe_pmc(session, query):
    cursor = '*'
    records = []
    seen = set()
    count = None
    while True:
        result = get(session, EPMC, dict(query=query, resultType='core',
                     format='json', pageSize=1000, cursorMark=cursor))
        if count is None:
            count = result['hitCount']
        elif result['hitCount'] != count:
            raise RuntimeError('Index hitCount changed during retrieval; retry')
        for record in result['resultList']['result']:
            key = (record['source'], record['id'])
            if key in seen:
                raise RuntimeError(f'Duplicate Europe PMC record: {key}')
            seen.add(key)
            records.append(record)
        print(f'Europe PMC: {len(records)}/{count}', flush=True)
        if len(records) == count:
            break
        next_cursor = result.get('nextCursorMark')
        if not next_cursor or next_cursor == cursor:
            raise RuntimeError('Pagination stopped before all records arrived')
        cursor = next_cursor
    return dict(searchedDate=datetime.date.today().isoformat(), query=query,
                hitCount=count, retrievedCount=len(records), results=records)


def validate(path):
    for name in ('literature-gene-title-abstract.json', 'literature-gene-unrestricted.json'):
        result = json.loads((path / name).read_text())
        rows = result['results']
        unique = len({(r['source'], r['id']) for r in rows})
        assert len(rows) == unique == result['hitCount'], name
        print(name, 'complete', len(rows), 'abstracts',
              sum(bool(r.get('abstractText')) for r in rows))
    result = json.loads((path / 'literature-pubmed-esearch.json').read_text())['esearchresult']
    assert len(set(result['idlist'])) == len(result['idlist']) == int(result['count'])
    print('PubMed complete', result['count'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True,
                        help='External output directory; refuses to overwrite snapshots')
    parser.add_argument('--date', default='2026-09-21', help='PubMed publication upper bound')
    parser.add_argument('--verify-only', action='store_true', help='Validate existing files without network')
    args = parser.parse_args()
    if args.verify_only:
        validate(args.out)
        return
    args.out.mkdir(parents=True, exist_ok=True)
    names = ('literature-gene-title-abstract.json', 'literature-gene-unrestricted.json',
             'literature-pubmed-esearch.json', 'literature-pubmed-additional.xml')
    if any((args.out / n).exists() for n in names):
        raise SystemExit('Choose a new --out directory to preserve prior snapshots')
    with requests.Session() as session:
        session.headers['User-Agent'] = 'BUB1B-literature-census/2026-09-21'
        narrow = europe_pmc(session, TITLE_QUERY)
        (args.out / names[0]).write_text(json.dumps(narrow, ensure_ascii=False, indent=2))
        broad = europe_pmc(session, GENE_QUERY)
        (args.out / names[1]).write_text(json.dumps(broad, ensure_ascii=False, indent=2))
        params = dict(db='pubmed', term=GENE_QUERY, retmode='json', retmax=9999,
                      datetype='pdat', mindate='1800/01/01', maxdate=args.date.replace('-', '/'))
        pubmed = get(session, NCBI + 'esearch.fcgi', params)
        (args.out / names[2]).write_text(json.dumps(pubmed, ensure_ascii=False, indent=2))
        known = {r['id'] for r in narrow['results'] if r['source'] == 'MED'}
        additional = [p for p in pubmed['esearchresult']['idlist'] if p not in known]
        if additional:
            time.sleep(0.4)
            data = get(session, NCBI + 'efetch.fcgi',
                       dict(db='pubmed', id=','.join(additional), retmode='xml'), False)
            (args.out / names[3]).write_bytes(data)
        print('Additional PubMed records', len(additional))
        (args.out / 'retrieval-recipe.json').write_text(json.dumps(dict(
            retrievedUTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),
            europePMC=dict(queries=[TITLE_QUERY, GENE_QUERY], resultType='core',
                           pageSize=1000, dateFilter=None), pubMed=params), indent=2))
    validate(args.out)


if __name__ == '__main__':
    main()
