#!/usr/bin/env python3
"""Display exact saved abstract strings in small, source-bound reading batches.

No HTML stripping, word substitution, summarization or biological selection is
performed. The caller must check the end-of-batch marker and avoid truncating
the output. Raw article text belongs in the external source cache, not in Git.
"""
import argparse
import hashlib
import json
from pathlib import Path


def sha(data):
    return hashlib.sha256(data).hexdigest()


def select(payload, identifiers):
    if not 1 <= len(identifiers) <= 8:
        raise ValueError('Request one to eight source:ID records per batch')
    if len(set(identifiers)) != len(identifiers):
        raise ValueError('Duplicate requested source:ID')
    wanted = set(identifiers)
    found = {}
    for record in payload['results']:
        key = record['source'] + ':' + record['id']
        if key in wanted:
            if key in found:
                raise ValueError('Duplicate source record: ' + key)
            found[key] = record
    if wanted != set(found):
        raise ValueError('Missing source records: ' + ', '.join(sorted(wanted - set(found))))
    result = []
    for key in identifiers:
        record = found[key]
        value = record.get('abstractText')
        if not isinstance(value, str) or not value:
            raise ValueError('No saved abstract: ' + key)
        result.append({
            'source_id': key, 'title': record['title'],
            'abstract_characters': len(value),
            'abstract_sha256': sha(value.encode()),
            'abstractText': value,
        })
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--corpus', type=Path, required=True)
    parser.add_argument('--expected-sha256')
    parser.add_argument('records', nargs='+', help='Source:ID, e.g. MED:27829147')
    args = parser.parse_args()
    raw = args.corpus.read_bytes()
    digest = sha(raw)
    if args.expected_sha256 and args.expected_sha256 != digest:
        raise ValueError('Source snapshot hash mismatch')
    rows = select(json.loads(raw), args.records)
    print(json.dumps({'corpus_sha256': digest, 'requested_records': len(rows)}))
    for row in rows:
        print(json.dumps(row, ensure_ascii=False))
    print(json.dumps({'end_of_batch': True, 'emitted_records': len(rows),
                      'total_abstract_characters': sum(r['abstract_characters'] for r in rows)}))


if __name__ == '__main__':
    main()
