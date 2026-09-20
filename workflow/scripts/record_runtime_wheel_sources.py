#!/usr/bin/env python3
"""Reconcile captured wheel identities with the public PyPI version metadata."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import datetime
import hashlib
import json
from pathlib import Path
from urllib.parse import urlsplit
from urllib.request import Request, urlopen


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--wheelhouse', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Existing source audit preserved')
    receipt = json.loads((args.wheelhouse / 'status.json').read_text())
    if receipt['status'] != 'complete':
        parser.error('Completed wheel preparation required')
    args.output.mkdir(parents=True)

    def inspect(wheel):
        url = 'https://pypi.org/pypi/{}/{}/json'.format(wheel['name'], wheel['version'])
        with urlopen(Request(url, headers={'User-Agent': 'MVA-runtime-reproduction/1.0'}), timeout=30) as stream:
            raw = stream.read(8_000_001)
        if len(raw) > 8_000_000:
            raise ValueError('Unexpectedly large package metadata response')
        data = json.loads(raw)
        metadata_file = wheel['name'] + '-' + wheel['version'] + '.pypi.json'
        (args.output / metadata_file).write_bytes(raw)
        locally_built = wheel['source'].startswith('pre-existing')
        if locally_built:
            origin = wheel['source_origin']
            expected = origin['archive_info']['hashes']['sha256']
            filename = origin['url'].rsplit('/', 1)[-1].removeprefix('local-archive:')
            matches = [row for row in data['urls'] if row['filename'] == filename and row['digests']['sha256'] == expected]
        else:
            expected = wheel['sha256']
            matches = [row for row in data['urls'] if row['filename'] == wheel['filename'] and row['digests']['sha256'] == expected]
        if len(matches) != 1:
            raise ValueError('Publisher identity did not match exactly once: ' + wheel['name'])
        artifact = matches[0]
        source = urlsplit(artifact['url'])
        if source.scheme != 'https' or source.hostname != 'files.pythonhosted.org' or source.query or source.username:
            raise ValueError('Unexpected publisher artifact origin')
        return {'name': wheel['name'], 'version': wheel['version'], 'local_wheel': wheel['filename'],
                'local_wheel_sha256': wheel['sha256'], 'locally_built_wheel': locally_built,
                'publisher_match_scope': 'original source archive identity' if locally_built else 'complete downloaded wheel bytes',
                'publisher_filename': artifact['filename'], 'publisher_url': artifact['url'],
                'publisher_sha256': expected, 'publisher_size': artifact['size'],
                'pypi_metadata_url': url, 'pypi_metadata_sha256': hashlib.sha256(raw).hexdigest(),
                'wheel_hash_is_publisher_wheel_hash': not locally_built}

    with ThreadPoolExecutor(max_workers=4) as pool:
        sources = list(pool.map(inspect, receipt['wheels']))
    result = {'verified_at': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'wheelhouse_receipt_sha256': hashlib.sha256((args.wheelhouse / 'status.json').read_bytes()).hexdigest(),
              'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'publisher_wheels': sum(not r['locally_built_wheel'] for r in sources),
              'locally_built_wheels_with_publisher_source_identity': sum(r['locally_built_wheel'] for r in sources),
              'sources': sources,
              'limits': ['Five cached builds retain locally measured wheel hashes, not publisher wheel signatures.',
                         'Rebuilding their source archives with a different compiler or build environment can change wheel hashes.',
                         'Publisher metadata equality records file identity; it is not a security or clinical validation.']}
    (args.output / 'summary.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ('publisher_wheels', 'locally_built_wheels_with_publisher_source_identity')}))


if __name__ == '__main__':
    main()
