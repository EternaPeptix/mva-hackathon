#!/usr/bin/env python3
"""Resume one download under an OS lock and verify it before publication."""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
import gzip


class DownloadError(Exception):
    pass


class PrivateRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        result = super().redirect_request(req, fp, code, msg, headers, newurl)
        if result is not None and urllib.parse.urlparse(req.full_url).netloc != urllib.parse.urlparse(newurl).netloc:
            result.remove_header('Authorization')
        return result


def digest(path, expected):
    algorithm = 'sha256' if len(expected) == 64 else 'sha1'
    h = hashlib.new(algorithm)
    with path.open('rb') as source:
        for block in iter(lambda: source.read(8 * 1024 * 1024), b''):
            h.update(block)
    if h.hexdigest().lower() != expected.lower():
        raise DownloadError('publisher checksum mismatch; existing bytes preserved')
    return h.hexdigest()


def verify(path, name, expected_size, expected_hash):
    if expected_size is not None and path.stat().st_size != expected_size:
        raise DownloadError('file size differs from publisher metadata')
    if expected_hash:
        return digest(path, expected_hash)
    if name.endswith('.zip'):
        with zipfile.ZipFile(path) as archive:
            if archive.testzip() is not None:
                raise DownloadError('archive CRC verification failed')
    elif name.endswith(('.gz', '.tgz')):
        with gzip.open(path, 'rb') as source:
            while source.read(8 * 1024 * 1024):
                pass
    return None


def fetch(url, dest, expected_hash=None, expected_size=None, retries=5):
    name = Path(urllib.parse.urlparse(url).path).name
    if not name or name in ('.', '..'):
        raise DownloadError('URL must name a file')
    if expected_hash and not re.fullmatch(r'(?:[0-9a-fA-F]{40}|[0-9a-fA-F]{64})', expected_hash):
        raise DownloadError('expected checksum must be SHA-1 or SHA-256')
    dest = Path(dest)
    dest.mkdir(parents=True, exist_ok=True)
    final = dest / name
    part = dest / (name + '.part')
    receipt = dest / (name + '.verified.json')
    lock = open(str(final) + '.lock', 'a+')
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        lock.close()
        raise DownloadError('another process holds the download lock')
    try:
        headers = {'User-Agent': 'mva-reproducible-fetch/1', 'Accept-Encoding': 'identity'}
        extra = os.environ.get('FETCH_HEADER', '')
        if extra:
            key, sep, value = extra.partition(':')
            if not sep:
                raise DownloadError('FETCH_HEADER must contain a header name and value')
            headers[key.strip()] = value.strip()
        elif urllib.parse.urlparse(url).hostname == 'huggingface.co':
            token = os.environ.get('HF_TOKEN')
            token_file = Path(os.environ.get('HF_TOKEN_FILE', str(Path.home() / '.cache/huggingface/token')))
            if not token and token_file.is_file():
                token = token_file.read_text().strip()
            if token:
                headers['Authorization'] = 'Bearer ' + token
        opener = urllib.request.build_opener(PrivateRedirect())
        metadata = {}
        try:
            with opener.open(urllib.request.Request(url, headers=headers, method='HEAD'), timeout=45) as response:
                metadata = {'content_length': response.headers.get('Content-Length'), 'etag': response.headers.get('ETag'), 'last_modified': response.headers.get('Last-Modified')}
                if expected_size is None and metadata['content_length']:
                    expected_size = int(metadata['content_length'])
        except urllib.error.HTTPError as error:
            if error.code not in (405, 501):
                raise DownloadError('metadata request returned HTTP %d' % error.code)
        if final.exists():
            print('[verify] existing ' + name, flush=True)
            actual_hash = verify(final, name, expected_size, expected_hash)
        else:
            for attempt in range(retries + 1):
                offset = part.stat().st_size if part.exists() else 0
                if expected_size is not None and offset > expected_size:
                    raise DownloadError('partial file exceeds publisher size; preserved for inspection')
                if expected_size is not None and offset == expected_size:
                    break
                request_headers = dict(headers)
                if offset:
                    request_headers['Range'] = 'bytes=%d-' % offset
                print('[get] %s resume=%d expected=%s attempt=%d' % (name, offset, expected_size, attempt + 1), flush=True)
                try:
                    with opener.open(urllib.request.Request(url, headers=request_headers), timeout=60) as response:
                        if offset:
                            content_range = response.headers.get('Content-Range', '')
                            match = re.fullmatch(r'bytes (\d+)-(\d+)/(\d+|\*)', content_range)
                            if response.status != 206 or not match or int(match.group(1)) != offset:
                                raise DownloadError('server did not honor requested range; partial preserved')
                            if expected_size is not None and match.group(3) != '*' and int(match.group(3)) != expected_size:
                                raise DownloadError('range total differs from publisher size')
                        elif response.status != 200:
                            raise DownloadError('unexpected initial HTTP response')
                        started = last = time.monotonic()
                        transferred = 0
                        with part.open('ab' if offset else 'wb') as target:
                            while True:
                                block = response.read(4 * 1024 * 1024)
                                if not block:
                                    break
                                if expected_size is not None and offset + transferred + len(block) > expected_size:
                                    raise DownloadError('server sent more data than publisher size')
                                target.write(block)
                                transferred += len(block)
                                if time.monotonic() - last >= 20:
                                    print('[progress] %s bytes=%d MB_per_second=%.2f' % (name, offset + transferred, transferred / (time.monotonic() - started) / 1e6), flush=True)
                                    last = time.monotonic()
                            target.flush()
                            os.fsync(target.fileno())
                    if expected_size is None or part.stat().st_size == expected_size:
                        break
                    raise OSError('response ended before the expected file size')
                except (urllib.error.URLError, TimeoutError, OSError) as error:
                    if isinstance(error, urllib.error.HTTPError) and error.code in (401, 403, 404):
                        raise DownloadError('download returned HTTP %d' % error.code)
                    if attempt == retries:
                        raise DownloadError('transfer failed after retries; partial preserved')
                    print('[retry] transient transfer error; retaining downloaded bytes', flush=True)
                    time.sleep(min(2 ** attempt, 20))
            print('[verify] downloaded ' + name, flush=True)
            actual_hash = verify(part, name, expected_size, expected_hash)
            os.replace(part, final)
        record = {'filename': name, 'verified_at': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), 'size': final.stat().st_size, 'expected_size': expected_size, 'publisher_checksum': expected_hash, 'verified_checksum': actual_hash, 'archive_integrity_checked': not bool(expected_hash) and name.endswith(('.zip', '.gz', '.tgz')), 'source': url, 'source_metadata': metadata}
        receipt_tmp = Path(str(receipt) + '.tmp')
        receipt_tmp.write_text(json.dumps(record, indent=2) + '\n')
        os.replace(receipt_tmp, receipt)
        print('[done] verified ' + name, flush=True)
        return final
    finally:
        lock.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('url')
    parser.add_argument('destination')
    parser.add_argument('checksum', nargs='?')
    parser.add_argument('size', nargs='?', type=int)
    args = parser.parse_args()
    try:
        fetch(args.url, args.destination, args.checksum, args.size)
    except DownloadError as error:
        print('[FAIL] ' + str(error), file=sys.stderr)
        return 1
    except Exception as error:
        print('[FAIL] %s; files preserved' % type(error).__name__, file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
