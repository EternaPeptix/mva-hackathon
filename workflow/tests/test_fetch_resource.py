"""Exercise interrupted/corrupt downloads using synthetic data only."""
from contextlib import contextmanager
import fcntl
import hashlib
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import importlib.util
import io
import json
import os
from pathlib import Path
import tempfile
import threading
import unittest
from unittest.mock import patch
import zipfile


SPEC = importlib.util.spec_from_file_location(
    'fetch_resource', Path(__file__).resolve().parents[1] / 'scripts/fetch_resource.py')
FETCH = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(FETCH)
DATA = b'synthetic sequence data, not a human genome\n' * 200
SHA = hashlib.sha256(DATA).hexdigest()


@contextmanager
def server(payload=DATA, ignore_range=False):
    requests = []

    class Handler(BaseHTTPRequestHandler):
        def log_message(self, *args):
            pass

        def respond(self, head=False):
            requests.append((self.command, self.path, dict(self.headers)))
            if self.path.startswith('/redirect/'):
                self.send_response(302)
                self.send_header('Location', 'http://localhost:%d/file.bin' % self.server.server_port)
                self.end_headers()
                return
            offset = 0
            if self.headers.get('Range') and not ignore_range:
                offset = int(self.headers['Range'].split('=')[1].split('-')[0])
                self.send_response(206)
                self.send_header('Content-Range', 'bytes %d-%d/%d' % (offset, len(payload) - 1, len(payload)))
            else:
                self.send_response(200)
            self.send_header('Content-Length', str(len(payload) - offset))
            self.end_headers()
            if not head:
                self.wfile.write(payload[offset:])

        def do_HEAD(self):
            self.respond(head=True)

        def do_GET(self):
            self.respond()

    httpd = ThreadingHTTPServer(('127.0.0.1', 0), Handler)
    worker = threading.Thread(target=httpd.serve_forever, daemon=True)
    worker.start()
    try:
        yield 'http://127.0.0.1:%d' % httpd.server_port, requests
    finally:
        httpd.shutdown()
        httpd.server_close()
        worker.join()


class FetchTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.dest = Path(self.tmp.name)

    def test_range_resume_verifies_checksum_and_publishes_receipt(self):
        (self.dest/'file.bin.part').write_bytes(DATA[:1234])
        with server() as (url, requests):
            output = FETCH.fetch(url + '/file.bin', self.dest, SHA, len(DATA))
        self.assertEqual(output.read_bytes(), DATA)
        self.assertFalse((self.dest/'file.bin.part').exists())
        self.assertEqual([headers.get('Range') for method, _, headers in requests if method == 'GET'], ['bytes=1234-'])
        receipt = json.loads((self.dest/'file.bin.verified.json').read_text())
        self.assertEqual(receipt['verified_checksum'], SHA)

    def test_same_size_bad_existing_file_is_preserved_and_rejected(self):
        bad = b'x' * len(DATA)
        (self.dest/'file.bin').write_bytes(bad)
        with server() as (url, requests), self.assertRaises(FETCH.DownloadError):
            FETCH.fetch(url + '/file.bin', self.dest, SHA, len(DATA))
        self.assertEqual((self.dest/'file.bin').read_bytes(), bad)
        self.assertFalse((self.dest/'file.bin.verified.json').exists())
        self.assertFalse(any(method == 'GET' for method, _, _ in requests))

    def test_bad_prefix_never_becomes_final(self):
        (self.dest/'file.bin.part').write_bytes(b'bad prefix')
        with server() as (url, _), self.assertRaises(FETCH.DownloadError):
            FETCH.fetch(url + '/file.bin', self.dest, SHA, len(DATA))
        self.assertFalse((self.dest/'file.bin').exists())
        self.assertTrue((self.dest/'file.bin.part').read_bytes().startswith(b'bad prefix'))

    def test_ignored_range_does_not_append_a_second_full_file(self):
        prefix = DATA[:1234]
        (self.dest/'file.bin.part').write_bytes(prefix)
        with server(ignore_range=True) as (url, _), self.assertRaises(FETCH.DownloadError):
            FETCH.fetch(url + '/file.bin', self.dest, SHA, len(DATA))
        self.assertEqual((self.dest/'file.bin.part').read_bytes(), prefix)
        self.assertFalse((self.dest/'file.bin').exists())

    def test_concurrent_writer_lock_is_rejected(self):
        with (self.dest/'file.bin.lock').open('a+') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            with self.assertRaisesRegex(FETCH.DownloadError, 'download lock'):
                FETCH.fetch('http://127.0.0.1:1/file.bin', self.dest, SHA, len(DATA))

    def test_authentication_is_removed_on_host_redirect(self):
        with server() as (url, requests), patch.dict(os.environ, {'FETCH_HEADER': 'Authorization: Bearer synthetic-test-credential'}):
            FETCH.fetch(url + '/redirect/file.bin', self.dest, SHA, len(DATA))
        for _, path, headers in requests:
            if path == '/file.bin':
                self.assertNotIn('Authorization', headers)

    def test_zip_crc_checked_when_publisher_has_no_checksum(self):
        data = io.BytesIO()
        with zipfile.ZipFile(data, 'w') as archive:
            archive.writestr('synthetic.txt', DATA)
        payload = data.getvalue()
        with server(payload=payload) as (url, _):
            FETCH.fetch(url + '/data.zip', self.dest, expected_size=len(payload))
        receipt = json.loads((self.dest/'data.zip.verified.json').read_text())
        self.assertTrue(receipt['archive_integrity_checked'])

    def test_truncated_zip_is_not_published(self):
        with server(payload=b'PK incomplete archive') as (url, _), self.assertRaises(zipfile.BadZipFile):
            FETCH.fetch(url + '/data.zip', self.dest)
        self.assertTrue((self.dest/'data.zip.part').exists())
        self.assertFalse((self.dest/'data.zip').exists())


if __name__ == '__main__':
    unittest.main()
