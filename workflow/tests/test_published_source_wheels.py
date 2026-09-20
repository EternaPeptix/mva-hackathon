"""Archive boundary checks for the source-build reproduction path."""
import hashlib
import io
from pathlib import Path
import sys
import tarfile
import tempfile
import unittest
import zipfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from rebuild_published_python_wheels import clean_source, wheel_info


class ArchiveChecks(unittest.TestCase):
    def test_top_level_metadata_with_bundled_dependency(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'example-1.0-py3-none-any.whl'
            with zipfile.ZipFile(path, 'w') as archive:
                archive.writestr('example-1.0.dist-info/METADATA', 'Name: example\nVersion: 1.0\n')
                archive.writestr('example-1.0.dist-info/WHEEL', 'Wheel-Version: 1.0\n')
                archive.writestr('example/_vendor/other-2.dist-info/METADATA', 'Name: other\nVersion: 2\n')
                archive.writestr('example/_vendor/other-2.dist-info/WHEEL', 'Wheel-Version: 1.0\n')
            self.assertEqual(wheel_info(path)['name'], 'example')
            self.assertEqual(wheel_info(path)['version'], '1.0')

    def test_two_top_level_metadata_roots_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'ambiguous.whl'
            with zipfile.ZipFile(path, 'w') as archive:
                for name in ('one', 'two'):
                    archive.writestr(name + '-1.dist-info/METADATA', 'Name: ' + name + '\nVersion: 1\n')
                    archive.writestr(name + '-1.dist-info/WHEEL', 'Wheel-Version: 1.0\n')
            with self.assertRaisesRegex(ValueError, 'metadata count'):
                wheel_info(path)

    def test_unsafe_wheel_path_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'unsafe.whl'
            with zipfile.ZipFile(path, 'w') as archive:
                archive.writestr('../escape', 'untrusted member')
            with self.assertRaisesRegex(ValueError, 'unsafe wheel member'):
                wheel_info(path)

    def test_prebuilt_members_omitted_only_from_build_copy(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); path = root / 'example.tar.gz'
            with tarfile.open(path, 'w:gz') as archive:
                for name, data in [('example/code.c', b'int main(void) { return 0; }'),
                                   ('example/tool', b'\x7fELFcompiled'),
                                   ('example/extension.so', b'placeholder prebuilt')]:
                    item = tarfile.TarInfo(name); item.size = len(data)
                    archive.addfile(item, io.BytesIO(data))
            before = hashlib.sha256(path.read_bytes()).hexdigest()
            source, record = clean_source(path, root / 'copy')
            self.assertEqual((source / 'code.c').read_bytes(), b'int main(void) { return 0; }')
            self.assertFalse((source / 'tool').exists())
            self.assertFalse((source / 'extension.so').exists())
            self.assertEqual(record['retained_source_files'], 1)
            self.assertEqual(len(record['omitted_prebuilt_members']), 2)
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), before)

    def test_source_links_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); path = root / 'linked.tar.gz'
            with tarfile.open(path, 'w:gz') as archive:
                item = tarfile.TarInfo('example/link'); item.type = tarfile.SYMTYPE; item.linkname = '/outside'
                archive.addfile(item)
            with self.assertRaisesRegex(ValueError, 'link or special'):
                clean_source(path, root / 'copy')


if __name__ == '__main__':
    unittest.main()
