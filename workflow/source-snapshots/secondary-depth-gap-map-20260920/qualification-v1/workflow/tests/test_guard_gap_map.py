import copy
import importlib.util
from pathlib import Path
import unittest

P=Path(__file__).resolve().parents[1]/'scripts/guard_gap_map.py'
spec=importlib.util.spec_from_file_location('gap_guard',P);G=importlib.util.module_from_spec(spec);spec.loader.exec_module(G)


class GuardControls(unittest.TestCase):
    def good(self):
        pins=G.sources()
        return dict(status='complete',fictional_only=True,tests_run=39,errors=0,failures=0,exit_code=0,skipped=[],source_files=pins),pins
    def test_wholly_passing_receipt(self):r,p=self.good();G.qualification_ok(r,p)
    def test_failed_receipt(self):
        r,p=self.good();r['failures']=1
        with self.assertRaises(ValueError):G.qualification_ok(r,p)
    def test_stale_source(self):
        r,p=self.good();r=copy.deepcopy(r);r['source_files'][G.SOURCES[0]]['sha256']='0'*64
        with self.assertRaises(ValueError):G.qualification_ok(r,p)
    def test_missing_source(self):
        r,p=self.good();r['source_files'].pop(G.SOURCES[0])
        with self.assertRaises(ValueError):G.qualification_ok(r,p)
    def test_boolean_count(self):
        r,p=self.good();r['tests_run']=True
        with self.assertRaises(ValueError):G.qualification_ok(r,p)
    def test_all_limit_boundaries(self):
        l=G.LIMITS
        self.assertIsNone(G.limit_reason(0,0,0,l['minimum_free_bytes'],l))
        for elapsed,rss,size,free,expected in [(l['wall_seconds'],0,0,l['minimum_free_bytes'],'wall_seconds'),
            (0,l['group_RSS_bytes'],0,l['minimum_free_bytes'],'group_RSS_bytes'),
            (0,0,l['output_bytes'],l['minimum_free_bytes'],'output_bytes'),
            (0,0,0,l['minimum_free_bytes']-1,'minimum_free_bytes')]:
            self.assertEqual(G.limit_reason(elapsed,rss,size,free,l),expected)
    def test_only_owned_group_rss(self):
        self.assertEqual(G.owned_rss([{'pgid':20,'RSS_bytes':5},{'pgid':21,'RSS_bytes':99}],20),5)


if __name__=='__main__':unittest.main()
