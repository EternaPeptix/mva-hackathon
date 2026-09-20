"""Fictional launcher refusals and owned-group isolation only."""
import importlib.util,json,os
from pathlib import Path
import subprocess,sys,tempfile,unittest
from unittest.mock import patch
s=importlib.util.spec_from_file_location('launcher',Path(__file__).with_name('guard-verification.py'));m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
class Controls(unittest.TestCase):
 def good(self):return dict(status='complete',tests_run=40,errors=0,failures=0,skipped=0)
 def test_success(self):m.complete_control(self.good(),40)
 def test_failed_or_skipped(self):
  for k in ['errors','failures','skipped']:
   q=self.good();q[k]=1
   with self.assertRaises(ValueError):m.complete_control(q,40)
 def test_wrong_count_or_bool(self):
  for n in [38,True]:
   q=self.good();q['tests_run']=n
   with self.assertRaises(ValueError):m.complete_control(q,40)
 def test_existing_outputs_refused(self):
  with tempfile.TemporaryDirectory() as d:
   out=Path(d)/'out'
   with patch.object(m,'OUTPUT',out):
    out.mkdir()
    with self.assertRaises(ValueError):m.fresh_output(out)
 def test_existing_launch_refused(self):
  with tempfile.TemporaryDirectory() as d:
   out=Path(d)/'out';out.with_name('out-launch').mkdir()
   with patch.object(m,'OUTPUT',out):
    with self.assertRaises(ValueError):m.fresh_output(out)
 def test_other_output_refused(self):
  with self.assertRaises(ValueError):m.fresh_output(Path('/tmp/unowned'))
 def test_all_finite_limits(self):
  v=[0,0,0,m.LIMITS['minimum_free_bytes']]
  self.assertIsNone(m.g.limit_reason(*v,m.LIMITS))
  for i,k in enumerate(['wall_seconds','group_RSS_bytes','output_bytes']):
   x=v.copy();x[i]=m.LIMITS[k];self.assertEqual(m.g.limit_reason(*x,m.LIMITS),k)
  v[3]-=1;self.assertEqual(m.g.limit_reason(*v,m.LIMITS),'minimum_free_bytes')
 def test_termination_only_own_session(self):
  unrelated=subprocess.Popen([sys.executable,'-c','import time;time.sleep(10)'],start_new_session=True)
  owned=subprocess.Popen([sys.executable,'-c','import time;time.sleep(10)'],start_new_session=True)
  try:
   self.assertNotEqual(owned.pid,os.getpgrp());m.g.stop_owned(owned)
   self.assertIsNotNone(owned.poll());self.assertIsNone(unrelated.poll())
  finally:
   if owned.poll() is None:m.g.stop_owned(owned)
   if unrelated.poll() is None:unrelated.terminate();unrelated.wait(timeout=5)
if __name__=='__main__':unittest.main()
