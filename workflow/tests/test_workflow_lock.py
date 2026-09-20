"""The command must retain its writer lock even if the supervising wrapper exits."""
import fcntl
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest


class WorkflowLockTests(unittest.TestCase):
    def test_child_retains_lock_if_wrapper_is_terminated(self):
        script = Path(__file__).resolve().parents[1] / 'scripts/with_lock.py'
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            lock, ready, done = root/'writer.lock', root/'ready', root/'done'
            child = 'from pathlib import Path; import time; p=Path(%r); p.write_text("ready"); time.sleep(1.0); Path(%r).write_text("done")' % (str(ready), str(done))
            wrapper = subprocess.Popen([sys.executable,str(script),str(lock),sys.executable,'-c',child])
            try:
                deadline = time.monotonic()+5
                while not ready.exists() and time.monotonic()<deadline:
                    time.sleep(0.01)
                self.assertTrue(ready.exists())
                wrapper.terminate()
                wrapper.wait(timeout=5)
                with lock.open('a+') as file:
                    with self.assertRaises(BlockingIOError):
                        fcntl.flock(file,fcntl.LOCK_EX|fcntl.LOCK_NB)
                deadline = time.monotonic()+5
                while not done.exists() and time.monotonic()<deadline:
                    time.sleep(0.02)
                self.assertTrue(done.exists())
                with lock.open('a+') as file:
                    deadline = time.monotonic()+2
                    while True:
                        try:
                            fcntl.flock(file,fcntl.LOCK_EX|fcntl.LOCK_NB)
                            break
                        except BlockingIOError:
                            if time.monotonic()>deadline:
                                raise
                            time.sleep(0.01)
            finally:
                if wrapper.poll() is None:
                    wrapper.terminate()
                    wrapper.wait(timeout=5)


if __name__ == '__main__':
    unittest.main()
