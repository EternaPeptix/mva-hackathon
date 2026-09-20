import importlib.util
from pathlib import Path
import unittest

spec=importlib.util.spec_from_file_location('resume_calling',Path(__file__).parents[1]/'scripts/resume_calling.py')
resume=importlib.util.module_from_spec(spec)
spec.loader.exec_module(resume)


class CallingDependencyTests(unittest.TestCase):
    def test_live_dependency_waits_and_missing_handle_never_restarts(self):
        self.assertEqual(resume.dependency_action({'status':'running'},True),'wait_for_existing_alignment')
        with self.assertRaises(RuntimeError): resume.dependency_action({'status':'running'},False)
        with self.assertRaises(RuntimeError): resume.dependency_action({'status':'failed'},True)

    def test_terminal_dependency_needs_output_provenance(self):
        with self.assertRaises(ValueError): resume.dependency_action({'status':'complete'},False)
        complete={'status':'complete','cram_path':'alignment.cram','cram_bytes':100,'cram_sha256':'abc'}
        self.assertEqual(resume.dependency_action(complete,False),'verify_completed_alignment')


if __name__=='__main__': unittest.main()
