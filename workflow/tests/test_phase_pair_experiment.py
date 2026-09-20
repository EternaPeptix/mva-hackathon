import importlib.util
from pathlib import Path
import unittest

spec = importlib.util.spec_from_file_location('phase_experiment', Path(__file__).parents[1]/'scripts/phase_pair_experiment.py')
experiment = importlib.util.module_from_spec(spec)
spec.loader.exec_module(experiment)


class PhasePairExperimentTests(unittest.TestCase):
    def test_only_target_gt_changes(self):
        line = '3\t10046723\t.\tAG\tA\t99\tPASS\t.\tGT:DP:PGT:PID\t0/1:89:0|1:10046720_C_T\n'
        changed, key = experiment.phase_line(line)
        self.assertEqual(key, ('3','10046723','AG','A'))
        self.assertEqual(changed, line.replace('0/1:89','0|1:89'))
        unrelated = line.replace('10046723','10046724')
        self.assertEqual(experiment.phase_line(unrelated), (unrelated,None))
        with self.assertRaises(ValueError):
            experiment.phase_line(line.replace('0|1:', '1|0:'))
        with self.assertRaises(ValueError):
            experiment.phase_line(line.replace('10046720_C_T', '10046719_A_G'))


if __name__ == '__main__':
    unittest.main()
