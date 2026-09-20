import importlib.util
from pathlib import Path
import unittest

spec = importlib.util.spec_from_file_location('ranked_phase', Path(__file__).parents[1]/'scripts/audit_ranked_phase.py')
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


class Sample(dict):
    phased = True


class RankedPhaseTests(unittest.TestCase):
    def test_explicit_blocks_and_opposite_alt_orientation(self):
        a = audit.explicit_phase({'GT': (0,1), 'PGT':'0|1', 'PID':'10_A_G'}, 1, '1')
        b = audit.explicit_phase({'GT': (1,0), 'PGT':'1|0', 'PID':'10_A_G'}, 1, '1')
        self.assertEqual(audit.pair_relation(a,a),'cis')
        self.assertEqual(audit.pair_relation(a,b),'trans')
        self.assertEqual(audit.pair_relation(a,('1','PID','20_A_G',1)),'unknown')
        self.assertEqual(audit.pair_relation(a,('2','PID','10_A_G',1)),'unknown')

    def test_invalid_unphased_and_homozygous_evidence_stays_unknown(self):
        self.assertIsNone(audit.explicit_phase({'GT':(0,1),'PGT':'0|2','PID':'x'},1,'1'))
        self.assertIsNone(audit.explicit_phase({'GT':(1,1),'PGT':'1|1','PID':'x'},1,'1'))
        self.assertIsNone(audit.explicit_phase({'GT':(0,1),'PGT':'0/1','PID':'x'},1,'1'))
        self.assertIsNone(audit.explicit_phase(Sample(GT=(0,1)),1,'1'))
        self.assertEqual(audit.explicit_phase(Sample(GT=(1,0),PS=123),1,'1'),('1','PS','123',0))
        self.assertEqual(audit.pair_relation(None,None),'unknown')


if __name__ == '__main__':
    unittest.main()
