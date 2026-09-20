import importlib.util
from pathlib import Path
import tempfile
import unittest

import pysam

spec = importlib.util.spec_from_file_location('phase_audit', Path(__file__).parents[1]/'scripts/audit_candidate_phase.py')
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


class CandidatePhaseTests(unittest.TestCase):
    def test_nonoverlapping_deletions_and_reference_check(self):
        result = audit.haplotypes('CAGCTAAGTT', 20, [(22, 'AG', 'A'), (25, 'TAAG', 'T')])
        self.assertEqual(result['CACTTT'], ['1|1'])
        self.assertEqual(result['CAGCTAAGTT'], ['0|0'])
        with self.assertRaises(ValueError):
            audit.haplotypes('CAGCTAAGTT', 20, [(22, 'TG', 'T'), (25, 'TAAG', 'T')])

    def test_fragment_counts_phase_overlap_duplicates_and_quality(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            reference = root/'ref.fa'
            reference.write_text('>1\n' + 'C'*20 + 'CAGCTAAGTT' + 'C'*100 + '\n')
            pysam.faidx(str(reference))
            header = {'HD': {'VN': '1.6'}, 'SQ': [{'SN': '1', 'LN': 130}], 'RG': [{'ID': 'lane1'}]}
            bam = root/'test.bam'
            # Exact cis, reference, and trans alternatives; duplicate/low-quality observations excluded.
            observations = [('ref', 'CAGCTAAGTT', '10M', 0, 40),
                            ('cis', 'CACTTT', '2M1D2M3D2M', 0, 40),
                            ('cis', 'CACTTT', '2M1D2M3D2M', 16, 40),
                            ('first', 'CACTAAGTT', '2M1D7M', 0, 40),
                            ('second', 'CAGCTTT', '5M3D2M', 0, 40),
                            ('dup', 'CACTTT', '2M1D2M3D2M', 1024, 40),
                            ('lowq', 'CACTTT', '2M1D2M3D2M', 0, 10)]
            with pysam.AlignmentFile(str(bam), 'wb', header=header) as stream:
                for name, sequence, cigar, flag, quality in observations:
                    read = pysam.AlignedSegment(stream.header)
                    read.query_name, read.query_sequence = name, sequence
                    read.flag, read.reference_id, read.reference_start = flag, 0, 20
                    read.mapping_quality, read.cigarstring = 60, cigar
                    read.query_qualities = [quality]*len(sequence)
                    read.set_tag('RG', 'lane1')
                    stream.write(read)
            pysam.index(str(bam))
            evidence = audit.count_haplotypes(bam, reference, '1', [(22,'AG','A'), (25,'TAAG','T')])
            self.assertEqual(evidence['counts'], {'0|0':1,'1|0':1,'0|1':1,'1|1':1,'other':0,'ambiguous':0,'mate_conflict':0})
            self.assertEqual(evidence['strand_counts']['1|1']['both'], 1)


if __name__ == '__main__':
    unittest.main()
