"""Check fragment counting and quality/duplicate filters with synthetic alignments."""
import importlib.util
from pathlib import Path
import tempfile
import unittest


class CandidateReadTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if importlib.util.find_spec('pysam') is None:
            raise unittest.SkipTest('requires the Studio pysam runtime')
        import pysam
        cls.pysam = pysam
        spec = importlib.util.spec_from_file_location('candidate_reads', Path(__file__).resolve().parents[1] / 'scripts/candidate_read_validation.py')
        cls.audit = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.audit)

    def test_counts_fragments_and_excludes_flagged_or_low_quality_reads(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            reference = root/'synthetic.fa'
            # Non-repetitive context keeps BAQ from converting this into an ambiguous alignment.
            sequence = ('ACGTGTCAGTGCATCGATGCTAGCATACGTCAGATCGTAGCTACGATCGTACGATCGACTGATCGTACGATC' * 4)
            reference.write_text('>1\n' + sequence + '\n')
            self.pysam.faidx(str(reference))
            position = 121
            ref = sequence[position-1]
            alt = next(base for base in 'ACGT' if base != ref)
            header = {'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'1','LN':len(sequence)}]}
            bam = root/'synthetic.bam'
            with self.pysam.AlignmentFile(str(bam), 'wb', header=header) as output:
                items = [('ref',ref,0,60,40),('alt-forward',alt,0,60,40),('alt-reverse',alt,16,60,40),
                         ('overlap',alt,0,60,40),('overlap',alt,128,60,40),
                         ('duplicate',alt,1024,60,40),('supplementary',alt,2048,60,40),
                         ('low-mapq',alt,0,10,40),('low-baseq',alt,0,60,10)]
                for name,base,flag,mapq,quality in items:
                    read = self.pysam.AlignedSegment()
                    read.query_name = name
                    bases = list(sequence[100:150]); bases[position-101] = base
                    read.query_sequence = ''.join(bases)
                    read.flag = 1 | 2 | (128 if flag & 128 else 64) | flag
                    read.reference_id = 0; read.reference_start = 100
                    read.mapping_quality = mapq; read.cigar = [(0,50)]
                    read.next_reference_id = 0; read.next_reference_start = 100
                    read.template_length = 50
                    read.query_qualities = [quality]*50
                    output.write(read)
            self.pysam.index(str(bam))
            result = self.audit.pileup_evidence(bam,reference,[('1',position,ref,alt)])[0]
            self.assertEqual(result['ref_fragments'],1)
            self.assertEqual(result['alt_fragments'],3)
            self.assertEqual(result['filtered_fragment_depth'],4)
            self.assertEqual(result['alt_reverse'],1)
            self.assertEqual(result['alt_forward'],2)


if __name__ == '__main__':
    unittest.main()
