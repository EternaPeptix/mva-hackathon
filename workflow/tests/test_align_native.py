"""Run the real native tools on synthetic reads, including partial-lane resumption."""
import gzip
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest


class NativeAlignmentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bin = Path(os.environ.get('MVA_BIO_BIN', str(Path(sys.executable).parent)))
        if not all((cls.bin / name).is_file() for name in ('bwa', 'samtools')):
            raise unittest.SkipTest('native tools are installed in the Studio runtime')
        import pysam
        cls.pysam = pysam
        cls.script = Path(__file__).resolve().parents[1] / 'scripts/01_align.sh'

    def test_partial_lanes_finalize_once_and_preserve_inputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rng = random.Random(9818)
            sequence = ''.join(rng.choice('ACGT') for _ in range(2500))
            reference = root/'synthetic.fa'
            reference.write_text('>1\n' + sequence + '\n')
            self.pysam.faidx(str(reference))
            subprocess.run([str(self.bin/'bwa'), 'index', str(reference)], check=True, capture_output=True)
            data, output = root/'data', root/'align'
            data.mkdir()

            def add_lane(lane):
                for mate in (1, 2):
                    with gzip.open(data/('SYNTH_L%03d_R%d_001.fastq.gz' % (lane, mate)), 'wt') as fastq:
                        for index in range(20):
                            start = 100 + index * 60
                            read = sequence[start:start+100] if mate == 1 else sequence[start+200:start+300].translate(str.maketrans('ACGT','TGCA'))[::-1]
                            fastq.write('@SIM:1:FLOWCELL:%d:1101:%d:1000 %d:N:0:ACGT\n%s\n+\n%s\n' % (lane, index, mate, read, 'I'*len(read)))

            env = dict(os.environ, PATH=str(self.bin)+os.pathsep+os.environ.get('PATH',''), MVA_PYTHON=sys.executable, SORT_MEM='32M')
            command = ['bash',str(self.script),'-r',str(reference),'-o',str(output),'-s','SYNTH','-f',str(data),'-t','2','-n','2']
            add_lane(1)
            first = subprocess.run(command, env=env, capture_output=True, text=True)
            self.assertEqual(first.returncode, 0, first.stdout + first.stderr)
            self.assertTrue((output/'SYNTH.L001.bam').is_file())
            self.assertFalse((output/'SYNTH.cram').exists())
            add_lane(2)
            second = subprocess.run(command, env=env, capture_output=True, text=True)
            self.assertEqual(second.returncode, 0, second.stdout + second.stderr)
            cram = output/'SYNTH.cram'
            self.assertTrue(cram.is_file())
            self.assertTrue((output/'SYNTH.cram.crai').is_file())
            self.assertTrue((output/'SYNTH.L001.bam').is_file())
            self.assertTrue((output/'SYNTH.L002.bam').is_file())
            self.assertEqual(len(list(output.glob('SYNTH.merged.*.bam'))), 1)
            with self.pysam.AlignmentFile(str(cram), reference_filename=str(reference)) as reads:
                self.assertEqual(sum(1 for _ in reads.fetch(until_eof=True)), 80)
            modified = cram.stat().st_mtime_ns
            repeat = subprocess.run(command, env=env, capture_output=True, text=True)
            self.assertEqual(repeat.returncode, 0, repeat.stdout + repeat.stderr)
            self.assertEqual(cram.stat().st_mtime_ns, modified)
            reference.touch()
            changed = subprocess.run(command, env=env, capture_output=True, text=True)
            self.assertNotEqual(changed.returncode, 0)
            self.assertEqual(cram.stat().st_mtime_ns, modified)


if __name__ == '__main__':
    unittest.main()
