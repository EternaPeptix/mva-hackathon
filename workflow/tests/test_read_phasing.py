"""Native phase truth and guards against invented linkage or changed genotypes."""
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest

import pysam

sys.path.insert(0, str(Path(__file__).parents[1] / 'scripts'))
import run_read_phasing as audit


HEADER = ('##fileformat=VCFv4.2\n##contig=<ID=1,length=1200>\n'
          '##FILTER=<ID=LowQual,Description="Low quality">\n'
          '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
          '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">\n'
          '##FORMAT=<ID=PGT,Number=1,Type=String,Description="Legacy phase">\n'
          '##FORMAT=<ID=PID,Number=1,Type=String,Description="Legacy block">\n'
          '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">\n'
          '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Genotype likelihoods">\n'
          '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTH\n')


class PhaseGuards(unittest.TestCase):
    def test_missing_blocks_and_legacy_fields_do_not_create_linkage(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'input.vcf'
            path.write_text(HEADER + '1\t100\t.\tA\tC\t99\tPASS\t.\tGT:PS:PGT:PID\t0|1:100:1|0:legacy\n'
                            '1\t180\t.\tA\tG\t99\tPASS\t.\tGT:PS:PGT:PID\t0|1:180:0|1:legacy\n'
                            '1\t250\t.\tA\tT\t99\tPASS\t.\tGT:PS\t0/1:100\n'
                            '1\t350\t.\tA\tT\t99\tPASS\t.\tGT:PS\t1|1:100\n'
                            '1\t450\t.\tA\tT\t99\tPASS\t.\tGT:PS\t1|.:100\n')
            with pysam.VariantFile(str(path)) as reader:
                records = list(reader)
            phases = [audit.ps_phase(r, 'SYNTH', 1) for r in records]
            self.assertEqual(phases[0], ('1', '100', 1))
            self.assertEqual(audit.relation(phases[0], phases[1]), 'unknown')
            self.assertEqual(phases[2:], [None, None, None])
            self.assertEqual(audit.relation(('1', '100', 0), ('2', '100', 1)), 'unknown')
            self.assertEqual(audit.relation(('1', '100', 0), ('1', '100', 1)), 'trans')

    def test_selection_clears_old_phase_without_changing_source_or_genotype(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ref = root / 'ref.fa'; ref.write_text('>1\n' + 'A' * 1200 + '\n'); pysam.faidx(str(ref))
            path = root / 'input.vcf'
            path.write_text(HEADER + '1\t100\t.\tA\tC\t99\tPASS\t.\tGT:PS:PGT:PID\t1|0:100:1|0:legacy\n'
                            '1\t180\t.\tA\tG\t99\t.\t.\tGT:PS:PGT:PID\t0|1:100:0|1:legacy\n'
                            '1\t250\t.\tA\tT\t10\tLowQual\t.\tGT\t0/1\n')
            source = Path(pysam.tabix_index(str(path), preset='vcf'))
            before = source.read_bytes(); output = root / 'selected.vcf'
            counts, signatures = audit.prepare_vcf(source, output, ref, 'SYNTH', {'chrom': '1', 'start': 1, 'end': 1200})
            self.assertEqual(counts['retained_records'], 2)
            self.assertEqual(counts['excluded_filter_records'], 1)
            self.assertEqual(before, source.read_bytes())
            with pysam.VariantFile(str(output)) as reader:
                records = list(reader)
            self.assertEqual(signatures, [audit.genotype_signature(r, 'SYNTH') for r in records])
            for r in records:
                self.assertFalse(r.samples['SYNTH'].phased)
                self.assertFalse(set(audit.PHASE_FIELDS) & set(r.format))
            with self.assertRaisesRegex(ValueError, 'configured sample'):
                audit.prepare_vcf(source, root / 'bad.vcf', ref, 'WRONG', {'chrom': '1', 'start': 1, 'end': 1200})

    def test_genotype_or_variant_changes_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); first = root / 'before.vcf'; second = root / 'after.vcf'
            first.write_text(HEADER + '1\t100\t.\tA\tC\t99\tPASS\t.\tGT\t0/1\n')
            for changed in ['1\t100\t.\tA\tC\t99\tPASS\t.\tGT:PS\t1|1:100\n',
                            '1\t100\t.\tA\tG\t99\tPASS\t.\tGT:PS\t0|1:100\n']:
                second.write_text(HEADER + changed)
                with self.assertRaisesRegex(ValueError, 'changed variant records'):
                    audit.summarize(first, second, 'SYNTH', [])

    def test_corrupted_text_and_changed_likelihoods_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); first = root / 'before.vcf'; second = root / 'after.vcf'
            row = '1\t100\t.\tA\tC\t99\tPASS\t.\tGT:PL\t0/1:100,0,100\n'
            first.write_text(HEADER + row)
            second.write_text(HEADER + row.replace('100,0,100', '100,0,99'))
            with self.assertRaisesRegex(ValueError, 'non-phase sample evidence'):
                audit.summarize(first, second, 'SYNTH', [])
            second.write_bytes((HEADER + row).encode() + b'\x00')
            with self.assertRaisesRegex(ValueError, 'NUL byte'):
                audit.summarize(first, second, 'SYNTH', [])


class NativeReadPhasing(unittest.TestCase):
    def run_truth(self, relationship):
        runtime = Path(os.environ.get('MVA_PROJECT_ROOT', '/Users/jeweled/Documents/mva-hackathon'))
        tool = runtime / 'tools/whatshap-2.8/bin/whatshap'
        if not tool.exists():
            self.skipTest('requires the Studio WhatsHap runtime')
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); rng = random.Random(71439)
            sequence = ''.join(rng.choice('ACGT') for _ in range(1200))
            reference = root / 'ref.fa'; reference.write_text('>1\n' + sequence + '\n'); pysam.faidx(str(reference))
            positions = [100, 180, 850]
            alts = {p: next(b for b in 'ACGT' if b != sequence[p - 1]) for p in positions}
            source = root / 'source.vcf'
            source.write_text(HEADER + ''.join(
                f'1\t{p}\t.\t{sequence[p - 1]}\t{alts[p]}\t99\tPASS\t.\tGT:PGT:PID:AD:PL\t0|1:0|1:legacy:6,6:100,0,100\n'
                for p in positions))
            source = Path(pysam.tabix_index(str(source), preset='vcf'))
            vcf = root / 'input.vcf'
            audit.prepare_vcf(source, vcf, reference, 'SYNTH', {'chrom': '1', 'start': 1, 'end': 1200})
            bam = root / 'reads.bam'
            header = {'HD': {'VN': '1.6', 'SO': 'coordinate'}, 'SQ': [{'SN': '1', 'LN': 1200}],
                      'RG': [{'ID': 'L1', 'SM': 'SYNTH', 'PL': 'ILLUMINA'}]}
            with pysam.AlignmentFile(str(bam), 'wb', header=header) as writer:
                for group_start, group_sites in [(40, [100, 180]), (780, [850])]:
                    for i in range(12):
                        start = group_start + i; bases = list(sequence[start:start + 230])
                        for p in group_sites:
                            alt = (i % 2 == 1)
                            if relationship == 'trans' and p == 180:
                                alt = not alt
                            if alt:
                                bases[p - 1 - start] = alts[p]
                        read = pysam.AlignedSegment(writer.header)
                        read.query_name = f'SYNTH_{group_start}_{i}'; read.query_sequence = ''.join(bases)
                        read.query_qualities = [40] * len(bases); read.reference_id = 0; read.reference_start = start
                        read.cigarstring = f'{len(bases)}M'; read.flag = 0; read.mapping_quality = 60; read.set_tag('RG', 'L1')
                        writer.write(read)
            pysam.index(str(bam))
            out = root / 'phased.vcf'
            cmd = audit.phase_command(tool, vcf, bam, reference, 'SYNTH', '1', out, root / 'reads.tsv')
            process = subprocess.run(cmd, capture_output=True, text=True)
            self.assertEqual(process.returncode, 0, process.stderr)
            targets = [{'chrom': '1', 'pos': p, 'ref': sequence[p - 1], 'alt': alts[p]} for p in positions]
            result = audit.summarize(vcf, out, 'SYNTH', targets[:2])
            self.assertEqual(result['candidate_pair_relation'], relationship)
            self.assertTrue(result['records_and_unordered_genotypes_preserved'])
            disconnected = audit.summarize(vcf, out, 'SYNTH', [targets[0], targets[2]])
            self.assertEqual(disconnected['candidate_pair_relation'], 'unknown')
            with pysam.VariantFile(str(out)) as reader:
                for record in reader:
                    self.assertEqual(record.samples['SYNTH'].get('AD'), (6, 6))
                    self.assertEqual(record.samples['SYNTH'].get('PL'), (100, 0, 100))
            stats = subprocess.run([str(tool), 'stats', '--sample', 'SYNTH', '--chromosome', '1',
                                    '--tsv', str(root / 'stats.tsv'), '--block-list', str(root / 'blocks.tsv'),
                                    str(out)], capture_output=True, text=True)
            self.assertEqual(stats.returncode, 0, stats.stderr)

    def test_actual_tool_recovers_cis_but_not_disconnected_marker(self):
        self.run_truth('cis')

    def test_actual_tool_recovers_trans_but_not_disconnected_marker(self):
        self.run_truth('trans')


if __name__ == '__main__':
    unittest.main()
