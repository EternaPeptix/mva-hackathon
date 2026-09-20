"""Reference evidence must survive missing calls, alternate alleles and mate conflicts."""
import copy
import json
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest

import pysam

sys.path.insert(0, str(Path(__file__).parents[1] / 'scripts'))
import audit_targeted_snvs as audit


POLICY = {'minimum_caller_DP': 10, 'minimum_reference_GQ': 20, 'minimum_reference_fragments': 10,
          'required_modes': ['GATK_BP_gVCF', 'GATK_genotyped', 'bcftools_allsites', 'bcftools_targeted'],
          'GQ_required_modes': ['GATK_BP_gVCF'], 'target_PL_gap_required_modes': ['bcftools_targeted'],
          'minimum_target_reference_PL_gap': 20,
          'require_both_strands': True, 'meaning': 'Selected SNP evidence only.'}
FILTERS = {'minimum_mapping_quality': 20, 'minimum_base_quality': 20, 'exclude_flags': '0xF04'}


def target(pos, ref='A', alt='C'):
    return dict(key=f'1:{pos}:{ref}:{alt}', chrom='1', pos=pos, ref=ref, alt=alt)


def reference_reads():
    return {'counts': {'REF': 30, 'ALT': 0, 'other': 0, 'mate_conflict': 0, 'deletion': 0, 'refskip': 0},
            'strand_counts': {'REF': {'forward': 15, 'reverse': 15}}}


def homref():
    return {'exact_site_records': [{'state': 'hom_reference', 'DP': 30, 'GQ': 99, 'FILTER': [],
                                   'GT': [0, 0], 'target_ALT_index': 1, 'alts': ['C'], 'PL': [0, 90, 200]}],
            'spanning_or_different_ref_records': []}


def callers():
    return {name: homref() for name in POLICY['required_modes']}


class TargetedSnvTests(unittest.TestCase):
    def test_absent_records_never_become_reference_calls(self):
        missing = {'exact_site_records': [], 'spanning_or_different_ref_records': []}
        self.assertEqual(audit.reference_support(reference_reads(), {'caller': missing}, POLICY)['status'], 'reference_support_not_established')
        self.assertIn('GATK_genotyped_missing_mode', audit.reference_support(reference_reads(), {}, POLICY)['reasons'])
        self.assertEqual(audit.reference_support(reference_reads(), callers(), POLICY)['status'], 'reference_supported_in_retained_lanes')

    def test_quality_conflicts_and_overlapping_records_block_reference_support(self):
        for field, value in [('DP', None), ('GQ', 0), ('FILTER', ['LowQual']), ('state', 'no_call_or_non_diploid')]:
            rows = callers(); rows['GATK_BP_gVCF']['exact_site_records'][0][field] = value
            self.assertEqual(audit.reference_support(reference_reads(), rows, POLICY)['status'], 'reference_support_not_established')
        for label in ('ALT', 'other', 'mate_conflict', 'deletion', 'refskip'):
            reads = reference_reads(); reads['counts'][label] = 1
            self.assertIn('nonreference_or_conflicting_fragments', audit.reference_support(reads, callers(), POLICY)['reasons'])
        rows = callers(); rows['GATK_genotyped']['spanning_or_different_ref_records'] = [{'state': 'other_alt_or_representation'}]
        self.assertIn('GATK_genotyped_overlapping_representation_requires_review', audit.reference_support(reference_reads(), rows, POLICY)['reasons'])
        reads = reference_reads(); reads['strand_counts']['REF'] = {'forward': 30}
        self.assertIn('reference_not_observed_on_both_strands', audit.reference_support(reads, callers(), POLICY)['reasons'])

    def test_reference_likelihood_gap_is_not_fabricated_from_missing_or_wrong_alleles(self):
        rows = callers(); rows['GATK_genotyped']['exact_site_records'][0]['GQ'] = None
        rows['bcftools_targeted']['exact_site_records'][0]['GQ'] = None
        self.assertEqual(audit.reference_support(reference_reads(), rows, POLICY)['status'], 'reference_supported_in_retained_lanes')
        for field, value in [('PL', None), ('PL', [0, 0, 99]), ('PL', [0, None, 99]), ('PL', [0, 90, 99, 100]), ('target_ALT_index', None), ('GT', [0, 1])]:
            changed = copy.deepcopy(rows); changed['bcftools_targeted']['exact_site_records'][0][field] = value
            self.assertIn('bcftools_targeted_low_or_unknown_reference_PL_gap', audit.reference_support(reference_reads(), changed, POLICY)['reasons'])

    def test_vcf_retains_multiallelic_genotypes_missingness_and_spanning_deletions(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / 'calls.vcf'
            p.write_text('##fileformat=VCFv4.2\n##contig=<ID=1,length=300>\n'
                         '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                         '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
                         '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Quality">\n'
                         '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTH\n'
                         '1\t50\t.\tA\tC,G\t99\t.\t.\tGT:DP:GQ\t1/2:30:99\n'
                         '1\t60\t.\tA\tC\t.\t.\t.\tGT:DP:GQ\t./.:0:0\n'
                         '1\t69\t.\tAAA\tA\t99\t.\t.\tGT:DP:GQ\t0/1:30:99\n'
                         '1\t80\t.\tA\tC\t99\t.\t.\tDP:GQ\t30:99\n')
            gz = pysam.tabix_index(str(p), preset='vcf')
            result = audit.target_calls(gz, [target(x) for x in (50, 60, 70, 80, 100)], 'SYNTH')
            row = result['1:50:A:C']['exact_site_records'][0]
            self.assertEqual(row['GT'], [1, 2]); self.assertEqual(row['target_ALT_copies'], 1)
            self.assertEqual(result['1:60:A:C']['exact_site_records'][0]['state'], 'no_call_or_non_diploid')
            self.assertEqual(len(result['1:70:A:C']['spanning_or_different_ref_records']), 1)
            self.assertIsNone(result['1:80:A:C']['exact_site_records'][0]['target_ALT_copies'])
            self.assertEqual(result['1:80:A:C']['exact_site_records'][0]['state'], 'no_call_or_non_diploid')
            self.assertEqual(result['1:100:A:C']['exact_site_records'], [])

    def test_baq_fragment_counts_keep_read_groups_conflicts_and_no_coverage(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); rng = random.Random(741)
            seq = ''.join(rng.choice('ACGT') for _ in range(300))
            ref = root / 'ref.fa'; ref.write_text('>1\n' + seq + '\n'); pysam.faidx(str(ref))
            row = target(100, seq[99], next(b for b in 'ACGT' if b != seq[99]))
            bam = root / 'reads.bam'
            header = {'HD': {'VN': '1.6', 'SO': 'coordinate'}, 'SQ': [{'SN': '1', 'LN': 300}], 'RG': [{'ID': 'L1'}, {'ID': 'L2'}]}
            specs = [('same', False, 0, 'L1', 60), ('same', False, 16, 'L1', 60),
                     ('same', True, 0, 'L2', 60), ('conflict', False, 0, 'L1', 60),
                     ('conflict', True, 16, 'L1', 60), ('dup', True, 1024, 'L1', 60), ('lowmap', True, 0, 'L1', 0)]
            with pysam.AlignmentFile(str(bam), 'wb', header=header) as writer:
                for name, alternate, flag, group, mq in specs:
                    read = pysam.AlignedSegment(writer.header); bases = list(seq[60:140])
                    if alternate: bases[39] = row['alt']
                    read.query_name = name; read.query_sequence = ''.join(bases)
                    read.query_qualities = [40] * 80; read.reference_id = 0; read.reference_start = 60
                    read.cigarstring = '80M'; read.flag = flag; read.mapping_quality = mq; read.set_tag('RG', group)
                    writer.write(read)
            pysam.index(str(bam)); audit.validate_targets([row], ref)
            result = audit.fragment_evidence(bam, ref, row, FILTERS)
            self.assertEqual(result['counts']['REF'], 1)
            self.assertEqual(result['counts']['ALT'], 1)
            self.assertEqual(result['counts']['mate_conflict'], 1)
            empty = audit.fragment_evidence(bam, ref, target(250, seq[249], next(b for b in 'ACGT' if b != seq[249])), FILTERS)
            self.assertEqual(sum(empty['counts'].values()), 0)
            wrong = copy.deepcopy(row); wrong['ref'] = row['alt']; wrong['alt'] = row['ref']; wrong['key'] = f"1:100:{wrong['ref']}:{wrong['alt']}"
            with self.assertRaisesRegex(ValueError, 'reference mismatch'): audit.validate_targets([wrong], ref)
            with self.assertRaisesRegex(ValueError, 'duplicate target'): audit.validate_targets([row, row], ref)


class NativeReferenceEvidenceTest(unittest.TestCase):
    def test_real_callers_distinguish_homref_heterozygous_and_uncovered_sites(self):
        runtime = Path(os.environ.get('MVA_PROJECT_ROOT', '/Users/jeweled/Documents/mva-hackathon'))
        if not (runtime / 'tools/bio/bin/bcftools').exists():
            self.skipTest('native Studio runtime required')
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); (root / 'tools').symlink_to(runtime / 'tools', target_is_directory=True)
            (root / 'gatk-4.6.1.0').symlink_to(runtime / 'gatk-4.6.1.0', target_is_directory=True)
            refdir = root / 'ref/vendor'; refdir.mkdir(parents=True)
            rng = random.Random(741963); seq = ''.join(rng.choice('ACGT') for _ in range(3000))
            reference = refdir / 'GRCh38_nochr.fa'; reference.write_text('>1\n' + seq + '\n'); pysam.faidx(str(reference))
            subprocess.run([str(runtime / 'tools/bio/bin/samtools'), 'dict', '-o', str(refdir / 'GRCh38_nochr.dict'), str(reference)], capture_output=True, check=True)
            targets = [target(pos, seq[pos - 1], next(b for b in 'ACGT' if b != seq[pos - 1])) for pos in (500, 1000, 2000)]
            config = {'scope': 'Synthetic truth: reference, heterozygous and uncovered.', 'sample': 'SYNTH',
                      'reference_sha256': audit.fingerprint(reference)['sha256'], 'targets': targets,
                      'read_filters': FILTERS, 'reference_support_policy': POLICY}
            config_path = root / 'config.json'; config_path.write_text(json.dumps(config))
            slices = root / 'slices'; slices.mkdir()
            bam = slices / 'L001.bam'
            header = {'HD': {'VN': '1.6', 'SO': 'coordinate'}, 'SQ': [{'SN': '1', 'LN': 3000}],
                      'RG': [{'ID': 'SYNTH', 'SM': 'SYNTH', 'PL': 'ILLUMINA'}]}
            with pysam.AlignmentFile(str(bam), 'wb', header=header) as writer:
                for locus in targets[:2]:
                    for j in range(40):
                        start = locus['pos'] - 100 + j
                        bases = list(seq[start:start + 150])
                        if locus['pos'] == 1000 and j % 2: bases[locus['pos'] - 1 - start] = locus['alt']
                        read = pysam.AlignedSegment(writer.header); read.query_name = f"SYNTH{locus['pos']}_{j}"
                        read.query_sequence = ''.join(bases); read.query_qualities = [40] * 150
                        read.reference_id = 0; read.reference_start = start; read.cigarstring = '150M'
                        read.mapping_quality = 60; read.flag = 16 if j % 4 < 2 else 0; read.set_tag('RG', 'SYNTH'); writer.write(read)
            manifest = {'status': 'complete', 'config_sha256': audit.fingerprint(config_path)['sha256'],
                        'files': [{'slice_name': bam.name, 'slice_sha256': audit.fingerprint(bam)['sha256']}],
                        'regions_BED_0based_half_open': [['1', 200, 2200]]}
            (slices / 'slice-provenance.json').write_text(json.dumps(manifest))
            output = root / 'result'
            result = subprocess.run([sys.executable, audit.__file__, '--root', str(root), '--config', str(config_path),
                                     '--slices', str(slices), '--output', str(output)], capture_output=True, text=True)
            if result.returncode:
                self.fail(result.stderr + '\n'.join(p.read_text()[-1600:] for p in output.glob('*.log')))
            evidence = json.loads((output / 'evidence.json').read_text())
            self.assertEqual(evidence['reference_support'][targets[0]['key']]['status'], 'reference_supported_in_retained_lanes',
                json.dumps({'decision': evidence['reference_support'][targets[0]['key']],
                            'reads': evidence['reads'][targets[0]['key']],
                            'calls': {k: v[targets[0]['key']] for k, v in evidence['calls'].items()}}, indent=2))
            for row in targets[1:]:
                self.assertEqual(evidence['reference_support'][row['key']]['status'], 'reference_support_not_established')
            for name in ('GATK_genotyped', 'bcftools_allsites', 'bcftools_targeted'):
                records = evidence['calls'][name][targets[1]['key']]['exact_site_records']
                self.assertEqual(len(records), 1, name); self.assertEqual(records[0]['GT'], [0, 1], name)
                self.assertEqual(records[0]['target_ALT_copies'], 1, name)
            self.assertEqual(sum(evidence['reads'][targets[2]['key']]['counts'].values()), 0)


if __name__ == '__main__':
    unittest.main()
