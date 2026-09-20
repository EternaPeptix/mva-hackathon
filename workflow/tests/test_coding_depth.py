import copy
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

SCRIPTS = Path(__file__).resolve().parents[1] / 'scripts'
sys.path.insert(0, str(SCRIPTS))
from coding_depth_targets import prepare_targets, summarize_depth
from measure_coding_depth import check_quality_availability


GFF = '''##gff-version 3
1\tx\tgene\t1\t30\t.\t+\t.\tID=gene:A;Name=ALPHA;biotype=protein_coding
1\tx\tmRNA\t1\t30\t.\t+\t.\tID=transcript:A1;Parent=gene:A
1\tx\tCDS\t9\t18\t.\t+\t0\tParent=transcript:A1
1\tx\tCDS\t26\t29\t.\t+\t2\tParent=transcript:A1
1\tx\tmRNA\t1\t30\t.\t+\t.\tID=transcript:A2;Parent=gene:A
1\tx\tCDS\t9\t18\t.\t+\t0\tParent=transcript:A2
1\tx\tgene\t1\t35\t.\t-\t.\tID=gene:B;Name=BETA;biotype=protein_coding
1\tx\tmRNA\t1\t35\t.\t-\t.\tID=transcript:B1;Parent=gene:B
1\tx\tCDS\t13\t24\t.\t-\t0\tParent=transcript:B1
2\tx\tgene\t1\t20\t.\t+\t.\tID=gene:C;Name=GAMMA;biotype=protein_coding
2\tx\tmRNA\t1\t20\t.\t+\t.\tID=transcript:C1;Parent=gene:C
2\tx\tCDS\t1\t4\t.\t+\t0\tParent=transcript:C1
'''


def annotation(root, content=GFF):
    p = Path(root) / 'fixture.gff3.gz'
    with gzip.open(p, 'wt') as f:
        f.write(content)
    return p


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


class CodingTargetTests(unittest.TestCase):
    def design(self, content=GFF):
        with tempfile.TemporaryDirectory() as temp:
            return prepare_targets(annotation(temp, content), ['ALPHA', 'BETA', 'GAMMA'], {'1': 60, '2': 20})

    def test_isoforms_and_overlapping_genes_keep_exact_union_denominators(self):
        d = self.design()
        self.assertEqual(d['panel_union_bases'], 24)
        self.assertEqual(d['gene_union_bases_sum'], 30)
        self.assertEqual([(w['chrom'], w['start'], w['end'], w['gene_ids']) for w in d['windows']],
            [('1', 8, 12, ['gene:A']), ('1', 12, 18, ['gene:A', 'gene:B']),
             ('1', 18, 24, ['gene:B']), ('1', 25, 29, ['gene:A']), ('2', 0, 4, ['gene:C'])])
        self.assertEqual(d['genes'][0]['transcript_ids'], ['transcript:A1', 'transcript:A2'])
        self.assertEqual(d['genes'][1]['strand'], '-')

    def test_annotation_order_does_not_drop_cds_before_transcript(self):
        original = self.design()
        reversed_rows = '\n'.join(reversed(GFF.strip().splitlines()[1:])) + '\n'
        d = self.design(reversed_rows)
        self.assertEqual(d['windows'], original['windows'])
        self.assertEqual([g['CDS_union_bases'] for g in d['genes']], [14, 12, 4])

    def test_missing_or_ambiguous_gene_refused(self):
        for content in [GFF.replace('Name=ALPHA', 'Name=OTHER'),
                        GFF + '1\tx\tgene\t40\t50\t.\t+\t.\tID=gene:EXTRA;Name=ALPHA;biotype=protein_coding\n']:
            with self.subTest(content=content), self.assertRaisesRegex(ValueError, 'missing or ambiguous'):
                self.design(content)

    def test_child_contig_strand_and_boundaries_refused(self):
        original = '1\tx\tCDS\t9\t18\t.\t+\t0\tParent=transcript:A1'
        for replacement in [original.replace('1\tx', '2\tx'), original.replace('\t+\t', '\t-\t'),
                            original.replace('\t9\t18\t', '\t9\t31\t')]:
            with self.subTest(replacement=replacement), self.assertRaisesRegex(ValueError, 'geometry disagrees'):
                self.design(GFF.replace(original, replacement))

    def test_unresolved_parent_and_absent_cds_refused(self):
        with self.assertRaisesRegex(ValueError, 'unresolved CDS'):
            self.design(GFF.replace('Parent=transcript:C1', 'Parent=transcript:UNKNOWN'))
        with self.assertRaisesRegex(ValueError, 'no annotated CDS'):
            self.design('\n'.join(x for x in GFF.splitlines() if 'Parent=transcript:C1' not in x))

    def test_duplicate_name_and_reference_bounds_refused(self):
        with tempfile.TemporaryDirectory() as temp:
            p = annotation(temp)
            with self.assertRaisesRegex(ValueError, 'unique gene'):
                prepare_targets(p, ['ALPHA', 'ALPHA'], {'1': 60, '2': 20})
            with self.assertRaisesRegex(ValueError, 'outside reference'):
                prepare_targets(p, ['ALPHA'], {'1': 20, '2': 20})

    def test_zero_depth_and_noncanonical_reference_denominators(self):
        # Literal fake depth rows include a fully ambiguous, fully uncovered gene.
        d = self.design()
        rows = []
        for w in d['windows']:
            n = w['end'] - w['start']
            acgt = 0 if w['chrom'] == '2' else n
            rows.append(dict(w, span_depth={'histogram': {0: n}}, ACGT_depth={'histogram': {0: acgt}}))
        s = summarize_depth(rows, d)
        self.assertEqual(s['panel_union']['span']['zero_depth_bases'], 24)
        self.assertEqual(s['gene:C']['span']['bases'], 4)
        self.assertIsNone(s['gene:C']['ACGT']['mean'])
        self.assertEqual(s['gene:C']['non_ACGT_reference_bases'], 4)
        self.assertEqual(s['gene:A']['ACGT']['bases'], 14)
        self.assertEqual(s['gene:B']['ACGT']['bases'], 12)
        for bad in [rows[:-1], list(reversed(rows))]:
            with self.assertRaises(ValueError):
                summarize_depth(bad, d)

    def test_changed_membership_or_histogram_denominator_refused(self):
        d = self.design()
        rows = [dict(w, span_depth={'histogram': {0: w['end'] - w['start']}},
                     ACGT_depth={'histogram': {0: w['end'] - w['start']}}) for w in d['windows']]
        wrong = copy.deepcopy(rows)
        wrong[0]['gene_ids'] = ['gene:B']
        with self.assertRaisesRegex(ValueError, 'target identity'):
            summarize_depth(wrong, d)
        wrong = copy.deepcopy(rows)
        wrong[0]['span_depth']['histogram'][0] += 1
        with self.assertRaisesRegex(ValueError, 'denominator'):
            summarize_depth(wrong, d)


class NativeCodingDepthTests(unittest.TestCase):
    def test_native_unavailable_base_quality_is_not_a_quality_pass(self):
        import pysam
        samtools = Path(os.environ['MVA_SAMTOOLS']).resolve()
        header = '@HD\tVN:1.6\n@SQ\tSN:1\tLN:20\n'
        for qual, expected in [('*', 1), ('!', 0), ('I', 1)]:
            sam = header + 'read\t0\t1\t5\t60\t1M\t*\t0\t0\tA\t' + qual + '\n'
            result = subprocess.run([samtools, 'depth', '-q', '20', '-Q', '20', '-G', '0xF04', '-s', '-'],
                                    input=sam, capture_output=True, text=True, check=True)
            self.assertEqual(int(result.stdout.split()[2]), expected)
            read = pysam.AlignedSegment.fromstring(sam.splitlines()[-1], pysam.AlignmentHeader.from_text(header))
            if qual == '*':
                with self.assertRaisesRegex(ValueError, 'unknown base qualities'):
                    check_quality_availability(read)
            else:
                check_quality_availability(read)
        read.query_sequence = 'AAA'
        read.cigarstring = '3M'
        read.query_qualities = [40, 255, 40]
        with self.assertRaisesRegex(ValueError, 'unknown base qualities'):
            check_quality_availability(read)
        read.flag = 1024
        check_quality_availability(read)  # Excluded duplicate does not enter any depth arm.

    def test_cli_bam_cram_quality_overlap_zero_and_input_refusals(self):
        import pysam
        samtools = Path(os.environ['MVA_SAMTOOLS']).resolve()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            reference = root / 'reference.fa'
            reference.write_text('>1\n' + 'A' * 15 + 'N' + 'A' * 44 + '\n>2\n' + 'N' * 4 + 'A' * 16 + '\n')
            pysam.faidx(str(reference))
            gff = annotation(root)
            header = {'HD': {'VN': '1.6', 'SO': 'coordinate'}, 'SQ': [{'SN': '1', 'LN': 60}, {'SN': '2', 'LN': 20}],
                      'RG': [{'ID': 'rg' + str(i), 'SM': 'fictional'} for i in range(4)]}
            records = [('pair', 8, '12M', 12, 99, 60, 40, 12),
                       ('lowmap', 8, '10M', 10, 0, 0, 40, -1),
                       ('lowbase', 8, '10M', 10, 0, 60, 10, -1),
                       ('duplicate', 8, '10M', 10, 1024, 60, 40, -1),
                       ('supplementary', 8, '10M', 10, 2048, 60, 40, -1),
                       ('deletion', 8, '3M2D5M', 8, 0, 60, 40, -1),
                       ('skip', 8, '3M2N5M', 8, 0, 60, 40, -1),
                       ('pair', 12, '12M', 12, 147, 60, 40, 8)]
            bam = root / 'fictional.bam'
            with pysam.AlignmentFile(str(bam), 'wb', header=header) as out:
                for i, (name, start, cigar, length, flag, mq, bq, mate) in enumerate(records):
                    read = pysam.AlignedSegment(out.header)
                    read.query_name, read.query_sequence = name, 'A' * length
                    read.query_qualities = [bq] * length
                    read.reference_id, read.reference_start = 0, start
                    read.cigarstring, read.flag, read.mapping_quality = cigar, flag, mq
                    read.set_tag('RG', 'rg0' if name == 'pair' else 'rg' + str(i % 4))
                    if mate >= 0:
                        read.next_reference_id, read.next_reference_start = 0, mate
                        read.template_length = 16 if start == 8 else -16
                    out.write(read)
            pysam.index(str(bam))
            cram = root / 'fictional.cram'
            subprocess.run([samtools, 'view', '-C', '-T', reference, '-o', cram, bam], check=True, capture_output=True)
            subprocess.run([samtools, 'index', cram], check=True, capture_output=True)
            observed = []
            for alignment in (bam, cram):
                inputs = dict(alignment=alignment, alignment_index=Path(str(alignment) + ('.bai' if alignment == bam else '.crai')),
                              reference=reference, reference_index=Path(str(reference) + '.fai'), gff=gff, samtools=samtools)
                config = dict(schema_version=1, assembly='GRCh38', sample='fictional', expected_read_groups=4,
                    genes=['ALPHA', 'BETA', 'GAMMA'], inputs={k: {'path': str(p), 'sha256': sha(p)} for k, p in inputs.items()},
                    source_sha256={p.name: sha(p) for p in SCRIPTS.glob('*.py')},
                    depth_modes=[dict(name='MQ0_BQ1', minimum_MQ=0, minimum_BQ=1),
                                 dict(name='MQ20_BQ20', minimum_MQ=20, minimum_BQ=20),
                                 dict(name='MQ60_BQ20', minimum_MQ=60, minimum_BQ=20)])
                c = root / (alignment.suffix[1:] + '.json');c.write_text(json.dumps(config))
                output = root / (alignment.suffix[1:] + '-result')
                cmd = [sys.executable, '-B', str(SCRIPTS / 'measure_coding_depth.py'), '--config', str(c),
                       '--expected-config-sha256', sha(c), '--output', str(output)]
                completed = subprocess.run(cmd, text=True, capture_output=True)
                self.assertEqual(completed.returncode, 0, completed.stderr)
                d = json.loads((output / 'summary.json').read_text());observed.append(d['depth_modes'])
                self.assertEqual(d['panel_union_bases'], 24)
                for label, expected in [('MQ20_BQ20', [3,3,3,1,1,3,3,3,3,3] + [1]*6 + [0]*8),
                                        ('MQ60_BQ20', [3,3,3,1,1,3,3,3,3,3] + [1]*6 + [0]*8),
                                        ('MQ0_BQ1', [5,5,5,3,3,5,5,5,5,5] + [1]*6 + [0]*8)]:
                    with gzip.open(output / (label + '.depth.tsv.gz'), 'rt') as f:
                        rows = [x.split() for x in f]
                    self.assertEqual([int(x[2]) for x in rows], expected)
                    self.assertEqual(d['depth_modes'][label]['groups']['gene:C']['span']['zero_depth_bases'], 4)
                    self.assertIsNone(d['depth_modes'][label]['groups']['gene:C']['ACGT']['mean'])
                    self.assertEqual(d['depth_modes'][label]['groups']['panel_union']['non_ACGT_reference_bases'], 5)
                saved = sha(output / 'status.json')
                rerun = subprocess.run(cmd, text=True, capture_output=True)
                self.assertNotEqual(rerun.returncode, 0)
                self.assertEqual(sha(output / 'status.json'), saved)
                bad = copy.deepcopy(config);bad['inputs']['gff']['sha256'] = '0' * 64
                c.write_text(json.dumps(bad));refused = root / (alignment.suffix[1:] + '-refused')
                refusal = subprocess.run(cmd[:5] + ['--expected-config-sha256', sha(c), '--output', str(refused)], text=True, capture_output=True)
                self.assertNotEqual(refusal.returncode, 0)
                self.assertFalse(refused.exists())
            self.assertEqual(observed[0], observed[1])
            for name, qualities in [('absent', None), ('partial255', [40] * 5 + [255] + [40] * 6)]:
                unknown = root / ('unknown-' + name + '.bam')
                with pysam.AlignmentFile(str(bam)) as source, pysam.AlignmentFile(str(unknown), 'wb', template=source) as dest:
                    for i, read in enumerate(source):
                        if i == 0:
                            read.query_qualities = qualities
                        dest.write(read)
                pysam.index(str(unknown))
                bad = copy.deepcopy(config)
                bad['inputs']['alignment'] = {'path': str(unknown), 'sha256': sha(unknown)}
                bad['inputs']['alignment_index'] = {'path': str(unknown) + '.bai', 'sha256': sha(str(unknown) + '.bai')}
                c = root / (name + '.json');c.write_text(json.dumps(bad))
                output = root / (name + '-result')
                result = subprocess.run([sys.executable, '-B', str(SCRIPTS / 'measure_coding_depth.py'),
                    '--config', str(c), '--expected-config-sha256', sha(c), '--output', str(output)],
                    text=True, capture_output=True)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('unknown base qualities', result.stderr)
                self.assertEqual(json.loads((output / 'status.json').read_text())['status'], 'failed')
                self.assertFalse((output / 'summary.json').exists())


if __name__ == '__main__':
    unittest.main()
