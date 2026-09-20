"""Native fictional CRAM, raw-call, extraction and constructed-annotation controls.

Caller names are provenance labels on synthetic records. No actual callers or
case genomes are used, and constructed haplotypes never become genotype calls.
"""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import pysam


WORKFLOW = Path(__file__).resolve().parents[1]
SCRIPT = WORKFLOW / 'scripts/review_joint_haplotypes.py'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def pin(path):
    return {'path': str(path), 'bytes': path.stat().st_size, 'sha256': digest(path)}


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def snapshot(paths):
    return {str(path): (path.read_bytes(), path.stat().st_mtime_ns) for path in paths}


class FictionalJointFixture:
    def __init__(self, folder, runtime, difficult=False):
        self.folder, self.runtime = folder.resolve(), runtime.resolve()
        self.bin = self.runtime / 'tools/bio/bin'
        self.reference = self.folder / 'reference.fa'
        sequence = list('ATG' + 'A' * 291 + 'TAA' + 'AAA')
        sequence[100:108] = 'ACAACTGT'
        self.sequence = ''.join(sequence)
        self.reference.write_text('>1\n' + self.sequence + '\n')
        pysam.faidx(str(self.reference))
        self.gff = self.folder / 'fictional.gff3'
        self.gff.write_text(
            '##gff-version 3\n'
            '##sequence-region 1 1 300\n'
            '1\tensembl\tgene\t1\t297\t.\t+\t.\tID=gene:ENSG00000000001;Name=FICTIONAL;biotype=protein_coding\n'
            '1\tensembl\tmRNA\t1\t297\t.\t+\t.\tID=transcript:ENST00000000001;Parent=gene:ENSG00000000001;biotype=protein_coding\n'
            '1\tensembl\texon\t1\t297\t.\t+\t.\tID=exon:ENSE00000000001;Parent=transcript:ENST00000000001\n'
            '1\tensembl\tCDS\t1\t297\t.\t+\t0\tID=CDS:ENSP00000000001;Parent=transcript:ENST00000000001\n')
        self.edits = {
            'gatk': [(101, 'A', 'ACAGC'), (104, 'ACTGT', 'A')],
            'bcftools': [(104, 'A', 'G'), (106, 'T', 'C'), (107, 'G', 'A'), (108, 'T', 'A')],
        }
        if difficult:
            # The first deletion starts well before the widest window, but its
            # literal REF span overlaps it. A POS-only regional filter misses it.
            self.edits['gatk'].extend([(60, self.sequence[59:102], self.sequence[59]),
                                      (106, 'T', 'C'),
                                      (131, self.sequence[130:136], self.sequence[130])])
        self.raw, self.expected_context, sources = {}, {}, []
        for caller in ('gatk', 'bcftools'):
            seed = self.folder / (caller + '-seed.vcf')
            header = ('##fileformat=VCFv4.2\n##contig=<ID=1,length=300>\n'
                      '##FORMAT=<ID=GT,Number=1,Type=String,Description="Native genotype">\n'
                      '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Native allele depths">\n'
                      '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Native depth">\n'
                      '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Native genotype quality">\n'
                      '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Native likelihoods">\n'
                      '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTH\n')
            rows = []
            for position, reference, alternate in sorted(self.edits[caller] + [(200, 'A', 'C')]):
                fields, values = ('GT:AD:DP:GQ:PL', '0/1:12,3:15:35:99,0,99') if caller == 'gatk' else (
                    'GT:AD:DP:PL', '0/1:12,3:15:99,0,99')
                rows.append('1\t%d\t.\t%s\t%s\t40\t.\t.\t%s\t%s\n' %
                            (position, reference, alternate, fields, values))
            seed.write_text(header + ''.join(rows))
            raw = self.folder / ('gatk.raw.vcf.gz' if caller == 'gatk' else 'bcftools.raw.bcf')
            self.native('bcftools', ['view', '-Oz' if caller == 'gatk' else '-Ob', '-o', str(raw), str(seed)])
            self.raw[caller] = raw
            payload, count = hashlib.sha256(), 0
            with pysam.VariantFile(str(raw)) as reader:
                for ordinal, record in enumerate(reader, 1):
                    count += 1
                    payload.update(str(record).encode())
                    self.expected_context[(caller, ordinal)] = {
                        'rendered_record': str(record), 'FORMAT': json.loads(json.dumps(dict(record.samples['SYNTH']))),
                        'pos': record.pos, 'ref': record.ref, 'alts': list(record.alts),
                    }
            sources.append(dict(input=caller + '_raw', caller=caller, expected_records=count,
                                expected_rendered_record_payload_sha256=payload.hexdigest()))
        header = {'HD': {'VN': '1.6', 'SO': 'coordinate'}, 'SQ': [{'SN': '1', 'LN': 300}],
                  'RG': [{'ID': 'L1', 'SM': 'SYNTH', 'PL': 'ILLUMINA'}]}
        bam = self.folder / 'original.bam'
        altered = self.sequence[70:100] + 'ACAGCCAA' + self.sequence[108:140]
        reads = [dict(name='indel-template', sequence=altered, cigar='31M4I3M4D32M', flag=0, mq=60),
                 dict(name='snp-template', sequence=altered, cigar='70M', flag=65, mq=60),
                 dict(name='snp-template', sequence=altered, cigar='70M', flag=145, mq=60)]
        if difficult:
            first, second = list(altered), list(altered)
            first[31], second[31] = 'T', 'G'  # no native hypothesis edits this core base
            reads.extend([dict(name='unmatched-conflict', sequence=''.join(first), cigar='70M', flag=65, mq=60),
                          dict(name='unmatched-conflict', sequence=''.join(second), cigar='70M', flag=145, mq=60),
                          dict(name='unknown-map', sequence=altered, cigar='70M', flag=0, mq=255)])
        with pysam.AlignmentFile(str(bam), 'wb', header=header) as writer:
            for row in reads:
                record = pysam.AlignedSegment(writer.header)
                record.query_name, record.query_sequence = row['name'], row['sequence']
                record.query_qualities = [40] * len(row['sequence'])
                record.reference_id, record.reference_start = 0, 70
                record.cigarstring, record.flag, record.mapping_quality = row['cigar'], row['flag'], row['mq']
                record.set_tag('RG', 'L1')
                if record.is_paired:
                    record.next_reference_id, record.next_reference_start, record.template_length = 0, 70, 70
                writer.write(record)
            off = pysam.AlignedSegment(writer.header)
            off.query_name, off.query_sequence = 'outside-target', self.sequence[220:240]
            off.query_qualities = [40] * 20
            off.reference_id, off.reference_start, off.cigarstring, off.mapping_quality = 0, 220, '20M', 60
            off.set_tag('RG', 'L1')
            writer.write(off)
        self.alignment = self.folder / 'original.cram'
        self.native('samtools', ['view', '-C', '--output-fmt-option', 'version=3.0', '-T', str(self.reference),
                                 '-o', str(self.alignment), str(bam)])
        self.native('samtools', ['index', str(self.alignment)])
        self.expected_extracted_records = len(reads)
        inputs = {'alignment': self.alignment, 'alignment_index': Path(str(self.alignment) + '.crai'),
                  'reference': self.reference, 'reference_index': Path(str(self.reference) + '.fai'),
                  'gff': self.gff, 'samtools': self.bin / 'samtools', 'bcftools': self.bin / 'bcftools',
                  'gatk_raw': self.raw['gatk'], 'bcftools_raw': self.raw['bcftools']}
        self.config_data = dict(inputs={name: pin(path) for name, path in inputs.items()},
            source_sha256={name: digest(WORKFLOW / 'scripts' / name) for name in
                           ('review_joint_haplotypes.py', 'joint_haplotype_helpers.py', 'native_coding_records.py')},
            sample='SYNTH', expected_read_groups=1, caller_sources=sources,
            targets=[dict(id='FICTIONAL', chrom='1', core_start0=100, core_end0=108)],
            window_flanks_bp=[3, 10, 25], mapping_quality_thresholds=[0, 20, 60],
            base_quality_thresholds=[20, 30], maximum_hypothesis_states=256,
            maximum_haplotype_bases=100000, raw_base_qualities_no_BAQ=True)
        self.config = self.folder / 'config.json'
        write_json(self.config, self.config_data)
        self.input_paths = [path for name, path in inputs.items() if name not in ('samtools', 'bcftools')] + [self.config]

    def native(self, tool, arguments):
        subprocess.run([str(self.bin / tool), *arguments], capture_output=True, text=True, check=True)

    def execute(self, name):
        output = self.folder / name
        process = subprocess.run([sys.executable, str(SCRIPT), '--root', str(self.runtime),
                                  '--config', str(self.config), '--output', str(output)],
                                 capture_output=True, text=True)
        return process, output


class NativeJointWorkflowTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.runtime = Path(os.environ.get('MVA_PROJECT_ROOT', '/Users/jeweled/Documents/mva-hackathon'))
        if not all((cls.runtime / 'tools/bio/bin' / tool).is_file() for tool in ('samtools', 'bcftools')):
            raise unittest.SkipTest('installed Studio pysam, samtools and bcftools required')

    def completed(self, process, output):
        logs = '\n'.join(path.read_text(errors='replace')[-3000:] for path in output.glob('*.log'))
        self.assertEqual(process.returncode, 0, process.stdout + process.stderr + logs)
        summary = json.loads((output / 'summary.json').read_text())
        self.assertEqual(summary['status'], 'complete')
        self.assertEqual(json.loads((output / 'status.json').read_text())['summary_sha256'], digest(output / 'summary.json'))
        return summary, json.loads((output / 'evidence.json').read_text())

    def failed(self, process, output, reason):
        self.assertNotEqual(process.returncode, 0, process.stdout + process.stderr)
        status = json.loads((output / 'status.json').read_text())
        self.assertEqual(status['status'], 'failed')
        self.assertIn(reason, status['error'])
        self.assertFalse((output / 'summary.json').exists())

    def test_native_indel_and_four_snp_lineages_share_two_fragments_and_sites_only_annotations(self):
        with tempfile.TemporaryDirectory() as temporary:
            fixture = FictionalJointFixture(Path(temporary), self.runtime)
            before = snapshot(fixture.input_paths)
            process, output = fixture.execute('review')
            summary, evidence = self.completed(process, output)
            self.assertEqual(summary['extraction_payload_check']['records'], 3)
            self.assertEqual(summary['quality_arms'], 18)
            self.assertTrue(summary['annotation_check']['annotation_payload_preserved'])
            self.assertTrue(summary['annotation_check']['every_normalized_constructed_allele_sequence_verified'])
            contexts = json.loads((output / 'regional-native-context.json').read_text())
            self.assertEqual({row['caller']: row['records'] for row in contexts['source_census']}, {'gatk': 3, 'bcftools': 5})
            self.assertEqual(len(contexts['targets']['FICTIONAL']), 6)
            for row in contexts['targets']['FICTIONAL']:
                context = row['context']
                original = fixture.expected_context[(context['caller'], context['source_ordinal'])]
                for field in ('rendered_record', 'FORMAT', 'pos', 'ref', 'alts'):
                    self.assertEqual(context[field], original[field])
                self.assertEqual(context['FORMAT']['GT'], [0, 1])
                self.assertEqual(context['FORMAT']['PL'], [99, 0, 99])
            for window in evidence['windows']:
                expected = fixture.sequence[window['start0']:100] + 'ACAGCCAA' + fixture.sequence[108:window['end0']]
                equivalent = [evidence['hypotheses'][identifier] for identifier in window['hypothesis_ids']
                              if evidence['hypotheses'][identifier]['sequence'] == expected]
                self.assertEqual(len(equivalent), 1)
                hypothesis = equivalent[0]
                self.assertEqual(sorted((row['caller'], len(row['edits'])) for row in hypothesis['source_combinations']),
                                 [('bcftools', 4), ('gatk', 2)])
                self.assertEqual(hypothesis['annotation_status'], 'local_constructed_hypothesis_annotation')
                self.assertTrue(hypothesis['joint_consequences'])
                for arm in window['arms']:
                    checked = arm['evidence']
                    self.assertEqual(checked['counts']['accepted_observations'], 3)
                    self.assertEqual(checked['counts']['matched_fragments'], 2)
                    self.assertEqual(checked['fragment_sequence_counts'][expected], 2)
                    self.assertEqual(checked['counts']['conflicting_fragments'], 0)
                    pair = next(row for row in checked['fragments'] if row['QNAME'] == 'snp-template')
                    self.assertEqual(pair['accepted_observations'], 2)
                    self.assertEqual(pair['hypothesis_labels'], [hypothesis['id']])
            for name in ('constructed-hypotheses.vcf', 'constructed-normalized.bcf', 'constructed-annotated.bcf'):
                with pysam.VariantFile(str(output / name)) as reader:
                    self.assertEqual(list(reader.header.samples), [])
                    self.assertIn('Constructed sequence edits', str(reader.header))
                    self.assertTrue(all(len(record.samples) == 0 for record in reader))
            self.assertEqual(snapshot(fixture.input_paths), before)

    def test_unknown_mq_unmatched_mates_and_overlapping_or_crossing_edits_stay_explicit(self):
        with tempfile.TemporaryDirectory() as temporary:
            fixture = FictionalJointFixture(Path(temporary), self.runtime, difficult=True)
            before = snapshot(fixture.input_paths)
            process, output = fixture.execute('review-difficult')
            summary, evidence = self.completed(process, output)
            self.assertEqual(summary['extraction_payload_check']['records'], 6)
            contexts = json.loads((output / 'regional-native-context.json').read_text())['targets']['FICTIONAL']
            self.assertTrue(any(row['context']['pos'] == 60 for row in contexts))
            for window in evidence['windows']:
                omitted = window['omitted_overlapping_native_alleles']
                self.assertTrue(any(row['pos'] == 60 and row['reason'] == 'crosses_window_boundary' for row in omitted))
                self.assertGreater(window['enumeration_metadata']['counters']['conflicting_REF_span_combinations'], 0)
                if window['flank_bp'] == 25:
                    self.assertTrue(any(row['pos'] == 131 and row['reason'] == 'crosses_window_boundary' for row in omitted))
                for arm in window['arms']:
                    checked = arm['evidence']
                    self.assertEqual(checked['rejection_counts']['mapping_quality_unavailable_255'], 1)
                    self.assertEqual(checked['counts']['conflicting_fragments'], 1)
                    conflict = next(row for row in checked['fragments'] if row['QNAME'] == 'unmatched-conflict')
                    self.assertEqual(conflict['status'], 'mate_or_repeated_observation_conflict')
                    self.assertEqual(len(conflict['sequences']), 2)
                    observations = [row for row in checked['observations'] if row['QNAME'] == 'unmatched-conflict']
                    self.assertTrue(all(row['status'] == 'accepted' and not row['matched_hypothesis'] for row in observations))
                    self.assertTrue(checked['thresholds']['reject_MQ_255'])
                    self.assertFalse(checked['thresholds']['BAQ_applied'])
            self.assertEqual(snapshot(fixture.input_paths), before)

    def test_raw_identity_full_payload_and_read_source_hash_refusals(self):
        with tempfile.TemporaryDirectory() as temporary:
            fixture = FictionalJointFixture(Path(temporary), self.runtime)
            original = json.loads(json.dumps(fixture.config_data))
            mutations = [
                ('bad-raw-hash', lambda value: value['inputs']['gatk_raw'].update(sha256='0' * 64), 'input SHA differs'),
                ('bad-read-hash', lambda value: value['inputs']['alignment'].update(sha256='0' * 64), 'input SHA differs'),
                ('bad-raw-payload', lambda value: value['caller_sources'][0].update(expected_rendered_record_payload_sha256='0' * 64), 'native source rendered payload differs'),
                ('bad-raw-census', lambda value: value['caller_sources'][0].update(expected_records=999), 'native source census count differs'),
            ]
            for name, mutate, reason in mutations:
                with self.subTest(refusal=name):
                    config = json.loads(json.dumps(original))
                    mutate(config)
                    write_json(fixture.config, config)
                    before = snapshot(fixture.input_paths)
                    process, output = fixture.execute(name)
                    self.failed(process, output, reason)
                    self.assertEqual(snapshot(fixture.input_paths), before)


if __name__ == '__main__':
    unittest.main()
