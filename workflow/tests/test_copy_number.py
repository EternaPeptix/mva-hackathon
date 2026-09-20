"""Synthetic depth truth, imported-count guards, and final-alignment dependency checks."""
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
import pysam

sys.path.insert(0, str(Path(__file__).parents[1] / 'scripts'))
import run_copy_number as cnv


def make_read(header, name, start, flag=0, mq=60, cigar='150M', group='L1'):
    read = pysam.AlignedSegment(header)
    read.query_name = name
    read.query_sequence = 'A' * 150
    read.query_qualities = [40] * 150
    read.reference_id = 0
    read.reference_start = start
    read.cigarstring = cigar
    read.mapping_quality = mq
    read.flag = flag
    if group is not None:
        read.set_tag('RG', group)
    return read


def write_fixture(folder, scale=1, bin_sizes=None):
    """Aligned read-density simulation, not an end-to-end SV/mapping benchmark."""
    folder.mkdir()
    rng = np.random.default_rng(20260919)
    length = 2_000_000 * scale
    reference = folder / 'synthetic.fa'
    sequences = {}
    with reference.open('w') as stream:
        for chrom in ('1', '2'):
            sequence = rng.choice(np.frombuffer(b'ACGT', dtype='S1'), length).tobytes().decode()
            sequences[chrom] = sequence
            stream.write('>' + chrom + '\n' + sequence + '\n')
    pysam.faidx(str(reference))
    truth = [{'chrom': '1', 'start_0based': 300_000, 'end_exclusive': 400_000,
              'type': 'deletion', 'relative_depth': 0.5},
             {'chrom': '1', 'start_0based': 1_200_000, 'end_exclusive': 1_350_000,
              'type': 'duplication', 'relative_depth': 1.5}]
    for row in truth:
        row['start_0based'] *= scale
        row['end_exclusive'] *= scale
    header = {'HD': {'VN': '1.6', 'SO': 'coordinate'},
              'SQ': [{'SN': c, 'LN': length} for c in ('1', '2')],
              'RG': [{'ID': 'L1', 'SM': 'SYNTH', 'PL': 'ILLUMINA'}]}
    bam = folder / 'synthetic.bam'
    counts = Counter()
    with pysam.AlignmentFile(str(bam), 'wb', header=header) as writer:
        for rid, chrom in enumerate(('1', '2')):
            for begin in range(0, length, 100):
                scale = next((r['relative_depth'] for r in truth if r['chrom'] == chrom
                              and r['start_0based'] <= begin < r['end_exclusive']), 1.0)
                for j, center in enumerate(sorted(rng.integers(begin, begin + 100, rng.poisson(20 * scale)))):
                    start = min(length - 150, max(0, int(center) - 75))
                    read = make_read(writer.header, f'synth_{chrom}_{begin}_{j}', start)
                    read.reference_id = rid
                    read.query_sequence = sequences[chrom][start:start + 150]
                    read.query_qualities = [40] * 150
                    writer.write(read)
                    counts[chrom] += 1
    pysam.index(str(bam))
    config = {'scope': 'synthetic', 'sample': 'SYNTH',
              'reference_sha256': hashlib.sha256(reference.read_bytes()).hexdigest(),
              'chromosomes': ['1', '2'], 'bin_sizes': bin_sizes or [1000, 10000],
              'mask_policy': 'unmasked_discovery_with_q0_and_reference_gap_metrics'}
    (folder / 'truth.json').write_text(json.dumps({'seed': 20260919, 'truth': truth,
        'read_records': dict(counts), 'method': 'Poisson aligned read-midpoint density; no realignment, GC-bias, tumor or mosaic simulation.'}, indent=2) + '\n')
    (folder / 'config.json').write_text(json.dumps(config, indent=2) + '\n')
    return bam, reference, config, truth, dict(counts)


class InputGuards(unittest.TestCase):
    def make_edge_bam(self, root, missing_group=False):
        path = root / 'edges.bam'
        header = {'HD': {'VN': '1.6', 'SO': 'coordinate'}, 'SQ': [{'SN': '1', 'LN': 2000}],
                  'RG': [{'ID': 'L1', 'SM': 'SYNTH'}]}
        with pysam.AlignmentFile(str(path), 'wb', header=header) as writer:
            writer.write(make_read(writer.header, 'start_zero', 0, group=None if missing_group else 'L1'))
            writer.write(make_read(writer.header, 'mq0', 200, mq=0))
            writer.write(make_read(writer.header, 'mq255', 300, mq=255))
            for i, flag in enumerate((0x100, 0x200, 0x400, 0x800)):
                writer.write(make_read(writer.header, 'flag_' + str(flag), 400 + i * 100, flag=flag))
            writer.write(make_read(writer.header, 'reference_skip', 900, cigar='75M100N75M'))
        pysam.index(str(path))
        return path

    def test_unconditional_flags_mq_and_coordinate_zero(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self.make_edge_bam(Path(tmp))
            with pysam.AlignmentFile(str(path)) as reader:
                parity, positive, counts = cnv.count_chromosome(reader, '1', 2000, {'L1'})
            self.assertEqual(int(parity.sum()), 2)
            self.assertEqual(int(positive.sum()), 1)
            self.assertEqual(int(parity[0]), 1)
            self.assertEqual(counts['excluded_flags'], 4)
            self.assertEqual(counts['excluded_unknown_MQ_255'], 1)
            self.assertEqual(counts['excluded_reference_skip'], 1)
            self.assertEqual(counts['included_reference_start_zero'], 1)
            from cnvpytor.bam import Bam
            native = Bam(str(path), init_log=False)
            native_p, _, _ = native.read_chromosome('1'); native.file.close()
            # Pin the observed source-version discrepancy rather than silently inheriting it.
            self.assertEqual(int(native_p.sum()), 7)
            self.assertEqual(int(native_p[0]), 0)

    def test_missing_read_group_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self.make_edge_bam(Path(tmp), missing_group=True)
            with pysam.AlignmentFile(str(path)) as reader:
                with self.assertRaisesRegex(ValueError, 'read group'):
                    cnv.count_chromosome(reader, '1', 2000, {'L1'})

    def test_completion_requires_exact_finished_four_lane_receipt(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); cram = root / 'final.cram'; state = root / 'state.json'
            config = {'scope': 'verified_four_lane_autosomal_discovery'}
            expected = {'bytes': 100, 'sha256': 'expected'}
            good = {'status': 'complete', 'expected_lanes': 4, 'cram_path': str(cram),
                    'cram_bytes': 100, 'cram_sha256': 'expected'}
            with self.assertRaisesRegex(ValueError, 'receipt required'):
                cnv.verify_completion(config, cram, expected, None)
            for field, value in [('status', 'running'), ('expected_lanes', 2),
                                 ('cram_sha256', 'wrong'), ('cram_bytes', 99), ('cram_path', str(root/'other.cram'))]:
                state.write_text(json.dumps(dict(good, **{field: value})))
                with self.assertRaisesRegex(ValueError, 'does not match'):
                    cnv.verify_completion(config, cram, expected, state)
            state.write_text(json.dumps(good))
            self.assertEqual(cnv.verify_completion(config, cram, expected, state)['path'], str(state.resolve()))

    def test_clean_midpoint_counts_match_native_importer(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); path = root / 'clean.bam'
            header = {'HD': {'VN': '1.6', 'SO': 'coordinate'}, 'SQ': [{'SN': '1', 'LN': 2000}],
                      'RG': [{'ID': 'L1', 'SM': 'SYNTH'}]}
            with pysam.AlignmentFile(str(path), 'wb', header=header) as writer:
                for i, (start, mq) in enumerate([(20, 60), (99, 0), (100, 10), (190, 60), (1840, 60)]):
                    writer.write(make_read(writer.header, 'clean_' + str(i), start, mq=mq))
            pysam.index(str(path))
            with pysam.AlignmentFile(str(path)) as reader:
                parity, positive, _ = cnv.count_chromosome(reader, '1', 2000, {'L1'})
            from cnvpytor.bam import Bam
            native = Bam(str(path), init_log=False)
            native_p, native_u, _ = native.read_chromosome('1'); native.file.close()
            np.testing.assert_array_equal(parity, native_p)
            np.testing.assert_array_equal(positive, native_u)

    def test_scope_and_bins_do_not_silently_shrink(self):
        good = {'scope': 'synthetic', 'chromosomes': ['1', '2'], 'bin_sizes': [1000, 10000],
                'mask_policy': 'unmasked_discovery_with_q0_and_reference_gap_metrics'}
        cnv.validate_definition(good, {'1': 2000, '2': 2000})
        for changes in ({'scope': 'verified_four_lane_autosomal_discovery'}, {'bin_sizes': [1001]},
                        {'chromosomes': ['1', 'X']}, {'chromosomes': ['1', '1']}):
            with self.assertRaises(ValueError):
                cnv.validate_definition(dict(good, **changes), {'1': 2000, '2': 2000, 'X': 2000})


class NativeCopyNumber(unittest.TestCase):
    def check_truth(self, scale=1, bin_sizes=None):
        with tempfile.TemporaryDirectory() as tmp:
            parent = Path(os.environ.get('MVA_CNV_TEST_ARTIFACT_DIR', tmp))
            if scale != 1:
                parent = parent / 'coarse'
            parent.mkdir(parents=True, exist_ok=True)
            fixture = parent / 'fixture'
            bam, reference, config, truth, counts = write_fixture(fixture, scale, bin_sizes)
            result = cnv.run_analysis(bam, reference, config, parent / 'analysis')
            self.assertTrue(result['counts_round_trip_verified'])
            self.assertEqual({c: r['included_read_midpoints'] for c, r in result['read_midpoint_import'].items()}, counts)
            calls = json.loads((parent / 'analysis/calls.json').read_text())['calls']
            matched = []
            for size, rows in calls.items():
                for target in truth:
                    overlaps = []
                    for row in rows:
                        if row['type'] != target['type'] or row['chrom'] != target['chrom']:
                            continue
                        overlap = max(0, min(row['end_inclusive'], target['end_exclusive'])
                                      - max(row['start_1based'] - 1, target['start_0based']))
                        truth_span = target['end_exclusive'] - target['start_0based']
                        if overlap >= 0.8 * truth_span and overlap >= 0.8 * row['size_bp']:
                            overlaps.append(row)
                    self.assertTrue(overlaps, (size, target, rows))
                    row = min(overlaps, key=lambda r: abs(r['relative_depth'] - target['relative_depth']))
                    self.assertAlmostEqual(row['relative_depth'], target['relative_depth'], delta=0.15)
                    self.assertTrue(row['within_reference_bounds'])
                    matched.append({'bin_size': int(size), 'truth': target, 'matched_call': row})
                self.assertFalse([r for r in rows if r['chrom'] == '2' and r['size_bp'] >= 50_000], rows)
            cnv.save(parent / 'synthetic-validation.json', {
                'matched': matched, 'control_chromosome_calls_over_50kb': 0,
                'result_summary_sha256': cnv.fingerprint(parent / 'analysis/summary.json')['sha256'],
                'scope': 'Positive aligned-density fixtures and a copy-neutral control at configured bin sizes; not clinical sensitivity/specificity or end-to-end variant truth.'})
            for label, changes, message in [('reference', {'reference_sha256': 'wrong'}, 'reference checksum'),
                                             ('sample', {'sample': 'WRONG'}, 'sample metadata')]:
                bad_output = parent / ('rejected-' + label)
                with self.assertRaisesRegex(ValueError, message):
                    cnv.run_analysis(bam, reference, dict(config, **changes), bad_output)
                self.assertEqual(json.loads((bad_output / 'status.json').read_text())['status'], 'failed')

    def test_pipeline_recovers_density_deletion_duplication_and_control(self):
        self.check_truth()

    def test_coarse_bins_recover_larger_density_changes_and_control(self):
        self.check_truth(scale=5, bin_sizes=[100000])


if __name__ == '__main__':
    unittest.main()
