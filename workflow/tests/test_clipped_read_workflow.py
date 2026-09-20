"""Native fictional controls for clipped-read recovery, fitting and remapping.

Uses actual BWA indexes/remaps and samtools CRAM conversion. Caller labels refer
to synthetic raw records; no caller or case data is used. The native Studio
qualification must run this module with no skips.
"""
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import random
import subprocess
import tempfile
import unittest

try:
    import pysam
except ImportError:
    pysam = None


WORKFLOW = Path(__file__).resolve().parents[1]
SCRIPT = WORKFLOW / 'scripts/review_clipped_reads.py'
SOURCE_NAMES = ('review_clipped_reads.py', 'clipped_read_helpers.py',
                'joint_haplotype_helpers.py', 'review_joint_haplotypes.py',
                'native_coding_records.py')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def pin(path):
    return dict(path=str(path.resolve()), bytes=path.stat().st_size, sha256=sha(path))


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def snapshot(paths):
    return {str(path): (path.read_bytes(), path.stat().st_mtime_ns) for path in paths}


def rc(sequence):
    return sequence.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]


class FictionalClipFixture:
    """149-base reads with a real two-base insertion and a distant exact copy."""
    def __init__(self, directory, runtime, mode='normal', alignment_format='cram'):
        self.directory, self.runtime = directory.resolve(), runtime.resolve()
        self.bin = self.runtime / 'tools/bio/bin'
        self.mode = mode
        rng = random.Random(20260920)
        self.local = ''.join(rng.choice('ACGT') for _ in range(6000))
        distant = list(''.join(rng.choice('ACGT') for _ in range(6000)))
        # Local insertion is after reference base 500. A copy on another contig
        # makes full-query compatibility intentionally insufficient for origin.
        self.inserted = 'AC'
        distant[2000:2452] = self.local[350:500] + self.inserted + self.local[500:800]
        self.distant = ''.join(distant)
        self.reference = self.directory / 'fictional.fa'
        self.reference.write_text('>1\n' + self.local + '\n>2\n' + self.distant + '\n')
        pysam.faidx(str(self.reference))
        self.native('bwa', ['index', str(self.reference)])
        suffixes = ('.amb', '.ann', '.bwt', '.pac', '.sa')
        self.bwa_index_record = self.directory / 'bwa-index.json'
        write_json(self.bwa_index_record, dict(reference_sha256=sha(self.reference),
            index_sha256={suffix: sha(Path(str(self.reference) + suffix)) for suffix in suffixes}))

        first = self.local[451:500] + self.inserted + self.local[500:598]
        second = self.local[480:500] + self.inserted + self.local[500:627]
        assert len(first) == len(second) == 149
        self.expected = {}
        reads = []

        def add(name, end, start, sequence, cigar, flag, mate_start=-1, qualities=True):
            quality = [25 + (i % 16) for i in range(len(sequence))] if qualities else None
            row = dict(name=name, end=end, start=start, sequence=sequence, cigar=cigar,
                       flag=flag, mate_start=mate_start, qualities=quality)
            reads.append(row)
            self.expected[(name, end)] = dict(original_sequence=sequence,
                original_qualities=quality, forward_sequence=rc(sequence) if flag & 16 else sequence,
                forward_qualities=list(reversed(quality)) if flag & 16 and quality is not None else quality)
            return row

        # Both ends overlap the seed window; fetching reciprocal mates must not
        # double their canonical identities. Reverse-end FASTQ uses RC and reversed BQ.
        first_row = add('alt-pair', 1, 451, first, '49M100S', 99, 480)
        add('alt-pair', 2, 480, second, '20M2I127M', 147, 451)
        add('retrieved-pair', 1, 451, first, '49M100S', 99, 900)
        outside_mate = add('retrieved-pair', 2, 900, self.local[900:1049], '149M', 147, 451)
        add('reference-read', 0, 451, self.local[451:600], '149M', 0)
        # This seed ends before the insertion: both local sequence hypotheses
        # must tie at zero despite complete-query fitting.
        add('pre-insertion-tie', 0, 341, self.local[341:490], '149M', 0)
        if mode == 'incomplete':
            ambiguous_sequence = first[:60] + 'N' + first[61:]
            add('N-query', 0, 451, ambiguous_sequence, '49M100S', 0)
            add('missing-quality', 0, 451, self.local[451:600], '149M', 0, qualities=False)
        elif mode == 'duplicate_seed':
            reads.append(dict(first_row))
        elif mode == 'ambiguous_mate':
            reads.append(dict(outside_mate))
        elif mode == 'nonreciprocal_mate':
            outside_mate['mate_start'] = 452
        elif mode != 'normal':
            raise ValueError('unknown fixture mode')
        self.expected_seed_records = sum(row['start'] < 525 and
            row['start'] + (49 if row['cigar'] == '49M100S' else 147 if row['cigar'] == '20M2I127M' else 149) > 475
            for row in reads)
        header = dict(HD=dict(VN='1.6', SO='coordinate'),
            SQ=[dict(SN='1', LN=6000), dict(SN='2', LN=6000)],
            RG=[dict(ID='RG1', SM='FICTIONAL', PL='ILLUMINA')])
        self.original_bam = self.directory / 'original.bam'
        with pysam.AlignmentFile(str(self.original_bam), 'wb', header=header) as writer:
            for row in sorted(reads, key=lambda item: item['start']):
                read = pysam.AlignedSegment(writer.header)
                read.query_name, read.query_sequence = row['name'], row['sequence']
                read.query_qualities = row['qualities']
                read.reference_id, read.reference_start = 0, row['start']
                read.cigarstring, read.flag, read.mapping_quality = row['cigar'], row['flag'], 60
                read.set_tag('RG', 'RG1')
                if row['flag'] & 1:
                    read.next_reference_id, read.next_reference_start = 0, row['mate_start']
                    # TLEN is explicit original context; assertions do not infer topology from it.
                    read.template_length = (627 - 451) if row['name'] == 'alt-pair' else (1049 - 451)
                    if row['end'] == 2:
                        read.template_length *= -1
                writer.write(read)
        self.native('samtools', ['index', str(self.original_bam)])
        self.original_cram = self.directory / 'original.cram'
        self.native('samtools', ['view', '-C', '--output-fmt-option', 'version=3.0',
            '-T', str(self.reference), '-o', str(self.original_cram), str(self.original_bam)])
        self.native('samtools', ['index', str(self.original_cram)])
        self.alignment = self.original_cram if alignment_format == 'cram' else self.original_bam
        index_suffix = '.crai' if alignment_format == 'cram' else '.bai'
        self.seed = self.directory / 'seed.bam'
        # Produce the seed from the chosen original to preserve decoded CRAM tags/order.
        with pysam.AlignmentFile(str(self.alignment), reference_filename=str(self.reference)) as original, \
                pysam.AlignmentFile(str(self.seed), 'wb', header=original.header) as writer:
            for read in original.fetch('1', 475, 525):
                writer.write(read)
        self.native('samtools', ['index', str(self.seed)])

        inputs = dict(reference=self.reference, reference_index=Path(str(self.reference) + '.fai'),
            alignment=self.alignment, alignment_index=Path(str(self.alignment) + index_suffix),
            seed_bam=self.seed, seed_bam_index=Path(str(self.seed) + '.bai'),
            bwa=self.bin / 'bwa', samtools=self.bin / 'samtools', bwa_index_record=self.bwa_index_record)
        inputs.update({'bwa_index' + suffix: Path(str(self.reference) + suffix) for suffix in suffixes})
        self.raw_expected = {}
        sources = []
        for caller in ('gatk', 'bcftools'):
            path = self.directory / ('gatk.raw.vcf' if caller == 'gatk' else 'bcftools.raw.bcf')
            header = pysam.VariantHeader()
            for chrom in ('1', '2'):
                header.contigs.add(chrom, length=6000)
            for key, number, kind in [('GT', 1, 'String'), ('AD', 'R', 'Integer'),
                                       ('DP', 1, 'Integer'), ('PL', 'G', 'Integer')]:
                header.formats.add(key, number, kind, 'Fictional original ' + key)
            if caller == 'gatk':
                header.formats.add('GQ', 1, 'Integer', 'Fictional original quality')
            header.add_sample('FICTIONAL')
            with pysam.VariantFile(str(path), 'w' if caller == 'gatk' else 'wb', header=header) as writer:
                for position, ref, alt in [(500, self.local[499], self.local[499] + self.inserted),
                                           (1800, self.local[1799], next(b for b in 'ACGT' if b != self.local[1799]))]:
                    record = writer.new_record(contig='1', start=position-1, alleles=(ref, alt))
                    sample = record.samples['FICTIONAL']
                    sample['GT'], sample['AD'], sample['DP'], sample['PL'] = (0, 1), (12, 3), 15, (90, 0, 150)
                    if caller == 'gatk':
                        sample['GQ'] = 35
                    writer.write(record)
            payload = hashlib.sha256()
            with pysam.VariantFile(str(path)) as reader:
                count = 0
                for ordinal, record in enumerate(reader, 1):
                    count = ordinal
                    payload.update(str(record).encode())
                    self.raw_expected[(caller, ordinal)] = dict(rendered_record=str(record),
                        FORMAT=json.loads(json.dumps(dict(record.samples['FICTIONAL']))))
            inputs[caller + '_raw'] = path
            sources.append(dict(input=caller + '_raw', caller=caller, expected_records=count,
                expected_rendered_record_payload_sha256=payload.hexdigest()))
        self.config = self.directory / 'config.json'
        self.config_data = dict(sample='FICTIONAL', expected_read_groups=1,
            inputs={name: pin(path) for name, path in inputs.items()},
            source_sha256={name: sha(WORKFLOW/'scripts'/name) for name in SOURCE_NAMES},
            seed_window=dict(chrom='1', start0=475, end0=525),
            targets=[dict(id='FICTIONAL', chrom='1', core_start0=499, core_end0=500)],
            caller_sources=sources, window_flanks_bp=[160], bwa_threads=1,
            maximum_seed_records=64, maximum_canonical_ends=128,
            maximum_hypothesis_states=32, maximum_haplotype_bases=5000,
            maximum_fit_cells=100000, maximum_total_fit_cells=5000000, maximum_tracebacks=32)
        self.save_config()
        self.input_paths = sorted({path for path in self.directory.iterdir() if path.is_file()} | {self.config})

    def native(self, tool, arguments):
        subprocess.run([str(self.bin/tool), *arguments], text=True, capture_output=True, check=True)

    def save_config(self):
        write_json(self.config, self.config_data)

    def execute(self, name='review'):
        output = self.directory / name
        process = subprocess.run([str(self.bin/'python'), str(SCRIPT), '--root', str(self.runtime),
            '--config', str(self.config), '--output', str(output)], text=True, capture_output=True)
        return process, output


class NativeClippedReadWorkflowTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.runtime = Path(os.environ.get('MVA_PROJECT_ROOT', '/Users/jeweled/Documents/mva-hackathon'))
        if pysam is None or not all((cls.runtime/'tools/bio/bin'/tool).is_file()
                                   for tool in ('python', 'samtools', 'bwa')):
            raise unittest.SkipTest('native Studio tools/bio Python/pysam, samtools and BWA required')

    def complete(self, process, output):
        diagnostics = '\n'.join(path.read_text(errors='replace')[-2000:] for path in output.glob('*.log'))
        self.assertEqual(process.returncode, 0, process.stdout + process.stderr + diagnostics)
        summary = json.loads((output/'summary.json').read_text())
        self.assertEqual(summary['status'], 'complete')
        self.assertEqual(json.loads((output/'status.json').read_text())['summary_sha256'], sha(output/'summary.json'))
        for name, expected in summary['outputs'].items():
            self.assertEqual(sha(output/name), expected['sha256'])
        return (summary, json.loads((output/'original-ends.json').read_text()),
                json.loads((output/'local-fits.json').read_text()), json.loads((output/'remapping.json').read_text()))

    def failed(self, process, output, reason):
        self.assertNotEqual(process.returncode, 0, process.stdout + process.stderr)
        status = json.loads((output/'status.json').read_text())
        self.assertEqual(status['status'], 'failed')
        self.assertIn(reason, status['error'])
        self.assertFalse((output/'summary.json').exists())

    def test_real_insertion_full_queries_distant_copy_mate_dedup_and_complete_native_remaps(self):
        with tempfile.TemporaryDirectory() as directory:
            fixture = FictionalClipFixture(Path(directory), self.runtime)
            before = snapshot(fixture.input_paths)
            process, output = fixture.execute()
            summary, originals, local, remapping = self.complete(process, output)
            self.assertEqual(summary['recovery']['raw_seed_records'], 5)
            self.assertEqual(summary['recovery']['primary_seed_ends'], 5)
            self.assertEqual(summary['recovery']['canonical_ends'], 6)
            self.assertEqual(summary['recovery']['unique_fragments'], 4)
            rows = {(row['QNAME'], row['end_number']): row for row in originals['ends'].values()}
            self.assertEqual(set(rows), set(fixture.expected))
            for identity, expected in fixture.expected.items():
                row = rows[identity]
                self.assertEqual(row['original_alignment_query'], expected['original_sequence'])
                self.assertEqual(row['raw_base_qualities'], expected['original_qualities'])
                self.assertEqual(row['forward_query'], expected['forward_sequence'])
                self.assertEqual(row['forward_base_qualities'], expected['forward_qualities'])
            self.assertEqual(rows[('alt-pair', 1)]['roles'], ['mate', 'seed'])
            self.assertEqual(rows[('alt-pair', 2)]['roles'], ['mate', 'seed'])
            self.assertEqual(rows[('retrieved-pair', 2)]['roles'], ['mate'])
            self.assertEqual(Counter(r['status'] for r in originals['mate_retrieval']),
                             {'retrieved_and_reciprocal': 3, 'unpaired_or_mate_reported_unmapped': 2})
            clipped = rows[('alt-pair', 1)]
            self.assertEqual(clipped['CIGAR'], '49M100S')
            self.assertEqual(clipped['terminal_clips'], [dict(side='right', query_start0=49, query_end0=149,
                sequence=fixture.expected[('alt-pair', 1)]['original_sequence'][49:],
                raw_base_qualities=fixture.expected[('alt-pair', 1)]['original_qualities'][49:])])

            contexts = json.loads((output/'native-context.json').read_text())
            self.assertEqual({r['caller']: r['records'] for r in contexts['source_census']}, {'gatk': 2, 'bcftools': 2})
            self.assertEqual(len(contexts['targets']['FICTIONAL']), 2)
            for entry in contexts['targets']['FICTIONAL']:
                context = entry['context']
                for key, value in fixture.raw_expected[(context['caller'], context['source_ordinal'])].items():
                    self.assertEqual(context[key], value)
                self.assertEqual(context['source_ordinal'], 1)
                self.assertEqual(context['FORMAT']['PL'], [90, 0, 150])
                self.assertEqual(context['FORMAT']['GT'], [0, 1])
                self.assertEqual(context['FILTER'], [])
                if context['caller'] == 'bcftools':
                    self.assertNotIn('GQ', context['FORMAT'])
            self.assertEqual(len(local['panels']), 1)
            panel = local['panels'][0]
            expected_reference = fixture.local[339:660]
            expected_alternate = fixture.local[339:500] + 'AC' + fixture.local[500:660]
            self.assertEqual({h['sequence'] for h in panel['hypotheses'].values()},
                             {expected_reference, expected_alternate})
            alternate_id = next(i for i,h in panel['hypotheses'].items() if h['sequence'] == expected_alternate)
            alternate = panel['hypotheses'][alternate_id]
            self.assertEqual(sorted((lineage['caller'], len(lineage['edits'])) for lineage in alternate['source_combinations']),
                             [('bcftools', 1), ('gatk', 1)])
            comparisons = {r['canonical_end_id']: r for r in panel['end_comparisons']}
            for identity in [('alt-pair', 1), ('alt-pair', 2), ('retrieved-pair', 1)]:
                comparison = comparisons[rows[identity]['id']]
                self.assertEqual(comparison['minimum_cost'], 0)
                self.assertEqual(comparison['reference_cost'], 2)
                self.assertEqual(comparison['best_hypothesis_ids'], [alternate_id])
                self.assertFalse(comparison['quality_used_in_score'])
                best_fit = local['fits'][comparison['fits'][alternate_id]]['result']
                self.assertEqual(best_fit['query_length'], 149)
                self.assertTrue(all(a['query_start0'] == 0 and a['query_end0'] == 149 for a in best_fit['alignments']))
            reference_comparison = comparisons[rows[('reference-read', 0)]['id']]
            self.assertEqual(reference_comparison['reference_cost'], 0)
            self.assertEqual(reference_comparison['costs'][alternate_id], 2)
            tie = comparisons[rows[('pre-insertion-tie', 0)]['id']]
            self.assertEqual(set(tie['best_hypothesis_ids']), set(panel['hypotheses']))
            self.assertTrue(all(cost == 0 for cost in tie['costs'].values()))
            self.assertNotIn(rows[('retrieved-pair', 2)]['id'], comparisons)

            for label, expected_ends, expected_pairs in [('single', 6, None), ('paired', 4, 2)]:
                arm = remapping[label]
                self.assertEqual(arm['status'], 'complete')
                self.assertEqual((arm['exported_ends'], arm['exported_pairs']), (expected_ends, expected_pairs))
                self.assertTrue(arm['primary_query_payloads_verified'])
                identity_map = {(x['query'], x['end_number']): x['canonical_end_id'] for x in arm['export_identity']}
                primary_counts = Counter()
                sam = output/('single-whole-reference.sam' if label == 'single' else 'paired-whole-reference.sam')
                with pysam.AlignmentFile(str(sam), 'r') as reader:
                    records = list(reader)
                self.assertEqual(len(records), len(arm['reported_alignments']))
                for read in records:
                    number = 1 if read.is_read1 else 2 if read.is_read2 else 0
                    canonical_id = identity_map[(read.query_name, number)]
                    original = originals['ends'][canonical_id]
                    if not read.is_secondary and not read.is_supplementary:
                        primary_counts[canonical_id] += 1
                        self.assertEqual(read.get_forward_sequence(), original['forward_query'])
                        self.assertEqual(list(read.get_forward_qualities()), original['forward_base_qualities'])
                self.assertEqual(primary_counts, {identifier: 1 for identifier in identity_map.values()})
            # Exact distant copy is a planted competing locus. Do not require
            # a particular tie-breaking primary or exhaustive BWA alternatives.
            alt_alignments = [r for r in remapping['single']['reported_alignments']
                              if r['canonical_end_id'] == clipped['id']]
            self.assertTrue(any(r['chrom'] == '2' and r['NM'] == 0 and r['aligned_query_bases'] == 149
                                for r in alt_alignments), alt_alignments)
            reference_alignments = [r for r in remapping['single']['reported_alignments']
                                    if r['canonical_end_id'] == rows[('reference-read', 0)]['id']]
            self.assertTrue(any(r['chrom'] == '1' and r['NM'] == 0 and r['aligned_query_bases'] == 149
                                for r in reference_alignments), reference_alignments)
            self.assertEqual(snapshot(fixture.input_paths), before)

    def test_N_and_missing_quality_are_retained_with_explicit_different_eligibility(self):
        with tempfile.TemporaryDirectory() as directory:
            fixture = FictionalClipFixture(Path(directory), self.runtime, mode='incomplete', alignment_format='bam')
            before = snapshot(fixture.input_paths)
            summary, originals, local, remapping = self.complete(*fixture.execute())
            rows = {row['QNAME']: row for row in originals['ends'].values() if row['end_number'] == 0}
            ambiguous, missing = rows['N-query'], rows['missing-quality']
            self.assertIn('query_missing_or_nonconcrete_ACGT', ambiguous['local_fit_ineligibility'])
            self.assertIn('N', ambiguous['original_alignment_query'])
            self.assertIsNone(missing['raw_base_qualities'])
            self.assertIsNone(missing['forward_base_qualities'])
            self.assertEqual(missing['remap_ineligibility'], 'quality_missing_or_not_FASTQ_representable')
            comparisons = {row['canonical_end_id']: row for row in local['panels'][0]['end_comparisons']}
            self.assertEqual(comparisons[ambiguous['id']]['status'], 'ineligible')
            self.assertEqual(comparisons[ambiguous['id']]['fits'], {})
            # Unit-distance fitting explicitly ignores BQ; missing BQ is not imputed.
            self.assertEqual(comparisons[missing['id']]['status'], 'scored_descriptively')
            self.assertFalse(comparisons[missing['id']]['quality_used_in_score'])
            exported = {row['canonical_end_id'] for row in remapping['single']['export_identity']}
            self.assertNotIn(missing['id'], exported)
            self.assertIn(ambiguous['id'], exported)  # BWA may remap N-containing FASTQ; local DP does not.
            self.assertEqual(summary['recovery']['canonical_ends'], 8)
            self.assertEqual(remapping['single']['exported_ends'], 7)
            self.assertEqual(snapshot(fixture.input_paths), before)

    def test_ambiguous_or_nonreciprocal_mate_stays_unresolved_and_is_not_paired(self):
        for mode, status in [('ambiguous_mate', 'mate_lookup_missing_or_ambiguous'),
                             ('nonreciprocal_mate', 'mate_coordinates_not_reciprocal')]:
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as directory:
                fixture = FictionalClipFixture(Path(directory), self.runtime, mode=mode)
                before = snapshot(fixture.input_paths)
                summary, originals, _, remapping = self.complete(*fixture.execute())
                seed = next(row for row in originals['ends'].values()
                            if (row['QNAME'], row['end_number']) == ('retrieved-pair', 1))
                retrieval = next(row for row in originals['mate_retrieval'] if row['seed_id'] == seed['id'])
                self.assertEqual(retrieval['status'], status)
                if mode == 'ambiguous_mate':
                    self.assertEqual(retrieval['candidate_records'], 2)
                self.assertFalse(any((row['QNAME'],row['end_number']) == ('retrieved-pair',2)
                                     for row in originals['ends'].values()))
                self.assertEqual(summary['recovery']['canonical_ends'], 5)
                self.assertEqual(remapping['paired']['exported_pairs'], 1)
                self.assertEqual(remapping['paired']['exported_ends'], 2)
                self.assertNotIn(seed['id'], {x['canonical_end_id'] for x in remapping['paired']['export_identity']})
                self.assertEqual(snapshot(fixture.input_paths), before)

    def test_repeated_original_primary_end_identity_is_refused(self):
        with tempfile.TemporaryDirectory() as directory:
            fixture = FictionalClipFixture(Path(directory), self.runtime, mode='duplicate_seed')
            before = snapshot(fixture.input_paths)
            process, output = fixture.execute()
            self.failed(process, output, 'repeated original primary end identity')
            self.assertFalse((output/'remapping.json').exists())
            self.assertEqual(snapshot(fixture.input_paths), before)

    def test_bad_original_input_hash_is_refused_before_recovery_or_remapping(self):
        with tempfile.TemporaryDirectory() as directory:
            fixture = FictionalClipFixture(Path(directory), self.runtime)
            fixture.config_data['inputs']['alignment']['sha256'] = '0' * 64
            fixture.save_config()
            before = snapshot(fixture.input_paths)
            process, output = fixture.execute()
            self.failed(process, output, 'input SHA differs')
            self.assertFalse((output/'original-ends.json').exists())
            self.assertFalse((output/'commands.json').exists())
            self.assertEqual(snapshot(fixture.input_paths), before)


if __name__ == '__main__':
    unittest.main()
