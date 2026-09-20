"""Native CRAM compatibility, held-out variants and evidence-loss rejection."""
import argparse
import gzip
import importlib.util
import itertools
import json
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest

import pysam

SPEC = importlib.util.spec_from_file_location('compat', Path(__file__).parents[1] /
                                             'scripts/prepare_gatk_compatible_cram.py')
compat = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(compat)


def identity(read):
    return (read.query_name, read.flag, read.reference_id, read.reference_start,
            read.mapping_quality, read.cigarstring, read.next_reference_id,
            read.next_reference_start, read.template_length, read.query_sequence,
            list(read.query_qualities) if read.query_qualities is not None else None,
            sorted(read.get_tags()))


class NativeCompatibilityTests(unittest.TestCase):
    def test_cram30_preserves_evidence_and_two_caller_truth(self):
        root = Path(os.environ.get('MVA_PROJECT_ROOT', '/Users/jeweled/Documents/mva-hackathon'))
        samtools = root / 'tools/bio/bin/samtools'
        bwa = root / 'tools/bio/bin/bwa'
        java = root / 'tools/bio/lib/jvm/bin/java'
        jar = root / 'gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar'
        if not all(p.exists() for p in (samtools, bwa, java, jar)):
            self.skipTest('native Studio runtime required')

        def command(args, **kwargs):
            result = subprocess.run(list(map(str, args)), capture_output=True, **kwargs)
            self.assertEqual(result.returncode, 0, result.stderr[-2500:])
            return result

        with tempfile.TemporaryDirectory(prefix='mva-synthetic-cram-') as tmp:
            folder = Path(tmp).resolve()
            rng = random.Random(20260920)
            sequences = {c: ''.join(rng.choice('ACGT') for _ in range(3000)) for c in ('1', '2')}
            snv = 1000
            alt = 'A' if sequences['1'][snv] != 'A' else 'C'
            deletion = 1000
            while sequences['2'][deletion] == sequences['2'][deletion + 2]:
                deletion += 1
            mutated = {'1': sequences['1'][:snv] + alt + sequences['1'][snv + 1:],
                       '2': sequences['2'][:deletion + 1] + sequences['2'][deletion + 3:]}
            truth = {('1', snv + 1, sequences['1'][snv], alt),
                     ('2', deletion + 1, sequences['2'][deletion:deletion + 3], sequences['2'][deletion])}
            reference = folder / 'reference.fa'
            reference.write_text(''.join('>' + c + '\n' + seq + '\n' for c, seq in sequences.items()))
            pysam.faidx(str(reference))
            command([samtools, 'dict', '-o', reference.with_suffix('.dict'), reference])
            command([bwa, 'index', reference])
            fastqs = [folder / 'R1.fastq.gz', folder / 'R2.fastq.gz']
            for mate, path in enumerate(fastqs, 1):
                with gzip.open(path, 'wt') as stream:
                    for contig in sequences:
                        for allele in (0, 1):
                            genome = (sequences if allele == 0 else mutated)[contig]
                            for index in range(24):
                                start = 900 + index * 2 + (220 if mate == 2 else 0)
                                read = genome[start:start + 150]
                                if mate == 2:
                                    read = read.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
                                name = f'SYNTH:{contig}:{allele}:{index}'
                                stream.write('@' + name + f' {mate}:N:0:ACGT\n' + read + '\n+\n' + 'I' * len(read) + '\n')
            sam = folder / 'reads.sam'
            result = command([bwa, 'mem', '-t', '2', '-R', r'@RG\tID:SIM\tSM:SYNTH\tPL:ILLUMINA', reference, *fastqs])
            sam.write_bytes(result.stdout)
            unsorted = folder / 'decorated.bam'
            with pysam.AlignmentFile(str(sam)) as reader, \
                    pysam.AlignmentFile(str(unsorted), 'wb', header=reader.header) as writer:
                original = list(reader)
                for index, read in enumerate(original):
                    read.set_tag('ZZ', 'held-out-tag-' + str(index))
                    read.set_tag('ZI', index)
                    writer.write(read)
                # Cover duplicate, supplementary and secondary evidence as well.
                for index, flag in enumerate((0x400, 0x800, 0x100)):
                    read = pysam.AlignedSegment.fromstring(original[index].to_string(), reader.header)
                    read.query_name = 'SYNTH_EXTRA_' + str(index)
                    read.flag |= flag
                    writer.write(read)
            bam = folder / 'reads.bam'
            command([samtools, 'sort', '-o', bam, unsorted])
            command([samtools, 'index', bam])
            source = folder / 'source31.cram'
            command([samtools, 'view', '-T', reference, '-O', 'CRAM,version=3.1,store_md=1,store_nm=1',
                     '--write-index', '-o', source, bam])
            count = len(original) + 3
            failed = subprocess.run(list(map(str, [java, '-XX:ActiveProcessorCount=2', '-Xmx2g', '-jar', jar,
                                                   'CountReads', '-R', reference, '-I', source, '-L', '1:900-1300'])),
                                    capture_output=True, text=True)
            self.assertNotEqual(failed.returncode, 0)
            self.assertIn('CRAM version 3.1 is not supported', failed.stderr)
            args = argparse.Namespace(root=root, source=source, reference=reference,
                                      output=folder / 'compatible', sample='SYNTH', threads=2,
                                      expected_source_sha256=compat.digest(source),
                                      expected_reference_sha256=compat.digest(reference),
                                      expected_records=count, probe_region='1:900-1300')
            summary = compat.prepare(args)
            converted = Path(summary['cram_path'])
            self.assertEqual(summary['decoded_record_verification']['records_compared'], count)
            self.assertEqual(compat.cram_version(converted), [3, 0])
            with pysam.AlignmentFile(str(source), reference_filename=str(reference)) as before, \
                    pysam.AlignmentFile(str(converted), reference_filename=str(reference)) as after:
                for a, b in itertools.zip_longest(before, after):
                    self.assertIsNotNone(a)
                    self.assertIsNotNone(b)
                    self.assertEqual(identity(a), identity(b))
            with self.assertRaisesRegex(ValueError, 'existing compatibility output'):
                compat.prepare(args)
            check = [samtools, 'checksum', '-P', '-C', '-M', '-c', '-O', '-O', '-b', '0xffff',
                     '-f', '0', '-F', '0', '-z', 'off', '-t', '*', '-q', '-v', '-T',
                     '--input-fmt-option', 'reference=' + str(reference)]
            original_checksums = folder / 'original.checksums.txt'
            command(check + ['-o', original_checksums, source])
            # Real changed records must fail the same full-data verification gate.
            for field in ('quality', 'tag', 'duplicate_flag', 'order'):
                changed = folder / (field + '.bam')
                with pysam.AlignmentFile(str(converted), reference_filename=str(reference)) as reader:
                    reads = list(reader)
                    if field == 'quality':
                        qualities = reads[0].query_qualities
                        qualities[0] -= 1
                        reads[0].query_qualities = qualities
                    elif field == 'tag':
                        reads[0].set_tag('ZZ', 'deliberate-change')
                    elif field == 'duplicate_flag':
                        reads[0].flag ^= 0x400
                    else:
                        reads[0], reads[1] = reads[1], reads[0]
                    with pysam.AlignmentFile(str(changed), 'wb', header=reader.header) as writer:
                        for read in reads:
                            writer.write(read)
                receipt = folder / (field + '.checksums.txt')
                command(check + ['-o', receipt, changed])
                if field == 'duplicate_flag':
                    # Reproduce the observed samtools 1.24 checksum blind spot.
                    compat.compare_checksums(original_checksums, receipt, count)
                else:
                    with self.assertRaisesRegex(ValueError, 'checksum mismatch', msg=field):
                        compat.compare_checksums(original_checksums, receipt, count)
                direct = subprocess.run(list(map(str, [args.output / 'verify-alignment-records',
                                                       reference, source, changed, count, 2])),
                                        capture_output=True, text=True)
                self.assertNotEqual(direct.returncode, 0, field)
                self.assertIn('Verification failed', direct.stderr)
            output = folder / 'calls'
            command([sys.executable, Path(__file__).parents[1] / 'scripts/call_genome.py', '--root', root,
                     '--reference', reference, '--alignment', converted, '--sample', 'SYNTH',
                     '--output', output, '--workers', '2', '--shard-max-bases', '4000'])
            status = json.loads((output / 'status.json').read_text())
            self.assertEqual(status['status'], 'complete')
            for method, result in status['outputs'].items():
                with pysam.VariantFile(result['path']) as calls:
                    observed = {(row.chrom, row.pos, row.ref, row.alts[0]): row.samples['SYNTH']['GT']
                                for row in calls.fetch()}
                for site in truth:
                    self.assertEqual(tuple(sorted(observed[site])), (0, 1), (method, site))


if __name__ == '__main__':
    unittest.main()
