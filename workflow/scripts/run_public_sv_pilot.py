#!/usr/bin/env python3
"""Run a pinned public HG002 chr20 structural-caller pilot on Mac Studio.

This entry point only accepts the publisher-checksummed public BAM. It does not
relax the independent four-lane patient runner. The entire chr20 truth context is
retained, with the official chr20 benchmark BED as the primary evaluation region.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess
import time

import pysam

from build_structural_exclusions import (TEMPLATE_SHA, fingerprint, native_query_intervals,
                                        parse_exclusions, query_rows, save, stamp, write_bed)
from prepare_public_sv_benchmark import ALIGNMENT, REFERENCE, FILES, checksums
from run_structural_variants import NATIVE_OPTIONS, alignment_census, review_calls

DELLY_SHA = '09a83090cdb28337f7647b94ce04396373fe244321048ae22d74a773cf18bac0'
SAMTOOLS_SHA = '7518cbb230a35393d7d3a9f02fc8b4519e4a37b5a17fbed06b9752ca1d962a91'
TRUTH_REFERENCE_MD5 = '939ce19062d1462c09b88c55faca4d76'
CHROM = 'chr20'
SAMPLE = 'HG002'


def public_identity(folder):
    receipt = json.loads((folder / 'input-receipt.json').read_text())
    if receipt['scope'] != 'public_HG002_chr20_prealigned_Illumina_35x_pilot' or receipt['status'] != 'inputs_verified':
        raise ValueError('verified public input receipt required')
    identities = {}
    for url, expected_md5, size in FILES:
        p = folder / url.rsplit('/', 1)[1]; d = checksums(p)
        saved = next(row for row in receipt['files'] if row['url'] == url)
        if d['md5'] != expected_md5 or d['sha256'] != saved['sha256'] or (size is not None and d['bytes'] != size):
            raise ValueError('public input differs from publisher or saved checksum')
        identities[p.name] = d
    return identities


def pilot(root, output):
    if output.resolve().is_relative_to(Path(__file__).resolve().parents[2]):
        raise ValueError('raw public benchmark outputs must remain outside the repository')
    output.mkdir(parents=True, exist_ok=False); started = time.monotonic()
    state = {'status': 'running', 'scope': 'public_HG002_chr20_prealigned_caller_component', 'started_at': stamp()}
    save(output / 'status.json', state); commands = []
    def command(label, args):
        state['stage'] = label; save(output / 'status.json', state)
        log = output / ('step-%02d.log' % len(commands)); tick = time.monotonic()
        with log.open('w') as handle: result = subprocess.run(list(map(str, args)), stdout=handle, stderr=subprocess.STDOUT)
        commands.append({'stage': label, 'command': list(map(str, args)), 'exit_code': result.returncode,
                         'elapsed_seconds': time.monotonic() - tick, 'log_sha256': fingerprint(log)['sha256']})
        save(output / 'commands.json', commands)
        if result.returncode: raise RuntimeError(label + ' failed; log retained')
    try:
        inputs_dir = root / 'ref/giab-sv-20260919/inputs'
        inputs = public_identity(inputs_dir)
        reference = inputs_dir / REFERENCE; alignment = inputs_dir / ALIGNMENT
        truth_ref = root / 'ref/giab-sv-20260919/truth-reference/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz'
        truth_identity = checksums(truth_ref)
        if truth_identity['md5'] != TRUTH_REFERENCE_MD5: raise ValueError('benchmark reference checksum mismatch')
        tool_paths = {'delly': root/'tools/delly-2.6.0/bin/delly', 'samtools': root/'tools/bio/bin/samtools',
                      'native_parser': root/'results/delly-exclusion-parser-20260919/parser'}
        tools = {name: fingerprint(path) for name, path in tool_paths.items()}
        if tools['delly']['sha256'] != DELLY_SHA or tools['samtools']['sha256'] != SAMTOOLS_SHA:
            raise ValueError('unreviewed public pilot tool version')
        if tools['native_parser']['sha256'] != '3669695b1dab52341fc5f7bdce1b96585ee984520f4cfb351596aaf52dd7acfc':
            raise ValueError('unreviewed native exclusion parser')
        template = root/'tools/cache/delly-review-20260919/delly-2.6.0/excludeTemplates/human.hg38.excl.tsv'
        if fingerprint(template)['sha256'] != TEMPLATE_SHA: raise ValueError('template checksum mismatch')
        command('public alignment format check', [tool_paths['samtools'], 'quickcheck', '-v', alignment])
        with pysam.FastaFile(str(reference)) as ref, pysam.FastaFile(str(truth_ref)) as tref:
            dictionary = dict(zip(ref.references, ref.lengths)); seq = ref.fetch(CHROM).upper()
            if len(seq) != 64444167 or seq != tref.fetch(CHROM).upper():
                raise ValueError('chr20 alignment reference differs from benchmark reference')
            with pysam.FastaFile(str(root/'ref/vendor/GRCh38_nochr.fa')) as patient_ref:
                same_as_patient_reference = seq == patient_ref.fetch('20').upper()
            ref_check = {'chromosome': CHROM, 'bases_compared': len(seq),
                         'sha256_uppercase_sequence': hashlib.sha256(seq.encode()).hexdigest(),
                         'alignment_and_benchmark_chr20_identical': True,
                         'patient_reference_chr20_identical': same_as_patient_reference,
                         'note': 'Sequence equality on chr20 only. Whole reference dictionaries and masking policies differ.'}
            selected = []
            for line in template.read_text().splitlines():
                row = line.split('\t')
                if row[0] == CHROM: selected.append(row[:3])
            native_mask = output/'native-exclusions.tsv'
            native_mask.write_text(''.join('\t'.join(row)+'\n' for row in selected)
                                   + ''.join(chrom+'\n' for chrom in dictionary if chrom != CHROM))
            intervals, whole = parse_exclusions(native_mask, dictionary)
            expected_queries = query_rows(dictionary, intervals, whole)
            proc = subprocess.run([str(tool_paths['native_parser']), str(native_mask), str(reference)+'.fai'], capture_output=True, text=True)
            (output/'native-parser.stdout').write_text(proc.stdout); (output/'native-parser.stderr').write_text(proc.stderr)
            observed = [[f[0], int(f[1]), int(f[2])] for line in proc.stdout.splitlines() if (f := line.split('\t'))]
            if proc.returncode or observed != expected_queries: raise ValueError('public native exclusion validation failed')
            write_bed(output/'query-intervals.bed', expected_queries)
            bed = inputs_dir/'HG002_GRCh38_v5.0q_stvar.benchmark.bed'
            regions = []
            for line in bed.read_text().splitlines():
                row = line.split('\t')
                if row[0] == CHROM: regions.append([CHROM, int(row[1]), int(row[2])])
            if not regions or any(not 0 <= a < b <= len(seq) for _, a, b in regions): raise ValueError('invalid benchmark regions')
            write_bed(output/'official-chr20-benchmark.bed', regions)
            truth = inputs_dir/'HG002_GRCh38_v5.0q_stvar.vcf.gz'; truth_counts = Counter()
            with pysam.VariantFile(str(truth)) as reader, pysam.VariantFile(str(output/'truth-chr20.full-context.vcf.gz'), 'wz', header=reader.header) as writer:
                if list(reader.header.samples) != [SAMPLE]: raise ValueError('benchmark sample mismatch')
                for record in reader.fetch(CHROM):
                    if seq[record.start:record.start+len(record.ref)] != record.ref.upper():
                        raise ValueError('benchmark REF differs from verified reference')
                    writer.write(record); truth_counts['all_chr20_records'] += 1
                    truth_counts['records_with_star_alt'] += int('*' in (record.alts or ()))
            pysam.tabix_index(str(output/'truth-chr20.full-context.vcf.gz'), preset='vcf')
        save(output/'reference-and-truth-check.json', {**ref_check, 'truth': dict(truth_counts),
              'benchmark_regions': len(regions), 'benchmark_bases': sum(b-a for _,a,b in regions),
              'selected_native_template_rows': selected, 'native_query_intervals': expected_queries,
              'benchmark_policy': 'Use the full official chr20 BED; do not restrict the primary benchmark to caller query regions.'})
        with pysam.AlignmentFile(str(alignment)) as reader:
            header = reader.header.to_dict(); groups = header.get('RG', [])
            if list(dictionary.items()) != list(zip(reader.references, reader.lengths)):
                raise ValueError('public alignment/reference dictionary differs')
            if not reader.has_index() or header.get('HD', {}).get('SO') != 'coordinate': raise ValueError('indexed sorted BAM required')
            if len(groups) != 1 or groups[0].get('SM') != SAMPLE or groups[0].get('ID') != 'HG002.pcr-free':
                raise ValueError('unexpected public sample/read group')
            if not any(pg.get('ID') == 'MarkDuplicates' for pg in header.get('PG', [])):
                raise ValueError('publisher duplicate-marking provenance absent')
            if any(row.mapped and row.contig != CHROM for row in reader.get_index_statistics()):
                raise ValueError('unexpected mapped chromosome in public subset')
            state['stage'] = 'full public read census'; save(output/'status.json', state)
            tick = time.monotonic(); census = alignment_census(reader, {groups[0]['ID']})
            census['elapsed_seconds'] = time.monotonic()-tick; save(output/'alignment-census.json', census)
        if census['counts'].get('mapped_unknown_MQ_255') or not census['counts'].get('mapped_records'):
            raise ValueError('public MQ policy or empty alignment failed')
        command('native public structural discovery', [tool_paths['delly'], 'sr', '-h', '2', *NATIVE_OPTIONS,
                  '-x', native_mask, '-g', reference, '-o', output/'calls.bcf', alignment])
        reviewed = review_calls(output/'calls.bcf', output/'calls-review.jsonl', SAMPLE, dictionary, [CHROM], intervals, whole)
        represented = Counter(); omitted = Counter()
        with pysam.VariantFile(str(output/'calls.bcf')) as reader, pysam.VariantFile(str(output/'del-ins-dup.vcf.gz'), 'wz', header=reader.header) as writer:
            for record in reader:
                kind = record.info['SVTYPE']
                if kind in ('DEL', 'INS', 'DUP'): writer.write(record); represented[kind] += 1
                else: omitted[kind] += 1
        pysam.tabix_index(str(output/'del-ins-dup.vcf.gz'), preset='vcf')
        result = {'scope': state['scope'], 'completed_at': stamp(), 'elapsed_seconds': time.monotonic()-started,
                  'inputs': inputs, 'truth_reference': truth_identity, 'tools': tools,
                  'code': {name: fingerprint(Path(__file__).with_name(name)) for name in
                           ('run_public_sv_pilot.py','prepare_public_sv_benchmark.py','run_structural_variants.py','build_structural_exclusions.py')},
                  'reference_check': ref_check, 'census': census, 'native_calls': reviewed,
                  'comparison_input_types': dict(represented), 'unbenchmarked_native_types': dict(omitted),
                  'benchmark_status': 'native calls prepared; comparison and discordance curation are separate',
                  'limitations': ['Prealigned single chromosome, different mapper and duplicate marker versions from patient pipeline.',
                                  'Interchromosomal and balanced rearrangements are outside this pilot truth comparison.',
                                  'No sensitivity estimate or causal/clinical conclusion follows from discovery alone.']}
        save(output/'summary.json', result); state.update(status='complete', completed_at=result['completed_at']); save(output/'status.json', state)
        return result
    except BaseException as error:
        state.update(status='failed', error_type=type(error).__name__, error=str(error), failed_at=stamp())
        save(output/'status.json', state); raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', required=True, type=Path); parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args(); result = pilot(args.root, args.output)
    print(json.dumps({'scope': result['scope'], 'elapsed_seconds': result['elapsed_seconds'],
                      'native_calls': result['native_calls'], 'benchmark_status': result['benchmark_status']}))
