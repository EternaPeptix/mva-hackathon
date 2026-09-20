#!/usr/bin/env python3
"""Validate the two prior candidates using retained lane slices; publish only aggregate evidence."""
import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import statistics
import subprocess
import time

import pysam


SITES = [('15', 40209701, 'T', 'G'), ('15', 40220612, 'T', 'G')]
REGION = '15:40160000-40260000'


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as source:
        for block in iter(lambda: source.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def save(path, value):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    os.replace(temporary, path)


def pileup_evidence(bam, reference, sites=SITES):
    results = []
    with pysam.AlignmentFile(str(bam), 'rb') as alignment, pysam.FastaFile(str(reference)) as fasta:
        for contig, pos, ref, alt in sites:
            if fasta.fetch(contig, pos - 1, pos).upper() != ref:
                raise ValueError('candidate REF does not match the reference')
            observations = []
            for column in alignment.pileup(contig, pos - 1, pos, truncate=True, stepper='samtools',
                    fastafile=fasta, compute_baq=True, min_base_quality=20, min_mapping_quality=20,
                    flag_filter=0xF04, ignore_overlaps=True, ignore_orphans=True, max_depth=100000):
                if column.reference_pos != pos - 1:
                    continue
                for item in column.pileups:
                    if item.is_del or item.is_refskip or item.query_position is None:
                        continue
                    read = item.alignment
                    q = item.query_position
                    observations.append({'template': read.query_name, 'base': read.query_sequence[q],
                        'base_quality': read.query_qualities[q], 'mapping_quality': read.mapping_quality,
                        'reverse': read.is_reverse, 'edge_distance': min(q, read.query_length - 1 - q),
                        'read_start': read.reference_start})
            # Count one observation per fragment even if mates overlap the site.
            fragments = {}
            for item in observations:
                if item['template'] not in fragments or item['base_quality'] > fragments[item['template']]['base_quality']:
                    fragments[item['template']] = item
            values = list(fragments.values())
            counts = Counter(item['base'] for item in values)
            alt_values = [item for item in values if item['base'] == alt]
            informative = counts[ref] + counts[alt]
            results.append({'chromosome': contig, 'position': pos, 'ref': ref, 'alt': alt,
                'filtered_fragment_depth': len(values), 'ref_fragments': counts[ref], 'alt_fragments': counts[alt],
                'other_fragments': len(values) - informative,
                'alt_fraction': counts[alt] / informative if informative else None,
                'alt_forward': sum(not item['reverse'] for item in alt_values),
                'alt_reverse': sum(item['reverse'] for item in alt_values),
                'alt_distinct_starts': len({item['read_start'] for item in alt_values}),
                'alt_median_base_quality': statistics.median(item['base_quality'] for item in alt_values) if alt_values else None,
                'alt_median_mapping_quality': statistics.median(item['mapping_quality'] for item in alt_values) if alt_values else None,
                'alt_median_distance_from_read_end': statistics.median(item['edge_distance'] for item in alt_values) if alt_values else None})
    return results


def candidate_calls(path, sites=SITES):
    calls = []
    with pysam.VariantFile(str(path)) as vcf:
        for contig, pos, ref, alt in sites:
            matches = []
            for record in vcf.fetch(contig, pos - 1, pos):
                if record.pos != pos or record.ref != ref or alt not in (record.alts or ()):
                    continue
                sample = record.samples[0]
                matches.append({'genotype': list(sample.get('GT') or ()), 'phased': bool(sample.phased),
                    'depth': sample.get('DP'), 'allele_depths': list(sample.get('AD') or ()),
                    'genotype_quality': sample.get('GQ'), 'quality': record.qual, 'filters': list(record.filter.keys()),
                    'physical_phase_present': sample.get('PGT') not in (None, '.', ''),
                    'phase_set_present': sample.get('PID') not in (None, '.', '') or sample.get('PS') not in (None, '.', '')})
            calls.append({'chromosome': contig, 'position': pos, 'ref': ref, 'alt': alt, 'calls': matches})
    return calls


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root, output = args.root.resolve(), args.output.resolve()
    if output.exists():
        parser.error('existing validation preserved; choose a new output directory')
    output.mkdir(parents=True)
    os.nice(10)
    env = dict(os.environ, PATH=str(root/'tools/bio/bin') + os.pathsep + os.environ.get('PATH', ''))
    source_dir = root/'work/read-validation-20260919'
    manifest = json.loads((source_dir/'slice-provenance.json').read_text())
    reference = root/'ref/vendor/GRCh38_nochr.fa'
    state = {'status': 'running', 'started_at': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
             'scope': manifest['scope'], 'script_sha256': sha256(Path(__file__)), 'commands': []}
    status_path = output/'status.json'
    save(status_path, state)
    started = time.monotonic()
    try:
        def execute(command, label):
            state['current_stage'] = label
            state['commands'].append(command)
            save(status_path, state)
            with (output/(label+'.log')).open('w') as log:
                subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)

        slices = []
        with pysam.FastaFile(str(reference)) as fasta:
            lengths = dict(zip(fasta.references, fasta.lengths))
        for item in manifest['files']:
            path = source_dir/item['slice_name']
            if sha256(path) != item['slice_sha256']:
                raise ValueError('transferred slice checksum mismatch')
            with pysam.AlignmentFile(str(path), 'rb') as bam:
                if dict(zip(bam.references, bam.lengths)) != lengths:
                    raise ValueError('slice reference dictionary differs from FASTA')
            slices.append(path)
        execute(['samtools','merge','-@','1','-o',str(output/'merged.bam')] + [str(p) for p in slices], 'merge')
        execute(['samtools','collate','-@','1','-o',str(output/'collated.bam'),str(output/'merged.bam')], 'collate')
        execute(['samtools','fixmate','-@','1','-m',str(output/'collated.bam'),str(output/'fixmate.bam')], 'fixmate')
        execute(['samtools','sort','-@','1','-m','128M','-o',str(output/'sorted.bam'),str(output/'fixmate.bam')], 'sort')
        bam = output/'two-lane-markdup.bam'
        execute(['samtools','markdup','-@','1','-f',str(output/'duplicate-marking.txt'),str(output/'sorted.bam'),str(bam)], 'markdup')
        execute(['samtools','quickcheck','-v',str(bam)], 'bam-validation')
        execute(['samtools','index',str(bam)], 'bam-index')
        dictionary = reference.with_suffix('.dict')
        if not dictionary.exists():
            execute(['samtools','dict','-o',str(dictionary)+'.tmp',str(reference)], 'reference-dictionary')
            os.replace(str(dictionary)+'.tmp', dictionary)
        gatk_vcf = output/'gatk.vcf.gz'
        execute([str(root/'tools/bio/lib/jvm/bin/java'),'-XX:ActiveProcessorCount=4','-Xmx8g','-jar',
                 str(root/'gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar'),'HaplotypeCaller',
                 '-R',str(reference),'-I',str(bam),'-O',str(gatk_vcf),'-L',REGION,
                 '--pair-hmm-implementation','LOGLESS_CACHING'], 'gatk-haplotypecaller')
        # Preserve the pileup and raw calls, avoiding a pipeline that could obscure a failure.
        execute(['bcftools','mpileup','-Ob','-o',str(output/'pileup.bcf'),'-f',str(reference),'-r',REGION,
                 '-a','AD,DP,SP','-q','20','-Q','20','-d','100000','--ff','0xF04',str(bam)], 'bcftools-pileup')
        execute(['bcftools','call','-mv','-Ob','-o',str(output/'calls.bcf'),str(output/'pileup.bcf')], 'bcftools-call')
        bcf_vcf = output/'bcftools.vcf.gz'
        execute(['bcftools','norm','-f',str(reference),'-m','-any','-Oz','-o',str(bcf_vcf),str(output/'calls.bcf')], 'bcftools-normalize')
        execute(['bcftools','index','-t',str(bcf_vcf)], 'bcftools-index')
        result = {'scope': manifest['scope'], 'lanes': [item['lane'] for item in manifest['files']],
            'candidate_region': REGION, 'reference_sha256': sha256(reference),
            'slice_provenance': manifest, 'evidence_filters': {'minimum_base_quality':20,'minimum_mapping_quality':20,
                'duplicates_excluded':True,'secondary_supplementary_qcfail_unmapped_excluded':True,
                'overlapping_mates_counted_once':True,'baq_enabled':True},
            'pileup': pileup_evidence(bam, reference), 'gatk': candidate_calls(gatk_vcf),
            'bcftools': candidate_calls(bcf_vcf),
            'limitations': ['Only two of four lanes; retained historical alignments, not the new four-lane run.',
                'Variant calls are unfiltered discovery calls; clinical classification is a separate review.',
                'Read support does not establish phase, functional effect, or pathogenicity of the VUS.']}
        result['completed_at'] = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
        result['elapsed_seconds'] = round(time.monotonic()-started,2)
        save(output/'candidate-evidence.json', result)
        state.update(status='complete', completed_at=result['completed_at'], elapsed_seconds=result['elapsed_seconds'],
                     current_stage='complete', evidence_sha256=sha256(output/'candidate-evidence.json'))
        save(status_path, state)
        print('Two-lane candidate read validation completed; aggregate evidence saved.')
    except BaseException as error:
        state.update(status='failed', error_type=type(error).__name__)
        save(status_path, state)
        raise


if __name__ == '__main__':
    main()
