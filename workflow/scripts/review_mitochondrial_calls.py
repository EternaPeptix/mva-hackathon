#!/usr/bin/env python3
"""Reconcile mitochondrial allele identities without replacing original genotype/filter context."""
import argparse
from collections import Counter
import json
from pathlib import Path
import subprocess

import pysam

from prepare_gatk_compatible_cram import digest, save, stamp


def normalized_records(source, reference, bcftools, output_prefix):
    """Bounded mitochondrial extraction, contig aliasing, normalization and original-row links."""
    original, skipped = [], []
    with pysam.VariantFile(str(source)) as reader, pysam.FastaFile(str(reference)) as fasta:
        candidates = [c for c in reader.header.contigs if c in ('M', 'MT', 'chrM', 'chrMT')]
        if len(candidates) != 1 or fasta.references != ['chrM']:
            raise ValueError('require one mitochondrial source alias and a chrM reference')
        header = reader.header.copy()
        if 'MVA_RECORD' in header.info or 'MVA_ALT_INDEX' in header.info:
            raise ValueError('reserved normalization tags already exist')
        if 'chrM' not in header.contigs:
            header.contigs.add('chrM', length=fasta.get_reference_length('chrM'))
        header.info.add('MVA_RECORD', 1, 'Integer', 'Index of unchanged source-record context')
        header.info.add('MVA_ALT_INDEX', 'A', 'Integer', 'One-based ALT index in the original source record')
        extract = output_prefix.with_suffix('.input.vcf')
        with pysam.VariantFile(str(extract), 'w', header=header) as writer:
            for record in reader.fetch(candidates[0]):
                row = {'chrom': record.contig, 'pos': record.pos, 'id': record.id, 'ref': record.ref,
                       'alts': list(record.alts or []), 'filters': list(record.filter),
                       'info': dict(record.info),
                       'samples': {name: dict(sample) | {'phased': sample.phased} for name, sample in record.samples.items()}}
                original.append(row)
                if not record.alts or any(any(c not in 'ACGT' for c in allele) for allele in record.alleles):
                    skipped.append({'source_record_index': len(original)-1, 'reason': 'nonconcrete_or_symbolic_allele'})
                    continue
                if fasta.fetch('chrM', record.start, record.stop).upper() != record.ref:
                    raise ValueError(f'source reference mismatch: {source.name}:{record.pos}')
                copy = record.copy()
                copy.translate(header)
                copy.contig = 'chrM'
                copy.info['MVA_RECORD'] = len(original)-1
                copy.info['MVA_ALT_INDEX'] = tuple(range(1, len(copy.alts)+1))
                writer.write(copy)
    normalized = output_prefix.with_suffix('.normalized.vcf')
    command = list(map(str, [bcftools, 'norm', '-c', 'e', '-f', reference, '-m', '-any',
                            '--old-rec-tag', 'MVA_ORIG', '-Ov', '-o', normalized, extract]))
    with output_prefix.with_suffix('.normalization.log').open('w') as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
    rows = []
    with pysam.VariantFile(str(normalized)) as reader:
        for record in reader:
            if len(record.alts) != 1:
                raise ValueError('normalization did not produce an allele identity')
            record_index = record.info['MVA_RECORD']
            alt_index = record.info['MVA_ALT_INDEX'][0]
            if not 1 <= alt_index <= len(original[record_index]['alts']):
                raise ValueError('original ALT index lost')
            rows.append({'key': f'chrM:{record.pos}:{record.ref}:{record.alts[0]}',
                         'pos': record.pos, 'ref': record.ref, 'alt': record.alts[0],
                         'source_record_index': record_index, 'original_ALT_index': alt_index})
    return {'source_sha256': digest(source), 'original_records': original, 'normalized_identities': rows,
            'skipped_records': skipped, 'normalization_command': command}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('root', 'run', 'vendor-vcf', 'clinvar-vcf', 'output'):
        p.add_argument('--' + name, type=Path, required=True)
    args = p.parse_args()
    out = args.output.resolve()
    if out.exists():
        raise ValueError('existing mitochondrial review preserved')
    status = json.loads((args.run / 'status.json').read_text())
    if status['status'] != 'complete':
        raise ValueError('case calling is incomplete')
    for name, expected in status['output_sha256'].items():
        if digest(args.run / name) != expected:
            raise ValueError('completed case output changed: ' + name)
    if digest(args.vendor_vcf) != 'bec64c68d3ad10d05b8071d07cdf87d217013d66d3ba02e24b505035c9d42175':
        raise ValueError('unexpected vendor VCF')
    if digest(args.clinvar_vcf) != '47ea3626cd90154d82304afe3b8f717aab2557a1330f4ea2c809892386e08d49':
        raise ValueError('unexpected dated ClinVar source')
    out.mkdir(parents=True, exist_ok=False)
    reference = args.root / 'ref/mitochondrial-20260920/Homo_sapiens_assembly38.chrM.fasta'
    bcftools = args.root / 'tools/bio/bin/bcftools'
    sources = {}
    for label, source in [('case', args.run/'masked.assuming-no-contamination.vcf.gz'),
                          ('vendor', args.vendor_vcf), ('clinvar', args.clinvar_vcf)]:
        sources[label] = normalized_records(source, reference, bcftools, out/label)
        save(out/(label+'.records.json'), sources[label])
    indices = {}
    for label, source in sources.items():
        index = {}
        for row in source['normalized_identities']:
            index.setdefault(row['key'], []).append(row)
        indices[label] = index
    rows = []
    for key in sorted(set(indices['case']) | set(indices['vendor']), key=lambda k:(int(k.split(':')[1]), k)):
        row = {'key': key}
        for label in ('case', 'vendor', 'clinvar'):
            row[label] = [{'normalized_identity': x,
                           'original_record': sources[label]['original_records'][x['source_record_index']]}
                          for x in indices[label].get(key, [])]
        rows.append(row)
    result = {'completed_at': stamp(), 'scope': 'Initial mitochondrial discovery versus vendor identities and exact dated ClinVar alleles; no new clinical classification.',
        'case_status_sha256': digest(args.run/'status.json'),
        'source_and_tool_sha256': {str(q): digest(q) for q in [Path(__file__), bcftools, reference,
            args.vendor_vcf, args.clinvar_vcf, Path(__file__).with_name('prepare_gatk_compatible_cram.py')]},
        'counts': {label: {'original_records':len(data['original_records']),
            'normalized_identities':len(data['normalized_identities']), 'skipped_records':len(data['skipped_records'])}
            for label,data in sources.items()},
        'comparison': {'shared': len(set(indices['case']) & set(indices['vendor'])),
                       'case_only': len(set(indices['case'])-set(indices['vendor'])),
                       'vendor_only': len(set(indices['vendor'])-set(indices['case']))},
        'alleles': rows,
        'limits': ['Original multi-allele GT, AD, AF and both site/allele filter annotations remain authoritative; normalized identities are not genotypes.',
                   'Case FILTER status assumes zero contamination. No clinical-pass or negative result follows.',
                   'ClinVar record classifications, conflicts and containing-haplotype assertions remain separate raw annotations; exact matching is not causality.',
                   'Vendor and new callsets use the same child data and are not biological replication.',
                   'Nuclear-copy mapping, low-fraction evidence, tissue identity, contamination and clinical assessment remain open.']}
    save(out/'summary.json',result)
    receipt={'completed_at':stamp(),'status':'complete','source_sha256':digest(Path(__file__)),
             'output_sha256':{p.name:digest(p) for p in out.iterdir() if p.is_file()}}
    save(out/'status.json',receipt)
    print(json.dumps({'status':'complete','counts':result['counts'],'comparison':result['comparison']}))


if __name__ == '__main__':
    main()
