#!/usr/bin/env python3
"""Review explicit phase blocks for reported recessive contributors without changing ranks."""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import itertools
import json
from pathlib import Path
import re
import time

import pysam


def explicit_phase(sample, alt_index, contig):
    gt = sample.get('GT')
    if not gt or len(gt) != 2 or None in gt or gt.count(alt_index) != 1:
        return None
    pgt, pid = sample.get('PGT'), sample.get('PID')
    if isinstance(pgt, str) and re.fullmatch(r'\d+\|\d+', pgt) and pid not in (None, '.', ''):
        alleles = tuple(map(int, pgt.split('|')))
        if sorted(alleles) == sorted(gt):
            return (contig, 'PID', str(pid), alleles.index(alt_index))
    ps = sample.get('PS')
    if getattr(sample, 'phased', False) and ps not in (None, '.', ''):
        return (contig, 'PS', str(ps), gt.index(alt_index))
    return None


def pair_relation(first, second):
    if first is None or second is None or first[:3] != second[:3]:
        return 'unknown'
    return 'cis' if first[3] == second[3] else 'trans'


def audit(vcf_path, variants_path, sample_name):
    groups = defaultdict(list)
    with variants_path.open() as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            if row['MOI'] in ('AR', 'XR') and row['CONTRIBUTING_VARIANT'] == '1':
                groups[(row['#RANK'], row['GENE_SYMBOL'], row['MOI'])].append(row)
    result = []
    with pysam.VariantFile(str(vcf_path)) as vcf:
        if sample_name not in vcf.header.samples:
            raise ValueError('sample missing from VCF')
        for (rank, gene, moi), rows in groups.items():
            phases = []
            for row in rows:
                position = int(row['START'])
                matches = [r for r in vcf.fetch(row['CONTIG'], position-1, position)
                           if r.pos == position and r.ref == row['REF'] and row['ALT'] in r.alts]
                # Exomiser can normalize/decompose alleles. Unmatched representations stay unknown.
                phases.append(explicit_phase(matches[0].samples[sample_name], matches[0].alts.index(row['ALT'])+1,
                                             matches[0].chrom) if len(matches) == 1 else None)
            relations = Counter(pair_relation(a, b) for a, b in itertools.combinations(phases, 2))
            if len(phases) == 2:
                status = next(iter(relations))
            elif len(phases) < 2:
                status = 'single_contributor_not_a_pair'
            else:
                status = 'multiple_contributors_review_pairs'
            result.append({'rank': int(rank), 'gene': gene, 'inheritance_model': moi,
                           'contributing_variants': len(rows), 'phase_status': status,
                           'pair_relations': dict(relations)})
    return {'scope': 'Explicit phase-block review of reported recessive contributors; original ranks are unchanged.',
            'source_variant_table_sha256': hashlib.sha256(variants_path.read_bytes()).hexdigest(),
            'sample': sample_name, 'entries': sorted(result, key=lambda x:x['rank']),
            'phase_status_counts': dict(Counter(row['phase_status'] for row in result)),
            'completed_at': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
            'limitations': ['Unknown means no comparable valid explicit block; it does not mean trans.',
                            'Exact allele matching only; decomposed/normalized representations may remain unknown.',
                            'Vendor physical phase is not parental origin or independent read confirmation.',
                            'A single contributor requires separate genotype/model review.']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--vcf', type=Path, required=True)
    parser.add_argument('--variants', type=Path, required=True)
    parser.add_argument('--sample', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('existing phase review preserved')
    result = audit(args.vcf, args.variants, args.sample)
    result['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'phase_status_counts':result['phase_status_counts'],
                      'cis_entries':[r for r in result['entries'] if r['phase_status']=='cis']}))


if __name__ == '__main__':
    main()
