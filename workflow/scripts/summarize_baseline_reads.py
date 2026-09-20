#!/usr/bin/env python3
"""Validate a completed bounded read review and export candidate-level evidence."""
import argparse
from collections import Counter
import datetime
import json
from pathlib import Path

import pysam

from baseline_transcript_sensitivity import checked_file
from screen_secondary_findings import fingerprint


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--vendor-vcf', type=Path, required=True)
    parser.add_argument('--prior-bub1b-review', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('existing receipt preserved; choose a new output')
    state = json.loads((args.run / 'run-status.json').read_text())
    if state['status'] != 'complete':
        raise ValueError('read review did not complete')
    for name in ('evidence.json', 'markdup.bam', 'gatk.normalized.vcf.gz', 'bcftools.normalized.vcf.gz'):
        checked_file(args.run / name, state['output_hashes'][name])
    evidence = json.loads((args.run / 'evidence.json').read_text())
    previous = json.loads(args.prior_bub1b_review.read_text())
    repeated = [row for row in evidence['SNP_BAQ_pileup'] if row['chromosome'] == '15']
    if repeated != previous['pileup']:
        raise ValueError('same-read BUB1B control pileup differs; inspect before exporting')
    vendor = checked_file(args.vendor_vcf, 'bec64c68d3ad10d05b8071d07cdf87d217013d66d3ba02e24b505035c9d42175')
    context = {'selection': 'Post hoc SERPINA1 alignment/caller review prompted by clipped alternate-supporting reads.',
               'interval_1based': '14:94379350-94379749', 'variant_record_counts': {}, 'alternate_read_geometries': []}
    for label, path in [('vendor', args.vendor_vcf), ('gatk', args.run / 'gatk.normalized.vcf.gz'),
                        ('bcftools', args.run / 'bcftools.normalized.vcf.gz')]:
        with pysam.VariantFile(str(path)) as reader:
            rows = list(reader.fetch('14', 94379349, 94379749))
            context['variant_record_counts'][label] = {'records': len(rows),
                'records_with_non_SNV_alleles': sum(any(len(a) != 1 for a in row.alleles) for row in rows),
                'called_variant_positions': sorted({row.pos for row in rows if any(x not in (0, None) for x in (row.samples[0].get('GT') or ()))})}
    sites = [('14', 94379526, 'C', 'T'), ('14', 94379563, 'A', 'T'),
             ('14', 94379593, 'T', 'A'), ('14', 94379607, 'C', 'T')]
    geometries = Counter()
    indel_mapping = []
    with pysam.AlignmentFile(str(args.run / 'markdup.bam'), 'rb') as alignment:
        for read in alignment.fetch('14', sites[0][1] - 1, sites[-1][1]):
            if read.flag & 0xF04 or read.mapping_quality < 20 or read.query_qualities is None:
                continue
            positions = {r: q for q, r in read.get_aligned_pairs() if r in {s[1] - 1 for s in sites}}
            bits = []
            for _, pos, ref, alt in sites:
                q = positions.get(pos - 1)
                bits.append('-' if q is None else '?' if read.query_qualities[q] < 20 else
                            '0' if read.query_sequence[q] == ref else '1' if read.query_sequence[q] == alt else '?')
            if '1' in bits:
                geometries[('|'.join(bits), read.reference_start, read.cigarstring, read.mapping_quality,
                            read.get_tag('NM') if read.has_tag('NM') else None, read.is_reverse)] += 1
        for row in evidence['config']['targets']:
            if len(row['ref']) == len(row['alt']) == 1:
                continue
            counts, flags = Counter(), Counter()
            for read in alignment.fetch(row['chrom'], row['pos'] - 1, row['pos']):
                if read.flag & 0xF04:
                    flags[hex(read.flag & 0xF04)] += 1
                else:
                    counts[read.mapping_quality] += 1
            indel_mapping.append({'key': row['key'], 'nonexcluded_reads_by_mapping_quality': dict(counts),
                                  'excluded_reads_by_flag_mask': dict(flags),
                                  'scope': 'Reads overlapping the anchor before a base-quality or complete-window test; not independent fragments.'})
    context['alternate_read_geometries'] = [{'states_at_four_sites': key[0], 'reference_start_0based': key[1],
        'cigar': key[2], 'mapping_quality': key[3], 'NM': key[4], 'reverse': key[5], 'read_count': value}
        for key, value in geometries.items()]
    context['geometry_scope'] = 'Per-read summaries, flags 0xF04 excluded, MQ>=20, target bases require BQ>=20. No read identifiers or sequences exported.'
    result = {'created_at': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'run': state, 'evidence': evidence, 'SERPINA1_context': context, 'indel_mapping_context': indel_mapping,
        'same_read_control': {'prior_receipt': fingerprint(args.prior_bub1b_review),
            'BUB1B_BAQ_pileup_reproduced_exactly': True,
            'scope': 'Same historical read data and helper, with different extraction windows; a consistency check, not independent biological replication.'},
        'vendor_VCF': vendor, 'exporter': fingerprint(Path(__file__))}
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps({'output': str(args.output), 'sha256': fingerprint(args.output)['sha256'],
                      'SERPINA1_context_record_counts': context['variant_record_counts'], 'indel_mapping_context': indel_mapping}))


if __name__ == '__main__':
    main()
