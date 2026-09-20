#!/usr/bin/env python3
"""Complement one-read SNP phase with conservative, reciprocal proper-mate review."""
import argparse
from collections import Counter, defaultdict
import datetime
import itertools
import json
from pathlib import Path

import pysam

from screen_secondary_findings import fingerprint


def merge_mate_calls(records, site_count, maximum_template_length=2000):
    if len(records) != 2 or {r['mate_number'] for r in records} != {1, 2}:
        return None, 'missing_or_ambiguous_mate'
    first, second = records
    if (first['mate_start'] != second['start'] or second['mate_start'] != first['start']
            or first['template_length'] != -second['template_length']
            or not 0 < abs(first['template_length']) <= maximum_template_length):
        return None, 'inconsistent_mate_geometry'
    calls = [{r['calls'][i] for r in records if i in r['calls']} for i in range(site_count)]
    if any(not choices for choices in calls):
        return None, 'incomplete_target_coverage'
    if any(len(choices) > 1 for choices in calls):
        return 'mate_conflict', None
    labels = [next(iter(choices)) for choices in calls]
    return ('other' if '?' in labels else '|'.join(labels)), None


def review(bam, reference, sites):
    contig = sites[0]['chrom']
    if not 2 <= len(sites) <= 8 or any(s['chrom'] != contig or len(s['ref']) != 1 or len(s['alt']) != 1 for s in sites):
        raise ValueError('two to eight SNPs on the same contig required')
    positions = {s['pos'] - 1 for s in sites}
    with pysam.FastaFile(str(reference)) as fasta:
        if any(fasta.fetch(contig, s['pos'] - 1, s['pos']).upper() != s['ref'] for s in sites):
            raise ValueError('reference mismatch')
    fragments, read_exclusions = defaultdict(list), Counter()
    with pysam.AlignmentFile(str(bam), 'rb') as alignment:
        for read in alignment.fetch(contig, max(0, min(positions) - 2000), max(positions) + 2001):
            if read.flag & 0xF04 or read.mapping_quality < 20:
                read_exclusions['flags_or_mapping_quality'] += 1
                continue
            if (not read.is_paired or not read.is_proper_pair or read.mate_is_unmapped
                    or read.next_reference_id != read.reference_id or read.query_qualities is None
                    or any(op == 3 for op, length in read.cigartuples or [])):
                read_exclusions['pair_or_quality_or_reference_skip'] += 1
                continue
            mapping = {r: q for q, r in read.get_aligned_pairs() if r in positions}
            calls = {}
            for i, site in enumerate(sites):
                q = mapping.get(site['pos'] - 1)
                if q is not None and read.query_qualities[q] >= 20:
                    base = read.query_sequence[q].upper()
                    calls[i] = '0' if base == site['ref'] else '1' if base == site['alt'] else '?'
            key = (read.get_tag('RG') if read.has_tag('RG') else '', read.query_name)
            fragments[key].append({'mate_number': 1 if read.is_read1 and not read.is_read2 else 2 if read.is_read2 and not read.is_read1 else 0,
                'start': read.reference_start, 'mate_start': read.next_reference_start,
                'template_length': read.template_length, 'calls': calls})
    counts, exclusions = Counter(), Counter()
    for records in fragments.values():
        label, reason = merge_mate_calls(records, len(sites))
        if reason:
            exclusions[reason] += 1
        else:
            counts[label] += 1
    labels = ['|'.join(row) for row in itertools.product('01', repeat=len(sites))] + ['other', 'mate_conflict']
    return {'sites': sites, 'counts': {label: counts[label] for label in labels},
        'excluded_read_counts': dict(read_exclusions), 'excluded_template_counts': dict(exclusions),
        'filters': {'MQ': 20, 'BQ': 20, 'excluded_flags': '0xF04', 'proper_pair_required': True,
                    'reciprocal_mate_starts_and_TLEN_required': True, 'maximum_absolute_TLEN': 2000,
                    'both_mates_observed': True, 'fragment_key': 'RG and query name', 'BAQ_applied': False},
        'limitation': 'Aligned mate evidence does not validate locus assignment; missing/low-quality target bases do not become reference calls. Counts overlap the one-read audit and must not be added to it.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bam', type=Path, required=True)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('existing output preserved')
    inputs = {name: fingerprint(path) for name, path in {
        'bam': args.bam, 'bam_index': Path(str(args.bam) + '.bai'), 'reference': args.reference,
        'config': args.config, 'script': Path(__file__),
        'hash_helper': Path(__file__).with_name('screen_secondary_findings.py')}.items()}
    config = json.loads(args.config.read_text())
    if inputs['reference']['sha256'] != config['reference_sha256']:
        raise ValueError('reference checksum mismatch')
    by_key = {row['key']: row for row in config['targets']}
    sites = sorted({key for pair in config['snv_phase_pairs'] for key in pair}, key=lambda k: by_key[k]['pos'])
    result = {'created_at': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'inputs': inputs,
        'scope': 'Post hoc complementary paired-read review on the same two historical lanes; not an independent biological sample.',
        'pairs': [review(args.bam, args.reference, [by_key[k] for k in pair]) for pair in config['snv_phase_pairs']],
        'four_site_cluster': review(args.bam, args.reference, [by_key[k] for k in sites])}
    for row in inputs.values():
        stat = Path(row['path']).stat()
        if (stat.st_size, stat.st_mtime_ns) != (row['bytes'], row['mtime_ns']):
            raise ValueError('input changed during review')
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps({'four_site_cluster': result['four_site_cluster']['counts'],
                      'pairs': [{'keys': [s['key'] for s in row['sites']], 'counts': row['counts']} for row in result['pairs']]}))


if __name__ == '__main__':
    main()
