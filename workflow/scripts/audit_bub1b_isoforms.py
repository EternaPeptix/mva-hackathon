#!/usr/bin/env python3
"""Bounded public-annotation audit; transcript existence is not functional rescue.

All coordinates are one-based inclusive. This fixed BUB1B audit joins GFF
features by parent ID, independently of row order, and checks its MANE exons
against the earlier reference/CDS audit. It does not classify patient alleles.
"""
import argparse
from collections import Counter
import gzip
import hashlib
import json
import os
from pathlib import Path
import sys
import time
from urllib.parse import unquote


GFF_SHA256 = '0857407b4a8be8792089a69fd961737256fb7129b6b799f8a8920ebea4fb134f'
GENE = 'gene:ENSG00000156970'
MANE = 'ENST00000287598'
TARGETS = (40209701, 40220612)


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def features(lines, chrom='15', start=40150000, end=40250000):
    rows = []
    for line in lines:
        if line.startswith('##FASTA'):
            break
        if line.startswith('#') or not line.strip():
            continue
        fields = line.rstrip('\n').split('\t')
        if len(fields) != 9:
            raise ValueError('invalid GFF row')
        if fields[0] != chrom or int(fields[3]) > end or int(fields[4]) < start:
            continue
        attributes = {}
        for item in fields[8].split(';'):
            if item:
                key, value = item.split('=', 1)
                if key in attributes:
                    raise ValueError('duplicate GFF attribute')
                attributes[key] = unquote(value)
        rows.append({'chrom': fields[0], 'type': fields[2], 'start': int(fields[3]),
                     'end': int(fields[4]), 'strand': fields[6], 'attrs': attributes})
    return rows


def summarize(rows, gene=GENE, targets=TARGETS):
    genes = [r for r in rows if r['attrs'].get('ID') == gene]
    if len(genes) != 1 or genes[0]['strand'] != '+':
        raise ValueError('expected one positive-strand BUB1B gene')
    children = [r for r in rows if r['attrs'].get('Parent') == gene]
    transcripts = []
    for tr in children:
        a = tr['attrs']
        tid = a['transcript_id']
        if a['ID'] != 'transcript:' + tid:
            raise ValueError('transcript ID mismatch')
        parts = [r for r in rows if a['ID'] in r['attrs'].get('Parent', '').split(',')]
        exons = sorted((r['start'], r['end']) for r in parts if r['type'] == 'exon')
        cds = sorted((r['start'], r['end']) for r in parts if r['type'] == 'CDS')
        if not exons or any(s2 <= e1 for (_, e1), (s2, _) in zip(exons, exons[1:])):
            raise ValueError('missing/overlapping transcript exons')
        if any(r['strand'] != '+' or r['chrom'] != tr['chrom'] for r in parts):
            raise ValueError('inconsistent feature strand or chromosome')
        if any(not any(s <= c <= d <= e for s, e in exons) for c, d in cds):
            raise ValueError('CDS outside transcript exon')
        if any(s2 <= e1 for (_, e1), (s2, _) in zip(cds, cds[1:])):
            raise ValueError('overlapping CDS')
        loci = {}
        for pos in targets:
            ranks = [i + 1 for i, (s, e) in enumerate(exons) if s <= pos <= e]
            offset = 0
            coding_position = None
            for s, e in cds:
                if s <= pos <= e:
                    coding_position = offset + pos - s + 1
                offset += e - s + 1
            loci[str(pos)] = {'exon_rank': ranks[0] if ranks else None,
                             'in_exon': bool(ranks), 'in_CDS': coding_position is not None,
                             'offset_in_annotated_CDS': coding_position}
        transcripts.append({'id': tid, 'version': a.get('version'), 'name': a.get('Name'),
                            'biotype': a['biotype'], 'tags': a.get('tag', '').split(','),
                            'support_level': a.get('transcript_support_level'),
                            'exons': exons, 'CDS_intervals': cds,
                            'annotated_CDS_bases': sum(e-s+1 for s, e in cds), 'loci': loci})
    if len({t['id'] for t in transcripts}) != len(transcripts):
        raise ValueError('duplicate transcript ID')
    return {'gene': genes[0], 'transcripts': sorted(transcripts, key=lambda x: x['id']),
            'biotype_counts': dict(Counter(t['biotype'] for t in transcripts))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gff', type=Path, required=True)
    parser.add_argument('--transcript-audit', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    os.nice(10)
    started = time.monotonic()
    if args.output.exists():
        raise FileExistsError(args.output)
    if sha256(args.gff) != GFF_SHA256:
        raise ValueError('GFF checksum mismatch')
    previous = json.loads(args.transcript_audit.read_text())
    if previous['transcript'] != MANE + '.11' or not previous['full_cds_matches_spliced_reference']:
        raise ValueError('unexpected reference/CDS audit')
    with gzip.open(args.gff, 'rt') as stream:
        result = summarize(features(stream))
    mane = [t for t in result['transcripts'] if 'MANE_Select' in t['tags']]
    if len(mane) != 1 or mane[0]['id'] != MANE or mane[0]['version'] != '11':
        raise ValueError('MANE identity mismatch')
    m = mane[0]
    if m['exons'] != [(e['start'], e['end']) for e in previous['exons']]:
        raise ValueError('MANE exons disagree with independent exon request')
    if m['annotated_CDS_bases'] != previous['coding_length']:
        raise ValueError('MANE coding length mismatch')
    if [m['loci'][str(p)]['offset_in_annotated_CDS'] for p in TARGETS] != [2210, 3006]:
        raise ValueError('MANE target coordinate mismatch')
    coding = [t for t in result['transcripts'] if t['biotype'] == 'protein_coding']
    result['target_summary'] = {str(p): {
        'exact_protein_coding_total': len(coding),
        'exact_protein_coding_target_in_CDS': sum(t['loci'][str(p)]['in_CDS'] for t in coding),
        'omitting_CDS_target': [t['id'] for t in coding if not t['loci'][str(p)]['in_CDS']]}
        for p in TARGETS}
    result.update({'created_at': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
                   'command': [sys.executable] + sys.argv,
                   'inputs': {'gff': {'path': str(args.gff), 'sha256': GFF_SHA256},
                              'transcript_audit': {'path': str(args.transcript_audit),
                                                   'sha256': sha256(args.transcript_audit)},
                              'script_sha256': sha256(__file__)},
                   'MANE_exact_exon_and_CDS_agreement': True,
                   'elapsed_seconds': time.monotonic() - started,
                   'scope': 'Public annotated isoforms only. No expression weighting, tissue relevance, patient RNA, exon-skipping measurement, functional rescue, NMD measurement or variant classification is established.'})
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'transcripts': len(result['transcripts']), 'biotypes': result['biotype_counts'],
                      'targets': result['target_summary'], 'elapsed_seconds': result['elapsed_seconds']}))


if __name__ == '__main__':
    main()
