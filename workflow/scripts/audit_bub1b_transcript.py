#!/usr/bin/env python3
"""Reconcile public BUB1B exons/CDS with the retained reference; no patient RNA inference."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import pysam


def digest(path):
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-directory', required=True, type=Path)
    p.add_argument('--reference', required=True, type=Path)
    p.add_argument('--output', required=True, type=Path)
    args = p.parse_args()
    root = args.source_directory
    transcript_id = 'ENST00000287598'
    exon_data = json.loads((root / 'exons.json').read_text())
    cds_data = json.loads((root / 'cds.json').read_text())
    exons = sorted((e for e in exon_data if e['Parent'] == transcript_id), key=lambda e: e['rank'])
    if len(exons) < 2 or [e['rank'] for e in exons] != list(range(1, len(exons) + 1)):
        raise ValueError('missing or duplicate transcript exons')
    if any(e['strand'] != 1 or e['seq_region_name'] != '15' or
           e['assembly_name'] != 'GRCh38' for e in exons):
        raise ValueError('this focused audit requires positive-strand GRCh38 chromosome 15')
    if any(a['end'] >= b['start'] for a, b in zip(exons, exons[1:])):
        raise ValueError('overlapping or unsorted exons')
    if cds_data['id'] != transcript_id or cds_data['molecule'] != 'dna':
        raise ValueError('unexpected coding-sequence identity')
    cds = cds_data['seq'].upper()
    with pysam.FastaFile(str(args.reference)) as fa:
        contigs = set(fa.references) & {'15', 'chr15'}
        if len(contigs) != 1:
            raise ValueError('ambiguous chromosome alias')
        contig = next(iter(contigs))
        spliced = ''.join(fa.fetch(contig, e['start'] - 1, e['end']).upper() for e in exons)
        if spliced.count(cds) != 1:
            raise ValueError('public CDS does not match exactly one spliced reference sequence')
        cds_offset = spliced.index(cds)
        final_junction = sum(e['end'] - e['start'] + 1 for e in exons[:-1])
        findings = []
        for pos, expected_c, protein in [(40209701, 2210, 'p.Leu737Ter'),
                                         (40220612, 3006, 'p.Asn1002Lys')]:
            offset, match = 0, None
            for exon in exons:
                if exon['start'] <= pos <= exon['end']:
                    match = exon
                    break
                offset += exon['end'] - exon['start'] + 1
            if match is None:
                raise ValueError('candidate not in retained transcript')
            transcript_pos = offset + pos - match['start'] + 1
            coding_pos = transcript_pos - cds_offset
            if coding_pos != expected_c or fa.fetch(contig, pos - 1, pos).upper() != 'T':
                raise ValueError('candidate coding coordinate or reference mismatch')
            codon_start = ((coding_pos - 1) // 3) * 3
            reference_codon = cds[codon_start:codon_start + 3]
            alt_codon = list(reference_codon)
            alt_codon[(coding_pos - 1) % 3] = 'G'
            alt_codon = ''.join(alt_codon)
            expected_codons = ('TTA', 'TGA') if expected_c == 2210 else ('AAT', 'AAG')
            if (reference_codon, alt_codon) != expected_codons:
                raise ValueError('unexpected codon translation context')
            finding = {'chrom': 'chr15', 'position': pos, 'ref': 'T', 'alt': 'G',
                       'coding_coordinate': coding_pos, 'transcript_coordinate': transcript_pos,
                       'exon_rank': match['rank'], 'exon_count': len(exons),
                       'reference_codon': reference_codon, 'alternate_codon': alt_codon,
                       'protein_annotation': protein}
            if expected_c == 2210:
                ptc_end_transcript = cds_offset + codon_start + 3
                distance = final_junction - ptc_end_transcript
                finding.update(stop_codon_end_transcript_coordinate=ptc_end_transcript,
                    final_junction_after_transcript_coordinate=final_junction,
                    spliced_bases_from_stop_end_to_final_junction=distance,
                    positional_nmd_rule_predicted=match['rank'] < len(exons) and distance > 50)
            findings.append(finding)
    result = {'created_at': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
              'scope': 'Reference/transcript and positional NMD prediction only. No patient RNA, '
                       'protein, exon-skipping, rescue-isoform, expression or functional result is measured.',
              'transcript': transcript_id + '.' + str(cds_data['version']),
              'ensembl_releases': json.loads((root / 'info.json').read_text())['releases'],
              'spliced_length': len(spliced), 'coding_length': len(cds),
              'cds_start_transcript_coordinate': cds_offset + 1,
              'full_cds_matches_spliced_reference': True,
              'reference_sha256': digest(args.reference),
              'script_sha256': digest(Path(__file__)),
              'public_source_hashes': {name: digest(root / name) for name in
                  ('exons.json', 'cds.json', 'info.json', 'map-stop.json',
                   'fetch-manifest.json', 'additional-fetch-manifest.json')},
              'exons': [{k: e[k] for k in ('rank', 'start', 'end', 'exon_id', 'version')} for e in exons],
              'findings': findings}
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps({k: result[k] for k in ('transcript', 'spliced_length', 'coding_length', 'findings')}))


if __name__ == '__main__':
    main()
