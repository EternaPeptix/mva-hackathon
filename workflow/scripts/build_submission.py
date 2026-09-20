#!/usr/bin/env python3
"""Build a reference-checked Track 1 CSV from an explicit candidate review.

EPCR is a subjective causal probability supplied by the review, never an Exomiser
score or a self-score against an assumed answer. This validates format/provenance,
not pathogenicity, phase or the private challenge answer. Existing files are preserved.
"""
import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import re
import time

import pysam

FIELDS = ['proband_id', 'chrom_1', 'pos_1', 'ref_1', 'alt_1',
          'chrom_2', 'pos_2', 'ref_2', 'alt_2', 'epcr', 'finding_type', 'notes']
PROBAND_ID = 'PROBAND01'


def digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def required_text(value, label):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(label + ' must be explicit nonempty text')
    return value.strip()


def allele(value):
    if not isinstance(value, dict) or set(value) != {'chrom', 'pos', 'ref', 'alt'}:
        raise ValueError('each allele needs exactly chrom, pos, ref and alt')
    chrom = value['chrom']
    if not isinstance(chrom, str) or not re.fullmatch(r'chr(?:[1-9]|1[0-9]|2[0-2]|X|Y|M)', chrom):
        raise ValueError('use GRCh38 primary chromosome names with chr prefix')
    pos = value['pos']
    if type(pos) is not int or pos < 1:
        raise ValueError('positions must be positive one-based integers')
    ref, alt = value['ref'], value['alt']
    if any(not isinstance(a, str) or not re.fullmatch('[ACGT]+', a) for a in (ref, alt)) or ref == alt:
        raise ValueError('use distinct, explicit uppercase DNA REF and ALT alleles')
    return chrom, pos, ref, alt


def check_reference(variant, reference):
    chrom, pos, ref, _ = variant
    aliases = [chrom, chrom[3:]] + (['MT'] if chrom == 'chrM' else [])
    found = [name for name in aliases if name in reference.references]
    if len(found) != 1:
        raise ValueError('reference chromosome missing or ambiguous: ' + chrom)
    observed = reference.fetch(found[0], pos - 1, pos - 1 + len(ref)).upper()
    if observed != ref:
        raise ValueError('REF mismatch at %s:%s' % (chrom, pos))


def build_rows(review, reference, draft=False):
    if review.get('schema_version') != 1 or review.get('assembly') != 'GRCh38':
        raise ValueError('review must specify schema_version 1 and GRCh38')
    if review.get('proband_id') != PROBAND_ID:
        raise ValueError('unexpected proband identifier')
    if review.get('review_status') not in ('draft', 'reviewed'):
        raise ValueError('explicit draft or reviewed status required')
    if not draft and review['review_status'] != 'reviewed':
        raise ValueError('final export requires a completed candidate review')
    for key in ('team', 'approach'):
        if not isinstance(review.get(key), str) or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]*', review[key]):
            raise ValueError('team and approach must be filename-safe labels')
    if not draft:
        required_text(review.get('reviewed_by'), 'reviewed_by')
        required_text(review.get('reviewed_at'), 'reviewed_at')
    candidates = review.get('candidates')
    if not isinstance(candidates, list) or not 1 <= len(candidates) <= 10:
        raise ValueError('one to ten candidate rows required')
    rows, seen = [], set()
    for candidate in candidates:
        epcr = candidate.get('epcr')
        if type(epcr) not in (int, float) or not math.isfinite(epcr) or not 0 < epcr <= 1:
            raise ValueError('EPCR must be an explicit finite probability in (0,1]')
        required_text(candidate.get('epcr_rationale'), 'epcr_rationale')
        finding_type = candidate.get('finding_type')
        if finding_type not in ('primary', 'secondary'):
            raise ValueError('invalid finding_type')
        variants = candidate.get('variants')
        if not isinstance(variants, list) or len(variants) not in (1, 2):
            raise ValueError('each row requires one variant or a proposed pair')
        variants = [allele(v) for v in variants]
        if len(set(variants)) != len(variants):
            raise ValueError('a pair cannot repeat the same allele')
        identity = tuple(sorted(variants))
        if identity in seen:
            raise ValueError('duplicate candidate, including a reversed pair')
        seen.add(identity)
        for variant in variants:
            check_reference(variant, reference)
        notes = required_text(candidate.get('notes'), 'notes')
        if len(variants) == 2:
            if variants[0][0] != variants[1][0]:
                raise ValueError('a compound-het pair must be on the same chromosome')
            phase = candidate.get('phase')
            if phase not in ('unknown', 'trans'):
                raise ValueError('cis pairs are not compound-het candidates; declare unknown or evidenced trans')
            if phase == 'trans':
                required_text(candidate.get('phase_evidence'), 'phase_evidence')
            notes += (' Proposed compound-het pair; phase unknown.' if phase == 'unknown'
                      else ' In-trans phase recorded in the review evidence.')
        elif candidate.get('phase') not in (None, 'not_applicable'):
            raise ValueError('phase applies only to a pair')
        row = dict.fromkeys(FIELDS, '')
        row.update(proband_id=PROBAND_ID, epcr=epcr, finding_type=finding_type, notes=notes)
        for i, variant in enumerate(variants, 1):
            for key, value in zip(('chrom', 'pos', 'ref', 'alt'), variant):
                row['%s_%s' % (key, i)] = value
        rows.append(row)
    # Equal credences are legitimate. Keep their input order and never move an
    # incidental finding merely to manipulate the competition's rank metric.
    return sorted(rows, key=lambda row: -row['epcr'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--review', required=True, type=Path)
    parser.add_argument('--reference', required=True, type=Path)
    parser.add_argument('--out', required=True, type=Path)
    parser.add_argument('--draft', action='store_true')
    args = parser.parse_args()
    review = json.loads(args.review.read_text())
    with pysam.FastaFile(str(args.reference)) as reference:
        rows = build_rows(review, reference, draft=args.draft)
    observed_hash = digest(args.reference)
    if review.get('reference_sha256') != observed_hash:
        raise ValueError('review reference identity missing or mismatched')
    expected_name = review['team'] + '_' + review['approach'] + ('.draft.csv' if args.draft else '.csv')
    if args.out.name != expected_name:
        raise ValueError('required team/approach filename is ' + expected_name)
    metadata = Path(str(args.out) + '.provenance.json')
    if args.out.exists() or metadata.exists():
        raise FileExistsError('existing export preserved; select a new reviewed approach/version')
    stream = io.StringIO(newline='')
    writer = csv.DictWriter(stream, fieldnames=FIELDS)
    writer.writeheader()
    writer.writerows(rows)
    data = stream.getvalue()
    if list(csv.DictReader(io.StringIO(data)).fieldnames) != FIELDS:
        raise ValueError('CSV round-trip header mismatch')
    receipt = {'created_at': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
               'status': 'draft' if args.draft else 'reviewed_format_checked', 'rows': len(rows),
               'review_sha256': digest(args.review), 'reference_sha256': observed_hash,
               'script_sha256': digest(Path(__file__)),
               'csv_sha256': hashlib.sha256(data.encode()).hexdigest(),
               'scope': 'Format and reference checks only; no private truth or clinical validation.'}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('x', newline='') as output:
        output.write(data)
    with metadata.open('x') as output:
        output.write(json.dumps(receipt, indent=2) + '\n')
    print('Wrote %d reference-checked %s row(s); no clinical or leaderboard validation.' %
          (len(rows), receipt['status']))


if __name__ == '__main__':
    main()
