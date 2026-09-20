#!/usr/bin/env python3
"""Audit genotype and AF suitability for ROH without emitting patient variant rows."""
import argparse
from collections import Counter, defaultdict
import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import time

import pysam


def fingerprint(path):
    before = path.stat()
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(block)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError('input changed while hashing')
    return {'bytes': after.st_size, 'mtime_ns': after.st_mtime_ns, 'sha256': digest.hexdigest()}


def genotype_kind(gt):
    if not gt or any(x is None for x in gt):
        return 'missing_or_partial'
    if len(gt) != 2:
        return 'non_diploid'
    if gt == (0, 0):
        return 'homref'
    return 'homalt' if gt[0] == gt[1] else 'heterozygous'


def audit(vcf_path, sample_name):
    counts, gt_all, af_counts, an_counts = Counter(), Counter(), Counter(), Counter()
    by_chrom = defaultdict(Counter)
    with pysam.VariantFile(str(vcf_path)) as reader:
        if list(reader.header.samples) != [sample_name]:
            raise ValueError('expected exactly the named sample')
        headers = {'samples': list(reader.header.samples),
                   'INFO': {key: {'number': str(reader.header.info[key].number),
                                 'description': reader.header.info[key].description}
                            for key in ('AF', 'AC', 'AN', 'MLEAF') if key in reader.header.info},
                   'FORMAT': list(reader.header.formats)}
        for record in reader:
            counts['all_records'] += 1
            sample = record.samples[sample_name]
            gt = sample.get('GT')
            kind = genotype_kind(gt)
            gt_all[kind] += 1
            chrom = record.contig.removeprefix('chr')
            if chrom not in {str(x) for x in range(1, 23)}:
                counts['non_autosomal_records'] += 1
                continue
            by_chrom[chrom]['records'] += 1
            by_chrom[chrom]['genotype_' + kind] += 1
            if tuple(record.filter) != ('PASS',):
                continue
            counts['autosomal_PASS'] += 1
            if len(record.alleles) != 2 or any(a not in 'ACGT' or len(a) != 1 for a in record.alleles):
                continue
            if record.ref == record.alts[0]:
                raise ValueError('identical alleles at a biallelic SNP')
            counts['autosomal_PASS_biallelic_SNP'] += 1
            if kind in ('missing_or_partial', 'non_diploid') or any(x not in (0, 1) for x in gt):
                continue
            counts['autosomal_PASS_biallelic_SNP_called_diploid'] += 1
            af, ac, an = record.info.get('AF'), record.info.get('AC'), record.info.get('AN')
            an_counts[str(an)] += 1
            if af is None or len(af) != 1 or af[0] is None or not math.isfinite(af[0]):
                af_counts['AF_missing_or_invalid'] += 1
            else:
                value = af[0]
                af_counts['AF_present'] += 1
                af_counts['AF_equals_sample_alt_dosage_over_2'] += abs(value - sum(gt) / 2) < 1e-6
                af_counts['AF_zero_or_one'] += value in (0.0, 1.0)
                if an is not None and an > 0 and ac is not None and len(ac) == 1 and ac[0] is not None:
                    af_counts['AC_AN_present'] += 1
                    af_counts['AF_equals_AC_over_AN'] += abs(value - ac[0] / an) < 1e-6
            dp, gq, pl = sample.get('DP'), sample.get('GQ'), sample.get('PL')
            pl_valid = (pl is not None and len(pl) == 3
                        and all(type(x) is int and x >= 0 for x in pl))
            counts['PL_three_nonnegative_integers'] += pl_valid
            if pl_valid:
                called_index = sum(gt)
                counts['PL_called_genotype_is_a_minimum'] += pl[called_index] == min(pl)
            if dp is None or gq is None:
                counts['missing_DP_or_GQ'] += 1
            if dp is not None and gq is not None and 10 <= dp <= 100 and gq >= 20:
                counts['nominal_DP10_100_GQ20'] += 1
                by_chrom[chrom]['nominal_' + kind] += 1
                if pl_valid:
                    counts['nominal_with_valid_PL'] += 1
    return {'header': headers, 'counts': dict(counts), 'all_record_genotypes': dict(gt_all),
            'AF_audit_on_PASS_autosomal_called_biallelic_SNPs': dict(af_counts),
            'AN_distribution_on_same_sites': dict(an_counts),
            'chromosome_counts': dict(by_chrom),
            'nominal_filter_scope': 'Descriptive eligibility census, not a validated clinical threshold or selected population marker panel.',
            'scope': 'Input-suitability audit only. No ROH intervals, UPD/parental-origin inference, callability or patient-negative claim.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--vcf', type=Path, required=True)
    parser.add_argument('--expected-sha256', required=True)
    parser.add_argument('--sample', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('choose a new output path')
    os.nice(10)
    started = time.monotonic()
    identity = fingerprint(args.vcf)
    if identity['sha256'] != args.expected_sha256:
        raise ValueError('source VCF checksum mismatch')
    result = audit(args.vcf, args.sample)
    end = args.vcf.stat()
    if (end.st_size, end.st_mtime_ns) != (identity['bytes'], identity['mtime_ns']):
        raise ValueError('source changed during audit')
    result.update(completed_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  elapsed_seconds=time.monotonic() - started, source=identity,
                  script_sha256=fingerprint(Path(__file__))['sha256'], pysam=pysam.__version__)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: result[key] for key in ('completed_at','elapsed_seconds','counts',
          'all_record_genotypes','AF_audit_on_PASS_autosomal_called_biallelic_SNPs',
          'AN_distribution_on_same_sites','scope')}, indent=2))


if __name__ == '__main__':
    main()
