#!/usr/bin/env python3
"""Qualify a pinned mutserve 2.0.0-rc12 SNP subset without changing genotypes.

Native text output is required because its VCF writer can silently retain only
the first heteroplasmy at a position. Unsupported sites are recorded, never
split or coerced. This is an input contract, not clinical validation.
"""
from collections import defaultdict
import csv
import math
from pathlib import Path

import pysam

MUTSERVE_HOTSPOTS = frozenset([301,302,303,304,305,306,307,308,309,310,311,
                             315,316,3105,3106,3107,16182])


def qualify(vcf, table, reference, sample, output, mask_hotspots=False):
    output = Path(output)
    if output.exists():
        raise ValueError('existing qualified VCF preserved')
    rows = defaultdict(list)
    with Path(table).open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            if row['ID'] != sample:
                raise ValueError('unexpected caller sample')
            rows[row['Pos']].append(row)
    accepted, excluded, seen, expected = [], [], set(), {}
    with pysam.VariantFile(str(vcf)) as reader:
        if list(reader.header.samples) != [sample]:
            raise ValueError('unexpected VCF sample')
        header = reader.header.copy()
        for record in reader:
            key = str(record.pos)
            if record.contig != 'chrM' or key in seen:
                raise ValueError('unexpected contig or duplicate VCF position')
            seen.add(key)
            if not 1 <= record.pos <= len(reference):
                raise ValueError('position outside reference')
            if record.ref != reference[record.start:record.stop]:
                raise ValueError('VCF reference mismatch')
            source = rows.get(key, [])
            reason = None
            if not source:
                raise ValueError('VCF position absent from native table')
            if len(source) != 1:
                reason = 'multiple_native_rows_at_position'
            elif list(record.filter) != ['PASS'] or source[0]['Filter'] != 'PASS':
                reason = 'native_filter_not_PASS'
            elif len(record.alts or []) != 1 or any(len(a) != 1 or a not in 'ACGT' for a in record.alleles):
                reason = 'non_biallelic_concrete_SNP'
            elif source[0]['Type'] not in ['1','2']:
                reason = 'unsupported_native_type'
            elif mask_hotspots and record.pos in MUTSERVE_HOTSPOTS:
                reason = 'explicit_mutserve_source_hotspot_mask'
            if reason:
                excluded.append(dict(pos=record.pos, reason=reason, native_rows=source))
                continue
            row = source[0]
            if row['Ref'] != record.ref or row['Variant'] != record.alts[0]:
                raise ValueError('native text/VCF allele disagreement')
            values = record.samples[sample]
            gt = values.get('GT')
            if not gt or any(i is None or i not in [0,1] for i in gt):
                raise ValueError('unsupported genotype')
            dp = int(row['Coverage'])
            if dp <= 0 or values.get('DP') != dp:
                raise ValueError('native text/VCF depth disagreement')
            level = float(row['VariantLevel'])
            if not math.isfinite(level) or not 0 < level <= 1:
                raise ValueError('invalid native variant level')
            af = values.get('AF')
            af = list(af) if isinstance(af, tuple) else [af]
            af = [float(x) for x in af]
            if any(not math.isfinite(x) or not 0 <= x <= 1 for x in af) or abs(af[0]-level) > .0005:
                raise ValueError('native text/VCF frequency disagreement')
            exp = dict(type=int(row['Type']), coverage=dp, ref=record.ref,
                       variantBase=record.alts[0])
            if row['Type'] == '1':
                if gt != (1,):
                    raise ValueError('homoplasmic native genotype differs')
            else:
                major, minor = row['MajorBase'], row['MinorBase']
                major_level, minor_level = float(row['MajorLevel']), float(row['MinorLevel'])
                if {major,minor} != set(record.alleles):
                    raise ValueError('native allele pair differs')
                if not (0 < minor_level <= major_level < 1 and abs(major_level+minor_level-1) <= .0015):
                    excluded.append(dict(pos=record.pos, reason='noncomplementary_native_fractions', native_rows=source))
                    continue
                if tuple(record.alleles[i] for i in gt) != (major,minor):
                    raise ValueError('native genotype order differs from native table')
                exp.update(major=major,minor=minor,majorLevel=major_level,minorLevel=minor_level)
            expected[key] = exp
            accepted.append(record.copy())
        for pos, source in rows.items():
            if pos not in seen:
                if any(row['Filter']=='PASS' and row['Type'] in ['1','2'] for row in source):
                    raise ValueError('unaccounted PASS native row without VCF record')
                excluded.append(dict(pos=pos,reason='native_table_only_filtered_or_unsupported',native_rows=source))
    with pysam.VariantFile(str(output),'w',header=header) as dest:
        for record in accepted:
            dest.write(record)
    return dict(sample=sample,accepted_records=len(accepted),expected_import=expected,
                excluded_positions=excluded,native_table_positions=len(rows),native_vcf_records=len(seen),
                hotspot_mask=sorted(MUTSERVE_HOTSPOTS) if mask_hotspots else [],
                scope='Original VCF records retained without GT/AF recoding. Entire ambiguous/unsupported sites excluded with accounting.')


def verify_import(observed, qualification):
    samples = observed['samples']
    sample = qualification['sample']
    if set(samples) != {sample}:
        raise ValueError('parser sample accounting differs')
    imported = samples[sample]['variants']
    expected = qualification['expected_import']
    if set(imported) != set(expected):
        raise ValueError('parser position accounting differs')
    for pos, values in expected.items():
        for key, value in values.items():
            actual = imported[pos].get(key)
            if key in ['majorLevel','minorLevel']:
                if actual is None or abs(actual-value) > .0015:
                    raise ValueError('parser frequency differs at '+pos)
            elif actual != value:
                raise ValueError('parser allele/class/depth differs at '+pos+':'+key)
    if samples[sample]['amountVariants'] != len(expected):
        raise ValueError('parser variant counter differs')
    return dict(verified_records=len(expected),homoplasmic=sum(x['type']==1 for x in expected.values()),
                heteroplasmic=sum(x['type']==2 for x in expected.values()),
                scope='Native importer preserves the selected caller SNP identities, classes, depths and major/minor fractions.')
