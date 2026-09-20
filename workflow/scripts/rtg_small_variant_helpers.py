#!/usr/bin/env python3
"""Shared provenance and native comparison helpers for the public control study."""
from collections import Counter
import json
import math
from pathlib import Path
import subprocess
import time

import pysam

from build_structural_exclusions import fingerprint, save, stamp


def execute(args, folder, label):
    args = list(map(str, args)); start = time.monotonic()
    with (folder / (label + '.stdout')).open('w') as stdout, (folder / (label + '.stderr')).open('w') as stderr:
        result = subprocess.run(args, stdout=stdout, stderr=stderr)
    receipt = {'command': args, 'exit_code': result.returncode,
               'elapsed_seconds': time.monotonic() - start,
               'stdout': fingerprint(folder / (label + '.stdout')),
               'stderr': fingerprint(folder / (label + '.stderr'))}
    save(folder / (label + '.command.json'), receipt)
    if result.returncode:
        raise RuntimeError('native command failed; output retained: ' + label)
    return receipt


def runtime(root, config):
    prefix = root / 'tools/rtg-tools-3.13'
    files = {name: fingerprint(prefix / name) for name in config['runtime_sha256']}
    for name, expected in config['runtime_sha256'].items():
        if files[name]['sha256'] != expected:
            raise ValueError('RTG runtime changed: ' + name)
    # The vendor wrapper reads optional system and host overrides after rtg.cfg.
    host = subprocess.check_output(['hostname', '-s'], text=True).strip()
    for p in (Path('/etc/rtg.cfg'), prefix / ('rtg.' + host + '.cfg')):
        if p.exists(): raise ValueError('unreviewed RTG configuration override: ' + str(p))
    java = fingerprint(root / 'tools/bio/lib/jvm/bin/java')
    if java['sha256'] != 'ec64ac719290ad080c7f8a4f5d83402e4152f98f486fe14150cf734a46ebcbb1':
        raise ValueError('RTG Java runtime changed')
    files['java'] = java
    return prefix / 'rtg', files


def record_class(record):
    alts = record.alts or ()
    if not alts or any(not a or any(c not in 'ACGTNacgtn' for c in a) for a in alts):
        return 'non_sequence'
    if len(record.ref) == 1 and all(len(a) == 1 for a in alts): return 'snv'
    if all(len(a) == len(record.ref) for a in alts): return 'mnp'
    if all(len(a) != len(record.ref) for a in alts): return 'indel'
    return 'mixed'


def audit_vcf(path, reference, sample, chromosome, require_raw_qual=False):
    """Read the entire supplied chr20 context, checking every REF and accounting GT/FILTER."""
    counts = Counter(); filters = Counter(); genotypes = Counter(); classes = Counter()
    previous = -1; minimum = None; maximum = None
    with pysam.VariantFile(str(path)) as vcf, pysam.FastaFile(str(reference)) as ref:
        if list(vcf.header.samples) != [sample]: raise ValueError('unexpected sample header')
        if vcf.header.contigs[chromosome].length != ref.get_reference_length(chromosome):
            raise ValueError('VCF/reference contig length differs')
        for record in vcf:
            if record.contig != chromosome or record.start < previous:
                raise ValueError('unexpected contig or unsorted VCF')
            previous = record.start
            if ref.fetch(chromosome, record.start, record.start + len(record.ref)).upper() != record.ref.upper():
                raise ValueError('VCF REF mismatch')
            counts['records'] += 1; classes[record_class(record)] += 1
            filt = ';'.join(record.filter.keys()) or '.'; filters[filt] += 1
            gt = record.samples[sample].get('GT'); genotypes[str(gt)] += 1
            fully_called_diploid = gt is not None and len(gt) == 2 and all(a is not None for a in gt)
            nonref = fully_called_diploid and any(a > 0 for a in gt)
            if fully_called_diploid: counts['fully_called_diploid'] += 1
            if nonref: counts['nonreference_diploid'] += 1
            if nonref and record_class(record) != 'non_sequence' and filt in ('.', 'PASS'):
                counts['eligible_sequence_records_before_BED'] += 1
                quality = record.qual
                if require_raw_qual and (quality is None or not math.isfinite(quality) or quality < 0):
                    raise ValueError('fixed raw QUAL threshold would exclude an eligible call')
                if quality is not None and math.isfinite(quality):
                    minimum = quality if minimum is None else min(minimum, quality)
                    maximum = quality if maximum is None else max(maximum, quality)
    return {'counts': dict(counts), 'classes': dict(classes), 'filters': dict(filters),
            'genotypes': dict(genotypes), 'eligible_QUAL_min': minimum, 'eligible_QUAL_max': maximum,
            'all_record_REF_alleles_verified': True}


def comparison_command(rtg, config, baseline, calls, sdf, bed, output, chromosome):
    return [rtg, 'vcfeval', '-b', baseline, '-c', calls, '-t', sdf, '-o', output,
            '--region', chromosome, '--evaluation-regions', bed] + config['vcfeval_arguments']


def count_outputs(output):
    counts = {}; subsets = {}; files = {}
    for key, name in [('tp_baseline', 'tp-baseline'), ('tp_calls', 'tp'), ('fp', 'fp'), ('fn', 'fn')]:
        p = output / (name + '.vcf.gz'); files[key] = fingerprint(p); kinds = Counter(); n = 0
        with pysam.VariantFile(str(p)) as vcf:
            for record in vcf:
                n += 1; kinds[record_class(record)] += 1
        counts[key] = n; subsets[key] = dict(kinds)
    def metrics(x):
        recall_den = x['tp_baseline'] + x['fn']; precision_den = x['tp_calls'] + x['fp']
        recall = x['tp_baseline'] / recall_den if recall_den else None
        precision = x['tp_calls'] / precision_den if precision_den else None
        return dict(x, truth_records=recall_den, call_records=precision_den, recall=recall, precision=precision)
    result = {'counts': metrics(counts), 'record_class_counts': subsets,
              'by_record_class': {kind: metrics({k: subsets[k].get(kind, 0) for k in counts})
                                  for kind in ('snv', 'indel', 'mnp', 'mixed', 'non_sequence')},
              'files': files, 'native_summary': (output / 'summary.txt').read_text()}
    return result
