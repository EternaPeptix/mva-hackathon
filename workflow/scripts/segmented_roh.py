#!/usr/bin/env python3
"""Run independent native ROH models on blocks of already-qualified markers.

This component does not perform population-panel selection, patient genotyping,
read QC, ancestry matching, clinical classification or parental-origin inference.
Its input must already have independently sourced POP_AF and usable PL values.
"""
import argparse
from collections import Counter
import datetime
import hashlib
import json
import math
from pathlib import Path
import subprocess
import time

import pysam

POPULATION_TAGS = ('POP_AF', 'POP_AFR_AF', 'POP_AMR_AF', 'POP_EAS_AF', 'POP_EUR_AF', 'POP_SAS_AF')

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def inspect_marker(record, sample, af_tag='POP_AF'):
    if af_tag not in POPULATION_TAGS:
        raise ValueError('explicit independent population tag required')
    if (len(record.alleles) != 2 or record.ref == record.alts[0]
            or any(len(a) != 1 or a not in 'ACGT' for a in record.alleles)):
        raise ValueError('prequalified biallelic SNPs required')
    af = record.info.get(af_tag)
    if isinstance(af, tuple):
        af = af[0] if len(af) == 1 else None
    if af is None or not math.isfinite(af) or not 0 < af < 1:
        raise ValueError('explicit nonfixed '+af_tag+' required; no sample-AF fallback')
    pl = record.samples[sample].get('PL')
    if (pl is None or len(pl) != 3 or any(type(x) is not int or x < 0 for x in pl)
            or min(pl) != 0
            or len(set(pl)) == 1):
        raise ValueError('usable diploid PL required; missing/flat likelihoods must be reported upstream')
    # Missing GT is allowed if valid likelihoods exist; native ROH uses the PLs.
    return af


def run_segmented(source, sample, executable, output, max_gap_bp=250000, min_markers=20, af_tag='POP_AF'):
    if max_gap_bp < 1 or min_markers < 2:
        raise ValueError('positive gap threshold and at least two markers required')
    output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    before = source.stat(); input_sha = digest(source)
    runs, gaps, skipped, counts = [], [], [], Counter()
    with pysam.VariantFile(str(source)) as reader:
        if list(reader.header.samples) != [sample]:
            raise ValueError('exact single-sample input required')
        if af_tag not in POPULATION_TAGS or af_tag not in reader.header.info or 'PL' not in reader.header.formats:
            raise ValueError('explicit population AF and PL headers required')
        order = {chrom:i for i,chrom in enumerate(reader.header.contigs)}
        completed_chroms = set()
        previous = None
        segment = []

        def flush():
            if not segment:
                return
            meta = {'chrom':segment[0].contig, 'start_1based':segment[0].pos,
                    'end_inclusive':segment[-1].pos, 'markers':len(segment)}
            if len(segment) < min_markers:
                skipped.append(dict(meta,reason='insufficient markers for the configured analysis'))
                counts['skipped_markers'] += len(segment)
                segment.clear()
                return
            folder = output / ('segment-%04d' % len(runs)); folder.mkdir()
            vcf = folder/'input.vcf'; native = folder/'roh.txt'
            with pysam.VariantFile(str(vcf),'w',header=reader.header) as writer:
                for record in segment:
                    writer.write(record)
            command = [str(executable),'roh','--skip-indels','--samples',sample,
                       '--AF-tag',af_tag,'--hw-to-az','6.7e-8','--az-to-hw','5e-9',
                       '--output-type','sr','--output',str(native),str(vcf)]
            tick = time.monotonic()
            result = subprocess.run(command,capture_output=True,text=True)
            (folder/'stderr.txt').write_text(result.stderr)
            if result.returncode:
                raise RuntimeError('native ROH segment failed')
            regions, states, observed_positions = [], Counter(), []
            for line in native.read_text().splitlines():
                fields = line.split('\t')
                if fields[0] == 'ST':
                    if len(fields) != 6 or fields[1] != sample or fields[2] != meta['chrom'] or fields[4] not in ('0','1'):
                        raise ValueError('unexpected native site schema')
                    observed_positions.append(int(fields[3]));states[fields[4]] += 1
                elif fields[0] == 'RG':
                    if len(fields) != 8 or fields[1] != sample or fields[2] != meta['chrom']:
                        raise ValueError('unexpected native region schema')
                    region = {'chrom':fields[2],'start_1based':int(fields[3]),
                              'end_inclusive':int(fields[4]),'length_bp':int(fields[5]),
                              'markers':int(fields[6]),'mean_model_phred':float(fields[7])}
                    if not (meta['start_1based'] <= region['start_1based'] <= region['end_inclusive'] <= meta['end_inclusive']):
                        raise ValueError('native interval escapes independently analyzed segment')
                    if region['length_bp'] != region['end_inclusive']-region['start_1based']+1 or not math.isfinite(region['mean_model_phred']):
                        raise ValueError('invalid native region values')
                    regions.append(region)
            if observed_positions != [record.pos for record in segment]:
                raise ValueError('native HMM silently omitted or reordered qualified markers')
            runs.append(dict(meta,command=command,elapsed_seconds=time.monotonic()-tick,
                             regions=regions,site_states=dict(states),input_sha256=digest(vcf),
                             output_sha256=digest(native),stderr_sha256=digest(folder/'stderr.txt')))
            counts['analyzed_markers'] += len(segment)
            segment.clear()

        for record in reader:
            inspect_marker(record,sample,af_tag)
            position = (order[record.contig],record.pos)
            if previous and position <= previous:
                raise ValueError('strictly ordered unique marker positions required')
            if segment:
                if record.contig != segment[-1].contig:
                    completed_chroms.add(segment[-1].contig)
                    flush()
                elif record.pos-segment[-1].pos > max_gap_bp:
                    gaps.append({'chrom':record.contig,'left_marker_1based':segment[-1].pos,
                                 'right_marker_1based':record.pos,'marker_distance_bp':record.pos-segment[-1].pos})
                    flush()
            if record.contig in completed_chroms:
                raise ValueError('chromosome revisited after completion')
            segment.append(record.copy());previous=position;counts['input_markers'] += 1
        flush()
    if not counts['input_markers'] or counts['input_markers'] != counts['analyzed_markers']+counts['skipped_markers']:
        raise ValueError('marker accounting failed')
    after = source.stat()
    if (before.st_size,before.st_mtime_ns) != (after.st_size,after.st_mtime_ns):
        raise ValueError('input changed during analysis')
    result = {'completed_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'elapsed_seconds':time.monotonic()-started,'input_sha256':input_sha,
        'script_sha256':digest(Path(__file__)),'bcftools_sha256':digest(executable),
        'max_marker_gap_bp':max_gap_bp,'minimum_segment_markers':min_markers,'population_af_tag':af_tag,
        'status':'complete_with_analyzed_segments' if runs else 'no_analyzable_segments',
        'counts':dict(counts),'gap_boundaries':gaps,'skipped_segments':skipped,'segments':runs,
        'scope':'Native ROH models independently restarted at large marker gaps. Upstream frequency provenance and genotype/read quality must be established separately; not UPD or a clinical negative test.'}
    (output/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--vcf',type=Path,required=True)
    p.add_argument('--expected-sha256',required=True)
    p.add_argument('--sample',required=True)
    p.add_argument('--bcftools',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--max-gap-bp',type=int,default=250000)
    p.add_argument('--min-markers',type=int,default=20)
    p.add_argument('--af-tag',choices=POPULATION_TAGS,default='POP_AF')
    a = p.parse_args()
    if digest(a.vcf) != a.expected_sha256:
        raise ValueError('input checksum differs from supplied provenance')
    result = run_segmented(a.vcf,a.sample,a.bcftools,a.output,a.max_gap_bp,a.min_markers,a.af_tag)
    print(json.dumps({'segments':len(result['segments']),'counts':result['counts'],
                      'gap_boundaries':len(result['gap_boundaries'])}))


if __name__ == '__main__':
    main()
