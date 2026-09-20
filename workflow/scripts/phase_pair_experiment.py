#!/usr/bin/env python3
"""Post hoc compatibility experiment: expose one verified FANCD2 cis block in GT."""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import time

import pysam

TARGETS = {('3', '10046723', 'AG', 'A'), ('3', '10046725', 'TAAG', 'T')}
SOURCE_SHA256 = 'bec64c68d3ad10d05b8071d07cdf87d217013d66d3ba02e24b505035c9d42175'


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def phase_line(line):
    if line.startswith('#') or not line.startswith('3\t1004672'):
        return line, None
    fields = line.rstrip('\n').split('\t')
    key = (fields[0], fields[1], fields[3], fields[4])
    if key not in TARGETS:
        return line, None
    if len(fields) != 10:
        raise ValueError('expected exactly one sample')
    keys, values = fields[8].split(':'), fields[9].split(':')
    sample = dict(zip(keys, values))
    if sample.get('GT') != '0/1' or sample.get('PGT') != '0|1' or sample.get('PID') != '10046720_C_T':
        raise ValueError('candidate phase annotations differ from the reviewed evidence')
    values[keys.index('GT')] = sample['PGT']
    fields[9] = ':'.join(values)
    return '\t'.join(fields)+'\n', key


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    source, output = args.source.resolve(), args.output.resolve()
    temporary = Path(str(output)+'.partial')
    if output.exists() or temporary.exists():
        parser.error('existing output or partial experiment preserved')
    if digest(source) != SOURCE_SHA256:
        parser.error('source does not match publisher-verified VCF')
    output.parent.mkdir(parents=True, exist_ok=True)
    os.nice(10)
    changed, count = [], 0
    with gzip.open(source, 'rt') as reader, pysam.BGZFile(str(temporary), 'wb') as writer:
        for line in reader:
            if line.startswith('#CHROM') and line.rstrip('\n').split('\t')[9:] != ['WGS_EX2312012']:
                raise ValueError('unexpected sample header')
            new, key = phase_line(line)
            writer.write(new.encode())
            if key:
                changed.append(key)
            if not line.startswith('#'):
                count += 1
    if len(changed) != 2 or set(changed) != TARGETS or count != 5012204:
        raise ValueError('unexpected changed sites or record count; partial output preserved')
    # Validate the two changes and a readable coordinate index before publishing.
    pysam.tabix_index(str(temporary), preset='vcf', force=False)
    with pysam.VariantFile(str(temporary)) as vcf:
        found = [rec for rec in vcf.fetch('3', 10046722, 10046726)
                 if (rec.chrom, str(rec.pos), rec.ref, rec.alts[0]) in TARGETS]
        if len(found) != 2 or any(not rec.samples[0].phased or rec.samples[0]['GT'] != (0,1) for rec in found):
            raise ValueError('derived phased GT validation failed')
    os.replace(temporary, output)
    os.replace(str(temporary)+'.tbi', str(output)+'.tbi')
    result = {'scope':'Post hoc targeted phase-compatibility experiment; not a general phasing pipeline.',
              'source_sha256':SOURCE_SHA256, 'derived_sha256':digest(output),
              'records':count, 'changed_fields':[{'contig':c,'position':int(p),'ref':r,'alt':a,'GT_before':'0/1','GT_after':'0|1'} for c,p,r,a in changed],
              'all_other_decompressed_lines_and_fields_preserved':True,
              'script_sha256':digest(Path(__file__)),
              'created_at':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())}
    output.with_suffix(output.suffix+'.experiment.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
