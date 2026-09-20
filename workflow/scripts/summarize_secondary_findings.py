#!/usr/bin/env python3
"""Export a small evidence receipt from a completed local secondary triage.

The returned genotypes are limited to explicitly selected follow-up variants.
Absence from a variant-only VCF is never converted to a reference genotype.
Clinical interpretation belongs in the accompanying reviewed report.
"""
import argparse
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import pysam


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def extract_public_record(entry):
    path = Path(entry['path'])
    if sha256(path) != entry['sha256']:
        raise ValueError('public XML checksum mismatch')
    root = ET.parse(path).getroot()
    archive = root if root.tag == 'VariationArchive' else root.find('.//VariationArchive')
    if archive is None or archive.get('VariationID') != str(entry['variation_id']):
        raise ValueError('wrong ClinVar variation identity')
    germline = archive.find('./ClassifiedRecord/Classifications/GermlineClassification')
    return dict(entry, name=archive.get('VariationName'), type=archive.get('VariationType'),
        classification=germline.findtext('Description') if germline is not None else None,
        review_status=germline.findtext('ReviewStatus') if germline is not None else None,
        grch38_alleles=[dict(x.attrib) for x in archive.findall('./ClassifiedRecord//SequenceLocation')
                       if x.get('Assembly') == 'GRCh38' and x.get('positionVCF')],
        haplotype_components=[dict(x.attrib) for x in archive.findall('./ClassifiedRecord/Haplotype/SimpleAllele')])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run', type=Path, required=True)
    p.add_argument('--public-manifest', type=Path, required=True)
    p.add_argument('--check-variation', type=int, action='append', default=[])
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    if args.output.exists():
        p.error('existing output preserved; choose a new path')
    status = json.loads((args.run / 'run-status.json').read_text())
    if status['status'] != 'complete':
        raise ValueError('screen did not complete')
    if sha256(args.run / 'triage.json') != status['output_hashes']['triage.json']:
        raise ValueError('triage checksum mismatch')
    data = json.loads((args.run / 'triage.json').read_text())
    public = [extract_public_record(x) for x in json.loads(args.public_manifest.read_text())]
    by_id = {x['variation_id']: x for x in public}
    followup = []
    if args.check_variation:
        source = status['inputs']['vcf']
        if sha256(source['path']) != source['sha256']:
            raise ValueError('subject VCF checksum mismatch')
        with pysam.VariantFile(source['path']) as reader:
            for variant in args.check_variation:
                for allele in by_id[variant]['grch38_alleles']:
                    chrom, pos = allele['Chr'], int(allele['positionVCF'])
                    matches = []
                    for record in reader.fetch(chrom, pos - 1, pos):
                        if (record.pos == pos and record.ref == allele['referenceAlleleVCF']
                                and allele['alternateAlleleVCF'] in (record.alts or ())):
                            call = record.samples[data['sample']]
                            index = record.alts.index(allele['alternateAlleleVCF']) + 1
                            gt = call.get('GT') or ()
                            matches.append({'GT': list(gt), 'alt_copies': gt.count(index),
                                            'DP': call.get('DP'), 'GQ': call.get('GQ'),
                                            'FILTER': list(record.filter)})
                    followup.append({'variation_id': variant, 'GRCh38': allele,
                                     'exact_source_vcf_matches': matches,
                                     'absence_is_not_reference_or_callability_proof': True})
    result = {'screen': data, 'run': status, 'public_records': public,
              'followup_checks': followup, 'summarizer_sha256': sha256(__file__),
              'public_manifest_sha256': sha256(args.public_manifest),
              'interpretation': 'See the reviewed report; this receipt does not automatically classify a clinical finding.'}
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'output': str(args.output), 'screen_matches': len(data['matches']),
                      'public_records': len(public), 'followup_checks': followup}))


if __name__ == '__main__':
    main()
