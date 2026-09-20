#!/usr/bin/env python3
"""Preserve downloaded resources and adapt HGNC columns for Phenol 2.1.4.

The bundled reader uses columns 0, 1 and 18, not header names. Current HGNC
puts entrez_id at 17. Reorder columns in a derived copy; never drop duplicate
Ensembl IDs or edit biological identifiers to work around the failure.
"""
import argparse
import collections
import csv
import json
from pathlib import Path
import re
import shutil
import xml.etree.ElementTree as ET

from run_lirical import DATA_FILES, digest, stamp


def adapt_hgnc(text, quarantine_ambiguous=False):
    lines = text.splitlines()
    if not lines:
        raise ValueError('empty HGNC table')
    columns = lines[0].split('\t')
    if len(set(columns)) != len(columns) or columns[:2] != ['hgnc_id', 'symbol'] or 'entrez_id' not in columns:
        raise ValueError('unexpected HGNC header')
    if len(columns) < 19:
        raise ValueError('HGNC table has too few columns')
    adjusted = list(columns)
    adjusted.remove('entrez_id')
    adjusted.insert(18, 'entrez_id')
    order = [columns.index(name) for name in adjusted]
    output = ['\t'.join(adjusted)]
    parsed = [line.split('\t') for line in lines[1:]]
    if any(len(values) != len(columns) for values in parsed):
        raise ValueError('HGNC row length mismatch')
    entrez_index = columns.index('entrez_id')
    counts = collections.Counter(values[entrez_index] for values in parsed if values[entrez_index])
    ambiguous = {identifier for identifier, count in counts.items() if count > 1}
    if ambiguous and not quarantine_ambiguous:
        raise ValueError('duplicate Entrez identifier; explicit audited quarantine required')
    ids, hgnc_ids, rows, missing, quarantined = set(), set(), 0, 0, []
    for line in lines[1:]:
        values = line.split('\t')
        if len(values) != len(columns):
            raise ValueError('HGNC row length mismatch')
        original = dict(zip(columns, values))
        identifier = original['entrez_id']
        if not re.fullmatch(r'HGNC:\d+', original['hgnc_id']) or not original['symbol']:
            raise ValueError('invalid HGNC identifier or symbol')
        if original['hgnc_id'] in hgnc_ids:
            raise ValueError('duplicate HGNC identifier')
        hgnc_ids.add(original['hgnc_id'])
        if identifier:
            if not re.fullmatch(r'[1-9]\d*', identifier):
                raise ValueError('invalid or duplicate Entrez identifier')
            if identifier not in ambiguous:
                ids.add(identifier)
        else:
            missing += 1
        reordered = [values[index] for index in order]
        if dict(zip(adjusted, reordered)) != original:
            raise AssertionError('adaptation changed a field value')
        if identifier in ambiguous:
            quarantined.append(original)
        else:
            output.append('\t'.join(reordered))
        rows += 1
    return '\n'.join(output) + '\n', {
        'input_rows': rows, 'output_rows': rows - len(quarantined),
        'nonempty_unique_entrez_ids': len(ids), 'missing_entrez_ids': missing,
        'source_entrez_column_zero_based': columns.index('entrez_id'),
        'derived_entrez_column_zero_based': 18,
        'all_retained_header_named_field_values_preserved': True,
        'rows_added': 0, 'rows_quarantined': len(quarantined),
        'quarantined': quarantined}


def require_no_disease_associations(source, quarantined):
    """Refuse quarantine of any mapping linked in the pinned disease resources."""
    ids = {r['entrez_id'] for r in quarantined}
    symbols = {r['symbol'] for r in quarantined}
    hgnc = {r['hgnc_id'].split(':')[1] for r in quarantined}
    with (source / 'mim2gene_medgen').open() as stream:
        if any(len(row) > 1 and row[1] in ids for row in csv.reader(stream, delimiter='\t')):
            raise ValueError('ambiguous ID has a MedGen disease association; do not quarantine')
    for gene in ET.parse(source / 'en_product6.xml').iter('Gene'):
        references = [(ref.findtext('Source'), ref.findtext('Reference'))
                      for ref in gene.iter('ExternalReference')]
        if gene.findtext('Symbol') in symbols or any(
                value in ids or (system == 'HGNC' and value in hgnc) for system, value in references):
            raise ValueError('ambiguous mapping has an Orphanet disease association; do not quarantine')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--quarantine-ambiguous', action='store_true')
    args = parser.parse_args()
    source, output = args.source.resolve(), args.output.resolve()
    inputs = {name: digest(source / name) for name in DATA_FILES}
    converted, audit = adapt_hgnc((source / 'hgnc_complete_set.txt').read_text(), args.quarantine_ambiguous)
    require_no_disease_associations(source, audit['quarantined'])
    audit['quarantined_medgen_or_orphanet_associations'] = 0
    output.mkdir(parents=True, exist_ok=False)
    for name in DATA_FILES:
        if name == 'hgnc_complete_set.txt':
            (output / name).write_text(converted)
        else:
            shutil.copy2(source / name, output / name)
    if any(digest(source / name) != inputs[name] for name in DATA_FILES):
        raise ValueError('download resources changed during snapshot')
    outputs = {name: digest(output / name) for name in DATA_FILES}
    if any(inputs[name] != outputs[name] for name in DATA_FILES if name != 'hgnc_complete_set.txt'):
        raise ValueError('non-HGNC resource changed during copy')
    report = {'created_at': stamp(), 'script_sha256': digest(Path(__file__)),
              'source_directory': str(source), 'source_hashes': inputs,
              'derived_hashes': outputs, 'hgnc_audit': audit,
              'reader_source': 'https://raw.githubusercontent.com/monarch-initiative/phenol/v2.1.4/phenol-annotations/src/main/java/org/monarchinitiative/phenol/annotations/assoc/HGNCGeneIdentifierLoader.java',
              'purpose': 'Column compatibility plus explicit quarantine of ambiguous reference mappings with no associations in pinned MedGen/Orphanet. No biological field values edited; all exclusions retained in this manifest.'}
    (output / 'compatibility-manifest.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(audit))


if __name__ == '__main__':
    main()
