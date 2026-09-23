#!/usr/bin/env python3
"""Read pinned Supplement 3 cells and replay six descriptive RNA ratios.

No workbook edits, downloads, normalization, differential expression or P/FDR.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET
import zipfile

EXPECTED_SHA256 = '50759d7180027b5be0cb11e5c4fdb226b69c37c320ea92995fbc4f3ebbbb1f00'
NS = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
REL = '{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id'


def read_workbook(path):
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != EXPECTED_SHA256:
        raise ValueError('Supplement 3 hash mismatch')
    result = {}
    with zipfile.ZipFile(path) as archive:
        strings = [''.join(e.itertext()) for e in
                   ET.fromstring(archive.read('xl/sharedStrings.xml'))]
        rels = {e.get('Id'): e.get('Target') for e in
                ET.fromstring(archive.read('xl/_rels/workbook.xml.rels'))}
        for sheet in ET.fromstring(archive.read('xl/workbook.xml')).find('m:sheets', NS):
            target = rels[sheet.get(REL)]
            target = target.lstrip('/') if target.startswith('/') else 'xl/' + target
            headers, hits = {}, []
            for row in ET.fromstring(archive.read(target)).findall('.//m:sheetData/m:row', NS):
                cells = {}
                for cell in row:
                    value = cell.find('m:v', NS)
                    value = value.text if value is not None else ''
                    if cell.get('t') == 's':
                        value = strings[int(value)]
                    cells[cell.get('r')] = value
                if int(row.get('r')) <= 2:
                    headers.update(cells)
                if 'BUB1B' in cells.values():
                    hits.append(cells)
            result[sheet.get('name')] = {'headers': headers, 'exact_symbol_rows': hits}
    return digest, result


def replay(path):
    digest, sheets = read_workbook(path)
    rows = sheets['total (25737)']['exact_symbol_rows']
    assert len(rows) == 1
    row = rows[0]
    assert row['B7664'] == row['X7664'] == 'BUB1B'
    assert row['Q7664'] == row['Y7664'] == '701'
    assert row['AE7664'] == 'NM_001211'
    headers = sheets['total (25737)']['headers']
    conditions = ['CON', 'LPS', 'PDLSC-CM', 'LPS+PDLSC-CM']
    assert [headers[c+'2'] for c in 'IJKL'] == conditions
    assert [headers[c+'2'] for c in 'MNOP'] == conditions
    pairs = [('C', 0, 1), ('D', 0, 2), ('E', 0, 3),
             ('F', 1, 2), ('G', 1, 3), ('H', 2, 3)]
    values = [float(row[c+'7664']) for c in 'IJKL']
    ratios = []
    for column, denominator, numerator in pairs:
        calculated = 2 ** (values[numerator] - values[denominator])
        reported = float(row[column+'7664'])
        assert math.isclose(calculated, reported, rel_tol=1e-10)
        ratios.append({'source_cell': column+'7664',
                       'source_header': headers[column+'2'],
                       'numerator': conditions[numerator],
                       'denominator': conditions[denominator],
                       'reported': reported, 'calculated': calculated})
    filtered = sheets['LPS VS. LPS+PDLSC-CM (776)']['exact_symbol_rows']
    assert len(filtered) == 1 and filtered[0]['B239'] == 'BUB1B'
    assert filtered[0]['C239'] == row['G7664']
    assert not sheets['CON VS. LPS (252)']['exact_symbol_rows']
    assert not sheets['CON VS. PDLSC-CM (1326)']['exact_symbol_rows']
    return {'status': 'PASS', 'source_sha256': digest, 'source_sheet': 'total (25737)',
            'source_row': 7664, 'source_values_verbatim': {k: v for k, v in row.items() if k in [c+'7664' for c in 'ABCDEFGHIJKLMNOPQ'] + ['X7664', 'Y7664', 'AE7664', 'AI7664', 'AJ7664', 'AK7664']},
            'conditions': [{'condition': condition, 'log2_cell': a+'7664',
                            'normalized_log2': float(row[a+'7664']),
                            'raw_cell': b+'7664', 'raw_count': int(row[b+'7664'])}
                           for condition, a, b in zip(conditions, 'IJKL', 'MNOP')],
            'ratios': ratios, 'filtered_source_row': filtered[0],
            'other_comparison_sheets_bub1b_present': False,
            'replicate_map': None, 'gene_specific_p': None, 'gene_specific_fdr': None,
            'scope': 'Descriptive replay of supplied processed RNA values only; no protein, causality, sample independence or functional inference.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workbook', type=Path, required=True)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    result = json.dumps(replay(args.workbook), indent=2) + '\n'
    if args.out:
        args.out.write_text(result)
    print(result, end='')
    print('END TARGETED BUB1B REPLAY')
