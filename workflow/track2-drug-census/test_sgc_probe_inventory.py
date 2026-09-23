"""Regression checks for identity and coverage failures in the SGC screen."""
import csv
import html
import json
import tempfile
import unittest
from pathlib import Path

import sgc_probe_inventory as s


def table(source, rows):
    header = '<tr>' + ''.join('<th>' + h + '</th>' for h in s.SOURCES[source][1]) + '</tr>'
    body = ''
    for name, family, target, control in rows:
        cells = ['<a href="/probe">' + html.escape(name) + '</a>', '',
                 family, target, html.escape(control), '', '']
        body += '<tr>' + ''.join('<td>' + c + '</td>' for c in cells) + '</tr>'
    return ('<table>' + header + body + '</table>').encode()


class IdentityCoverageTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        (self.root / 'catalog.tsv').write_text('catalog_id\tpreferred_name\nC1\t(+)-JQ1\n')
        (self.root / 'aliases.tsv').write_text('catalog_id\talias\nC1\t(+)-JQ1\n')
        (self.root / 'drugcentral-status-unverified.tsv').write_text('drugcentral_id\tname\taliases\n')
        (self.root / 'decisions.json').write_text('{"decisions": []}')

    def run_screen(self, rows, records):
        for kind in s.SOURCES:
            (self.root / (kind + '.html')).write_bytes(table(kind, rows))
        (self.root / 'corpus.json').write_text(json.dumps({'results': records}))
        return s.build(self.root / 'all.html', self.root / 'donated.html',
                       self.root / 'corpus.json', self.root, self.root / 'out')

    def read(self, name):
        return s.read_tsv(self.root / 'out' / name)

    def test_negative_enantiomer_does_not_become_positive_probe(self):
        records = [{'source': 'MED', 'id': str(i), 'title': text, 'abstractText': ''}
                   for i, text in enumerate(['(-)-JQ1 and BubR1', '(+)-JQ1', 'JQ1'])]
        self.run_screen([('(+)-JQ1', 'BET', 'BRD4', '(-)-JQ1')], records)
        hits = {r['record_id']: r for r in self.read('mentions.tsv')}
        self.assertEqual(hits['0']['role'], 'negative_control')
        self.assertEqual(hits['1']['role'], 'probe')
        self.assertEqual(hits['2']['match_policy'], 'ambiguous_unspecified_stereochemistry')
        self.assertNotEqual(s.identity('(+)-JQ1'), s.identity('(-)-JQ1'))

    def test_code_boundary_and_gene_identity(self):
        self.assertIsNone(s.combined_pattern(['BAY-179']).search('bay-1797'))
        self.assertIsNone(s.GENE.search('bub1 kinase'))
        self.assertIsNotNone(s.GENE.search('bub1b'))
        self.assertIn('<0.05', s.clean('GSK343 p <0.05; intact following words'))
        self.assertIn('following words', s.clean('GSK343 p <0.05; intact following words'))

    def test_duplicate_annotation_rows_preserved_without_double_count(self):
        summary = self.run_screen([('BAY1816032', 'Kinase', 'BUB1', ''),
                                   ('BAY1816032', '2OG', 'BUB1', '')], [])
        self.assertEqual(summary['total_source_rows'], 4)
        self.assertEqual(summary['distinct_probe_labels'], 1)
        self.assertEqual(summary['source_rows_with_literal_BUB1B_target'], 0)
        self.assertEqual(summary['multiple_family_labels'], ['BAY1816032'])
        self.assertEqual(self.read('mentions.tsv'), [])
        self.assertEqual(self.read('controls.tsv'), [])

    def test_contextual_control_is_not_global_compound(self):
        summary = self.run_screen([('LP99', 'BET', 'BRD7', '2S,3R-enantiomer')],
                                 [{'source': 'MED', 'id': '1', 'title': '2S,3R-enantiomer', 'abstractText': ''}])
        self.assertEqual(summary['contextual_control_descriptions_not_screened'], 1)
        self.assertEqual(summary['mention_rows'], 0)

    def test_changed_header_and_dangling_alias_fail_explicitly(self):
        with self.assertRaises(ValueError):
            s.parse_source(table('all', [('X', '', '', '')]).replace(b'Target Protein', b'Changed'), 'all')
        (self.root / 'aliases.tsv').write_text('catalog_id\talias\nMISSING\tX\n')
        with self.assertRaisesRegex(ValueError, 'missing catalog ID'):
            self.run_screen([('X', '', '', '')], [])


if __name__ == '__main__':
    unittest.main()
