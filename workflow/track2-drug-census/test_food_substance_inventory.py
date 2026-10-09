import html
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import food_substance_inventory as f
from food_substance_inventory import HEADERS, parse_source


def page(names, declared=None):
    header = '<tr>' + ''.join('<th>' + html.escape(x) + '</th>' for x in HEADERS) + '</tr>'
    rows = ''.join('<tr><td></td><td><a href="?set=FoodSubstances&amp;id=ROW' + str(i) + '">' +
                   html.escape(name) + '</a></td><td>Source effect</td><td>189.1</td></tr>'
                   for i, name in enumerate(names))
    # FDA's fixed-column header duplicate must not double the rows.
    return (f'Records Found: {len(names) if declared is None else declared}' +
            '<table>' + header + rows + '</table><table>' + header + '</table>').encode()


class FoodSourceTests(unittest.TestCase):
    def test_status_is_retained_while_only_known_suffix_is_removed(self):
        names = ['COUMARIN--PROHIBITED', 'COBALT SULFATE--PROHIBITED WITH EXCEPTIONS',
                 'ACETAMIDE--NLFG', 'FD&C RED NO. 1--DELISTED',
                 'POLYMER--DVB-(DEG-DIVINYL ETHER)', 'L-ARGININE']
        rows, _ = parse_source(page(names))
        self.assertEqual([r['substance_as_reported'] for r in rows], names)
        self.assertEqual(rows[0]['retrieval_name'], 'COUMARIN')
        self.assertEqual(rows[1]['explicit_name_status'], 'PROHIBITED WITH EXCEPTIONS')
        self.assertEqual(rows[2]['explicit_name_status'], 'NLFG')
        self.assertEqual(rows[3]['explicit_name_status'], 'DELISTED')
        self.assertEqual(rows[4]['retrieval_name'], names[4])
        self.assertEqual(rows[5]['current_permission'], 'NOT_DETERMINED_BY_INVENTORY')

    def test_partial_page_or_changed_header_fails(self):
        with self.assertRaisesRegex(ValueError, 'header/count mismatch'):
            parse_source(page(['ONE'], declared=50))
        with self.assertRaisesRegex(ValueError, 'header/count mismatch'):
            parse_source(page(['ONE']).replace(b'21 CFR*', b'Changed'))

    def test_duplicate_labels_and_blank_identifiers_are_not_merged(self):
        rows, _ = parse_source(page(['MIXTURE', 'MIXTURE']))
        self.assertEqual(len(rows), 2)
        self.assertEqual([r['cas_or_other_id_as_reported'] for r in rows], ['', ''])
        self.assertNotEqual(rows[0]['source_record_id'], rows[1]['source_record_id'])
        self.assertTrue(all(r['BUB1B_efficacy'] == 'UNKNOWN_NOT_ASSESSED_BY_INVENTORY' for r in rows))

    def test_hyphenated_amino_acids_and_common_classes_are_ambiguous(self):
        for name in ['L-CYSTEINE', 'L-ARGININE', 'L-GLUTAMINE', 'CASEIN', 'COLLAGEN']:
            self.assertEqual(f.retrieval_policy(name), 'COMMON_OR_SHORT_NAME_CONTEXT_AMBIGUOUS')
        self.assertEqual(f.retrieval_policy('HEXYL ISOTHIOCYANATE'), 'LITERAL_NAME_CONTEXT_UNASSESSED')

    def test_input_change_failure_preserves_existing_output_package(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / 'food.html').write_bytes(page(['L-ARGININE']))
            (root / 'corpus.json').write_text('{"results": []}')
            (root / 'catalog.tsv').write_text('catalog_id\tpreferred_name\n')
            (root / 'aliases.tsv').write_text('catalog_id\talias\n')
            (root / 'drugcentral-status-unverified.tsv').write_text('drugcentral_id\tname\taliases\n')
            (root / 'decisions.json').write_text('{"decisions": []}')
            out = root / 'output'
            out.mkdir()
            previous = {n: ('PREVIOUS ' + n).encode() for n in ['substances.tsv', 'mentions.tsv', 'summary.json']}
            for n, b in previous.items():
                (out / n).write_bytes(b)
            real_write = f.write_tsv
            def mutate_input(*args):
                real_write(*args)
                (root / 'decisions.json').write_text('{"decisions": [], "changed": true}')
            with patch.object(f, 'write_tsv', side_effect=mutate_input):
                with self.assertRaisesRegex(ValueError, 'Input changed'):
                    f.build(root / 'food.html', root / 'corpus.json', root, out)
            self.assertEqual({n: (out / n).read_bytes() for n in previous}, previous)
            self.assertEqual(list(root.glob('.fda-food-*')), [])


if __name__ == '__main__':
    unittest.main()
