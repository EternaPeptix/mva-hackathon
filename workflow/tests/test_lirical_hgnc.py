import importlib.util
from pathlib import Path
import sys
import unittest

SCRIPTS = Path(__file__).parents[1] / 'scripts'
sys.path.insert(0, str(SCRIPTS))
from prepare_lirical_data import adapt_hgnc, require_no_disease_associations
import tempfile


class HgncAdapterTests(unittest.TestCase):
    def fixture(self):
        columns = ['hgnc_id', 'symbol'] + ['field%d' % i for i in range(2, 17)] + ['entrez_id', 'ensembl_gene_id', 'final_field']
        row = ['HGNC:1149', 'BUB1B'] + ['value%d' % i for i in range(2, 17)] + ['701', 'ENSG00000156970', '']
        return columns, row

    def test_column_order_without_biological_edits(self):
        columns, row = self.fixture()
        text = '\t'.join(columns) + '\n' + '\t'.join(row) + '\n'
        result, audit = adapt_hgnc(text)
        h, r = [line.split('\t') for line in result.splitlines()]
        self.assertEqual(r[18], '701')
        self.assertEqual(r[1], 'BUB1B')
        self.assertEqual(dict(zip(h, r)), dict(zip(columns, row)))
        self.assertEqual(audit['rows_added'], 0)
        self.assertEqual(audit['rows_quarantined'], 0)
        self.assertEqual(adapt_hgnc(result)[0], result)

    def test_true_duplicates_or_nonnumeric_ids_fail_instead_of_being_dropped(self):
        columns, row = self.fixture()
        second = list(row)
        second[0] = 'HGNC:2'
        second[1] = 'OTHER'
        text = '\n'.join('\t'.join(r) for r in (columns, row, second)) + '\n'
        with self.assertRaisesRegex(ValueError, 'duplicate Entrez'):
            adapt_hgnc(text)
        row[17] = 'ENSG00000156970'
        with self.assertRaisesRegex(ValueError, 'invalid or duplicate Entrez'):
            adapt_hgnc('\t'.join(columns) + '\n' + '\t'.join(row) + '\n')

    def test_explicit_quarantine_preserves_evidence_and_refuses_disease_association(self):
        columns, row = self.fixture()
        second = list(row)
        second[0] = 'HGNC:2'
        second[1] = 'OTHER'
        text = '\n'.join('\t'.join(r) for r in (columns, row, second)) + '\n'
        result, audit = adapt_hgnc(text, quarantine_ambiguous=True)
        self.assertEqual(audit['rows_quarantined'], 2)
        self.assertEqual(audit['quarantined'][0]['symbol'], 'BUB1B')
        self.assertEqual(len(result.splitlines()), 1)
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory)
            (source / 'mim2gene_medgen').write_text('123456\t701\tphenotype\n')
            (source / 'en_product6.xml').write_text('<Root/>')
            with self.assertRaisesRegex(ValueError, 'MedGen disease association'):
                require_no_disease_associations(source, audit['quarantined'])
            (source / 'mim2gene_medgen').write_text('#header\n')
            (source / 'en_product6.xml').write_text('<Root><Gene><Symbol>BUB1B</Symbol></Gene></Root>')
            with self.assertRaisesRegex(ValueError, 'Orphanet disease association'):
                require_no_disease_associations(source, audit['quarantined'])


if __name__ == '__main__':
    unittest.main()
