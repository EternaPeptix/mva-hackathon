"""Guard against the source-loss mechanism found during annotation audit."""
import hashlib
import unittest
from read_abstracts import select


class ExactAbstractTests(unittest.TestCase):
    def test_preserves_inequalities_markup_and_named_intervention(self):
        text = ('<h4>Results</h4>p < 0.05; rhGAA rescued accumulation; '
                'output > 2; APC/C<sup>Cdc20</sup>; &lt; 1.')
        payload = {'results': [{'source': 'TEST', 'id': '1', 'title': 'Synthetic fixture',
                                'abstractText': text}]}
        row = select(payload, ['TEST:1'])[0]
        self.assertEqual(row['abstractText'], text)
        self.assertEqual(row['abstract_sha256'], hashlib.sha256(text.encode()).hexdigest())
        self.assertEqual(row['abstract_characters'], len(text))

    def test_missing_duplicate_and_empty_inputs_fail(self):
        record = {'source': 'TEST', 'id': '1', 'title': 'Synthetic fixture', 'abstractText': 'x'}
        for payload, requested in [({'results': []}, ['TEST:1']),
                                   ({'results': [record, record]}, ['TEST:1']),
                                   ({'results': [record]}, ['TEST:1', 'TEST:1']),
                                   ({'results': [dict(record, abstractText='')]}, ['TEST:1']),
                                   ({'results': [record]}, [])]:
            with self.subTest(payload=payload, requested=requested):
                with self.assertRaises(ValueError):
                    select(payload, requested)


if __name__ == '__main__':
    unittest.main()
