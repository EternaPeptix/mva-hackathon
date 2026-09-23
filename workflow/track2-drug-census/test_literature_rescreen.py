"""Regressions for the actual inequality-loss bug and unsafe source-note reuse."""
import unittest
import tempfile
from pathlib import Path
from literature_rescreen import clean, normalize, matches, make_trie, role_index, digest, triage, assert_unchanged


class RescreenTests(unittest.TestCase):
    def test_inequality_prose_survives(self):
        raw = 'P < 0.05; eplerenone reduced cycling. <b>BUB1B</b> was measured; n > 5.'
        self.assertIn('eplerenone', normalize(raw))
        self.assertIn('BUB1B', clean(raw))
        self.assertIn('< 0.05', clean(raw))
        trie = make_trie({'eplerenone': {}, 'bub1b': {}})
        self.assertEqual(matches(normalize(raw), trie), {'eplerenone', 'bub1b'})

    def test_encoded_formatting_and_comparisons(self):
        raw = 'P &lt; 0.05; &lt;italic&gt;selumetinib&lt;/italic&gt;, survival &gt; 50%.'
        self.assertEqual(normalize(raw), 'p 0 05 selumetinib survival 50')
        self.assertIn('< 0.05', clean(raw))
        self.assertNotIn('italic', clean(raw))

    def test_alias_boundaries(self):
        trie = make_trie({'metformin': {}, 'mitomycin c': {}})
        self.assertEqual(matches(normalize('metforminate; mitomycin'), trie), set())
        self.assertEqual(matches(normalize('metformin and mitomycin-C'), trie), {'metformin', 'mitomycin c'})

    def test_wrong_source_notes_fail(self):
        note = dict(source='MED', record_id='1', abstract_sha256=digest(b'old text'), recovered_alias_roles=[],
                    annotation_status='source_revalidated', review_depth='complete_saved_abstract')
        review = dict(status='source_role_review_complete', corpus_sha256='same', records=[note])
        with self.assertRaisesRegex(ValueError, 'abstract changed'):
            role_index(review, {('MED', '1'): {'abstractText': 'new text'}}, 'same')
        with self.assertRaisesRegex(ValueError, 'another corpus'):
            role_index(review, {}, 'different')

    def test_rejected_unread_and_duplicate_notes_fail(self):
        note = dict(source='MED', record_id='1', abstract_sha256=digest(b'text'), recovered_alias_roles=[],
                    annotation_status='source_revalidated', review_depth='complete_saved_abstract')
        records = {('MED', '1'): {'abstractText': 'text'}}
        review = dict(status='source_role_review_complete', corpus_sha256='same', records=[note])
        for field, invalid in [('annotation_status', 'rejected'), ('review_depth', 'not_read')]:
            changed = dict(review, records=[dict(note, **{field: invalid})])
            with self.assertRaisesRegex(ValueError, 'not a complete revalidated reading'):
                role_index(changed, records, 'same')
        with self.assertRaisesRegex(ValueError, 'Duplicate source-role note'):
            role_index(dict(review, records=[note, note]), records, 'same')

    def test_concurrent_input_change_fails(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/'policy.tsv'
            path.write_bytes(b'read input')
            snapshot = {path: path.read_bytes()}
            assert_unchanged(snapshot)
            path.write_bytes(b'later input')
            with self.assertRaisesRegex(ValueError, 'Input changed during replay'):
                assert_unchanged(snapshot)

    def test_drug_specific_recovered_role(self):
        row = dict(catalog_id='X', drug_name='EOSIN', source='MED', record_id='1', pmid='1',
                   doi='', title='Study', matched_normalized_aliases=['eosin'],
                   source_url='https://europepmc.org/article/MED/1', raw_abstract_sha256='hash',
                   direct_gene_name_in_title_or_abstract=True)
        note = dict(uncertainty='No eosin treatment benefit.', followup_reason='No therapeutic exposure.', requires_primary_followup=False)
        roles = {('X', 'MED', '1'): [(note, dict(normalized_alias='eosin', role='assay_reagent',
                                                       role_evidence='Histological stain.'))]}
        result, _ = triage(row, {'abstractText': 'BUB1B and eosin staining improved.'}, roles)
        self.assertEqual(result['triage_category'], 'assay_reagent')
        self.assertEqual(result['source_revalidated_aliases'], ['eosin'])
        self.assertEqual(result['actual_BUB1B_variant_muscle_rescue'], 'NOT_ESTABLISHED_BY_THIS_ROW')
        roles[('X', 'MED', '1')][0][1]['normalized_alias'] = 'other drug'
        with self.assertRaisesRegex(ValueError, 'alias missing'):
            triage(row, {'abstractText': 'BUB1B and eosin'}, roles)


if __name__ == '__main__':
    unittest.main()
