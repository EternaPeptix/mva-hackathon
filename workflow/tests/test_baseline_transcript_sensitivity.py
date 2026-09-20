import copy
import io
from pathlib import Path
import sys
import unittest

SCRIPTS = Path(__file__).resolve().parents[1] / 'scripts'
sys.path.insert(0, str(SCRIPTS))
import baseline_transcript_sensitivity as sensitivity
from rare_biallelic_baseline import models_at_frequency
from test_rare_biallelic_baseline import variant


class TranscriptSensitivityTests(unittest.TestCase):
    def test_exact_tag_and_transcript_identity(self):
        data = ('##gff-version 3\n'
                '1\tensembl\tmRNA\t1\t30\t.\t+\t.\tID=transcript:TX1;Parent=gene:G1;transcript_id=TX1;biotype=protein_coding;tag=MANE_Select,Ensembl_canonical;version=2\n'
                '1\tensembl\tmRNA\t1\t30\t.\t+\t.\tID=transcript:TX2;Parent=gene:G1;transcript_id=TX2;biotype=protein_coding;tag=MANE_Select_like\n'
                '1\tensembl\tmRNA\t1\t30\t.\t+\t.\tID=transcript:TXLOF;Parent=gene:G2;transcript_id=TXLOF;biotype=protein_coding_LoF;tag=Ensembl_canonical\n'
                '1\tensembl\tmRNA\t1\t30\t.\t+\t.\tID=transcript:TX3;Parent=gene:G1;transcript_id=TX3;biotype=nonsense_mediated_decay;tag=MANE_Select\n')
        transcripts, selected = sensitivity.tagged_transcripts(io.StringIO(data))
        self.assertEqual(set(transcripts), {'TX1', 'TX2', 'TXLOF'})
        self.assertEqual(transcripts['TX1']['version'], '2')
        self.assertEqual(selected['MANE_Select'], {'TX1'})
        self.assertEqual(selected['Ensembl_canonical'], {'TX1', 'TXLOF'})
        self.assertEqual(selected['GFF_exact_protein_coding'], {'TX1', 'TX2'})
        with self.assertRaisesRegex(ValueError, 'ambiguous'):
            sensitivity.tagged_transcripts(io.StringIO(data.replace('ID=transcript:TX1', 'ID=transcript:WRONG')))

    def test_restriction_keeps_missingness_and_phase_without_mutating_baseline(self):
        a = variant(1, {'ALTERNATIVE': 2, 'SELECT': 1}, af=None)
        b = variant(10, {'ALTERNATIVE': 2, 'SELECT': 1})
        source = [a, b]
        snapshot = copy.deepcopy(source)
        primary, _ = models_at_frequency(source, 0.001)
        filtered = sensitivity.restrict_transcripts(source, {'SELECT'})
        models, _ = models_at_frequency(filtered, 0.001)
        self.assertEqual(source, snapshot)
        self.assertEqual(models[0]['frequency_group'], 'frequency_unresolved')
        self.assertEqual(models[0]['phase'], 'unknown')
        self.assertEqual(models[0]['tier'], 3)
        self.assertEqual(sensitivity.verify_subset(primary, models)['tier_transitions'], {'1->3': 1})
        self.assertEqual(sensitivity.restrict_transcripts(source, {'ABSENT'}), [])
        models[0]['phase'] = 'trans_from_explicit_source_block'
        with self.assertRaisesRegex(ValueError, 'changed phase'):
            sensitivity.verify_subset(primary, models)

    def test_cannot_create_pair_from_different_selected_transcripts_or_promote_tier(self):
        a = variant(1, {'ALL_SHARED': 1, 'SELECT_A': 1})
        b = variant(10, {'ALL_SHARED': 1, 'SELECT_B': 1})
        primary, _ = models_at_frequency([a, b], 0.001)
        models, _ = models_at_frequency(sensitivity.restrict_transcripts([a, b], {'SELECT_A', 'SELECT_B'}), 0.001)
        self.assertEqual(models, [])
        self.assertEqual(sensitivity.verify_subset(primary, models)['removed_models'], 1)
        promoted = copy.deepcopy(primary)
        promoted[0]['tier'] = 1
        with self.assertRaisesRegex(ValueError, 'improved'):
            sensitivity.verify_subset(primary, promoted)
        promoted[0]['gene'] = 'DIFFERENT'
        with self.assertRaisesRegex(ValueError, 'new model'):
            sensitivity.verify_subset(primary, promoted)


if __name__ == '__main__':
    unittest.main()
