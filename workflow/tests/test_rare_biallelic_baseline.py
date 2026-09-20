import importlib.util
import json
from pathlib import Path
import sys
import unittest

SCRIPTS = Path(__file__).resolve().parents[1] / 'scripts'
sys.path.insert(0, str(SCRIPTS))
spec = importlib.util.spec_from_file_location('baseline', SCRIPTS / 'rare_biallelic_baseline.py')
baseline = importlib.util.module_from_spec(spec)
spec.loader.exec_module(baseline)
CONFIG = json.loads((SCRIPTS.parent / 'configs/rare-biallelic-baseline.json').read_text())


def variant(pos, transcripts=None, af=0.0001, gt=None, pgt=None, pid=None, index=1, alts=None):
    gt = [0, 1] if gt is None else gt
    alt = (alts or ['G'])[index - 1]
    allele = {'chrom': '1', 'pos': pos, 'ref': 'A', 'alt': alt}
    return {'key': '1:%d:A:%s' % (pos, alt), 'allele': allele,
            'source_call': {'source_allele': allele, 'source_all_alts': alts or ['G'],
                'source_alt_index': index, 'GT_values': gt, 'PGT': pgt, 'PID': pid,
                'PS': None, 'GT_phased': False, 'alt_copies': gt.count(index), 'called_ploidy': len(gt)},
            'genes': {'GENE': {key: {'group': value} for key, value in (transcripts or {'TX': 1}).items()}},
            'frequency': {'status': 'unresolved' if af is None else 'observed', 'max_AF_fraction': af}}


class RareBiallelicBaselineTests(unittest.TestCase):
    def test_missing_frequency_is_not_zero_and_count_units_are_fractions(self):
        empty = baseline.frequency_state({'allele_found': False, 'frequencies': []})
        self.assertEqual(empty['status'], 'unresolved')
        self.assertIsNone(empty['max_AF_fraction'])
        data = {'allele_found': True, 'frequencies': [
            {'source': 'GNOMAD_E_NFE', 'ac': 1, 'an': 1000},
            {'source': 'GNOMAD_G_AFR', 'ac': 4, 'an': 2000},
            {'source': 'GNOMAD_E_OTH', 'ac': 0, 'an': 0}]}
        self.assertEqual(baseline.frequency_state(data)['max_AF_fraction'], 0.002)
        data['frequencies'][0]['ac'] = 1001
        with self.assertRaisesRegex(ValueError, 'invalid public'):
            baseline.frequency_state(data)

    def test_only_shared_transcript_effects_define_tier(self):
        a = variant(1, {'TX1': 2, 'TX2': 1})
        b = variant(10, {'TX1': 1, 'TX2': 2})
        models, _ = baseline.models_at_frequency([a, b], 0.001)
        self.assertEqual(len(models), 1)
        self.assertEqual(models[0]['tier'], 2)
        b['genes'] = {'GENE': {'UNRELATED_TX': {'group': 2}}}
        models, counts = baseline.models_at_frequency([a, b], 0.001)
        self.assertEqual(models, [])
        self.assertEqual(counts['heterozygous_pair_without_shared_transcript'], 1)

    def test_cis_exclusion_unknown_retention_and_block_validation(self):
        a = variant(1, pgt='0|1', pid='block')
        b = variant(10, pgt='0|1', pid='block')
        models, counts = baseline.models_at_frequency([a, b], 0.001)
        self.assertFalse(models)
        self.assertEqual(counts['explicit_cis_pairs_excluded'], 1)
        b['source_call']['PID'] = 'different_block'
        models, _ = baseline.models_at_frequency([a, b], 0.001)
        self.assertEqual(models[0]['phase'], 'unknown')
        b['source_call'].update(PID='block', PGT='1|0')
        self.assertEqual(baseline.pair_phase(a, b), 'trans_from_explicit_source_block')
        b['source_call']['PGT'] = '1|1'
        self.assertEqual(baseline.pair_phase(a, b), 'unknown')

    def test_diploid_multiallelic_trans_is_not_two_inferred_cis_hits(self):
        a = variant(1, gt=[1, 2], index=1, alts=['C', 'G'])
        b = variant(1, gt=[1, 2], index=2, alts=['C', 'G'])
        self.assertEqual(baseline.pair_phase(a, b), 'trans_from_diploid_multiallelic_GT')
        models, counts = baseline.models_at_frequency([a, b], 0.001)
        self.assertEqual(len(models), 1)
        self.assertEqual(models[0]['tier'], 3)
        self.assertNotIn('overlapping_pair_requiring_review', counts)

    def test_homozygous_and_unknown_frequency_remain_distinct(self):
        a = variant(1, {'TX': 2}, gt=[1, 1])
        b = variant(10, {'TX': 1}, af=None)
        c = variant(20, {'TX': 2}, af=0.0001)
        models, _ = baseline.models_at_frequency([a, b, c], 0.001)
        summary = baseline.summarize_models(models, {})
        self.assertEqual(summary['groups']['observed_rare']['best_gene_tier_counts'], {1: 1})
        self.assertEqual(summary['groups']['frequency_unresolved']['best_gene_tier_counts'], {2: 1})
        # No gene is promoted by treating absent population data as observed AF=0.
        self.assertEqual(summary['groups']['observed_rare']['models'], 1)

    def test_csq_noncoding_and_haplotype_pointer_are_not_protein_evidence(self):
        entries = ['stop_gained|GENE|TX|protein_coding|+|10Q>10*|1A>T',
                   'missense|GENE|TX|protein_coding|+|10Q>10H|1A>T',
                   'missense&NMD_transcript|GENE|TXNMD|NMD|+|10Q>10H|1A>T',
                   'intron|GENE|TX|protein_coding']
        self.assertEqual(set(baseline.effects_by_gene(entries, CONFIG)['GENE']), {'TX'})
        self.assertEqual(baseline.effects_by_gene(entries, CONFIG)['GENE']['TX']['group'], 2)
        with self.assertRaisesRegex(ValueError, 'haplotype pointer'):
            baseline.effects_by_gene(['@123'], CONFIG)


if __name__ == '__main__':
    unittest.main()
