"""ROH boundary/input safeguards with native execution where it affects inference."""
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest

import pysam

sys.path.insert(0,str(Path(__file__).parents[1]/'scripts'))
from segmented_roh import inspect_marker,run_segmented
from probe_roh_native import fixture


class SegmentedROH(unittest.TestCase):
    def setUp(self):
        self.bcftools=Path(os.environ.get('MVA_BCFTOOLS','/Users/jeweled/Documents/mva-hackathon/tools/bio/bin/bcftools'))

    def test_native_gap_restarts_and_accounts_for_every_marker(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=root/'input.vcf';fixture(source,gap=True)
            result=run_segmented(source,'SYNTH',self.bcftools,root/'output')
            self.assertEqual(result['counts']['input_markers'],3600)
            self.assertEqual(result['counts']['analyzed_markers'],3600)
            self.assertEqual(len(result['gap_boundaries']),1)
            self.assertEqual(len(result['segments']),3)
            regions=[r for segment in result['segments'] for r in segment['regions']]
            self.assertEqual(len([r for r in regions if r['chrom']=='1' and r['length_bp']>=2500000]),2)
            self.assertFalse([r for r in regions if r['chrom']=='1' and r['start_1based']<=8000000 and r['end_inclusive']>=12000000])
            self.assertFalse([r for r in regions if r['chrom']=='2'])

    def test_insufficient_segments_are_explicitly_unanalyzed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=root/'input.vcf';fixture(source)
            result=run_segmented(source,'SYNTH',self.bcftools,root/'output',min_markers=2001)
            self.assertEqual(result['status'],'no_analyzable_segments')
            self.assertEqual(result['counts']['skipped_markers'],4000)
            self.assertEqual(len(result['skipped_segments']),2)

    def test_missing_population_frequency_and_unusable_likelihoods_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=root/'input.vcf';fixture(source)
            with pysam.VariantFile(str(source)) as reader:
                row=next(reader).copy()
            inspect_marker(row,'SYNTH')
            for af in (0.0,1.0,float('nan'),None):
                bad=row.copy();bad.info['POP_AF']=af
                with self.assertRaisesRegex(ValueError,'POP_AF'):
                    inspect_marker(bad,'SYNTH')
            for pl in ((0,0,0),(1,2,3),(0,None,60)):
                bad=row.copy();bad.samples['SYNTH']['PL']=pl
                with self.assertRaisesRegex(ValueError,'PL required'):
                    inspect_marker(bad,'SYNTH')
            bad=row.copy();bad.samples['SYNTH']['GT']=(None,None)
            inspect_marker(bad,'SYNTH')  # likelihood-based inference is explicit

    def test_duplicate_marker_positions_fail_before_native_inference(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);original=root/'original.vcf';fixture(original)
            source=root/'duplicate.vcf'
            with pysam.VariantFile(str(original)) as reader,pysam.VariantFile(str(source),'w',header=reader.header) as writer:
                first=next(reader);writer.write(first);writer.write(first)
            with self.assertRaisesRegex(ValueError,'unique marker'):
                run_segmented(source,'SYNTH',self.bcftools,root/'output')


if __name__=='__main__':
    unittest.main()
