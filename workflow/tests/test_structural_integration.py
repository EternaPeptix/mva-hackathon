"""Evidence-preservation and provenance failures in the guarded structural runner."""
import copy
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import pysam

sys.path.insert(0, str(Path(__file__).parents[1]/'scripts'))
from build_structural_exclusions import fingerprint, parse_exclusions, save
from probe_structural_integration import fixture_config
from run_structural_variants import (PATIENT_SCOPE, alignment_census, load_mask,
                                     review_calls, run_analysis, validate_config, verify_completion)


class StructuralIntegration(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name)
        self.project=Path(os.environ.get('MVA_STUDIO_ROOT','/Users/jeweled/Documents/mva-hackathon'))
        self.experiment=self.project/'results/structural-synthetic-20260919-v2'
        self.reference=self.experiment/'fixture/reference.fa'
        self.delly=self.project/'tools/delly-2.6.0/bin/delly';self.sam=self.project/'tools/bio/bin/samtools'
        self.config=fixture_config(self.experiment,self.root/'mask','control',self.delly,self.sam)

    def tearDown(self):self.temp.cleanup()

    def run_fixture(self, config=None, alignment=None, output='analysis'):
        return run_analysis(alignment or self.experiment/'aligned.bam',self.reference,self.root/'mask',
                            config or self.config,self.delly,self.sam,self.root/output)

    def test_wrong_pins_and_sample_fail_before_discovery(self):
        for field,bad in [('reference_sha256','0'*64),('reference_index_sha256','0'*64),
                          ('delly_sha256','0'*64),('samtools_sha256','0'*64),('exclusion_summary_sha256','0'*64)]:
            with self.subTest(field=field):
                config=copy.deepcopy(self.config);config[field]=bad
                with self.assertRaises(ValueError):self.run_fixture(config,output=field)
                self.assertFalse((self.root/field/'calls.bcf').exists())
                self.assertEqual(json.loads((self.root/field/'status.json').read_text())['status'],'failed')
        config=copy.deepcopy(self.config);config['sample']='WRONG'
        with self.assertRaisesRegex(ValueError,'reserved synthetic'):self.run_fixture(config,output='sample')

    def test_interval_sidecar_cannot_be_forged_by_only_updating_its_digest(self):
        folder=self.root/'mask';path=folder/'query-intervals.bed'
        path.write_text(path.read_text().replace('1\t0\t1000000','1\t0\t999999'))
        summary=json.loads((folder/'summary.json').read_text());summary['outputs'][path.name]=fingerprint(path);save(folder/'summary.json',summary)
        config=copy.deepcopy(self.config);config['exclusion_summary_sha256']=fingerprint(folder/'summary.json')['sha256']
        with self.assertRaisesRegex(ValueError,'native parser semantics'):self.run_fixture(config)
        self.assertFalse((self.root/'analysis/calls.bcf').exists())

    def test_unknown_contigs_and_out_of_bounds_exclusions_are_errors(self):
        for text in ['chr1\t0\t10\n','1\t-1\t10\n','1\t0\t1001\n','1\t5\t5\n','1\t0\t10\textra\n']:
            with self.subTest(text=text):
                p=self.root/'bad.tsv';p.write_text(text)
                with self.assertRaises(ValueError):parse_exclusions(p,{'1':1000})

    def test_patient_requires_all_primary_nuclear_chromosomes_and_completed_cram(self):
        config=json.loads(Path(__file__).parents[1].joinpath('configs/structural-discovery.json').read_text())
        narrowed=copy.deepcopy(config);narrowed['chromosomes']=list(map(str,range(1,23)))
        with self.assertRaisesRegex(ValueError,'all 24'):validate_config(narrowed)
        alignment=self.root/'full.cram';identity={'bytes':123,'sha256':'1'*64};receipt=self.root/'completion.json'
        with self.assertRaisesRegex(ValueError,'receipt required'):verify_completion(config,alignment,identity,None)
        valid={'status':'complete','expected_lanes':4,'cram_path':str(alignment),'cram_bytes':123,'cram_sha256':'1'*64}
        save(receipt,valid);verify_completion(config,alignment,identity,receipt)
        for field,bad in [('status','running'),('expected_lanes',2),('cram_path',str(self.root/'old.cram')),
                          ('cram_bytes',122),('cram_sha256','2'*64)]:
            with self.subTest(field=field):
                value=dict(valid);value[field]=bad;save(receipt,value)
                with self.assertRaisesRegex(ValueError,'does not match'):verify_completion(config,alignment,identity,receipt)
        save(receipt,{'status':'running','expected_lanes':4})
        with self.assertRaisesRegex(ValueError,'not complete'):
            run_analysis(self.root/'nonexistent.cram',self.reference,self.root/'mask',config,
                         self.delly,self.sam,self.root/'partial',receipt)
        self.assertFalse((self.root/'partial/calls.bcf').exists())

    def test_patient_output_inside_repository_is_rejected_before_writing(self):
        config=json.loads(Path(__file__).parents[1].joinpath('configs/structural-discovery.json').read_text())
        target=Path(__file__).parents[2]/'results/forbidden-patient-test-output'
        self.assertFalse(target.exists())
        with self.assertRaisesRegex(ValueError,'outside Git'):
            run_analysis(self.root/'nonexistent.cram',self.reference,self.root/'mask',config,self.delly,self.sam,target)
        self.assertFalse(target.exists())

    def quality_fixture(self):
        with pysam.FastaFile(str(self.reference)) as genome:
            dictionary=dict(zip(genome.references,genome.lengths))
            header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':c,'LN':n} for c,n in dictionary.items()],
                    'RG':[{'ID':'SYNTH.1','SM':'SYNTH_SV','PL':'ILLUMINA'}]}
            source=self.root/'quality-source.bam'
            with pysam.AlignmentFile(str(source),'wb',header=header) as writer:
                for i,(flag,mq) in enumerate([(0,255),(0x100,60),(0x800,60)]):
                    read=pysam.AlignedSegment();read.query_name='SYNTH_QUALITY_'+str(i);read.flag=flag
                    read.reference_id=0;read.reference_start=1000+100*i;read.mapping_quality=mq;read.cigar=[(0,50)]
                    read.query_sequence=genome.fetch('1',read.reference_start,read.reference_start+50)
                    read.query_qualities=pysam.qualitystring_to_array('I'*50);read.set_tag('RG','SYNTH.1');writer.write(read)
        path=self.root/'quality.bam';subprocess.run([str(self.sam),'markdup',str(source),str(path)],check=True,capture_output=True)
        pysam.index(str(path));return path

    def test_unknown_mapping_quality_blocks_discovery_and_supplementary_records_remain(self):
        path=self.quality_fixture()
        with self.assertRaisesRegex(ValueError,'unknown MQ 255'):self.run_fixture(alignment=path)
        census=json.loads((self.root/'analysis/alignment-census.json').read_text())
        self.assertEqual(census['counts']['alignment_records'],3)
        self.assertEqual(census['counts']['mapped_unknown_MQ_255'],1)
        self.assertEqual(census['counts']['secondary_records'],1)
        self.assertEqual(census['counts']['supplementary_records'],1)
        self.assertFalse((self.root/'analysis/calls.bcf').exists())

    def test_missing_read_group_and_unsorted_records_are_rejected(self):
        path=self.quality_fixture()
        with pysam.AlignmentFile(str(path)) as source:reads=list(source)
        class Stream(list):references=('1','2','3')
        with self.assertRaisesRegex(ValueError,'coordinate sorted'):alignment_census(Stream(list(reversed(reads))),{'SYNTH.1'})
        reads[0].set_tag('RG',None)
        with self.assertRaisesRegex(ValueError,'undeclared read group'):alignment_census(Stream(reads),{'SYNTH.1'})

    def test_wrong_alignment_sample_and_missing_duplicate_provenance_fail_early(self):
        original=self.quality_fixture()
        with pysam.AlignmentFile(str(original)) as source:
            header=source.header.to_dict();reads=list(source)
        for mode in ('wrong_sample','missing_markdup'):
            h=copy.deepcopy(header)
            if mode=='wrong_sample':h['RG'][0]['SM']='WRONG'
            else:h.pop('PG',None)
            path=self.root/(mode+'.bam')
            with pysam.AlignmentFile(str(path),'wb',header=h) as writer:
                for read in reads:writer.write(read)
            pysam.index(str(path))
            expected='sample/read-group' if mode=='wrong_sample' else 'duplicate-marking provenance'
            with self.assertRaisesRegex(ValueError,expected):self.run_fixture(alignment=path,output=mode)
            self.assertFalse((self.root/mode/'alignment-census.json').exists())
            self.assertFalse((self.root/mode/'calls.bcf').exists())

    def test_sex_genotype_and_endpoint_outside_query_are_retained_for_review(self):
        h=pysam.VariantHeader();h.contigs.add('X',length=1000)
        h.info.add('SVTYPE',1,'String','Synthetic type');h.formats.add('GT',1,'String','Genotype');h.add_sample('SYNTH_SV')
        path=self.root/'calls.vcf'
        with pysam.VariantFile(str(path),'w',header=h) as writer:
            row=writer.new_record(contig='X',start=149,stop=250,id='SYNTH_DEL',alleles=('A','<DEL>'))
            row.info['SVTYPE']='DEL';row.samples['SYNTH_SV']['GT']=(0,1);row.samples['SYNTH_SV'].phased=True;writer.write(row)
        output=self.root/'review.jsonl';result=review_calls(path,output,'SYNTH_SV',{'X':1000},['X'],{'X':[[100,200]]},set())
        row=json.loads(output.read_text())
        self.assertTrue(row['genotype_phased']);self.assertTrue(row['native_diploid_sex_genotype_uninterpreted'])
        self.assertEqual(row['endpoint_reference_bases_in_query_intervals'],[False,True])
        self.assertEqual(result['counts']['records_with_endpoint_base_outside_query_intervals'],1)
        self.assertEqual(row['filter'],[])

    def test_empty_native_header_is_explicit_but_missing_or_wrong_samples_with_calls_fail(self):
        h=pysam.VariantHeader();h.contigs.add('1',length=1000);h.info.add('SVTYPE',1,'String','Synthetic type')
        path=self.root/'empty.vcf'
        with pysam.VariantFile(str(path),'w',header=h):pass
        result=review_calls(path,self.root/'empty.jsonl','SYNTH_SV',{'1':1000},['1'],{'1':[]},set())
        self.assertEqual(result['counts'],{});self.assertEqual(result['native_header_samples'],[])
        self.assertEqual(result['native_sample_column_status'],'absent_empty_callset')
        path=self.root/'sites-only.vcf'
        with pysam.VariantFile(str(path),'w',header=h) as writer:
            row=writer.new_record(contig='1',start=99,stop=200,id='SYNTH',alleles=('A','<DEL>'));row.info['SVTYPE']='DEL';writer.write(row)
        with self.assertRaisesRegex(ValueError,'nonempty native callset has no sample'):
            review_calls(path,self.root/'invalid.jsonl','SYNTH_SV',{'1':1000},['1'],{'1':[]},set())
        h.add_sample('WRONG');path=self.root/'wrong.vcf'
        with pysam.VariantFile(str(path),'w',header=h):pass
        with self.assertRaisesRegex(ValueError,'sample identity mismatch'):
            review_calls(path,self.root/'wrong.jsonl','SYNTH_SV',{'1':1000},['1'],{'1':[]},set())


if __name__=='__main__':unittest.main()
