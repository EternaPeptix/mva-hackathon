"""Fictional controls only. Real normalization/integration require pysam and bcftools."""
import copy
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import sys
import tempfile
from types import ModuleType, SimpleNamespace
import unittest

try:
    import pysam
    REAL_PYSAM = hasattr(pysam, 'VariantFile')
except ImportError:
    pysam = ModuleType('pysam')
    sys.modules['pysam'] = pysam
    REAL_PYSAM = False
SCRIPT = Path(__file__).resolve().parents[1] / 'scripts/review_native_secondary.py'
spec = importlib.util.spec_from_file_location('native_secondary_under_test', SCRIPT)
a = importlib.util.module_from_spec(spec); spec.loader.exec_module(a)
sf, native = a.load_helpers(a.Pins())
BCF = shutil.which('bcftools')
NATIVE = REAL_PYSAM and BCF is not None


def pin(path):
    path = Path(path).resolve()
    return dict(path=str(path), bytes=path.stat().st_size, sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def dump(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, sort_keys=True) + '\n')
    return pin(path)


class PureControls(unittest.TestCase):
    def test_absent_GT_stays_unknown(self):
        s = a.genotype_shape({'FORMAT': {}}, 2)
        self.assertEqual((s['state'], s['encoded_GT_slots'], s['observed_GT_ploidy']), ('absent', 0, None))

    def test_literal_dot_diploid_context_stays_unknown(self):
        s = a.genotype_shape({'FORMAT': {'GT': [None]}}, 2)
        self.assertEqual((s['state'], s['encoded_GT_slots'], s['observed_GT_ploidy']), ('fully_missing', 1, None))

    def test_diploid_missing_tuple_stays_unknown(self):
        s = a.genotype_shape({'FORMAT': {'GT': [None, None]}}, 2)
        self.assertEqual((s['state'], s['encoded_GT_slots'], s['observed_GT_ploidy']), ('fully_missing', 2, None))

    def test_haploid_dot_stays_unknown(self):
        s = a.genotype_shape({'FORMAT': {'GT': [None]}}, 1)
        self.assertIsNone(s['observed_GT_ploidy'])

    def test_partial_diploid_is_informative_shape(self):
        s = a.genotype_shape({'FORMAT': {'GT': [None, 1]}}, 2)
        self.assertEqual((s['state'], s['observed_GT_ploidy']), ('partial', 2))
        with self.assertRaises(ValueError):
            a.genotype_shape({'FORMAT': {'GT': [None, 1]}}, 1)

    def test_informative_haploid_on_autosome_refused(self):
        with self.assertRaisesRegex(ValueError, 'ploidy'):
            a.genotype_shape({'FORMAT': {'GT': [1]}}, 2)

    def test_nonfinite_duplicate_json_keys_refused(self):
        for text in ['{"a":NaN}', '{"a":Infinity}', '{"a":1,"a":2}']:
            with self.subTest(text=text), self.assertRaises(ValueError): a.parse(text)

    def test_typed_equality_does_not_merge_bool_or_float(self):
        for value in [True, 1.0]:
            with self.subTest(value=value), self.assertRaises(ValueError): a.same(value, 1, 'typed')

    def test_wrong_path_same_bytes_refused(self):
        with tempfile.TemporaryDirectory() as d:
            p, q = Path(d)/'p', Path(d)/'q';p.write_text('same');q.write_text('same')
            with self.assertRaisesRegex(ValueError, 'path'):
                a.Pins().get(p, pin(q))

    def test_stale_hash_refused(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'p';p.write_text('first');expected=pin(p);p.write_text('changed')
            with self.assertRaisesRegex(ValueError, 'SHA256'): a.Pins().get(p, expected)

    def test_duplicate_proof_path_refused(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'p';p.write_text('same')
            with self.assertRaisesRegex(ValueError, 'duplicate proof'):a.bind(p,[pin(p),pin(p)],a.Pins())

    def test_unready_config_refused_without_output(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);p=root/'draft.json';proof=dump(p,dict(schema_version=1,ready=False))
            with self.assertRaises(ValueError):a.run(p,proof['sha256'],root/'out')
            self.assertFalse((root/'out').exists())

    def test_existing_output_refused(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);out=root/'out';out.mkdir();(out/'sentinel').write_text('keep')
            with self.assertRaisesRegex(ValueError,'existing output'):a.run(root/'notread','bad',out)
            self.assertEqual((out/'sentinel').read_text(),'keep')

    def test_clinvar_strata_and_non_germline_remain_distinct(self):
        examples=[({'CLNSIG':'Pathogenic','CLNREVSTAT':'criteria_provided,_multiple_submitters,_no_conflicts'},'aggregate_plp_2plus_stars'),
            ({'CLNSIG':'Conflicting_classifications_of_pathogenicity','CLNSIGCONF':'Pathogenic_(1)|Uncertain_significance_(2)','CLNREVSTAT':'criteria_provided,_conflicting_classifications'},'conflicting_or_mixed_classification_review'),
            ({'CLNSIGINCL':'123:Pathogenic','CLNREVSTAT':'no_classification_for_the_single_variant'},'haplotype_or_genotype_only_review')]
        for info,want in examples:self.assertEqual(sf.germline_stratum(info)[0],want)
        self.assertIsNone(sf.germline_stratum({'ONC':'Oncogenic','CLNSIG':'Uncertain_significance'}))
        with self.assertRaises(ValueError):sf.germline_stratum({'CLNSIG':'Pathogenic','CLNREVSTAT':'unknown'})

    def test_literal_edit_equivalence_refuses_netlength_only(self):
        fasta=SimpleNamespace(fetch=lambda c,s,e:'ACGTACGT'[s:e])
        norm=SimpleNamespace(contig='1',pos=2,ref='C',alts=('A',))
        with self.assertRaises(ValueError):a.literal_equivalence(fasta,norm,dict(chrom='1',pos=2,ref='C',alt='G'))

    def test_duplicate_full_native_evidence_detects_QUAL_phase_and_custom_FORMAT(self):
        context=dict(source_record_id='first',source_ordinal=1,source_sha256='a'*64,caller='gatk',sample='S',
            source_namespace='fixture',chrom='1',pos=1,ref='C',alts=['A'],FORMAT={'GT':[0,1],'PS':123,'XX':'original'},
            GT_phased=True,FILTER=[],QUAL=60.0,rendered_record='literal original context')
        row=dict(caller='gatk',normalized_allele=dict(chrom='1',pos=1,ref='C',alt='A'),
            source_call=dict(source_record_id='first',source_alt_index=1,native_context=context),
            quality=dict(called_ALT_present=True,fully_called_ALT_present=True),clinvar=[{'stratum':'aggregate_plp_1star'}])
        duplicate=copy.deepcopy(row);duplicate['source_call']['source_record_id']='second'
        duplicate['source_call']['native_context'].update(source_record_id='second',source_ordinal=2)
        groups,counts=a.summarize([row,duplicate]);self.assertFalse(groups[0]['differing_native_evidence'])
        self.assertEqual(counts['gatk']['unique_called_ALT_identities'],1)
        for label in ['QUAL','PS','XX']:
            changed=copy.deepcopy(duplicate)
            if label=='QUAL':changed['source_call']['native_context']['QUAL']=61.0
            else:changed['source_call']['native_context']['FORMAT'][label]='changed'
            with self.subTest(label=label):
                groups,after=a.summarize([row,changed]);self.assertTrue(groups[0]['differing_native_evidence'])
                self.assertEqual(after,counts)


    def test_frozen_header_date_filename_and_snapshot_are_distinct_strict_fields(self):
        config=dict(expected_genes=84,public_snapshot=dict(clinvar_fileDate='2026-09-13',clinvar_snapshot_label='20260913',clinvar_resource_filename='clinvar_20260913.vcf.gz'),inputs={k:dict(path='/fictional/'+('clinvar_20260913.vcf.gz' if k=='clinvar' else k),sha256=v) for k,v in a.FROZEN_PUBLIC.items()})
        a.frozen_public_config(config)
        for field,value in [('clinvar_fileDate','20260913'),('clinvar_fileDate','2026-09-14'),('clinvar_snapshot_label','2026-09-13'),('clinvar_resource_filename','clinvar_latest.vcf.gz')]:
            changed=copy.deepcopy(config);changed['public_snapshot'][field]=value
            with self.subTest(field=field,value=value),self.assertRaises(ValueError):a.frozen_public_config(changed)
        changed=copy.deepcopy(config);changed['inputs']['clinvar']['sha256']='0'*64
        with self.assertRaises(ValueError):a.frozen_public_config(changed)


def write_vcf(path, records, chroms, include_gq=True, public=False, public_clnsig_number='.', public_date='synthetic', public_clnalle_header=True):
    if public:
        h='##fileformat=VCFv4.2\n##reference=GRCh38\n##fileDate='+public_date+'\n'
        h+=''.join('##INFO=<ID='+k+',Number=.,Type=String,Description="fictional">\n' for k in ['GENEINFO','CLNSIG','CLNSIGCONF','CLNSIGINCL','CLNREVSTAT','CLNALLE','EXTRA'] if public_clnalle_header or k!='CLNALLE')
        h=h.replace('ID=CLNSIG,Number=.', 'ID=CLNSIG,Number='+public_clnsig_number)
        h+='#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    else:
        h='##fileformat=VCFv4.2\n'+''.join('##contig=<ID='+c+',length=100>\n' for c in chroms)
        for k,number,typ in [('GT','1','String'),('DP','1','Integer'),('AD','R','Integer'),('PL','G','Integer'),('XX','1','String')]+([('GQ','1','Integer')] if include_gq else []):
            h+='##FORMAT=<ID='+k+',Number='+number+',Type='+typ+',Description="fictional">\n'
        h+='##FILTER=<ID=Low,Description="fictional">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n'
    raw=path.with_suffix('');raw.write_text(h+''.join(r+'\n' for r in records))
    pysam.tabix_compress(str(raw),str(path),force=True);pysam.tabix_index(str(path),preset='vcf',force=True)
    return pin(path)


def fixture(root):
    sequence='CAAAAAGTAC'+'GTACGTACGT'*9
    y=list(sequence);y[2:5]='NNN';y[90:94]='NNNN';y=''.join(y)
    ref=root/'reference.fa';ref.write_text(''.join('>'+c+'\n'+(y if c=='Y' else sequence)+'\n' for c in a.AUTOSOMES+['X','Y']));pysam.faidx(str(ref))
    panel=root/'panel.json';dump(panel,dict(total=1,rows=[dict(symbol='G1')]))
    public=root/'clinvar.vcf.gz'
    records=['1\t1\t101\tC\tA\t.\t.\tGENEINFO=G1:1;CLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts;EXTRA=full_original',
        '1\t1\t102\tC\tA\t.\t.\tGENEINFO=G1:1;CLNSIG=Likely_pathogenic;CLNREVSTAT=criteria_provided,_single_submitter',
        '1\t1\t103\tCA\tC\t.\t.\tGENEINFO=G1:1;CLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_single_submitter',
        '1\t7\t104\tG\tA,T\t.\t.\tGENEINFO=G1:1;CLNSIG=Pathogenic;CLNALLE=1;CLNREVSTAT=criteria_provided,_single_submitter']
    records += ['1\t'+str(pos)+'\t'+str(100+pos)+'\t'+refbase+'\t'+alt+'\t.\t.\tGENEINFO=G1:1;CLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_single_submitter' for pos,refbase,alt in [(8,'T','C'),(9,'A','G'),(10,'C','T'),(11,'G','A'),(12,'T','A')]]
    records += [c+'\t25\t200\tG\tA\t.\t.\tGENEINFO=G1:1;CLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_single_submitter' for c in ['X','Y']]
    write_vcf(public,records,[],public=True)
    raw_records=['1\t1\t.\tC\tA\t60\t.\t.\tGT:DP:AD:PL:XX\t0/1:20:10,10:99,0,99:custom',
        '1\t1\t.\tC\tA\t60\t.\t.\tGT:DP:AD:PL:XX\t0/1:20:10,10:99,0,99:custom',
        '1\t4\t.\tAA\tA\t60\tPASS\t.\tGT:DP:AD:PL:XX\t0/1:20:10,10:99,0,99:custom',
        '1\t7\t.\tG\tA,T\t60\tLow\t.\tGT:DP:AD:PL:XX\t1/2:20:0,10,10:99,99,99,99,0,99:custom',
        '1\t8\t.\tT\tC\t60\t.\t.\tGT:DP\t./1:20',
        '1\t9\t.\tA\tG\t60\t.\t.\tGT:DP\t.:20',
        '1\t10\t.\tC\tT\t60\t.\t.\tGT:DP\t./.:20',
        '1\t11\t.\tG\tA\t60\t.\t.\tDP\t20',
        '1\t12\t.\tT\tA\t60\t.\t.\tGT:DP:GQ\t0/0:0:0']
    handoff=dict(status='complete',sample='S',selected_autosomes=a.AUTOSOMES,pending_autosomes=[],all_autosomes_complete=True,inputs={},shards=[])
    manifest=[]
    for index,c in enumerate(a.AUTOSOMES):
        shard=dict(shard_id=index,contigs=[[c,100]],methods={})
        for caller in a.CALLERS:
            path=root/(caller+'-'+c+'.vcf.gz');p=write_vcf(path,raw_records if c=='1' else [],[c]);manifest.append(p)
            with pysam.VariantFile(str(path)) as reader:rows=[str(r) for r in reader]
            shard['methods'][caller]=dict(replayed_split_payload_identical=True,original_raw_pin=p,
                original_raw=dict(counts={'records':len(rows)},rendered_record_payload_sha256=hashlib.sha256(''.join(rows).encode()).hexdigest()))
        handoff['shards'].append(shard)
    xy=root/'fixed-xy';xy.mkdir()
    outputs={};xy_counts={}
    for caller in a.CALLERS:
        path=xy/(caller+'.original.vcf.gz')
        write_vcf(path,['X\t25\t.\tG\tA\t60\t.\t.\tGT:DP:AD:PL\t1:20:0,20:99,0','Y\t25\t.\tG\tA\t60\t.\t.\tGT:DP\t.:20'],['X','Y'])
        for f in [path,Path(str(path)+'.tbi')]:outputs[f.name]=pin(f)['sha256']
        xy_counts['fixed_gatk' if caller=='gatk' else 'bcftools']={'records':2}
    model=dict(sample='S',reference_sha256=pin(ref)['sha256'],regions=[dict(chrom='X',length=100,intervals_0based_half_open=[[0,10,1],[10,20,2],[20,80,1],[80,90,2],[90,100,1]]),dict(chrom='Y',length=100,intervals_0based_half_open=[[0,100,1]])],required_N_masks=[dict(chrom='Y',start=2,end=5),dict(chrom='Y',start=90,end=94)])
    statuspin=dump(xy/'status.json',dict(status='complete',config=model,output_sha256=outputs))
    reviewpin=dump(root/'xy-review/summary.json',dict(status='complete',verified_at='2026-09-20T00:00:00+00:00',script_sha256='a'*64,fixed_outputs_verified=len(outputs),original_files_preserved=1,bcftools_original_reused_byte_for_byte=True,fixed_run_status_sha256=statuspin['sha256'],raw_record_accounting=xy_counts))
    proofpin=dump(root/'xy-review/verification.json',dict(verified_at='2026-09-20T00:00:01+00:00',script_sha256='a'*64,fixed_outputs=outputs,files=[reviewpin]))
    inputs=dict(reference=pin(ref),reference_index=pin(str(ref)+'.fai'),bcftools=pin(BCF),clinvar=pin(public),clinvar_index=pin(str(public)+'.tbi'),panel=pin(panel),fixed_xy_status=statuspin,fixed_xy_review=reviewpin,fixed_xy_verification=proofpin)
    handoff['inputs']={k:inputs[k] for k in ['reference','reference_index','bcftools']}
    inputs['handoff']=dump(root/'handoff.json',handoff);manifest.append(inputs['handoff'])
    inputs['autosome_verification']=dump(root/'native-verification.json',dict(status='complete',scope={k:handoff[k] for k in ['selected_autosomes','pending_autosomes','all_autosomes_complete']},files=manifest))
    cfg=dict(schema_version=1,ready=True,fixture_only=True,root=str(root),sample='S',expected_genes=1,quality_policy='descriptive_native_axes_no_eligibility_gate',public_snapshot={'clinvar_fileDate':'synthetic','panel':'fictional1gene'},inputs=inputs)
    config=root/'config.json';dump(config,cfg)
    return cfg,config


@unittest.skipUnless(NATIVE,'requires pysam and bcftools; full native qualification must have no skips')
class NativeControls(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name).resolve();self.cfg,self.config=fixture(self.root)
    def tearDown(self):self.tmp.cleanup()
    def sources(self):
        with pysam.FastaFile(self.cfg['inputs']['reference']['path']) as fasta:
            return a.prepare_sources(self.cfg,self.cfg['inputs'],a.Pins(),fasta)
    def alter_json(self,label,change):
        p=Path(self.cfg['inputs'][label]['path']);d=json.loads(p.read_text());change(d);self.cfg['inputs'][label]=dump(p,d)
    def test_complete_adapter_and_source_preservation(self):
        before={str(p):pin(p)['sha256'] for p in self.root.rglob('*') if p.is_file()}
        out=self.root/'out';result=a.run(self.config,pin(self.config)['sha256'],out)
        self.assertEqual(result['source_files'],46)
        for c in a.CALLERS:
            count=result['counts_by_caller'][c]
            self.assertEqual(count['exact_matching_occurrences'],12)
            self.assertEqual(count['unique_exact_matching_identities'],11)
            self.assertEqual(count['unique_called_ALT_identities'],6)
            self.assertEqual(count['unique_fully_called_ALT_identities'],5)
            self.assertEqual(count['unique_identities_without_called_ALT'],5)
            self.assertEqual(count['duplicate_normalized_identities'],1)
        with gzip.open(out/'exact-matching-occurrences.jsonl.gz','rt') as f:rows=[json.loads(s) for s in f]
        ordinary=[r for r in rows if r['caller']=='gatk' and r['normalized_allele']==dict(chrom='1',pos=1,ref='C',alt='A')]
        self.assertEqual(len(ordinary),2);self.assertNotEqual(ordinary[0]['source_call']['source_record_id'],ordinary[1]['source_call']['source_record_id'])
        self.assertEqual(ordinary[0]['source_call']['native_context']['FORMAT']['XX'],'custom')
        self.assertEqual(ordinary[0]['source_call']['PL'],[99,0,99])
        self.assertEqual(ordinary[0]['quality']['filter_state'],'unfiltered')
        self.assertEqual(ordinary[0]['quality']['GQ']['FORMAT_presence'],'absent')
        self.assertEqual(len(ordinary[0]['clinvar']),2)
        self.assertEqual(ordinary[0]['clinvar'][0]['source']['INFO']['EXTRA'],['full_original'])
        self.assertIn('EXTRA=full_original',ordinary[0]['clinvar'][0]['source']['rendered_record'])
        multi=[r for r in rows if r['caller']=='gatk' and r['source_call']['source_allele']['pos']==7]
        self.assertEqual([r['source_call']['source_alt_index'] for r in multi],[1,2])
        self.assertTrue(all(r['source_call']['GT_values']==[1,2] for r in multi))
        self.assertTrue(all(c['stratum']=='multiallelic_assertion_mapping_review' for r in multi for c in r['clinvar']))
        self.assertTrue(all(c['assertion_mapping']['original_CLNALLE']==['1'] for r in multi for c in r['clinvar']))
        self.assertEqual(result['public_counts']['selected_multiallelic_source_records'],1)
        self.assertTrue(all(r['quality']['genotype']['observed_GT_ploidy'] is None for r in rows if r['normalized_allele']['pos'] in [9,10,11]))
        y=[r for r in rows if r['normalized_allele']['chrom']=='Y'];self.assertTrue(all(r['quality']['genotype']['observed_GT_ploidy'] is None for r in y))
        self.assertEqual(before,{name:pin(name)['sha256'] for name in before})
    def test_allele_specific_public_assertion_retains_benign_split_as_ambiguous(self):
        path=Path(self.cfg['inputs']['clinvar']['path'])
        with pysam.VariantFile(str(path)) as reader:
            records=[str(r).rstrip('\n') for r in reader]
        records=[r.replace('CLNSIG=Pathogenic;', 'CLNSIG=Pathogenic,Benign;') if r.split('\t')[2]=='104' else r for r in records]
        self.cfg['inputs']['clinvar']=write_vcf(path,records,[],public=True,public_clnsig_number='A')
        self.cfg['inputs']['clinvar_index']=pin(str(path)+'.tbi');dump(self.config,self.cfg)
        out=self.root/'allele-specific-public';a.run(self.config,pin(self.config)['sha256'],out)
        with gzip.open(out/'exact-matching-occurrences.jsonl.gz','rt') as stream:
            rows=[json.loads(line) for line in stream]
        contexts=[c for r in rows if r['caller']=='gatk' and r['normalized_allele']['pos']==7 for c in r['clinvar']]
        self.assertEqual(len(contexts),2)
        self.assertEqual([c['source_alt_index'] for c in contexts],[1,2])
        self.assertTrue(all(c['stratum']=='multiallelic_assertion_mapping_review' for c in contexts))
        self.assertTrue(all(c['source_record_stratum']=='conflicting_or_mixed_classification_review' for c in contexts))
        self.assertTrue(all(c['source']['INFO']['CLNSIG']==['Pathogenic','Benign'] for c in contexts))
        self.assertEqual([c['classification'] for c in contexts],['Pathogenic','Benign'])
        self.assertEqual([c['normalized_record_stratum'] for c in contexts],['aggregate_plp_1star',None])
        self.assertEqual([c['normalized_record_stars'] for c in contexts],[1,None])

    def test_actual_literal_header_date_and_wrong_compact_format(self):
        path=Path(self.cfg['inputs']['clinvar']['path'])
        with pysam.VariantFile(str(path)) as reader:records=[str(r).rstrip('\n') for r in reader]
        self.cfg['inputs']['clinvar']=write_vcf(path,records,[],public=True,public_date='2026-09-13')
        self.cfg['inputs']['clinvar_index']=pin(str(path)+'.tbi')
        self.cfg['public_snapshot']['clinvar_fileDate']='2026-09-13';dump(self.config,self.cfg)
        result=a.run(self.config,pin(self.config)['sha256'],self.root/'literal-date')
        self.assertEqual(result['public_snapshot']['clinvar_fileDate'],'2026-09-13')
        self.cfg['public_snapshot']['clinvar_fileDate']='20260913';dump(self.config,self.cfg)
        with self.assertRaisesRegex(ValueError,'dated ClinVar header'):
            a.run(self.config,pin(self.config)['sha256'],self.root/'wrong-date')

    def test_public_CLNALLE_undeclared_in_header_remains_unknown(self):
        path=Path(self.cfg['inputs']['clinvar']['path'])
        with pysam.VariantFile(str(path)) as reader:records=[str(r).rstrip('\n').replace(';CLNALLE=1','') for r in reader]
        self.cfg['inputs']['clinvar']=write_vcf(path,records,[],public=True,public_clnalle_header=False)
        self.cfg['inputs']['clinvar_index']=pin(str(path)+'.tbi');dump(self.config,self.cfg)
        with pysam.VariantFile(str(path)) as reader:self.assertNotIn('CLNALLE',reader.header.info)
        out=self.root/'undeclared-public-optional';a.run(self.config,pin(self.config)['sha256'],out)
        with gzip.open(out/'exact-matching-occurrences.jsonl.gz','rt') as stream:rows=[json.loads(line) for line in stream]
        contexts=[c for row in rows for c in row['clinvar']]
        self.assertTrue(contexts)
        self.assertTrue(all(c['assertion_mapping']['normalized_CLNALLE'] is None and c['assertion_mapping']['original_CLNALLE'] is None for c in contexts))
        multi=[c for c in contexts if c['public_original_multiallelic']]
        self.assertTrue(multi);self.assertTrue(all(c['stratum']=='multiallelic_assertion_mapping_review' for c in multi))

    def test_complete_zero_record_chr22_valid(self):
        sources,_=self.sources();self.assertEqual(len(sources),46)
        self.assertTrue(all(s['expected_census']['counts']['records']==0 for s in sources if s['contigs']==['22']))
    def test_missing_chr22_source_refused_even_full_header(self):
        self.alter_json('handoff',lambda d:d['shards'].pop())
        with self.assertRaises(ValueError):self.sources()
    def test_wrong_same_byte_handoff_path_refused(self):
        original=Path(self.cfg['inputs']['handoff']['path']);other=self.root/'other-handoff.json';other.write_bytes(original.read_bytes());self.cfg['inputs']['handoff']=pin(other)
        with self.assertRaisesRegex(ValueError,'actual path absent'):self.sources()
    def test_stale_source_hash_refused(self):
        path=self.root/'gatk-22.vcf.gz';path.write_bytes(path.read_bytes()+b'changed')
        with self.assertRaisesRegex(ValueError,'SHA256'):self.sources()
    def test_fixed_status_link_refused(self):
        self.alter_json('fixed_xy_review',lambda d:d.update(fixed_run_status_sha256='0'*64))
        with self.assertRaisesRegex(ValueError,'binding'):self.sources()
    def test_ploidy_gap_and_mask_refused(self):
        with pysam.FastaFile(self.cfg['inputs']['reference']['path']) as fasta:
            model=json.loads(Path(self.cfg['inputs']['fixed_xy_status']['path']).read_text())['config']
            bad=copy.deepcopy(model);bad['regions'][0]['intervals_0based_half_open'][1][0]=11
            with self.assertRaises(ValueError):a.validate_regions(bad,fasta)
            bad=copy.deepcopy(model);bad['required_N_masks'][0]['start']=1
            with self.assertRaisesRegex(ValueError,'mask'):a.validate_regions(bad,fasta)
    def test_xy_crossing_native_span_refused(self):
        _,regions=self.sources();source=dict(kind='fixed_xy',contigs=['X','Y'])
        record=SimpleNamespace(contig='X',start=8,stop=12)
        with self.assertRaisesRegex(ValueError,'crosses'):a.expected_ploidy(record,source,regions)
    def test_native_unobserved_match_is_not_called_ALT(self):
        sources,regions=self.sources();source=next(s for s in sources if s['id']=='gatk-1')
        with pysam.FastaFile(self.cfg['inputs']['reference']['path']) as fasta:
            rows,_=a.scan_source(source,'S',regions,{('1',12,'T','A'):[{'stratum':'aggregate_plp_1star'}]},self.cfg['inputs'],self.root/'one',sf,native,fasta)
        groups,counts=a.summarize(rows)
        self.assertEqual(counts['gatk']['unique_exact_matching_identities'],1)
        self.assertEqual(counts['gatk']['unique_called_ALT_identities'],0)
        self.assertEqual(rows[0]['quality']['GQ']['state'],'explicit_zero')
        self.assertFalse(groups[0]['called_ALT_present'])
    def test_fixed_verification_without_success_evidence_refused(self):
        self.alter_json('fixed_xy_verification',lambda d:d.pop('verified_at'))
        with self.assertRaises(KeyError):self.sources()
    def test_fixed_verification_wrong_source_digest_refused(self):
        self.alter_json('fixed_xy_verification',lambda d:d.update(script_sha256='b'*64))
        with self.assertRaisesRegex(ValueError,'source identity'):self.sources()
    def test_fixed_unverified_output_count_refused(self):
        self.alter_json('fixed_xy_review',lambda d:d.update(fixed_outputs_verified=0))
        with self.assertRaisesRegex(ValueError,'outputs verified'):self.sources()
    def custom_scan(self,records,keys):
        path=self.root/'custom.vcf.gz';write_vcf(path,records,['1'])
        with pysam.VariantFile(str(path)) as reader:rendered=[str(r) for r in reader]
        source=dict(caller='gatk',kind='autosome',contigs=['1'],pin=pin(path),id='custom',expected_census=dict(counts={'records':len(rendered)},rendered_record_payload_sha256=hashlib.sha256(''.join(rendered).encode()).hexdigest()))
        _,regions=self.sources()
        with pysam.FastaFile(self.cfg['inputs']['reference']['path']) as fasta:
            return a.scan_source(source,'S',regions,{k:[{'stratum':'aggregate_plp_1star'}] for k in keys},self.cfg['inputs'],self.root/'custom-out',sf,native,fasta)[0]
    def test_concrete_match_preserves_symbolic_other_ALT_and_full_PL(self):
        rows=self.custom_scan(['1\t13\t.\tA\tC,*\t60\t.\t.\tGT:DP:AD:PL:XX\t1/2:20:0,10,10:99,99,99,99,0,99:with_symbolic'],[('1',13,'A','C')])
        self.assertEqual(len(rows),1)
        call=rows[0]['source_call'];self.assertEqual(call['source_all_alts'],['C','*'])
        self.assertEqual(call['GT_values'],[1,2]);self.assertEqual(call['PL'],[99,99,99,99,0,99])
        self.assertEqual(call['native_context']['FORMAT']['XX'],'with_symbolic')
        self.assertTrue(rows[0]['quality']['nonconcrete_in_GT'])
        self.assertTrue(rows[0]['quality']['called_ALT_present'])
    def test_present_missing_GQ_remains_distinct_from_absence_and_zero(self):
        rows=self.custom_scan(['1\t14\t.\tC\tG\t60\t.\t.\tGT:DP:GQ\t0/1:20:.'],[('1',14,'C','G')])
        self.assertEqual(rows[0]['quality']['GQ']['FORMAT_presence'],'present_missing')
        self.assertIsNone(rows[0]['quality']['GQ']['value'])
        roundtrip=json.loads(json.dumps(rows))
        self.assertEqual(roundtrip[0]['quality']['GQ'],rows[0]['quality']['GQ'])
    def test_nonmatching_scope_has_explicit_zero_counts(self):
        groups,counts=a.summarize([]);self.assertEqual(groups,[])
        self.assertEqual(counts['gatk']['unique_exact_matching_identities'],0)
        self.assertEqual(counts['gatk']['unique_called_ALT_identities'],0)


if __name__=='__main__':unittest.main(verbosity=2)
