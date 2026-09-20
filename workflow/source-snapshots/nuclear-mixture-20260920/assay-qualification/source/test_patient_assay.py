"""Fictional contract controls; native integrations preserve all generated files."""
import copy,gzip,hashlib,importlib.util,json,os,shutil,subprocess,sys,tempfile,unittest
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import patient_common as c
def load(name,file):
    spec=importlib.util.spec_from_file_location(name,HERE/file);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
a=load('independent_assay','audit-patient-output.py');r=load('assay_producer','run-patient-assay.py')
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write(path,value):Path(path).write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n')
def row(pos,depth,af='0.1'):return dict(contig='fictional1',position=str(pos),ref_count=str(depth),alt_count='0',other_alt_count='0',allele_frequency=af)
class Contract(unittest.TestCase):
    def setUp(self):self.tmp=tempfile.TemporaryDirectory();self.base=Path(self.tmp.name)
    def tearDown(self):self.tmp.cleanup()
    def table(self,rows,columns=None):
        p=self.base/'table.tsv';p.write_text('#<METADATA>SAMPLE=FICTIONAL_TEST\n'+'\t'.join(columns or c.FIELDS)+'\n'+'\n'.join(rows)+'\n');return p
    def test_strict_json(self):
        for parse in (c.parse,a.parse):
            for s in ('{"a":1,"a":2}','{"x":NaN}','{"x":Infinity}'):
                with self.assertRaises(ValueError):parse(s)
    def test_typed_equality(self):
        with self.assertRaises(ValueError):a.same({'count':True},{'count':1},'type')
    def test_extra_column(self):
        p=self.table(['fictional1\t1\t12\t0\t0\t0.1\textra'])
        for f in (lambda:c.table(p,'FICTIONAL_TEST'),lambda:a.read_rows(p,c.FIELDS,'FICTIONAL_TEST')):
            with self.assertRaisesRegex(ValueError,'malformed'):f()
    def test_missing_column(self):
        p=self.table(['fictional1\t1\t12\t0\t0'])
        for f in (lambda:c.table(p,'FICTIONAL_TEST'),lambda:a.read_rows(p,c.FIELDS,'FICTIONAL_TEST')):
            with self.assertRaisesRegex(ValueError,'malformed'):f()
    def test_wrong_column_header(self):
        p=self.table([],c.FIELDS+['unexpected'])
        with self.assertRaises(ValueError):c.table(p,'FICTIONAL_TEST')
        with self.assertRaises(ValueError):a.read_rows(p,c.FIELDS,'FICTIONAL_TEST')
    def test_metadata_sample(self):
        p=self.table([])
        with self.assertRaises(ValueError):a.read_rows(p,c.FIELDS,'FICTIONAL_OTHER')
    def test_af_equal_different_literal(self):
        d,p,_,_=a.reconcile([row(1,12,'0.1')],{('fictional1',1):('0.1000',.1)});self.assertEqual(d['coverage_retained'],1);self.assertEqual((p[0]['panel_AF_literal'],p[0]['table_AF_literal']),('0.1000','0.1'))
    def test_af_changed(self):
        with self.assertRaisesRegex(ValueError,'AF binary64'):a.reconcile([row(1,12,'0.10000000000000002')],{('fictional1',1):('0.1',.1)})
    def test_af_nonfinite(self):
        with self.assertRaises(ValueError):a.reconcile([row(1,12,'nan')],{('fictional1',1):('0.1',.1)})
    def test_extra_marker(self):
        with self.assertRaises(ValueError):a.reconcile([row(2,12)],{('fictional1',1):('0.1',.1)})
    def test_duplicate_marker(self):
        with self.assertRaises(ValueError):a.reconcile([row(1,12),row(1,12)],{('fictional1',1):('0.1',.1)})
    def test_bad_counts(self):
        for count in ('-1','1.2','NaN','true'):
            with self.assertRaises(ValueError):a.reconcile([row(1,count)],{('fictional1',1):('0.1',.1)})
    def test_missing_zero_low_distinct(self):
        markers={('fictional1',i):('0.1',.1) for i in range(1,5)};d,_,m,_=a.reconcile([row(2,0),row(3,10),row(4,11)],markers)
        self.assertEqual([x['category'] for x in m],['no_emission','emitted_zero','emitted_1_to_10','coverage_retained']);self.assertEqual(d['covered_mean'],11)
    def test_strict_low_boundary(self):
        markers={('fictional1',i):('0.1',.1) for i in range(1,4)};d,_,m,_=a.reconcile([row(1,11),row(2,22),row(3,33)],markers)
        self.assertEqual(d['covered_median'],22);self.assertEqual(m[0]['category'],'coverage_excluded')
    def test_strict_high_boundary(self):
        markers={('fictional1',i):('0.1',.1) for i in range(1,5)};d,_,m,_=a.reconcile([row(i,n) for i,n in enumerate([11,11,11,99],1)],markers)
        self.assertEqual(d['covered_mean'],33);self.assertEqual(m[-1]['category'],'coverage_excluded')
    def test_covered_statistics_ignore_low(self):
        d,*_=a.reconcile([row(1,0),row(2,10),row(3,20)],{('fictional1',i):('0.1',.1) for i in range(1,4)});self.assertEqual((d['covered_mean'],d['covered_median']),(20,20))
    def test_all_missing(self):
        d,_,m,_=a.reconcile([],{('fictional1',1):('0.1',.1)});self.assertIsNone(d['covered_mean']);self.assertEqual(d['coverage_retained'],0);self.assertIsNone(m[0]['total_read_observations'])
    def test_raw_fallback(self):
        result=a.interpret([dict(sample='FICTIONAL_TEST',contamination='0',error='1')],'FICTIONAL_TEST');self.assertIsNone(result['interpreted_estimate']);self.assertEqual(result['raw_error'],1)
    def test_finite_zero_is_not_absence(self):
        result=a.interpret([dict(sample='FICTIONAL_TEST',contamination='0',error='.001')],'FICTIONAL_TEST');self.assertEqual(result['interpreted_estimate'],0);self.assertFalse(result['zero_means_absence']);self.assertIsNone(result['clinical_pass'])
    def test_invalid_model_values(self):
        for val,err in [('nan','1'),('1.1','.1'),('-.1','.1'),('0','-1'),('0','inf')]:
            with self.assertRaises(ValueError):a.interpret([dict(sample='FICTIONAL_TEST',contamination=val,error=err)],'FICTIONAL_TEST')
    def test_duplicate_model_sample(self):
        item=dict(sample='FICTIONAL_TEST',contamination='0',error='.1')
        with self.assertRaises(ValueError):a.interpret([item,item],'FICTIONAL_TEST')
    def test_input_hash(self):
        p=self.base/'pin';p.write_text('original')
        with self.assertRaises(ValueError):c.fingerprint(p,'0'*64)
    def test_false_config(self):
        p=self.base/'config';write(p,dict(fixture_only=True,ready=False))
        with self.assertRaises(ValueError):c.config(p,sha(p))
    def test_empty_panel_wrong_denominator(self):
        p=self.base/'panel.gz'
        with gzip.open(p,'wt') as f:f.write('##fileformat=VCFv4.2\n')
        with self.assertRaises(ValueError):a.marker_map(p,1,['fictional1'])

class Native(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import pysam
        cls.pysam=pysam;cls.base=Path(os.environ['MVA_ASSAY_TEST_OUTPUT']).resolve();cls.base.mkdir();cls.fixture=load('qualified_fictional_fixture','fictional_controls.py');cls.f=cls.fixture.count_fixture(cls.base/'constructor');cls.counter=0
        cls.runtime={k:dict(path=str(c.ROOT/v),sha256=c.RUNTIME[k]) for k,v in {'java':'tools/bio/lib/jvm/bin/java','gatk':'gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar'}.items()}
        for pin in cls.runtime.values():c.fingerprint(pin['path'],pin['sha256'])
        cls.source={name:sha(HERE/name) for name in c.SOURCE_NAMES}
        cls.cfgs={scale:cls.make_fixture(scale) for scale in [1,12]}
    @classmethod
    def make_fixture(cls,scale):
        p=cls.base/('scale-'+str(scale));p.mkdir();ref=p/'fictional.fa';shutil.copyfile(cls.f['reference'],ref)
        for suffix in ['.fai']:shutil.copyfile(cls.f['reference']+suffix,str(ref)+suffix)
        shutil.copyfile(Path(cls.f['reference']).with_suffix('.dict'),ref.with_suffix('.dict'))
        panel=p/'panel.vcf';lines=[]
        with gzip.open(cls.f['sites'],'rt') as source:
            for line in source:
                if line.startswith('#') or .01<float(line.rstrip().split('\t')[-1].split('=')[1])<.2:lines.append(line)
        panel.write_text(''.join(lines));cls.pysam.tabix_compress(str(panel),str(panel)+'.gz');cls.pysam.tabix_index(str(panel)+'.gz',preset='vcf')
        cram=p/'fictional.cram'
        with cls.pysam.AlignmentFile(cls.f['bam'],'rb') as original:
            h=original.header;records=[]
            for read in original:
                for i in range(scale):
                    record=copy.copy(read);record.query_name+='_'+str(i);records.append(record)
        records.sort(key=lambda x:(x.reference_id,x.reference_start,x.query_name,x.flag))
        with cls.pysam.AlignmentFile(str(cram),'wc',header=h,reference_filename=str(ref),format_options=[b'version=3.0']) as out:
            for rec in records:out.write(rec)
        cls.pysam.index(str(cram))
        with cls.pysam.AlignmentFile(str(cram),'rc',reference_filename=str(ref),format_options=[b'decode_md=0']) as ins:decoded=[cls.fixture.decoded_record(x) for x in ins]
        expected=[cls.fixture.decoded_record(x) for x in records];cls.fixture.require(decoded==expected,'fictional record preservation')
        write(p/'payload.json',dict(expected=expected,decoded=decoded,equal=True))
        paths={'alignment':cram,'alignment_index':str(cram)+'.crai','reference':ref,'reference_index':str(ref)+'.fai','reference_dictionary':ref.with_suffix('.dict'),'panel':str(panel)+'.gz','panel_index':str(panel)+'.gz.tbi'}
        cfg=dict(ready=True,fixture_only=True,patient_assay_authorized=False,fixture_root=str(p),sample=cls.fixture.SAMPLE,expected_read_groups=2,expected_panel_markers=21,contigs=['fictional1'],limits=c.LIMITS,timeouts={'pileup':3600,'estimator':600},inputs={k:dict(path=str(v),sha256=sha(v)) for k,v in paths.items()},runtime=cls.runtime,source_sha256=cls.source,output=str(p/'assay'),launch=str(p/'launch'))
        write(p/'config.json',cfg);return cfg
    def test_header_and_dictionary_valid(self):
        self.assertEqual(len(r.header_check(self.cfgs[1])['RG']),2)
    def test_wrong_sample_refused(self):
        cfg=copy.deepcopy(self.cfgs[1]);cfg['sample']='FICTIONAL_OTHER'
        with self.assertRaisesRegex(ValueError,'sample'):r.header_check(cfg)
    def test_read_group_count_refused(self):
        cfg=copy.deepcopy(self.cfgs[1]);cfg['expected_read_groups']=1
        with self.assertRaisesRegex(ValueError,'group'):r.header_check(cfg)
    def test_dictionary_mismatch_refused(self):
        cfg=copy.deepcopy(self.cfgs[1]);bad=self.base/'wrong.dict';bad.write_text('@SQ\tSN:fictional1\tLN:23999\n');cfg['inputs']['reference_dictionary']['path']=str(bad)
        with self.assertRaisesRegex(ValueError,'dictionary'):r.header_check(cfg)
    def test_missing_source_pin_refused(self):
        cfg=copy.deepcopy(self.cfgs[1]);cfg['source_sha256'].pop('audit-patient-output.py');path=self.base/'missing-source.json';write(path,cfg)
        with self.assertRaisesRegex(ValueError,'source inventory'):c.config(path,sha(path))
    def test_bad_source_pin_refused(self):
        cfg=copy.deepcopy(self.cfgs[1]);cfg['source_sha256']['audit-patient-output.py']='0'*64;path=self.base/'bad-source.json';write(path,cfg)
        with self.assertRaisesRegex(ValueError,'SHA differs'):c.config(path,sha(path))
    def test_runtime_pin_refused(self):
        cfg=copy.deepcopy(self.cfgs[1]);cfg['runtime']['gatk']['sha256']='0'*64;path=self.base/'bad-runtime.json';write(path,cfg)
        with self.assertRaisesRegex(ValueError,'runtime'):c.config(path,sha(path))
    def test_boolean_resource_refused(self):
        cfg=copy.deepcopy(self.cfgs[1]);cfg['limits']['total_wall_seconds']=True;path=self.base/'boolean-resource.json';write(path,cfg)
        with self.assertRaisesRegex(ValueError,'resource'):c.config(path,sha(path))
    def execute(self,scale):
        cfg=self.cfgs[scale];path=Path(cfg['fixture_root'])/'config.json';cmd=[sys.executable,'-B',str(HERE/'guard-patient-assay.py'),'--config',str(path),'--expected-config-sha256',sha(path)]
        with (Path(cfg['fixture_root'])/'outer.stdout.log').open('xb') as stdout,(Path(cfg['fixture_root'])/'outer.stderr.log').open('xb') as stderr:completed=subprocess.run(cmd,stdout=stdout,stderr=stderr,timeout=120)
        self.assertEqual(completed.returncode,0,(Path(cfg['fixture_root'])/'outer.stderr.log').read_text());self.assertEqual(json.loads((Path(cfg['launch'])/'status.json').read_text())['status'],'complete')
        result=json.loads((Path(cfg['output'])/'verification/verification.json').read_text());self.assertEqual(result['status'],'complete');self.assertEqual(result['denominators']['panel_sites'],21);self.assertEqual(result['denominators']['emitted'],20);self.assertEqual(result['denominators']['no_emission'],1)
        producer=json.loads((Path(cfg['output'])/'producer/summary.json').read_text());self.assertEqual(producer['estimator_executed'],scale==12)
        expected={x['position']:{k:x[k]*scale for k in ['ref_count','alt_count','other_alt_count']} for x in self.f['expected_counts']}
        rows=c.table(Path(cfg['output'])/'producer/pileups.tsv',self.fixture.SAMPLE)
        self.assertEqual({int(x['position']):{k:int(x[k]) for k in ['ref_count','alt_count','other_alt_count']} for x in rows},expected)
        for v in cfg['inputs'].values():self.assertEqual(sha(v['path']),v['sha256'])
        return result
    def test_native_covered_full_wrapper(self):
        result=self.execute(12);self.assertEqual(result['denominators']['coverage_retained'],20);self.assertIn(result['model_result']['classification'],['limited_model_output_available','insufficient_information'])
    def test_native_zero_eligible_full_wrapper(self):
        result=self.execute(1);self.assertEqual(result['denominators']['coverage_retained'],0);self.assertIsNone(result['model_result']['interpreted_estimate']);self.assertEqual(result['model_result']['classification'],'insufficient_information')

if __name__=='__main__':unittest.main(verbosity=2)
