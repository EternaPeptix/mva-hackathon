"""Fictional public transformation and independent reconciliation controls."""
import copy,gzip,hashlib,importlib.util,json,os,tempfile,unittest
from contextlib import nullcontext
from pathlib import Path
from unittest.mock import patch

HERE=Path(__file__).resolve().parent
def module(name):
    spec=importlib.util.spec_from_file_location(name.replace('-','_'),HERE/name);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
v=module('verify-public-panel.py');g=module('guard-public-panel.py')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
class Controls(unittest.TestCase):
    def fixture(self):
        source=[];ledger=[];derived=[]
        for i,literal in enumerate(['0.001','0.01','0.0100001','0.1000','0.199999','0.2','0.21','1']):
            line=f'1\t{i+1}\trs{i}\tA\tC\t42\tPASS\tPOP_AF={literal};ORIGINAL=yes\n';source.append(line)
            reason=['below_lower_endpoint','at_lower_endpoint','included','included','included','at_upper_endpoint','above_upper_endpoint','above_upper_endpoint'][i]
            digest=hashlib.sha256(line.encode()).hexdigest();ledger.append(json.dumps(dict(source_sha256='a'*64,source_ordinal=i+1,source_record_sha256=digest,chrom='1',pos=i+1,ref='A',alt='C',population_ALT_frequency_literal=literal,selection=reason)))
            if reason=='included':derived.append(line.rstrip('\n')+f';AF={literal};NMP_ORD={i+1};NMP_SHA256={digest}\n')
        return source,ledger,derived
    def compare(self,data):return v.compare_records(*data,lambda c,a,b:'A',8,'a'*64)
    def test_exact_original_fields_literals_and_strict_endpoints(self):
        r=self.compare(self.fixture());self.assertEqual(r['counts'],dict(source_records=8,below_lower_endpoint=1,at_lower_endpoint=1,included=3,at_upper_endpoint=1,above_upper_endpoint=2));self.assertEqual(r['included_by_chromosome'],{'1':3})
    def test_selected_omission_extra_and_reorder_refused(self):
        for kind in ['omit','extra','reorder']:
            s,l,d=self.fixture()
            if kind=='omit':d.pop()
            elif kind=='extra':d.append(d[-1])
            else:d.reverse()
            with self.subTest(kind=kind),self.assertRaises(ValueError):self.compare((s,l,d))
    def test_original_info_orientation_literal_mutations_refused(self):
        for old,new in [('ORIGINAL=yes','ORIGINAL=no'),('AF=0.1000;NMP','AF=0.1;NMP'),('\tA\tC\t','\tC\tA\t')]:
            s,l,d=self.fixture();d=[x.replace(old,new) for x in d]
            with self.subTest(old=old),self.assertRaises(ValueError):self.compare((s,l,d))
    def test_full_ledger_mutations_refused(self):
        for key,value in [('source_ordinal',True),('source_record_sha256','b'*64),('source_sha256','b'*64),('selection','included'),('population_ALT_frequency_literal','0.1')]:
            s,l,d=self.fixture();r=json.loads(l[0]);r[key]=value;l[0]=json.dumps(r)
            with self.subTest(key=key),self.assertRaises(ValueError):self.compare((s,l,d))
    def test_ledger_omission_and_extra_refused(self):
        for kind in ['omit','extra']:
            s,l,d=self.fixture()
            if kind=='omit':l.pop()
            else:l.append(l[-1])
            with self.subTest(kind=kind),self.assertRaises((ValueError,StopIteration)):self.compare((s,l,d))
    def test_reference_mismatch_nonconcrete_and_duplicate_refused(self):
        s,l,d=self.fixture()
        with self.assertRaisesRegex(ValueError,'REF mismatch'):v.compare_records(s,l,d,lambda c,a,b:'G',8,'a'*64)
        for changed in [s[0].replace('\tA\tC\t','\tN\tC\t'),s[0]]:
            bad=copy.deepcopy(s);bad[1]=changed
            with self.assertRaises(ValueError):self.compare((bad,l,d))
    def test_duplicate_json_and_nonfinite_refused(self):
        for text in ['{"x":NaN}','{"x":1,"x":2}']:
            with self.assertRaises(ValueError):v.json_value(text)
    def authorization(self,root):
        plan=root/'plan.json';plan.write_text(json.dumps(dict(case_execution_authorized=False,ready_for_native_pilot=False)))
        a=dict(schema_version=1,public_transform_authorized=True,patient_assay_authorized=False,plan=dict(path=str(plan),sha256=sha(plan)),source_sha256={name:sha(HERE/name) for name in ['guard-public-panel.py','guard-fictional-pilot.py','prepare-public-panel.py','verify-public-panel.py']},output=str(g.OUT),launch=str(g.LAUNCH),limits=g.g.LIMITS)
        return a
    def test_authorization_scope_and_typed_flags(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);a=self.authorization(root);p=root/'authorization.json'
            with patch.object(g,'PLAN_SHA',a['plan']['sha256']):
                for key,value in [('public_transform_authorized',False),('public_transform_authorized',1),('patient_assay_authorized',True),('schema_version',True),('output','/other')]:
                    bad=copy.deepcopy(a);bad[key]=value;p.write_text(json.dumps(bad))
                    with self.subTest(key=key),self.assertRaises(ValueError):g.bind(p,sha(p))
                p.write_text(json.dumps(a));self.assertEqual(g.bind(p,sha(p)),a)
    def test_source_and_config_hash_refusals(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);a=self.authorization(root);p=root/'authorization.json';p.write_text(json.dumps(a))
            with self.assertRaises(ValueError):g.bind(p,'0'*64)
            with patch.object(g,'PLAN_SHA',a['plan']['sha256']):
                a['source_sha256']['verify-public-panel.py']='0'*64;p.write_text(json.dumps(a))
                with self.assertRaises(ValueError):g.bind(p,sha(p))
    def test_worker_requires_owned_parent(self):
        with patch.dict(g.os.environ,{},clear=True),self.assertRaisesRegex(ValueError,'owned guard'):g.worker({})
    def test_native_full_length_fictional_transform_and_reconciliation(self):
        try:import pysam
        except ImportError:self.skipTest('native pysam required')
        producer=module('prepare-public-panel.py')
        context=nullcontext(os.environ['MVA_PUBLIC_FIXTURE_ROOT']) if os.environ.get('MVA_PUBLIC_FIXTURE_ROOT') else tempfile.TemporaryDirectory()
        with context as d:
            root=Path(d);root.mkdir(parents=True,exist_ok=True);public=root/'ref';public.mkdir();ref=public/'fictional.fa';n=224487;ref.write_text('>1\n'+'A'*(n*2)+'\n');pysam.faidx(str(ref));raw=public/'sites.vcf'
            with raw.open('w') as f:
                f.write(f'##fileformat=VCFv4.2\n##contig=<ID=1,length={n*2}>\n##INFO=<ID=POP_AF,Number=1,Type=Float,Description="Fictional">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
                for i in range(n):f.write(f'1\t{2*i+1}\trs{i}\tA\tC\t.\tPASS\tPOP_AF='+['0.001','0.01','0.0100001','0.1000','0.199999','0.2','0.21','1'][i%8]+'\n')
            panel=Path(str(raw)+'.gz');pysam.tabix_compress(str(raw),str(panel));pysam.tabix_index(str(panel),preset='vcf')
            plan=dict(root=str(root),scope='public_marker_derivation_and_fictional_controls_only',case_execution_authorized=False,public_panel=dict(path=str(panel.relative_to(root)),sha256=sha(panel),index_path=str(panel.relative_to(root))+'.tbi',index_sha256=sha(Path(str(panel)+'.tbi')),minimum_exclusive='.01',maximum_exclusive='.2',frequency_source='POP_AF',frequency_target='AF',expected_records=n,source_ordinal_definition='fictional fixture ordinal'),reference=dict(path=str(ref.relative_to(root)),sha256=sha(ref),index_path=str(ref.relative_to(root))+'.fai',index_sha256=sha(Path(str(ref)+'.fai'))))
            cp=root/'plan.json';cp.write_text(json.dumps(plan));out=root/'derived'
            with patch.multiple(producer,EXPECTED_PANEL=sha(panel),EXPECTED_REFERENCE=sha(ref)):result=producer.prepare(cp,sha(cp),out)
            with pysam.FastaFile(str(ref)) as fasta,gzip.open(panel,'rt',newline='') as source,gzip.open(out/'source-ordinal-ledger.jsonl.gz','rt') as ledger,gzip.open(out/'panel.vcf.gz','rt',newline='') as derived:
                checked=v.compare_records((x for x in source if not x.startswith('#')),ledger,(x for x in derived if not x.startswith('#')),fasta.fetch,n,sha(panel))
            self.assertEqual(checked['counts'],result['counts']);self.assertEqual(checked['counts']['included'],84183)
            with pysam.TabixFile(str(out/'panel.vcf.gz')) as index:
                h=hashlib.sha256();count=0
                for line in index.fetch('1'):h.update((line+'\n').encode());count+=1
                self.assertEqual(count,84183);self.assertEqual(h.hexdigest(),checked['indexed_record_sha256']['1'])
            self.assertEqual(sha(panel),plan['public_panel']['sha256']);self.assertEqual(sha(ref),plan['reference']['sha256'])
if __name__=='__main__':unittest.main(verbosity=2)
