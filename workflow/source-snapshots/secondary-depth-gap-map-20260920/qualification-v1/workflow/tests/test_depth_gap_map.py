"""Fictional saved-depth controls; no patient/reference resource is opened."""
from collections import Counter
import copy
import gzip
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

SOURCE = Path(__file__).resolve().parents[1]/'scripts/map_depth_gaps.py'
spec = importlib.util.spec_from_file_location('mapper', SOURCE)
M = importlib.util.module_from_spec(spec); spec.loader.exec_module(M)


GFF = '''##gff-version 3
1\ttest\tgene\t1\t12\t.\t+\t.\tID=gene:A;Name=A;biotype=protein_coding
1\ttest\tmRNA\t1\t12\t.\t+\t.\tID=transcript:A1;Parent=gene:A
1\ttest\tCDS\t1\t4\t.\t+\t0\tID=cds:a1;Parent=transcript:A1
1\ttest\tCDS\t9\t12\t.\t+\t2\tID=cds:a2;Parent=transcript:A1
1\ttest\tmRNA\t1\t12\t.\t+\t.\tID=transcript:A2;Parent=gene:A
1\ttest\tCDS\t3\t6\t.\t+\t1\tID=cds:a3;Parent=transcript:A2
1\ttest\tgene\t5\t14\t.\t-\t.\tID=gene:B;Name=B;biotype=protein_coding
1\ttest\tmRNA\t5\t14\t.\t-\t.\tID=transcript:B1;Parent=gene:B
1\ttest\tCDS\t5\t10\t.\t-\t0\tID=cds:b1;Parent=transcript:B1
X\ttest\tgene\t20\t21\t.\t-\t.\tID=gene:C;Name=C;biotype=protein_coding
X\ttest\tmRNA\t20\t21\t.\t-\t.\tID=transcript:C1;Parent=gene:C
X\ttest\tCDS\t20\t21\t.\t-\t0\tID=cds:c1;Parent=transcript:C1
'''


def write(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False)+'\n')


def fixture(base, first=None, second=None):
    base=base.resolve()
    base.mkdir(parents=True)
    gff=base/'frozen.gff3.gz'
    with gzip.open(gff,'wt') as f:f.write(GFF)
    fai=base/'reference.fa.fai';fai.write_text('1\t100\t0\t100\t101\nX\t100\t0\t100\t101\n')
    design=M.prepare_targets(gff,['A','B','C'],{'1':100,'X':100})
    first=first if first is not None else [0,1,19,20,30,0,4,20,19,20,0,30,20,0]
    second=second if second is not None else [0,0,19,19,30,0,0,20,0,19,0,30,20,0]
    require_size=sum(w['end']-w['start'] for w in design['windows'])
    assert len(first)==len(second)==require_size
    files={'gff':gff,'fai':fai}
    results={}
    for mode,values,mode_value in zip(M.MODES,(first,second),M.MODE_VALUES):
        depth=base/(mode+'.depth.tsv.gz');files[mode+'_depth']=depth
        i=0;windows=[]
        with gzip.open(depth,'wt') as f:
            for w in design['windows']:
                h=Counter()
                for pos in range(w['start'],w['end']):
                    d=values[i];i+=1;f.write(f"{w['chrom']}\t{pos+1}\t{d}\n");h[d]+=1
                s=dict(bases=sum(h.values()),sum_depth=sum(d*n for d,n in h.items()),
                       mean=sum(d*n for d,n in h.items())/sum(h.values()),histogram=dict(h),
                       bases_at_least={str(t):sum(n for d,n in h.items() if d>=t) for t in (1,10,20,30)})
                windows.append(dict(w,span_depth=s,ACGT_depth=copy.deepcopy(s)))
        from coding_depth_targets import summarize_depth
        groups=summarize_depth(windows,design)
        result=dict(mode=mode_value,position_count=require_size,groups=groups,windows=windows)
        p=base/(mode+'.json');write(p,result);files[mode+'_summary']=p
        results[mode]={k:result[k] for k in ('mode','position_count','groups')}
    annotation_inputs={'gff':M.fingerprint(gff),'reference_index':M.fingerprint(fai)}
    public=copy.deepcopy(design)
    for w in public['windows']:w.update(reference_ACGT_bases=w['end']-w['start'],reference_non_ACGT_bases=0)
    for g in public['genes']:g['original_gene_feature']={'fictional':True}
    public.update(coordinates='0-based half-open',reads_inspected=False,inputs=annotation_inputs)
    files['public_design']=base/'reference-design.json';write(files['public_design'],public)
    public_summary={'status':'complete','panel_non_ACGT_bases':0,'panel_ACGT_bases':require_size}
    files['public_summary']=base/'target-summary.json';write(files['public_summary'],public_summary)
    files['public_status']=base/'public-status.json'
    write(files['public_status'],dict(status='complete',reads_inspected=False,inputs=annotation_inputs,
          outputs={files[k].name:M.fingerprint(files[k]) for k in ('public_design','public_summary')}))
    selection=copy.deepcopy(design);selection.update(coordinates='0-based half-open; all CDS phases retained, not translated')
    files['selection']=base/'selection.json';write(files['selection'],selection)
    summary=dict(status='complete',depth_modes=results,inputs=annotation_inputs,
                 genes=[{k:g[k] for k in ('id','name','chrom','strand','CDS_union_bases','transcript_ids')} for g in design['genes']],
                 panel_union_bases=require_size,gene_union_bases_sum=design['gene_union_bases_sum'])
    files['summary']=base/'summary.json';write(files['summary'],summary)
    files['status']=base/'status.json'
    write(files['status'],dict(status='complete',summary_sha256=M.fingerprint(files['summary'])['sha256'],
          outputs={p.name:M.fingerprint(p) for k,p in files.items() if k in ('selection','summary') or k.startswith(M.MODES)}))
    files['verification']=base/'verification.json'
    write(files['verification'],dict(status='complete',all_verified_inputs_stable=True,public_target_reconciled=True,
          producer_status=M.fingerprint(files['status']),producer_summary=M.fingerprint(files['summary']),
          verified_files={str(p):M.fingerprint(p) for p in files.values() if p!=files['verification']}))
    config=dict(schema_version=1,ready=True,assembly='GRCh38',scope='84-gene-all-transcript-CDS-saved-depth',
                depth_modes=list(M.MODE_VALUES),expected_scope=dict(genes=3,transcripts=4,windows=len(design['windows']),
                  panel_union_bases=require_size,gene_union_bases_sum=design['gene_union_bases_sum']),
                source_sha256={r:M.fingerprint(SOURCE.parent/r)['sha256'] for r in M.SOURCE_FILES},
                inputs={k:M.fingerprint(p) for k,p in files.items()})
    cfg=base/'config.json';write(cfg,config)
    return cfg,config,files,design


class MappingControls(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.base=Path(self.temp.name)
        self.cfg,self.config,self.files,self.design=fixture(self.base/'fixture')
    def tearDown(self):self.temp.cleanup()
    def bound(self):return M.bind(self.cfg,M.fingerprint(self.cfg)['sha256'])
    def calculate(self,first=None,second=None):
        if first is not None:
            _,_,self.files,self.design=fixture(self.base/'other',first,second)
        rows=[];out,hist,windows=M.calculate(self.design,[self.files[m+'_depth'] for m in M.MODES],rows.append)
        return out,rows,hist,windows
    def mutate_json(self,role,action):
        x=M.read_json(self.files[role]);action(x);write(self.files[role],x)
        self.config['inputs'][role]=M.fingerprint(self.files[role]);write(self.cfg,self.config)
    def test_full_binding_geometry_and_saved_reconciliation(self):
        c,p,d=self.bound();design=M.geometry(c,p,d);out,rows,h,w=self.calculate();M.reconcile(design,d,h,w)
        self.assertEqual(out['panel_union'][M.MODES[0]]['CDS_union_bases'],14)
    def test_transcript_overlap_not_panel_multiplication(self):
        out,rows,_,_=self.calculate();by={x['transcript_id']:x for x in out['transcripts']}
        self.assertEqual([by[t]['modes'][M.MODES[0]]['CDS_union_bases'] for t in ('transcript:A1','transcript:A2','transcript:B1','transcript:C1')],[8,4,6,2])
        self.assertEqual(sum(g['modes'][M.MODES[0]]['CDS_union_bases'] for g in out['genes']),18)
        self.assertEqual(out['panel_union'][M.MODES[0]]['CDS_union_bases'],14)
    def test_exact_zero_and_threshold_19_20(self):
        out,rows,_,_=self.calculate();self.assertEqual(out['panel_union'][M.MODES[0]]['zero_depth_bases'],4)
        self.assertEqual(out['panel_union'][M.MODES[0]]['below20_bases'],8)
        self.assertEqual(out['panel_union'][M.MODES[1]]['below20_bases'],10)
    def test_first_last_and_negative_strand_coordinates(self):
        out,rows,_,_=self.calculate();self.assertEqual(rows[0]['start'],0)
        self.assertEqual((rows[-1]['chrom'],rows[-1]['start'],rows[-1]['end']),('X',20,21))
        b=next(t for t in out['transcripts'] if t['transcript_id']=='transcript:B1');self.assertEqual(b['strand'],'-')
    def test_no_gaps_is_explicit_empty(self):
        out,rows,_,_=self.calculate([20]*14,[20]*14);self.assertEqual(rows,[]);self.assertEqual(out['gap_intervals'],0)
    def test_all_zero(self):
        out,rows,_,_=self.calculate([0]*14,[0]*14);self.assertEqual(out['panel_union'][M.MODES[1]]['zero_depth_bases'],14)
    def test_transitions_partition(self):
        out,rows,_,_=self.calculate();self.assertEqual(sum(x['bases'] for x in out['paired_categories']),14)
        self.assertEqual(out['MQ60_depth_greater_than_MQ20_bases'],0)
    def test_strict_increase_is_reported_not_silently_recoded(self):
        out,rows,_,_=self.calculate([19]*14,[20]*14);self.assertEqual(out['MQ60_depth_greater_than_MQ20_bases'],14)
    def test_adjacent_cds_membership_is_not_merged(self):
        out,rows,_,_=self.calculate([0]*14,[0]*14)
        self.assertTrue(any(r['end']==2 for r in rows));self.assertTrue(any(r['start']==2 for r in rows))
    def test_intron_and_chromosome_are_not_bridged(self):
        out,rows,_,_=self.calculate([0]*14,[0]*14)
        self.assertFalse(any(r['chrom']=='1' and r['end']>12 for r in rows));self.assertTrue(all(r['start']>=19 for r in rows if r['chrom']=='X'))
    def test_duplicate_cds_line_multiplicity_keeps_union(self):
        self.design['genes'][0]['cds_rows'].append(dict(self.design['genes'][0]['cds_rows'][0],ordinal=99))
        out,rows,_,_=self.calculate();self.assertEqual(out['panel_union'][M.MODES[0]]['CDS_union_bases'],14)
        self.assertTrue(any(m['gff_line_ordinal']==99 for r in rows for m in r['source_CDS_members']))
    def bad_depth(self,transform):
        p=self.files[M.MODES[0]+'_depth']
        with gzip.open(p,'rt') as f:rows=f.readlines()
        rows=transform(rows)
        with gzip.open(p,'wt') as f:f.writelines(rows)
        with self.assertRaises(ValueError):self.calculate()
    def test_missing_depth_refused(self):self.bad_depth(lambda r:r[:-1])
    def test_duplicate_depth_refused(self):self.bad_depth(lambda r:[r[0]]+r)
    def test_extra_depth_refused(self):self.bad_depth(lambda r:r+[r[-1]])
    def test_order_refused(self):self.bad_depth(lambda r:[r[1],r[0]]+r[2:])
    def test_negative_depth_refused(self):self.bad_depth(lambda r:[r[0].rsplit('\t',1)[0]+'\t-1\n']+r[1:])
    def test_float_depth_refused(self):self.bad_depth(lambda r:[r[0].rsplit('\t',1)[0]+'\t1.0\n']+r[1:])
    def test_coordinate_zero_refused(self):self.bad_depth(lambda r:['1\t0\t5\n']+r[1:])
    def test_config_sha_refused(self):
        with self.assertRaises(ValueError):M.bind(self.cfg,'0'*64)
    def test_input_hash_refused(self):
        self.config['inputs']['gff']['sha256']='0'*64;write(self.cfg,self.config)
        with self.assertRaises(ValueError):self.bound()
    def test_source_hash_refused(self):
        self.config['source_sha256']['map_depth_gaps.py']='0'*64;write(self.cfg,self.config)
        with self.assertRaises(ValueError):self.bound()
    def test_wrong_path_same_bytes_receipt_refused(self):
        self.mutate_json('verification',lambda x:x['producer_summary'].update(path='/fictional/wrong'))
        with self.assertRaises(ValueError):self.bound()
    def test_failed_independent_receipt_refused(self):
        self.mutate_json('verification',lambda x:x.update(status='failed'))
        with self.assertRaises(ValueError):self.bound()
    def test_missing_independent_output_pin_refused(self):
        self.mutate_json('verification',lambda x:x['verified_files'].pop(str(self.files[M.MODES[0]+'_depth'])))
        with self.assertRaises(KeyError):self.bound()
    def test_unapproved_config_refused(self):
        self.config['ready']=False;write(self.cfg,self.config)
        with self.assertRaises(ValueError):self.bound()
    def test_alignment_input_role_refused(self):
        self.config['inputs']['alignment']=self.config['inputs']['gff'];write(self.cfg,self.config)
        with self.assertRaises(ValueError):self.bound()
    def test_unknown_reference_bases_refused(self):
        c,p,d=self.bound();d['public_design']['windows'][0]['reference_non_ACGT_bases']=1
        with self.assertRaises(ValueError):M.geometry(c,p,d)
    def test_altered_public_gene_geometry_refused(self):
        c,p,d=self.bound();d['public_design']['genes'][0]['cds_rows'][0]['phase']=2
        with self.assertRaises(ValueError):M.geometry(c,p,d)
    def test_bool_count_refused(self):
        c,p,d=self.bound();c['expected_scope']['genes']=True
        with self.assertRaises(ValueError):M.geometry(c,p,d)
    def test_saved_histogram_mutation_refused(self):
        c,p,d=self.bound();out,rows,h,w=self.calculate();d[M.MODES[0]+'_summary']['windows'][0]['span_depth']['histogram']['0']=99
        with self.assertRaises(ValueError):M.reconcile(self.design,d,h,w)
    def test_nonfinite_and_duplicate_json_refused(self):
        for text in ('{"x":NaN}','{"x":Infinity}','{"x":1e999}','{"x":1,"x":2}'):
            p=self.base/'bad.json';p.write_text(text)
            with self.subTest(text=text),self.assertRaises(ValueError):M.read_json(p)
    def test_cli_end_to_end_and_existing_output_refused(self):
        out=self.base/'result';cmd=[sys.executable,str(SOURCE),'--config',str(self.cfg),
            '--expected-config-sha256',M.fingerprint(self.cfg)['sha256'],'--output',str(out)]
        result=subprocess.run(cmd,capture_output=True,text=True);self.assertEqual(result.returncode,0,result.stderr)
        state=M.read_json(out/'status.json');self.assertEqual(state['status'],'complete');self.assertIs(state['reads_opened'],False)
        original=(out/'summary.json').read_bytes();self.assertNotEqual(subprocess.run(cmd,capture_output=True).returncode,0)
        self.assertEqual((out/'summary.json').read_bytes(),original)


if __name__=='__main__':unittest.main()
