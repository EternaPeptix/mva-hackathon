#!/usr/bin/env python3
"""Independent saved public VCF/ledger readback; never import producer computation."""
import argparse,collections,datetime,hashlib,json,math
from pathlib import Path
import independent_vcf_reader_v2 as reader
OUTCOMES=('matched_by_both','DELLY_only','SvABA_only','missed_by_both')
TYPES=('DEL','INS')
BINS=('50-99','100-499','500-999','1000-9999','10000-50000')

def require(v,label):
 if not v:raise ValueError(label)
def load(p):return json.loads(Path(p).read_text())
def pin(p):
 p=Path(p).resolve();a=p.stat();h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 z=p.stat();require((a.st_size,a.st_mtime_ns)==(z.st_size,z.st_mtime_ns),'file changed during hash')
 return dict(path=str(p),bytes=z.st_size,mtime_ns=z.st_mtime_ns,sha256=h.hexdigest())
def bind(p,expect):
 q=pin(p)
 for k in ('bytes','sha256'):require(q[k]==expect[k],'fresh pin differs '+str(p))
 if 'path' in expect:require(q['path']==str(Path(expect['path']).resolve()),'path pin differs')
 return q
def eq(a,b,label):
 require(json.dumps(a,sort_keys=True,allow_nan=False)==json.dumps(b,sort_keys=True,allow_nan=False),label)
def gt(r):
 text=r['format_values']['GT'];phase='|' in text
 require(not(phase and '/' in text),'mixed GT separators')
 vals=tuple(None if x=='.' else int(x) for x in text.replace('|','/').split('/'))
 return tuple(sorted(vals,key=lambda x:-1 if x is None else x))
def shape(r):
 kind=r['info']['SVTYPE'];size=abs(int(r['info']['SVLEN']))
 require(kind in TYPES and 50<=size<=50000 and len(r['alts'])==1,'unsupported fixed selected truth')
 label=next(name for end,name in zip((100,500,1000,10000,50001),BINS) if size<end)
 return kind,size,label

def independent_tables(rows):
 def margin(group):
  c=collections.Counter(r['outcome'] for r in group);n=len(group)
  return dict(truth_records=n,outcomes={x:c[x] for x in OUTCOMES},shared_miss_fraction=c[OUTCOMES[3]]/n if n else None,union_matched_fraction=(n-c[OUTCOMES[3]])/n if n else None)
 def state(field):
  return ('not_encoded',None) if not field['present'] else ('explicit_missing',None) if field['value'] is None else ('value',field['value'])
 result={'total':margin(rows),'class_size':[],'TRF_RM_by_class_size':[],'RM_class_family_by_class_size':[],'REMAP_by_class_size':[],'raw_LCR_by_class_size_and_outcome':[]}
 for kind in TYPES:
  for size in BINS:
   subset=[r for r in rows if (r['type'],r['size_bin'])==(kind,size)];prefix=dict(type=kind,size_bin=size)
   result['class_size'].append({**prefix,**margin(subset)})
   cells={(t,m):[] for t in ('flag_present','flag_not_encoded') for m in ('value_present','fields_explicit_missing','not_encoded')}
   for r in subset:
    a=r['annotations'];t='flag_present' if a['TRF']['present'] else 'flag_not_encoded';fields=[a[x] for x in ('RM_score','RM_repeat','RM_clsfam')]
    m='value_present' if any(x['value'] is not None for x in fields) else 'fields_explicit_missing' if any(x['present'] for x in fields) else 'not_encoded'
    cells[(t,m)].append(r)
   for (t,m),group in cells.items():result['TRF_RM_by_class_size'].append({**prefix,'TRF_annotation_state':t,'RM_annotation_state':m,**margin(group)})
   for tag,outname in [('RM_clsfam','RM_class_family_by_class_size'),('REMAP','REMAP_by_class_size')]:
    groups={('not_encoded',None):[],('explicit_missing',None):[]}
    for r in subset:groups.setdefault(state(r['annotations'][tag]),[]).append(r)
    for (st,value),group in sorted(groups.items(),key=lambda x:str(x[0])):result[outname].append({**prefix,'state':st,'literal_value':value,**margin(group)})
   for outcome in OUTCOMES:
    group=[r for r in subset if r['outcome']==outcome];a=[r['annotations']['LCR'] for r in group]
    vals=sorted(x['value'] for x in a if x['value'] is not None);n=len(vals)
    median=None if not n else vals[n//2] if n%2 else (vals[n//2-1]+vals[n//2])/2
    result['raw_LCR_by_class_size_and_outcome'].append({**prefix,'outcome':outcome,'truth_records':len(group),'not_encoded':sum(not x['present'] for x in a),'explicit_missing':sum(x['present'] and x['value'] is None for x in a),'finite_values':n,'minimum':vals[0] if n else None,'median':median,'maximum':vals[-1] if n else None})
 return result


def run(args):
 out=Path(args.output).resolve();require(not out.exists(),'Existing verification output preserved')
 cfgpath=Path(args.config);cp=pin(cfgpath);require(cp['sha256']==args.expected_config_sha256,'config SHA differs')
 cfg=load(cfgpath);require(cfg['ready'] is True,'unapproved config');run=Path(args.run).resolve();launch=Path(args.launch).resolve()
 producer=load(run/'status.json');guard=load(launch/'status.json')
 require(producer['status']==guard['status']=='complete' and guard['exit_code']==0 and guard['all_pins_stable'] is True,'producer/guard incomplete')
 require(set(producer['outputs'])=={'ledger.jsonl','summary.json','headers.json'},'producer payload membership')
 pins={}
 def keep(p,expect=None):
  q=bind(p,expect) if expect else pin(p);pins[q['path']]=q;return q
 keep(cfgpath,cp);keep(__file__);keep(reader.__file__)
 keep(guard['guard']['path'],guard['guard']);keep(run/'status.json',guard['producer_status']);keep(run/'summary.json',guard['producer_summary']);keep(launch/'status.json')
 for role,row in cfg['inputs'].items():keep(row['path'],row)
 for name,row in producer['outputs'].items():
  require(Path(name).name==name,'output path escape');keep(run/name,row)
 plan_path=guard['plan']['path'];keep(plan_path,guard['plan']);plan=load(plan_path)
 require(Path(plan['output']).resolve()==run and Path(plan['launch']).resolve()==launch,'guard plan paths differ')
 for row in plan['pins']:keep(row['path'],row)
 summary=load(run/'summary.json');headers=load(run/'headers.json');rows=[json.loads(x) for x in (run/'ledger.jsonl').read_text().splitlines()]
 keep(summary['source']['path'],summary['source']);eq(summary['config'],{k:cp[k] for k in ('path','bytes','sha256')},'config summary pin')
 for role,p in cfg['inputs'].items():eq(summary['inputs'][role],{k:p[k] for k in ('path','bytes','sha256')},'summary input bind')
 parsed={name:reader.read_vcf(cfg['inputs'][role]['path']) for name,role in [('DELLY','delly_final_truth'),('SvABA','svaba_final_truth')]}
 keyed={}
 for name,data in parsed.items():
  keyed[name]={reader.identity(r):r for r in data['rows']};require(len(keyed[name])==len(data['rows']),'duplicate final full identity')
 keys=set(keyed['DELLY']);require(keys==set(keyed['SvABA']),'caller truth universes differ')
 original=reader.read_vcf(cfg['inputs']['original_truth_chr20_no_star']['path'],keys);parsed['original']=original
 require(original['total_records']==cfg['expected']['original_full_context_records'],'original count')
 require(len(original['rows'])==len(keys)==len(rows)==cfg['expected']['truth_records']==producer['truth_records'],'all ledger denominators')
 originals={reader.identity(r):r for r in original['rows']};require(len(originals)==len(keys),'duplicate original identity')
 for name,data in parsed.items():
  require(data['samples']==['HG002'],'public sample identity')
  eq(headers[name]['headers'],data['headers'],'literal headers '+name)
  eq(headers[name]['info_definitions'],{k:{x:v[x] for x in ('ID','Number','Type')} for k,v in data['info_definitions'].items()},'declared INFO headers')
 outcome_map={('TP','TP'):OUTCOMES[0],('TP','FN'):OUTCOMES[1],('FN','TP'):OUTCOMES[2],('FN','FN'):OUTCOMES[3]}
 bed=[]
 for line in Path(cfg['inputs']['official_chr20_benchmark']['path']).read_text().splitlines():
  chrom,a,b=line.split('\t');a,b=int(a),int(b);require(chrom=='chr20' and a<b and a>=0,'BED syntax');bed.append((a,b))
 expected=[]
 for row,orig in zip(rows,original['rows']):
  identity=reader.identity(orig);other=[keyed[x][identity] for x in ('DELLY','SvABA')]
  typ,size,binname=shape(orig);ann={k:{x:v[x] for x in ('present','value','raw')} for k,v in orig['annotations'].items()}
  require(gt(orig)==gt(other[0])==gt(other[1]),'genotype allele-content conservation')
  for r in other:
   require(shape(r)==(typ,size,binname),'type/size conservation')
   eq({k:(v['present'],v['value']) for k,v in orig['annotations'].items()},{k:(v['present'],v['value']) for k,v in r['annotations'].items()},'typed annotation conservation')
  labels=tuple(r['format_values']['BD'] for r in other);require(labels in outcome_map,'final decision domain')
  exp={'original_ordinal':orig['ordinal'],'chrom':orig['chrom'],'pos_1based':orig['pos'],'id':orig['id'] if orig['id'] is not None else '.',
       'ref':orig['ref'],'alts':orig['alts'],'allele_sha256':hashlib.sha256(json.dumps([orig['ref'],orig['alts']],separators=(',',':')).encode()).hexdigest(),
       'ref_span_0based_half_open':[orig['pos']-1,orig['pos']-1+len(orig['ref'])],'type':typ,'size':size,'size_bin':binname,'original_raw_record':orig['record_literal'],'original_GT':orig['format_values']['GT'],
       'REF_span_contained_in_one_benchmark_interval':any(a<=orig['pos']-1 and orig['pos']-1+len(orig['ref'])<=b for a,b in bed),'GT_encoding_matches_original':{name:r['format_values']['GT']==orig['format_values']['GT'] for name,r in zip(('DELLY','SvABA'),other)},
       'annotations':ann,'outcome':outcome_map[labels],'callers':{name:{'refined_baseline_ordinal':r['ordinal'],'decision':r['format_values']['BD'],'GT':r['format_values']['GT'],'raw_info':r['info'],'raw_sample':r['format_values']} for name,r in zip(('DELLY','SvABA'),other)}}
  eq(row,exp,'full independently reconstructed ledger row');expected.append(exp)
 tables=independent_tables(expected)
 for name,data in tables.items():eq(summary[name],data,'independent aggregate '+name)
 eq(tables['total']['outcomes'],cfg['expected']['outcomes'],'expected four outcomes')
 for caller in ('DELLY','SvABA'):
  counts=collections.Counter(('INS_like' if r['type']=='INS' else 'DEL')+':'+r['size_bin']+':'+r['callers'][caller]['decision'] for r in expected)
  eq(dict(counts),cfg['expected']['per_caller_by_size'][caller],'independent published class/size margin')
 cls=collections.Counter(r['type'] for r in expected);eq(dict(cls),cfg['expected']['classes'],'class margin')
 for kind in TYPES:
  counts=collections.Counter(r['outcome'] for r in expected if r['type']==kind)
  eq({x:counts[x] for x in OUTCOMES},cfg['expected']['outcomes_by_class'][kind],'class four-outcome margin')
 bed=[]
 for line in Path(cfg['inputs']['official_chr20_benchmark']['path']).read_text().splitlines():
  chrom,a,b=line.split('\t');a,b=int(a),int(b);require(chrom=='chr20' and a<b and a>=0,'BED syntax');bed.append((a,b))
 require(bed==sorted(bed) and all(bed[i][1]<=bed[i+1][0] for i in range(len(bed)-1)),'BED ordering/nonoverlap')
 require(len(bed)==cfg['expected']['benchmark_intervals'] and sum(b-a for a,b in bed)==cfg['expected']['benchmark_bases'],'official BED scope')
 for p,prior in list(pins.items()):require(pin(p)==prior,'input/output/source unstable')
 receipt={'status':'complete','verified_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'verification_kind':'Separate stdlib VCF-text parsing and full row/aggregate reconstruction; no producer computation imports','truth_records':len(rows),'original_context_records':original['total_records'],'all_ledger_rows_exact':True,'all_aggregate_tables_exact':True,'totals':tables['total'],'class_counts':dict(cls),'benchmark_intervals':len(bed),'benchmark_bases':sum(b-a for a,b in bed),'all_verified_files_stable':True,'verified_files':list(pins.values()),'limitations':['No original caller/Truvari benchmark or scientific sensitivity validation repeated.','Genotype allele content compared under historical order-independent policy; raw GT/phase encodings checked and retained individually, not asserted unchanged.','Literal existing annotations only; LCR remains continuous and missing fields do not imply repeat absence.','Historical labels are preserved; no nearby-call/filter-cause attribution or clinical conclusion.']}
 out.mkdir(parents=True,exist_ok=False);encoded=json.dumps(receipt,indent=2,sort_keys=True,allow_nan=False)+'\n'
 for name in ('verification.json','summary.json'):(out/name).write_text(encoded)
 (out/'status.json').write_text(json.dumps({'status':'complete','verification_sha256':pin(out/'verification.json')['sha256'],'outputs':{name:pin(out/name) for name in ('verification.json','summary.json')}},indent=2)+'\n')
 print(json.dumps({'status':'complete','truth_records':len(rows),'verification_sha256':pin(out/'verification.json')['sha256']}))

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__)
 for name in ('config','expected-config-sha256','run','launch','output'):p.add_argument('--'+name,required=True)
 run(p.parse_args())
