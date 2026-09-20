from pathlib import Path
from collections import defaultdict
import argparse,datetime,gzip,hashlib,json,re
parser=argparse.ArgumentParser(description='Reconcile saved depth report tables; no alignment access. Output must be new.')
parser.add_argument('--base',type=Path,default=Path('/Users/jeweled/mva-hackathon/artifacts/secondary-depth-gap-map-20260920'))
parser.add_argument('--original-depth',type=Path,default=Path('/Users/jeweled/mva-hackathon/artifacts/secondary-panel-depth-preparation-20260920/case-results/review-v2'))
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
BASE=args.base.resolve()
OUT=BASE/'independent-result-review'
assert not args.output.exists(), 'Output must be new'
MODES=('MQ20_BQ20','MQ60_BQ20')
def pin(p):return {'path':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def eq(a,b,label):
 if json.dumps(a,sort_keys=True)!=json.dumps(b,sort_keys=True):raise ValueError(label)
def load(p):return json.loads(p.read_text())
s=load(BASE/'map-v1/summary.json');v=load(OUT/'verification.json');loc=load(BASE/'gap-locations.json')
assert v['status']=='complete' and v['all_verified_files_stable'] is True
verified_increases=v['numerically_increased_MQ60_bases']
verified_crossings=v['upward_category_crossings']
assert type(verified_increases) is int and verified_increases>=0
assert type(verified_crossings) is int and verified_crossings>=0
for file in ['summary.json','status.json','gap-intervals.jsonl.gz']:
 expected=next(x for x in v['verified_files'] if x['path'].endswith('/map-v1/'+file));p=pin(BASE/'map-v1'/file);assert p['sha256']==expected['sha256'] and p['bytes']==expected['bytes']
with gzip.open(BASE/'map-v1/gap-intervals.jsonl.gz','rt') as f:atoms=[json.loads(x) for x in f]
def expected_runs(identity,kind,index,criterion):
 selected=[]
 for ordinal,a in enumerate(atoms,1):
  if identity is not None and identity not in a[kind]:continue
  cat=a['categories'][index]
  if not (cat=='zero' if criterion=='zero' else cat!='at-least-20'):continue
  selected.append((ordinal,a))
 runs=[]
 for ordinal,a in selected:
  if runs and runs[-1]['chrom']==a['chrom'] and a['start']<=runs[-1]['end']:
   r=runs[-1];r['end']=max(r['end'],a['end']);r['atomic_row_numbers'].append(ordinal)
   r['transcript_ids_overlapping_any_part']=sorted(set(r['transcript_ids_overlapping_any_part'])|set(a['transcript_ids']))
  else:runs.append(dict(chrom=a['chrom'],start=a['start'],end=a['end'],atomic_row_numbers=[ordinal],transcript_ids_overlapping_any_part=sorted(a['transcript_ids'])))
 return runs
expected_appendix=[];gene_count=0;transcript_count=0
for kind,key,idkey in [('genes','gene_ids','gene_id'),('transcripts','transcript_ids','transcript_id')]:
 assert len(loc[kind])==len(s[kind])
 for original,derivative in zip(s[kind],loc[kind]):
  expected_meta={k:val for k,val in original.items() if k!='cds_rows'};eq(expected_meta,{k:val for k,val in derivative.items() if k!='intervals'},kind+' metadata')
  for i,m in enumerate(MODES):
   for criterion in ['below20','zero']:
    expected=expected_runs(original[idkey],key,i,criterion);eq(expected,derivative['intervals'][m][criterion],kind+' exact merged runs')
    measured=sum(x['end']-x['start'] for x in expected);assert measured==original['modes'][m][criterion+'_bases' if criterion=='below20' else 'zero_depth_bases']
    if kind=='genes':
     for x in expected:expected_appendix.append(f"| {original['name']} | {m} | {criterion} | {x['chrom']} | {x['start']} | {x['end']} | {x['end']-x['start']} | {','.join(map(str,x['atomic_row_numbers']))} |")
  if kind=='genes':gene_count+=1
  else:transcript_count+=1
for i,m in enumerate(MODES):
 for criterion in ['below20','zero']:eq(expected_runs(None,None,i,criterion),loc['panel']['intervals'][m][criterion],'panel merged runs')
eq(loc['panel']['modes'],s['panel_union'],'panel metrics');eq(loc['paired_categories'],s['paired_categories'],'paired categories');eq(loc['strict_depth_numeric_increase_bases'],verified_increases,'strict increases');eq(s['MQ60_depth_greater_than_MQ20_bases'],verified_increases,'independent increase count')
appendix_lines=[l for l in (BASE/'GAP-LOCATIONS.md').read_text().splitlines() if re.match(r'^\| [A-Z0-9]+ \| MQ(?:20|60)_BQ20 \| (?:zero|below20) \|',l)]
# Human table uses a space after each comma; compare fields after whitespace removal only for row IDs.
def fields(line):return [v.strip().replace(', ', ',') for v in line.strip('|').split('|')]
eq([fields(x) for x in expected_appendix],[fields(x) for x in appendix_lines],'all appendix row coordinates/counts')
main=(BASE/'RESULT.md').read_text();expected_table=[]
for g in loc['genes']:
 if not any(g['modes'][m]['below20_bases'] for m in MODES):continue
 txs=[t for t in s['transcripts'] if t['gene_id']==g['gene_id']]
 count=[sum(t['modes'][m]['below20_bases']>0 for t in txs) for m in MODES]
 longest=min(g['intervals'][MODES[1]]['below20'],key=lambda x:(-(x['end']-x['start']),x['chrom'],x['start'],x['end']))
 expected_table.append(f"| {g['name']} | {g['modes'][MODES[0]]['CDS_union_bases']} | {g['modes'][MODES[0]]['below20_bases']} | {g['modes'][MODES[1]]['below20_bases']} | {g['modes'][MODES[1]]['zero_depth_bases']} | {count[0]}/{count[1]}/{len(txs)} | {longest['chrom']}:{longest['start']}–{longest['end']} |")
actual=[l for l in main.splitlines() if re.match(r'^\| [A-Z0-9]+ \| [0-9]+ \| [0-9]+ \| [0-9]+ \| [0-9]+ \| [0-9]+/[0-9]+/[0-9]+ \|',l)]
eq(expected_table,actual,f'all {len(expected_table)} main gene rows')
runs={m:{c:sum(len(g['intervals'][m][c]) for g in loc['genes']) for c in ('below20','zero')} for m in MODES}
eq(runs,load(BASE/'location-derivation.json')['gene_interval_counts'],'gene interval counts')
links=re.findall(r'\]\(([^)]+)\)',main)
for t in links:assert t.startswith('http') or (BASE/t).resolve().exists()
# Inspect only the original saved command metadata, never its referenced alignment.
original=args.original_depth.resolve()
cmd_path=original/'commands.json';cmd_pin=pin(cmd_path);prior=load(original/'status.json')['outputs']['commands.json'];assert cmd_pin['sha256']==prior['sha256'] and cmd_pin['bytes']==prior['bytes']
commands=[r for r in load(cmd_path) if r['label'] in MODES]
assert len(commands)==2
for row,mq in zip(commands,[20,60]):
 c=row['command'];assert row['exit_code']==0 and c[1]=='depth' and c[c.index('-q')+1]=='20' and c[c.index('-Q')+1]==str(mq) and c[c.index('-G')+1]=='0xF04' and '-s' in c and '-aa' in c and '-J' not in c
receipt={'status':'complete','reviewed_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source':pin(Path(__file__)),'inputs':[pin(BASE/n) for n in ['RESULT.md','GAP-LOCATIONS.md','gap-locations.json','location-derivation.json','map-v1/summary.json','map-v1/gap-intervals.jsonl.gz']]+[pin(OUT/'verification.json'),cmd_pin], 'counts':{'all_genes':gene_count,'all_transcripts':transcript_count,'main_gene_rows':len(actual),'location_appendix_rows':len(appendix_lines),'gene_runs':runs,'report_links':len(links)},'findings':[], 'interpretation_review':'Not assessed by this executable. Only the enumerated saved tables, links, and command metadata are reconciled; scientific interpretation requires separate review.', 'original_samtools_command_options':[['depth','-aa','-b','<regions.bed>','-q','20','-Q',str(q),'-G','0xF04','-s','<selected.bam>'] for q in [20,60]],'numeric_increase_note':{'saved_data_count':verified_increases,'upward_category_crossings':verified_crossings,'official_documentation':'https://www.htslib.org/doc/samtools-depth.html','documentation_scope':'Current official manual (samtools1.24) documents -q base-quality, -Q mapping-quality, -s first-read overlap suppression and no BAQ.','interpretation':'Mate-overlap suppression and quality filtering mean the result is not simply a per-base count of independently selected reads. This is an interaction to investigate, not a demonstrated explanation of the observed positions; no alignments or depth rerun were used.'}}
with args.output.open('x') as f:json.dump(receipt,f,indent=2);f.write('\n')
print(json.dumps(receipt['counts']))
