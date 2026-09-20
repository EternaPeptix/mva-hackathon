from pathlib import Path
import collections,datetime,gzip,hashlib,json,urllib.parse
root=Path('/Users/jeweled/Documents/mva-hackathon');case=root/'results/y-call-differences-20260920';out=root/'results/y-call-differences-context-20260920.json'
assert not out.exists()
annotation=root/'ref/baseline-20260919/Homo_sapiens.GRCh38.116.chr.gff3.gz'
h=hashlib.sha256(annotation.read_bytes()).hexdigest();assert h=='0857407b4a8be8792089a69fd961737256fb7129b6b799f8a8920ebea4fb134f'
config=json.loads((case/'evidence.json').read_text())['config'];targets=config['targets'];positions=sorted({(t['chrom'],t['pos']) for t in targets});overlaps={f'{c}:{p}':[] for c,p in positions};nearby={f'{c}:{p}':[] for c,p in positions};types=collections.Counter()
with gzip.open(annotation,'rt') as f:
 for line in f:
  if line.startswith('##FASTA'):break
  if line.startswith('#'):continue
  fields=line.rstrip('\n').split('\t')
  if len(fields)!=9 or fields[0] not in {c for c,p in positions}:continue
  c,source,kind,left,right,score,strand,phase,attributes=fields;a,b=int(left),int(right)
  types[kind]+=1
  relevant=[(chrom,pos) for chrom,pos in positions if chrom==c and a<=pos<=b]
  genes_near=[(chrom,pos) for chrom,pos in positions if chrom==c and kind in ('gene','pseudogene','ncRNA_gene') and a-100000<=pos<=b+100000]
  if not relevant and not genes_near:continue
  attrs={urllib.parse.unquote(k):urllib.parse.unquote(v) for k,v in (item.split('=',1) for item in attributes.split(';') if '=' in item)}
  row={'type':kind,'source':source,'start_1based':a,'end_1based':b,'strand':strand,'attributes':attrs}
  for chrom,pos in relevant:overlaps[f'{chrom}:{pos}'].append(row)
  for chrom,pos in genes_near:
   distance=0 if a<=pos<=b else a-pos if pos<a else pos-b
   nearby[f'{chrom}:{pos}'].append(dict(row,distance_bp=distance))
result={'recorded_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'annotation_path':str(annotation),'annotation_sha256':h,'source_url':'https://ftp.ensembl.org/pub/release-116/gff3/homo_sapiens/Homo_sapiens.GRCh38.116.chr.gff3.gz','config_sha256':hashlib.sha256((root/'results/y-call-differences-config-20260920.json').read_bytes()).hexdigest(),'coordinate_model':'GFF3 one-based inclusive interval overlap; ±100kb gene-span context is descriptive, not a regulatory assignment','overlaps':overlaps,'nearby_gene_spans':nearby,'chromosome_feature_type_counts':dict(types),'scope':'Annotation overlap only. No consequence prediction, pathogenicity classification, rare-variant claim or causal nomination. Lack of annotated overlap is not absence of function.'}
out.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'path':str(out),'sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'overlapping_feature_types':{k:dict(collections.Counter(x['type'] for x in v)) for k,v in overlaps.items()},'nearby_genes':{k:[{'name':x['attributes'].get('Name'),'id':x['attributes'].get('ID'),'distance':x['distance_bp'],'type':x['type']} for x in v] for k,v in nearby.items()}}))
