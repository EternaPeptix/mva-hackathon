#!/usr/bin/env python3
"""Replay a pinned, descriptive GSE132397 Bub1b subset; no download or DEG inference.

Requires the original three gz inputs. Reports all 3 exact GPL1261 features x
8 GSM samples, native detection calls/P, arithmetic means, sample SD and ratios.
Detection P is a within-sample call statistic, not a between-group DE P or FDR.
"""
import argparse,csv,gzip,hashlib,io,json,statistics
from decimal import Decimal
from pathlib import Path

PINNED={
 'matrix':'15e2053fe4e5c250b2ac1338517cfcc114b276bfd3bd67ce04654a7a617df5f2',
 'soft':'f6efca0ff9c1776aabe64e434645cf87c86355537882a8035b15d6d5ea5dd5f4',
 'annotation':'394c103cf4ca99500d06369e4fb2b9b9bffc2c10fd3cbeb0b0d64c51c3fc95f8'}
SAMPLES=[f'GSM{i}' for i in range(3861604,3861612)]
PROBES={'1416961_at','1447362_at','1447363_s_at'}
TITLES=[f'P35-45 3DIV ONCulture {g}_{i}' for g in ['Control','Lithium'] for i in range(1,5)]

def need(ok,message):
 if not ok:raise ValueError(message)

def read_pinned(path,label):
 raw=Path(path).read_bytes();digest=hashlib.sha256(raw).hexdigest()
 need(digest==PINNED[label],f'{label}: SHA256 differs from pinned official input')
 return gzip.decompress(raw).decode('utf-8')

def summarize(matrix,soft,annotation):
 ann=read_pinned(annotation,'annotation');mt=read_pinned(matrix,'matrix');st=read_pinned(soft,'soft')
 lines=ann.splitlines();start=next(i for i,l in enumerate(lines) if l.startswith('ID\t'))
 table=lines[start:];end=next((i for i,l in enumerate(table) if l.startswith('!platform_table_end')),len(table));selected={}
 for row in csv.DictReader(table[:end],delimiter='\t'):
  if row.get('Gene symbol','').strip()=='Bub1b':
   need(row['Gene ID'].strip()=='12236','Bub1b symbol has unexpected GeneID')
   probe=row['ID'];need(probe not in selected,'Duplicate exact-symbol annotation feature')
   selected[probe]={'gene_symbol':'Bub1b','entrez_gene_id':'12236','GenBank_accession':row['GenBank Accession']}
 need(set(selected)==PROBES,'Exact-symbol annotation does not yield the three pinned Bub1b probes')
 parsed=list(csv.reader(io.StringIO(mt),delimiter='\t'))
 titles=next(r[1:] for r in parsed if r and r[0]=='!Sample_title')
 accessions=next(r[1:] for r in parsed if r and r[0]=='!Sample_geo_accession')
 need(accessions==SAMPLES and titles==TITLES,'Unexpected matrix sample order/titles')
 begin=next(i for i,r in enumerate(parsed) if r and r[0]=='ID_REF');need(parsed[begin][1:]==SAMPLES,'Matrix table header differs')
 matrix_rows={}
 for row in parsed[begin+1:]:
  if not row:continue
  if row[0].startswith('!series_matrix_table_end'):break
  if row[0] in PROBES:
   need(row[0] not in matrix_rows and len(row)==9,'Duplicate/malformed matrix probe row')
   matrix_rows[row[0]]=row[1:]
 need(set(matrix_rows)==PROBES,'Matrix target coverage differs')
 soft_rows={};seen=[];definitions={}
 for block in st.split('^SAMPLE = ')[1:]:
  sample=block.split('\n',1)[0].strip()
  if sample not in SAMPLES:continue
  need(sample not in seen,'Duplicate SOFT sample');seen.append(sample)
  prefix,body=block.split('!sample_table_begin\n',1)
  title=next(l.split(' = ',1)[1] for l in prefix.splitlines() if l.startswith('!Sample_title = '))
  need(title==TITLES[SAMPLES.index(sample)],'SOFT title differs from matrix')
  for line in prefix.splitlines():
   if line.startswith('#VALUE =') or line.startswith('#ABS_CALL =') or line.startswith('!Sample_data_processing ='):
    k,v=line.split(' =',1);v=v.strip()
    need(k not in definitions or definitions[k]==v,'Inconsistent native field definitions');definitions[k]=v
  content=body.split('!sample_table_end',1)[0]
  reader=csv.DictReader(io.StringIO(content),delimiter='\t')
  need(reader.fieldnames==['ID_REF','VALUE','ABS_CALL','DETECTION P-VALUE'],'Unexpected SOFT table schema')
  sample_probes=set()
  for row in reader:
   probe=row['ID_REF']
   if probe not in PROBES:continue
   need(probe not in sample_probes,'Duplicate SOFT target row');sample_probes.add(probe)
   need(row['ABS_CALL'] in ['P','A','M'],'Unknown detection-call code')
   need(Decimal(row['VALUE'])==Decimal(matrix_rows[probe][SAMPLES.index(sample)]),'Matrix/SOFT target value mismatch')
   need(Decimal(0)<=Decimal(row['DETECTION P-VALUE'])<=Decimal(1),'Invalid detection P')
   soft_rows[(probe,sample)]=row
  need(sample_probes==PROBES,'Incomplete SOFT sample target coverage')
 need(seen==SAMPLES and len(soft_rows)==24,'SOFT eight-sample/24-cell coverage differs')
 need(definitions.get('#VALUE')=='MAS5.0 signal','Native signal unit differs')
 rows=[]
 for probe in sorted(PROBES):
  values=[Decimal(v) for v in matrix_rows[probe]];cm=sum(values[:4])/4;lm=sum(values[4:])/4
  calls=[soft_rows[(probe,s)]['ABS_CALL'] for s in SAMPLES]
  rows.append({'probe':probe,**selected[probe],
   'control_values':[float(v) for v in values[:4]],'lithium_values':[float(v) for v in values[4:]],
   'native_value_text':matrix_rows[probe],
   'control_mean':float(cm),'lithium_mean':float(lm),
   'control_sample_SD':statistics.stdev(float(v) for v in values[:4]),
   'lithium_sample_SD':statistics.stdev(float(v) for v in values[4:]),
   'ratio_lithium_mean_to_control_mean':float(lm/cm),
   'all_lithium_values_below_all_control_values':max(values[4:])<min(values[:4]),
   'control_detection_calls':calls[:4],'lithium_detection_calls':calls[4:],
   'control_detection_p_values':[soft_rows[(probe,s)]['DETECTION P-VALUE'] for s in SAMPLES[:4]],
   'lithium_detection_p_values':[soft_rows[(probe,s)]['DETECTION P-VALUE'] for s in SAMPLES[4:]],
   'gene_specific_differential_expression_p':None,'gene_specific_fdr':None})
 return {'dataset':'GSE132397','pmid':'31477687','platform':'GPL1261','gene':'Bub1b','entrez_gene_id':'12236',
  'source_sha256':PINNED,'controls':SAMPLES[:4],'lithium':SAMPLES[4:],'sample_titles':TITLES,
  'unit':'Deposited processed MAS5.0 signal, not counts/FPKM/TPM or protein',
  'native_field_definitions':definitions,'rows':rows,'matrix_soft_values_matched':24,
  'normalization_limit':'No raw CEL/CHP processing or paper normalization/filter pipeline reproduced.',
  'inference_limit':'Detection P is within-sample detection, not between-group DE significance. The all-absent 1447362_at feature is background-limited. Lower deposited probe signals are not a validated RNA-abundance decrease, protein/function result or three independent studies.',
  'biological_unit_limit':'Four arrays/group; paper/deposit mouse strain, pooling/donor mapping and cRNA-input descriptions differ. Arrays/probes are not automatically independent animals.',
  'inferential_analysis_performed':False,'new_reading_credit':0}

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--matrix',required=True,type=Path);p.add_argument('--soft',required=True,type=Path);p.add_argument('--annotation',required=True,type=Path);p.add_argument('--output',type=Path)
 a=p.parse_args()
 try:result=summarize(a.matrix,a.soft,a.annotation)
 except (ValueError,KeyError,IndexError,StopIteration) as exc:p.exit(1,f'ERROR: {exc}\n')
 out=json.dumps(result,ensure_ascii=False,indent=2)+'\n'
 if a.output:a.output.write_text(out)
 else:print(out,end='')
if __name__=='__main__':main()
