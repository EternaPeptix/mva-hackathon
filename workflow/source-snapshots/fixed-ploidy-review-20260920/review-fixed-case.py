from pathlib import Path
from collections import Counter
import datetime,hashlib,json,sys
import pysam
root=Path('/Users/jeweled/Documents/mva-hackathon');repo=root/'repo'
old=root/'calls/ploidy-aware-sex-20260920';new=root/'calls/fixed-ploidy-sex-20260920'
out=root/'results/fixed-ploidy-review-20260920'
def digest(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def save(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
status=json.loads((new/'status.json').read_text());assert status['status']=='complete'
for name,h in status['output_sha256'].items():assert digest(new/name)==h,name
prior=json.loads((repo/'results/ploidy-aware-sex-original-call-audit-20260920.json').read_text())
for name,row in prior['preserved_original_files'].items():assert digest(old/name)==row['sha256'],name
assert digest(new/'bcftools.original.vcf.gz')==digest(old/'bcftools.original.vcf.gz')
out.mkdir(exist_ok=False)
sample=status['config']['sample'];regions=status['config']['regions']
sys.path.insert(0,str(repo/'workflow/scripts'))
from call_ploidy_regions import inspect_calls
counts={}
def original(path):
 rows={};shapes=Counter();called=Counter()
 with pysam.VariantFile(str(path)) as stream:
  assert list(stream.header.samples)==[sample]
  for r in stream:
   key=(r.contig,r.pos,r.ref,tuple(r.alts or []));assert key not in rows,key
   values=r.samples[sample];gt=values.get('GT')
   assert gt and None not in gt
   row={'chrom':r.contig,'pos':r.pos,'ref':r.ref,'alts':list(r.alts or []),'filters':list(r.filter),'GT':list(gt),'phased':values.phased,'DP':values.get('DP'),'AD':values.get('AD'),'GQ':values.get('GQ'),'PL':values.get('PL')}
   rows[key]=row;shapes[str(len(gt))]+=1
   called['retained_ALT_entries']+=len(r.alts or [])
   called['original_GT_positive_ALT_entries']+=sum(i in gt for i in range(1,len(r.alts or [])+1))
 return rows,{'records':len(rows),'GT_lengths':dict(shapes),'allele_entry_counts':dict(called)}
a,counts['original_gatk']=original(old/'gatk.original.vcf.gz')
b,counts['fixed_gatk']=original(new/'gatk.original.vcf.gz')
c,counts['bcftools']=original(new/'bcftools.original.vcf.gz')
shared=a.keys()&b.keys();only_old=a.keys()-b.keys();only_new=b.keys()-a.keys()
fields=['GT','phased','DP','AD','GQ','PL','filters'];changed=[];field_counts=Counter()
for k in sorted(shared):
 differing=[field for field in fields if a[k][field]!=b[k][field]]
 if differing:
  field_counts.update(differing);changed.append({'changed_fields':differing,'original':a[k],'fixed':b[k]})
assert len(shared)+len(only_old)==len(a) and len(shared)+len(only_new)==len(b)
save(out/'raw-record-differences.json',{'shared_changes':changed,'original_only':[a[k] for k in sorted(only_old)],'fixed_only':[b[k] for k in sorted(only_new)],'scope':'Exact original-record key comparison. Raw representation changes and absence require interpretation; absence never establishes reference genotype.'})
def identities(path):
 keys=set();n=0;symbolic=0
 with pysam.VariantFile(str(path)) as stream:
  for r in stream:
   assert len(r.alts or [])==1
   n+=1
   if any(any(c not in 'ACGT' for c in allele) for allele in r.alleles):symbolic+=1;continue
   keys.add((r.contig,r.pos,r.ref,r.alts[0]))
 return keys,{'normalized_rows':n,'unique_concrete_identities':len(keys),'nonconcrete_rows':symbolic}
g,gc=identities(new/'gatk.normalized.vcf.gz');v,vc=identities(new/'bcftools.normalized.vcf.gz')
with (out/'normalized-identity-presence.tsv').open('w') as f:
 f.write('chrom\tpos\tref\talt\tgatk_retained\tbcftools_retained\n')
 for key in sorted(g|v):f.write('\t'.join(map(str,(*key,int(key in g),int(key in v))))+'\n')
site=('X',2781514,'C',('A',));assert site in a and site in b
assert a[site]['GT']==[1,1] and b[site]['GT']==[1]
result={'verified_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'host':'512S2','status':'complete','fixed_run_status_sha256':digest(new/'status.json'),'original_audit_sha256':digest(repo/'results/ploidy-aware-sex-original-call-audit-20260920.json'),'script_sha256':digest(Path(__file__)),'original_files_preserved':len(prior['preserved_original_files']),'fixed_outputs_verified':len(status['output_sha256']),'bcftools_original_reused_byte_for_byte':True,'raw_record_accounting':counts,'fixed_ploidy_accounting':{label:inspect_calls(new/(label+'.original.vcf.gz'),regions,sample) for label in ['gatk','bcftools']},'old_vs_fixed_gatk':{'shared_exact_record_keys':len(shared),'original_only_record_keys':len(only_old),'fixed_only_record_keys':len(only_new),'changed_shared_records':len(changed),'changed_field_counts':dict(field_counts)},'known_boundary_site':{'original':a[site],'fixed':b[site]},'retained_normalized_identity_presence':{'gatk':gc,'bcftools':vc,'shared':len(g&v),'gatk_only':len(g-v),'bcftools_only':len(v-g)},'limitations':['Full raw original genotype and ALT contexts remain authoritative.','Normalized identity presence is an inventory that may include uncalled retained ALT entries; it is not genotype concordance, precision, sensitivity or causal evidence.','Old/new record or identity absence is not reference genotype or callability.','Computational ploidy assumptions do not measure a karyotype or exclude mosaicism.','All differences require biological/quality interpretation before candidate integration.']}
save(out/'summary.json',result)
files=[{'path':str(p.relative_to(root)),'bytes':p.stat().st_size,'sha256':digest(p)} for p in sorted(out.iterdir()) if p.is_file()]
save(out/'verification.json',{'verified_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'files':files,'fixed_outputs':status['output_sha256'],'script_sha256':digest(Path(__file__))})
print(json.dumps(result))
