from pathlib import Path
import json,sys,subprocess,hashlib,datetime
root=Path('/Users/jeweled/Documents/mva-hackathon');out=root/'results/mutserve-contract-controls-20260920'
sys.path.insert(0,str(root/'repo/workflow/scripts'))
from prepare_haplocheck_input import qualify,verify_import
x=json.loads((out/'observed.json').read_text());jvm=root/'tools/bio/lib/jvm/bin';jar=root/'tools/haplocheck-1.3.3/haplocheck.jar'
assert x['read_contract']==dict(primary=1,secondary=1,supplementary=1,qc_failed=1,duplicate=0,mq255=1,mq29=0,as29=0,missing_as=1,bq19=0,length25=0,overlapping_pair=2)
v={r['pos']:r for r in x['vcf']};parsed=x['parser']['samples']['fictional.bam']['variants']
assert set(v)=={309,1000,2000,3000,4000,5000,6000,7000,8000}
assert v[1000]['samples']['fictional.bam']['GT']==[1]
assert v[5000]['samples']['fictional.bam']['GT']==[1]
assert parsed['2000']['minorLevel']==.25 and parsed['3000']['majorLevel']==.75
assert parsed['4000']['minorLevel']==.01
assert parsed['6000']['majorLevel']==.75 and parsed['6000']['minorLevel']==.25
assert parsed['7000']['majorLevel']==.8
assert sum(r['Pos']=='7000' for r in x['caller_table'])==2 and len(v[7000]['alts'])==1
assert v[309]['filters']==['PASS']
seq=''.join(l for l in (out/'reference.fasta').read_text().splitlines() if not l.startswith('>'))
results={}
for label,mask in [('supported',False),('hotspots-excluded',True)]:
    q=qualify(out/'calls.vcf.gz',out/'calls.txt',seq,'fictional.bam',out/(label+'.vcf'),mask)
    assert q['accepted_records']==(6 if mask else 7)
    dropped={str(r['pos']):r['reason'] for r in q['excluded_positions']}
    assert dropped['7000']=='multiple_native_rows_at_position' and dropped['6000']=='non_biallelic_concrete_SNP'
    if mask: assert dropped['309']=='explicit_mutserve_source_hotspot_mask'
    cmd=[str(jvm/'java'),'-XX:ActiveProcessorCount=2','-Xmx2g','-cp',str(jar)+':'+str(out),'HaplocheckParserProbe',str(out/(label+'.vcf')),str(out/(label+'-parser.json'))]
    with (out/(label+'-parser.log')).open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
    result=verify_import(json.loads((out/(label+'-parser.json')).read_text()),q)
    assert result['homoplasmic']==2 and result['heteroplasmic']==(4 if mask else 5)
    results[label]={'qualification':q,'import_verification':result,'command':cmd}
files={str(p.relative_to(root)):dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(out.iterdir()) if p.is_file()}
for p in [Path(__file__),root/'repo/workflow/scripts/prepare_haplocheck_input.py']:
    files[str(p.relative_to(root))]=dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
receipt=dict(verified_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),host='512S2',status='qualified_with_explicit_exclusions',fictional_pileup_scenarios=12,read_filter_scenarios=12,selected_native_imports_verified=13,results=results,files=files,limitations=['Error-free directly aligned reads test caller/parser behavior; they do not establish alignment accuracy or a clinical detection limit.','Read counter accepts secondary/supplementary/QC-failed/MQ255/missing-AS reads and counts overlapping mates twice. Explicit input selection is required; depths remain read observations.','Native VCF drops one allele at the three-allele control; all duplicate native-table sites are excluded.','Native heteroplasmy at source-blacklisted position 309 passes; a separate explicit hotspot exclusion is required.','Native 99.5 percent alternate control is classed homoplasmic under the 1 percent caller setting; it is not molecular homogeneity.'])
(out/'qualification-verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:v for k,v in receipt.items() if k not in ['files','results']}))
