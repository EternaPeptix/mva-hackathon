from pathlib import Path
import csv, datetime, hashlib, json, os, shutil, subprocess
import pysam

ROOT=Path('/Users/jeweled/Documents/mva-hackathon')
SRC=ROOT/'methods/mutserve-qualification-20260920'
OUT=ROOT/'results/mutserve-contract-controls-20260920'
JAR=ROOT/'tools/mutserve-2.0.0-rc12/mutserve.jar'
HAPLO=ROOT/'tools/haplocheck-1.3.3/haplocheck.jar'
JAVA=ROOT/'tools/bio/lib/jvm/bin/java'
JAVAC=JAVA.with_name('javac')
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,x): Path(p).write_text(json.dumps(x,indent=2)+'\n')
def stamp(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
assert sha(JAR)=='bb741acf57ed23388082260358aa44794b00006506297b9b115fcc3d58096387'
assert sha(HAPLO)=='2b30f573f98395593076d32d9e377818c3feb4505a08a53151a39656cf075782'
OUT.mkdir(exist_ok=False)
os.nice(10)
commands=[]
def run(cmd,label):
    row=dict(label=label,command=list(map(str,cmd)),started_at=stamp())
    with (OUT/(label+'.log')).open('w') as h:
        r=subprocess.run(row['command'],stdout=h,stderr=subprocess.STDOUT,cwd=OUT)
    row.update(exit_code=r.returncode,finished_at=stamp());commands.append(row);save(OUT/'commands.json',commands)
    if r.returncode: raise RuntimeError(label)
ref=OUT/'reference.fasta'
shutil.copy2(ROOT/'tools/mutserve-2.0.0-rc12/rCRS.fasta',ref)
seq=''.join(x for x in ref.read_text().splitlines() if not x.startswith('>'))
assert len(seq)==16569 and hashlib.sha256(seq.encode()).hexdigest()=='f156ff3f65bbcc80c7ebb9936dceb96b1477b4f8f535c4e1dbe7baea225cbc66'
pysam.faidx(str(ref))
cases=[(309,'blacklist_heteroplasmy',[.75,.25,0]),(1000,'homoplasmic',[0,1,0]),
       (2000,'minor_alt',[.75,.25,0]),(3000,'major_alt',[.25,.75,0]),
       (4000,'one_percent',[.99,.01,0]),(5000,'near_homoplasmic',[.005,.995,0]),
       (6000,'two_nonreference',[0,.75,.25]),(7000,'three_alleles',[.7,.2,.1]),
       (8000,'balanced',[.5,.5,0]),(9000,'reference_only',[1,0,0]),
       (10000,'low_quality_alt',[0,1,0]),(11000,'one_strand_alt',[.75,.25,0])]
header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'chrM','LN':16569}],'RG':[{'ID':'fictional','SM':'fictional'}]}
truth=[]
with pysam.AlignmentFile(str(OUT/'fictional.bam'),'wb',header=header) as dest:
    for pos,label,fractions in cases:
        bases=[seq[pos-1]]+[x for x in 'ACGT' if x!=seq[pos-1]][:2]
        counts={}
        for strand in [0,16]:
            fs=[1,0,0] if label=='one_strand_alt' and strand==16 else fractions
            counts[str(strand)]={base:round(frac*1000) for base,frac in zip(bases,fs)}
            for base,count in counts[str(strand)].items():
                for n in range(count):
                    read=pysam.AlignedSegment(dest.header);read.query_name=f'{label}_{strand}_{base}_{n}'
                    read.flag=strand;read.reference_id=0;read.reference_start=pos-51;read.mapping_quality=60
                    read.cigartuples=[(0,100)];s=list(seq[pos-51:pos+49]);s[50]=base;read.query_sequence=''.join(s)
                    read.query_qualities=[11 if label=='low_quality_alt' else 40]*100
                    read.set_tag('AS',100);read.set_tag('RG','fictional');dest.write(read)
        truth.append(dict(pos=pos,label=label,ref=bases[0],counts_by_flag=counts))
pysam.index(str(OUT/'fictional.bam'))
save(OUT/'truth.json',truth)
run([JAVA,'-XX:ActiveProcessorCount=2','-Xmx2g','-jar',JAR,'call','--level','0.01','--reference',ref,
     '--mapQ','30','--baseQ','20','--alignQ','30','--deletions','--write-raw','--contig-name','chrM',
     '--output',OUT/'calls.vcf.gz','--no-ansi','--threads','2',OUT/'fictional.bam'],'mutserve')
for filename in ['MutserveReadProbe.java','HaplocheckParserProbe.java']: shutil.copy2(SRC/filename,OUT/filename)
run([JAVAC,'-cp',JAR,OUT/'MutserveReadProbe.java'],'compile-read-probe')
run([JAVA,'-XX:ActiveProcessorCount=2','-Xmx2g','-cp',str(JAR)+':'+str(OUT),'MutserveReadProbe',ref,OUT/'read-contract.json'],'read-probe')
run([JAVAC,'-cp',HAPLO,OUT/'HaplocheckParserProbe.java'],'compile-parser-probe')
run([JAVA,'-XX:ActiveProcessorCount=2','-Xmx2g','-cp',str(HAPLO)+':'+str(OUT),'HaplocheckParserProbe',OUT/'calls.vcf.gz',OUT/'parser-observed.json'],'parser-probe')
with pysam.VariantFile(str(OUT/'calls.vcf.gz')) as handle:
    records=[dict(pos=x.pos,ref=x.ref,alts=x.alts,filters=list(x.filter),samples={n:dict(s) for n,s in x.samples.items()}) for x in handle]
with (OUT/'calls.txt').open() as h: table=list(csv.DictReader(h,delimiter='\t'))
save(OUT/'observed.json',dict(completed_at=stamp(),truth=truth,vcf=records,caller_table=table,
    read_contract=json.loads((OUT/'read-contract.json').read_text()),
    parser=json.loads((OUT/'parser-observed.json').read_text()),
    scope='Fictional error-free pileup and native input-format controls, not mapping validation or clinical sensitivity. No patient data.'))
print(json.dumps(dict(completed_at=stamp(),records=len(records),table_rows=len(table),read_contract=json.loads((OUT/'read-contract.json').read_text()))))
