#!/usr/bin/env python3
"""Generate finite fictional inputs only. Does not access case files or run tools."""
from decimal import Decimal,ROUND_HALF_UP
from fractions import Fraction
import hashlib,json,random
from pathlib import Path

SAMPLE='FICTIONAL_NUCLEAR_MIXTURE'
TABLE_COLUMNS=['contig','position','ref_count','alt_count','other_alt_count','allele_frequency']
SCENARIOS=[dict(name='clean_d100',depth=100,alpha='0',kind='diploid'),dict(name='clean_d400',depth=400,alpha='0',kind='diploid')]+[dict(name='mixture_'+str(a)+'_d'+str(d),depth=d,alpha=str(a),kind='diploid') for a in ['0.02','0.10'] for d in [100,400]]+[dict(name='same_donor_010_d100',depth=100,alpha='0.10',kind='same_donor'),dict(name='imbalance_one_third',depth=300,alpha='0',kind='imbalance1'),dict(name='imbalance_two_thirds',depth=300,alpha='0',kind='imbalance2'),dict(name='regional_copy_depth',depth=300,alpha='0',kind='cn_region'),dict(name='regional_LOH',depth=100,alpha='0',kind='loh'),dict(name='empty',depth=100,alpha='0',kind='empty'),dict(name='all_heterozygous',depth=100,alpha='0',kind='uninformative')]

def require(v,m):
    if not v:raise ValueError(m)

def save(path,data):
    with Path(path).open('x') as f:f.write(json.dumps(data,indent=2,sort_keys=True,allow_nan=False)+'\n')

def round_half_up(value):
    return int((Decimal(value.numerator)/Decimal(value.denominator)).quantize(Decimal('1'),rounding=ROUND_HALF_UP))

def genotype_grid(seed):
    donor=[0]*81+[1]*18+[2]
    pairs=[(a,b) for a in donor for b in donor]
    random.Random(seed).shuffle(pairs)
    return pairs

def table_rows(spec,seed=20260920):
    pairs=genotype_grid(seed)
    if spec['kind']=='empty':return [],[]
    rows=[];truth=[]
    for i,(recipient,donor) in enumerate(pairs):
        kind=spec['kind'];depth=spec['depth'];alpha=Fraction(spec['alpha']);effective_donor=recipient if kind=='same_donor' else donor
        fraction=(1-alpha)*Fraction(recipient,2)+alpha*Fraction(effective_donor,2);copy_number=2
        if kind in ['imbalance1','imbalance2'] and recipient==1:fraction=Fraction(1 if kind=='imbalance1' else 2,3)
        if kind=='cn_region' and i<5000:
            copy_number=3;depth=450
            if recipient==1:fraction=Fraction(1 if i%2 else 2,3)
        if kind=='loh' and i<5000 and recipient==1:fraction=Fraction(i%2,1)
        if kind=='uninformative':fraction=Fraction(1,2)
        exact=fraction*depth;alt=round_half_up(exact);pos=1000+1000*i
        rows.append(dict(contig='fictional1',position=pos,ref_count=depth-alt,alt_count=alt,other_alt_count=0,allele_frequency='0.1'))
        truth.append(dict(position=pos,recipient_ALT_copies=recipient,distinct_donor_ALT_copies=donor,effective_donor_ALT_copies=effective_donor,declared_mixing_fraction=spec['alpha'],constructed_ALT_fraction=str(fraction),copy_number_pattern=copy_number,depth=depth,rounded=exact.denominator!=1,kind=kind))
    return rows,truth

def write_tables(output):
    output=Path(output);output.mkdir()
    manifest=[]
    for spec in SCENARIOS:
        path=output/(spec['name']+'.tsv');rows,truth=table_rows(spec)
        with path.open('x') as f:
            f.write('#<METADATA>SAMPLE='+SAMPLE+'\n');f.write('\t'.join(TABLE_COLUMNS)+'\n')
            for row in rows:f.write('\t'.join(str(row[k]) for k in TABLE_COLUMNS)+'\n')
        truth_path=output/(spec['name']+'-truth.json');save(truth_path,truth)
        manifest.append(spec|dict(path=str(path),truth_path=str(truth_path),markers=len(rows),rounded_sites=sum(t['rounded'] for t in truth),scope='deterministic estimator-table control; not aligned donor reads'))
    save(output/'manifest.json',manifest);return manifest

def validate_site(sequence,pos,ref,alt):
    require(type(pos) is int and 1<=pos<=len(sequence),'site outside reference');require(sequence[pos-1]==ref and len(ref)==len(alt)==1 and ref in 'ACGT' and alt in 'ACGT' and ref!=alt,'invalid reference/ALT SNP')

def decoded_record(r):
    return dict(QNAME=r.query_name,FLAG=r.flag,reference_id=r.reference_id,start0=r.reference_start,MQ=r.mapping_quality,CIGAR=r.cigartuples,mate_reference_id=r.next_reference_id,mate_start0=r.next_reference_start,TLEN=r.template_length,SEQ=r.query_sequence,qualities=None if r.query_qualities is None else list(r.query_qualities),tags=sorted(r.get_tags(with_value_type=True)))

def count_fixture(output):
    import pysam
    from array import array
    output=Path(output);output.mkdir()
    ref=output/'fictional.fa';sequence=('ACGT'*6000)[:24000]
    with ref.open('x') as f:
        for contig,seq in [('fictional1',sequence),('fictional2',sequence)]:
            f.write('>'+contig+'\n')
            for i in range(0,len(seq),60):f.write(seq[i:i+60]+'\n')
    pysam.faidx(str(ref))
    with ref.with_suffix('.dict').open('x') as f:
        f.write('@HD\tVN:1.6\tSO:unsorted\n')
        for contig in ['fictional1','fictional2']:f.write('@SQ\tSN:'+contig+'\tLN:24000\tM5:'+hashlib.md5(sequence.encode()).hexdigest()+'\n')
    header=pysam.AlignmentHeader.from_dict(dict(HD=dict(VN='1.6',SO='coordinate'),SQ=[dict(SN=c,LN=24000) for c in ['fictional1','fictional2']],RG=[dict(ID=r,SM=SAMPLE,PL='ILLUMINA') for r in ['FICTIONAL_RG1','FICTIONAL_RG2']]))
    features=[('concrete_bases',3),('BQ0',1),('BQ30',1),('MQ49',0),('MQ50',1),('MQ255',0),('duplicate',0),('secondary',0),('supplementary',0),('QC_fail',0),('unmapped',0),('missing_BQ',0),('partial_BQ255',1),('unpaired',1),('mate_other_contig',0),('mate_unmapped',1),('overlap_agree',2),('overlap_disagree',1),('deletion_at_site',0),('skipped_reference_CIGAR',0),('zero_site',0),('AF_lower_endpoint',0),('AF_upper_endpoint',0)]
    rows=[];sites=[];expected={};ledger=[]
    def read(name,pos,base,flag=0,mq=60,bq=30,rg='FICTIONAL_RG1',cigar=None):
        r=pysam.AlignedSegment(header);r.query_name='FICTIONAL_'+name;r.flag=flag;r.reference_id=0;r.reference_start=pos-75;r.mapping_quality=mq;r.cigartuples=cigar or [(0,149)]
        # Regular count controls share a149-base query. Deletion/skip controls use148.
        length=148 if cigar else 149;query=list(sequence[r.reference_start:r.reference_start+length]);query[74]=base;r.query_sequence=''.join(query);r.query_qualities=array('B',[bq]*length);r.set_tag('RG',rg)
        if flag&1:r.next_reference_id=0;r.next_reference_start=r.reference_start;r.template_length=-149 if flag&128 else 149
        if flag&4:r.mapping_quality=0;r.cigartuples=None
        return r
    for i,(name,alt_expected) in enumerate(features):
        pos=501+900*i;refbase=sequence[pos-1];alt=next(v for v in 'ACGT' if v!=refbase);third=next(v for v in 'ACGT' if v not in [refbase,alt]);af='0.01' if name=='AF_lower_endpoint' else '0.2' if name=='AF_upper_endpoint' else '0.1'
        validate_site(sequence,pos,refbase,alt);sites.append((pos,name,refbase,alt,af));group=[]
        if name!='zero_site':group.append(read(name+'_baseline',pos,refbase))
        if name=='concrete_bases':
            group += [read(name+'_alt',pos,alt),read(name+'_other',pos,third),read(name+'_N',pos,'N')];alt_expected=1
        elif name in ['overlap_agree','overlap_disagree']:
            group += [read(name+'_pair',pos,alt,flag=99),read(name+'_pair',pos,alt if name=='overlap_agree' else refbase,flag=147,rg='FICTIONAL_RG1')]
        elif name in ['deletion_at_site','skipped_reference_CIGAR']:
            group.append(read(name,pos,alt,cigar=[(0,74),(2 if name=='deletion_at_site' else 3,1),(0,74)]))
        elif name!='zero_site':
            flag={'duplicate':1024,'secondary':256,'supplementary':2048,'QC_fail':512,'unmapped':4,'mate_other_contig':65,'mate_unmapped':73}.get(name,0)
            r=read(name,pos,alt,flag=flag,mq={'MQ49':49,'MQ50':50,'MQ255':255}.get(name,60),bq=0 if name=='BQ0' else 30,rg='FICTIONAL_RG2' if name=='BQ30' else 'FICTIONAL_RG1')
            if name=='missing_BQ':r.query_qualities=None
            if name=='partial_BQ255':q=r.query_qualities;q[74]=255;r.query_qualities=q
            if name=='mate_other_contig':r.next_reference_id=1
            group.append(r)
        for r in group:rows.append(r);ledger.append(dict(feature=name,position=pos,decoded_record=decoded_record(r)))
        if name not in ['zero_site','AF_lower_endpoint','AF_upper_endpoint']:
            expected[pos]=dict(contig='fictional1',position=pos,ref_count=2 if name=='overlap_disagree' else 1,alt_count=alt_expected,other_alt_count=1 if name=='concrete_bases' else 0,allele_frequency='0.1')
    rows.sort(key=lambda r:(r.reference_id,r.reference_start,r.query_name,r.flag))
    bam=output/'fictional.bam';cram=output/'fictional.cram'
    for path,mode in [(bam,'wb'),(cram,'wc')]:
        options={'format_options':[b'version=3.0']} if mode=='wc' else {}
        with pysam.AlignmentFile(str(path),mode,header=header,reference_filename=str(ref),**options) as writer:
            for r in rows:writer.write(r)
        pysam.index(str(path))
    with cram.open('rb') as stream:cram_magic=stream.read(6)
    require(cram_magic==b'CRAM\x03\x00','fictional CRAM must be version3.0')
    # Compare full original payload, permitting no missing-quality recoding or mate loss.
    with pysam.AlignmentFile(str(bam),'rb') as a,pysam.AlignmentFile(str(cram),'rc',reference_filename=str(ref),format_options=[b'decode_md=0']) as b:
        bam_records=[decoded_record(r) for r in a];cram_records=[decoded_record(r) for r in b]
    save(output/'decoded-payloads.json',dict(bam=bam_records,cram=cram_records,equal=bam_records==cram_records))
    vcf=output/'sites.vcf'
    with vcf.open('x') as f:
        f.write('##fileformat=VCFv4.2\n##contig=<ID=fictional1,length=24000>\n##contig=<ID=fictional2,length=24000>\n##INFO=<ID=AF,Number=A,Type=Float,Description="Fictional population AF">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        for pos,name,refbase,alt,af in sites:f.write(f'fictional1\t{pos}\t{name}\t{refbase}\t{alt}\t.\t.\tAF={af}\n')
    pysam.tabix_compress(str(vcf),str(vcf)+'.gz');pysam.tabix_index(str(vcf)+'.gz',preset='vcf')
    missing=output/'missing-AF.vcf';missing.write_text('\n'.join(line if line.startswith('#') else '\t'.join(line.split('\t')[:7]+['.']) for line in vcf.read_text().splitlines() if not line.startswith('##INFO='))+'\n')
    pysam.tabix_compress(str(missing),str(missing)+'.gz');pysam.tabix_index(str(missing)+'.gz',preset='vcf')
    save(output/'expected-counts.json',list(expected.values()));save(output/'read-construction.json',ledger)
    result=dict(reference=str(ref),bam=str(bam),cram=str(cram),cram_version='3.0',cram_magic_hex=cram_magic.hex(),sites=str(vcf)+'.gz',missing_AF_sites=str(missing)+'.gz',expected_counts=list(expected.values()),decoded_payload_equal=bam_records==cram_records,records=len(rows),source_scope='Generated fictional reads only; no donor BAM mixture')
    save(output/'fixture.json',result);return result
