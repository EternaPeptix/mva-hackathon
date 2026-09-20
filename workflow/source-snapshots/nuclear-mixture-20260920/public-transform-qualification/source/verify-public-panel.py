#!/usr/bin/env python3
"""Independent public-panel original/ordinal/literal/index reconciliation."""
import argparse,datetime,gzip,hashlib,json
from collections import Counter,defaultdict
from decimal import Decimal
from pathlib import Path

PUBLIC={'panel':'18b0420da6c4e8b24dda55bf4ca4a936a37deb6b60f1e906ca696a2f30b988b3','panel_index':'b9d9411a21a1690e56329c9923781f96618182f7c1b769652cc631e1583b05d7','reference':'187917c00fa7373a054f9768009d19e4f6e5fc87d3a6c125ce131164ca69fcd6','reference_index':'d81f672f40cfabf2510f37e3f790f2c79665023f4ac950bb9cc6496a40cebf32'}
def require(v,m):
    if not v:raise ValueError(m)
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()
def json_value(text):
    def pairs(items):
        d={}
        for k,v in items:require(k not in d,'duplicate JSON key');d[k]=v
        return d
    def bad(v):raise ValueError('nonfinite JSON '+v)
    return json.loads(text,object_pairs_hook=pairs,parse_constant=bad)
def compare_records(original,ledger,derived,fetch,expected_records,source_sha):
    selected=iter(derived);records=iter(ledger);counts=Counter();bychrom=Counter();seen=set();digests=defaultdict(hashlib.sha256)
    for ordinal,line in enumerate(original,1):
        f=line.rstrip('\r\n').split('\t');require(len(f)==8,'source columns');chrom,pos,_,ref,alt,_,_,info=f;pos=int(pos)
        require(chrom in set(map(str,range(1,23))) and pos>0,'autosomal source coordinate')
        require(len(ref)==len(alt)==1 and ref in 'ACGT' and alt in 'ACGT' and ref!=alt,'concrete biallelic SNP')
        require((chrom,pos) not in seen,'source duplicate');seen.add((chrom,pos));require(fetch(chrom,pos-1,pos).upper()==ref,'original REF mismatch')
        attrs={}
        for x in info.split(';'):
            k,sep,v=x.partition('=');require(k not in attrs,'source duplicate INFO');attrs[k]=v if sep else None
        literal=attrs['POP_AF'];value=Decimal(literal);require(value.is_finite() and 0<=value<=1,'source frequency invalid')
        reason='included' if Decimal('.01')<value<Decimal('.2') else 'at_lower_endpoint' if value==Decimal('.01') else 'at_upper_endpoint' if value==Decimal('.2') else 'below_lower_endpoint' if value<Decimal('.01') else 'above_upper_endpoint'
        digest=hashlib.sha256(line.encode()).hexdigest();record=json_value(next(records));want=dict(source_sha256=source_sha,source_ordinal=ordinal,source_record_sha256=digest,chrom=chrom,pos=pos,ref=ref,alt=alt,population_ALT_frequency_literal=literal,selection=reason)
        require(type(record.get('source_ordinal')) is int and type(record.get('pos')) is int and record==want,'full source ledger differs')
        counts['source_records']+=1;counts[reason]+=1
        if reason=='included':
            expected='\t'.join(f[:7]+[info+';AF='+literal+';NMP_ORD='+str(ordinal)+';NMP_SHA256='+digest])+'\n'
            require(next(selected,None)==expected,'selected identity/orientation/original literal differs');bychrom[chrom]+=1;digests[chrom].update(expected.encode())
    require(counts['source_records']==expected_records and counts['included']>0,'source or selected count differs')
    require(next(selected,None) is None and next(records,None) is None,'extra derived or ledger rows')
    return dict(counts=dict(counts),included_by_chromosome=dict(bychrom),indexed_record_sha256={k:v.hexdigest() for k,v in digests.items()})
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',required=True,type=Path);p.add_argument('--expected-plan-sha256',required=True);p.add_argument('--run',required=True,type=Path);p.add_argument('--output',required=True,type=Path);a=p.parse_args()
    require(not a.output.exists(),'verification output exists');require(sha(a.plan)==a.expected_plan_sha256,'plan SHA');plan=json_value(a.plan.read_text());require(plan['case_execution_authorized'] is False,'case forbidden');root=Path(plan['root']).resolve()
    require(a.run.resolve()==root/'results/nuclear-mixture-pilot-20260920/public-panel-v1/derivative','derivative scope')
    require(a.output.resolve()==root/'results/nuclear-mixture-pilot-20260920/public-panel-v1/verification','verification scope')
    summary=json_value((a.run/'summary.json').read_text());require(summary['status']=='complete','producer incomplete');require(set(summary['inputs'])==set(PUBLIC),'public source membership')
    verified=[]
    def pin(path,expected=None):
        path=Path(path);before=path.stat();digest=sha(path);after=path.stat();require((before.st_size,before.st_mtime_ns)==(after.st_size,after.st_mtime_ns),'input unstable')
        if expected is not None:require(digest==expected,'file digest differs '+str(path))
        v=dict(path=str(path),bytes=after.st_size,mtime_ns=after.st_mtime_ns,sha256=digest);verified.append(v);return v
    for key,expected in PUBLIC.items():
        src=summary['inputs'][key];require(Path(src['path']).resolve().is_relative_to(root/'ref'),'nonpublic input');require(src['sha256']==expected,'source digest declaration');pin(src['path'],expected)
    require(summary['plan']['sha256']==a.expected_plan_sha256,'producer plan binding')
    require(set(summary['outputs'])=={'panel.vcf','panel.vcf.gz','panel.vcf.gz.tbi','source-ordinal-ledger.jsonl.gz'},'output membership')
    for name,value in summary['outputs'].items():
        observed=pin(a.run/name,value['sha256']);require(observed['bytes']==value['bytes'],'output bytes')
    import pysam
    with pysam.FastaFile(summary['inputs']['reference']['path']) as ref,gzip.open(summary['inputs']['panel']['path'],'rt',newline='') as original,gzip.open(a.run/'panel.vcf.gz','rt',newline='') as derived,gzip.open(a.run/'source-ordinal-ledger.jsonl.gz','rt') as ledger:
        result=compare_records((x for x in original if not x.startswith('#')),ledger,(x for x in derived if not x.startswith('#')),ref.fetch,224487,PUBLIC['panel'])
    require(result['counts']==summary['counts'] and result['included_by_chromosome']==summary['included_by_chromosome'],'producer totals differ')
    h=hashlib.sha256()
    with gzip.open(a.run/'panel.vcf.gz','rb') as stream:
        for chunk in iter(lambda:stream.read(8*1024*1024),b''):h.update(chunk)
    require(h.hexdigest()==sha(a.run/'panel.vcf'),'compressed/raw payload differs')
    with pysam.TabixFile(str(a.run/'panel.vcf.gz'),index=str(a.run/'panel.vcf.gz.tbi')) as index:
        require(set(index.contigs)==set(result['included_by_chromosome']),'index contig scope')
        for chrom,count in result['included_by_chromosome'].items():
            h=hashlib.sha256();n=0
            for line in index.fetch(chrom):h.update((line+'\n').encode());n+=1
            require(n==count and h.hexdigest()==result['indexed_record_sha256'][chrom],'indexed full payload differs')
    for src in verified:
        current=Path(src['path']).stat();require((current.st_size,current.st_mtime_ns)==(src['bytes'],src['mtime_ns']),'verified file changed')
    pin(a.plan,a.expected_plan_sha256);pin(a.run/'summary.json');pin(Path(__file__))
    result.update(status='complete',completed_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),verified_files=verified,scope='Independent original/ordinal/REF/ALT/INFO-literal/selection/index reconciliation; same pysam/HTSlib I/O engine, no producer import. No patient input or estimator execution.')
    a.output.mkdir(parents=True)
    with (a.output/'verification.json').open('x') as f:f.write(json.dumps(result,indent=2,sort_keys=True)+'\n')
if __name__=='__main__':main()
