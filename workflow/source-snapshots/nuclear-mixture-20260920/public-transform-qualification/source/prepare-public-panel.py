#!/usr/bin/env python3
"""Public reference-only POP_AF -> AF derivative; never reads case data."""
import argparse
from collections import Counter
from decimal import Decimal, InvalidOperation
import datetime
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import time

import pysam

EXPECTED_PANEL='18b0420da6c4e8b24dda55bf4ca4a936a37deb6b60f1e906ca696a2f30b988b3'
EXPECTED_REFERENCE='187917c00fa7373a054f9768009d19e4f6e5fc87d3a6c125ce131164ca69fcd6'


def require(value,message):
    if not value:raise ValueError(message)


def fingerprint(path,expected=None):
    path=Path(path).resolve();old=path.stat();h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    new=path.stat();require((old.st_size,old.st_mtime_ns)==(new.st_size,new.st_mtime_ns),'input changed while hashing')
    pin=dict(path=str(path),bytes=new.st_size,mtime_ns=new.st_mtime_ns,sha256=h.hexdigest())
    if expected is not None:require(pin['sha256']==expected,'input hash mismatch')
    return pin


def prepare(plan_path,plan_sha,output):
    require(not output.exists(),'existing output preserved')
    plan_pin=fingerprint(plan_path,plan_sha);plan=json.loads(plan_path.read_text())
    require(plan['scope']=='public_marker_derivation_and_fictional_controls_only','wrong scope')
    require(plan['case_execution_authorized'] is False,'case access forbidden')
    root=Path(plan['root']).resolve();panel=plan['public_panel'];reference=plan['reference']
    require(panel['sha256']==EXPECTED_PANEL and reference['sha256']==EXPECTED_REFERENCE,'only declared public resources permitted')
    sources={}
    for name,path,digest in [('panel',panel['path'],panel['sha256']),('panel_index',panel['index_path'],panel['index_sha256']),('reference',reference['path'],reference['sha256']),('reference_index',reference['index_path'],reference['index_sha256'])]:
        p=(root/path).resolve();require(p.is_relative_to(root/'ref'),'public reference directory required');sources[name]=fingerprint(p,digest)
    lo,hi=Decimal(panel['minimum_exclusive']),Decimal(panel['maximum_exclusive'])
    require((lo,hi)==(Decimal('.01'),Decimal('.2')),'unreviewed AF limits')
    require(panel['frequency_source']=='POP_AF' and panel['frequency_target']=='AF','wrong AF translation')
    output.mkdir(parents=True);started=time.monotonic();counts=Counter();seen=set();bychrom=Counter()
    raw=output/'panel.vcf';ledger=output/'source-ordinal-ledger.jsonl.gz';source_ordinal=0
    try:
        with pysam.FastaFile(sources['reference']['path']) as fasta,gzip.open(sources['panel']['path'],'rt',newline='') as source,raw.open('x') as writer,gzip.open(ledger,'wt') as records:
            have_header=False
            for line in source:
                if line.startswith('##'):
                    require(not line.startswith('##INFO=<ID=AF,') and not line.startswith('##INFO=<ID=NMP_'),'reserved derivative INFO field')
                    writer.write(line);continue
                if line.startswith('#CHROM'):
                    require(not have_header,'duplicate VCF header');require(len(line.rstrip('\r\n').split('\t'))==8,'public sites-only input required')
                    writer.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Exact original public POP_AF for the unchanged ALT allele">\n')
                    writer.write('##INFO=<ID=NMP_ORD,Number=1,Type=Integer,Description="One-based ordinal in frozen ROH public panel">\n')
                    writer.write('##INFO=<ID=NMP_SHA256,Number=1,Type=String,Description="SHA256 of exact original panel record bytes including newline">\n')
                    writer.write(line);have_header=True;continue
                require(have_header and not line.startswith('#'),'missing or malformed VCF header')
                source_ordinal+=1;fields=line.rstrip('\r\n').split('\t');require(len(fields)==8,'unexpected genotype/sample columns')
                chrom,pos,identity,ref,alt,qual,filt,info=fields;position=int(pos)
                require(chrom in set(map(str,range(1,23))) and position>0,'nonautosomal marker')
                require(len(ref)==len(alt)==1 and ref in 'ACGT' and alt in 'ACGT' and ref!=alt,'not concrete biallelic SNP')
                require((chrom,position) not in seen,'duplicate source locus');seen.add((chrom,position))
                require(fasta.fetch(chrom,position-1,position).upper()==ref,'source REF does not match exact reference')
                attrs={}
                for field in info.split(';'):
                    k,sep,v=field.partition('=');require(k not in attrs,'duplicate INFO field');attrs[k]=v if sep else None
                require('POP_AF' in attrs and attrs['POP_AF'] is not None,'POP_AF missing')
                value=Decimal(attrs['POP_AF']);require(value.is_finite() and 0<=value<=1,'invalid population ALT frequency')
                reason='included' if lo<value<hi else 'at_lower_endpoint' if value==lo else 'at_upper_endpoint' if value==hi else 'below_lower_endpoint' if value<lo else 'above_upper_endpoint'
                digest=hashlib.sha256(line.encode()).hexdigest();counts[reason]+=1;counts['source_records']+=1
                records.write(json.dumps(dict(source_sha256=EXPECTED_PANEL,source_ordinal=source_ordinal,source_record_sha256=digest,chrom=chrom,pos=position,ref=ref,alt=alt,population_ALT_frequency_literal=attrs['POP_AF'],selection=reason),sort_keys=True)+'\n')
                if reason=='included':
                    fields[7]=info+';AF='+attrs['POP_AF']+';NMP_ORD='+str(source_ordinal)+';NMP_SHA256='+digest
                    writer.write('\t'.join(fields)+'\n');bychrom[chrom]+=1
            require(have_header,'missing VCF header')
        require(source_ordinal==panel['expected_records']==224487,'source record total differs')
        require(counts['included']>0,'empty eligible public panel')
        compressed=output/'panel.vcf.gz';pysam.tabix_compress(str(raw),str(compressed),force=False);pysam.tabix_index(str(compressed),preset='vcf',force=False)
        # Independently reread source, ledger and derivative in lockstep; compare exact identities/literals.
        got=0;replayed=0
        with gzip.open(sources['panel']['path'],'rt',newline='') as original,gzip.open(compressed,'rt',newline='') as derived,gzip.open(ledger,'rt') as ledger_stream:
            selected=(line for line in derived if not line.startswith('#'))
            for line in original:
                if line.startswith('#'):continue
                replayed+=1;cols=line.rstrip('\r\n').split('\t');info=dict(field.split('=',1) for field in cols[7].split(';') if '=' in field)
                record=json.loads(next(ledger_stream));digest=hashlib.sha256(line.encode()).hexdigest()
                require(record['source_ordinal']==replayed and record['source_record_sha256']==digest and record['source_sha256']==EXPECTED_PANEL,'ledger original identity differs')
                require([record[k] for k in ['chrom','pos','ref','alt']]==[cols[0],int(cols[1]),cols[3],cols[4]],'ledger allele differs')
                require(record['population_ALT_frequency_literal']==info['POP_AF'],'ledger frequency literal differs')
                eligible=Decimal('.01')<Decimal(info['POP_AF'])<Decimal('.2')
                require((record['selection']=='included')==eligible,'ledger selected membership differs')
                if not eligible:continue
                emitted=next(selected,None);require(emitted is not None,'selected source omitted')
                want=cols[:7]+[cols[7]+';AF='+info['POP_AF']+';NMP_ORD='+str(replayed)+';NMP_SHA256='+digest]
                require(emitted.rstrip('\r\n').split('\t')==want,'selected original fields/INFO/orientation/literal differ');got+=1
            require(next(selected,None) is None and next(ledger_stream,None) is None,'extra derivative or ledger record')
        require(got==counts['included'] and replayed==counts['source_records'],'serialized derivative/source count differs')
        for pin in sources.values():
            st=Path(pin['path']).stat();require((st.st_size,st.st_mtime_ns)==(pin['bytes'],pin['mtime_ns']),'public input changed')
        result=dict(status='complete',completed_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),elapsed_seconds=time.monotonic()-started,scope=plan['scope'],source_ordinal_definition=panel['source_ordinal_definition'],strict_ALT_AF_bounds=[str(lo),str(hi)],counts=dict(counts),included_by_chromosome=dict(bychrom),inputs=sources,plan=plan_pin,source=fingerprint(__file__),outputs={p.name:fingerprint(p) for p in sorted(output.iterdir()) if p.is_file()},runtime=dict(pysam=pysam.__version__),limitations=['Public-only resource transformation; no case genotype or alignment accessed.','ROH panel ordinal is not the original full IGSR VCF ordinal.','No estimator, ancestry inference, linkage-independence proof or clinical performance qualification.'])
        (output/'summary.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
        return result
    except BaseException as e:
        (output/'failure.json').write_text(json.dumps(dict(status='failed',error_type=type(e).__name__,error=str(e)),indent=2)+'\n')
        raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);p.add_argument('--expected-plan-sha256',required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();os.nice(10);r=prepare(a.plan.resolve(),a.expected_plan_sha256,a.output.resolve());print(json.dumps(r['counts'],sort_keys=True))
