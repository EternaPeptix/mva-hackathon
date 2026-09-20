#!/usr/bin/env python3
"""Independent saved-pileup/model-output audit; no producer/helper imports."""
import argparse,csv,datetime,gzip,hashlib,io,json,math,statistics
from collections import Counter
from pathlib import Path
def require(v,m):
    if not v:raise ValueError(m)
def same(a,b,m):require(json.dumps(a,sort_keys=True,allow_nan=False)==json.dumps(b,sort_keys=True,allow_nan=False),m)
def parse(text):
    def pairs(items):
        d={}
        for k,v in items:require(k not in d,'duplicate JSON key');d[k]=v
        return d
    def bad(v):raise ValueError('nonfinite JSON '+v)
    return json.loads(text,object_pairs_hook=pairs,parse_constant=bad)
def fingerprint(path,expected=None):
    p=Path(path).resolve();before=p.stat();h=hashlib.sha256()
    with p.open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
    after=p.stat();require((before.st_size,before.st_mtime_ns)==(after.st_size,after.st_mtime_ns),'input changed during hash')
    if expected is not None:require(h.hexdigest()==expected,'hash differs '+str(p))
    return dict(path=str(p),bytes=after.st_size,mtime_ns=after.st_mtime_ns,sha256=h.hexdigest())
def read_rows(path,columns,sample=None):
    lines=Path(path).read_text().splitlines()
    if sample is not None:require([x for x in lines if x.startswith('#<METADATA>SAMPLE=')]==['#<METADATA>SAMPLE='+sample],'sample metadata differs')
    data=csv.DictReader(io.StringIO('\n'.join(x for x in lines if not x.startswith('#'))),delimiter='\t');require(data.fieldnames==columns,'column schema differs');rows=list(data)
    require(all(set(row)==set(columns) and all(type(v) is str and v!='' for v in row.values()) for row in rows),'malformed table row width/value');return rows
def marker_map(path,count,contigs):
    result={}
    with gzip.open(path,'rt') as f:
        for line in f:
            if line.startswith('#'):continue
            fields=line.rstrip('\n').split('\t');require(len(fields)==8,'panel genotype/sample columns');chrom,pos,_,ref,alt,_,_,info=fields;key=(chrom,int(pos));attrs={}
            for item in info.split(';'):
                name,sep,value=item.partition('=');require(name not in attrs,'duplicate panel INFO');attrs[name]=value if sep else None
            literal=attrs['AF'];number=float(literal)
            require(key not in result and key[0] in contigs and key[1]>0 and len(ref)==len(alt)==1 and ref in 'ACGT' and alt in 'ACGT' and ref!=alt and math.isfinite(number) and .01<number<.2,'invalid fixed panel')
            result[key]=(literal,number)
    require(len(result)==count and {x[0] for x in result}==set(contigs),'fixed marker membership differs');return result
def reconcile(rows,markers):
    by_key={};payload=[];depths=[]
    for row in rows:
        key=(row['contig'],int(row['position']));require(key in markers and key not in by_key,'extra/duplicate emitted marker')
        counts=[int(row[k]) for k in ('ref_count','alt_count','other_alt_count')];require(all(n>=0 for n in counts),'invalid native count')
        af=float(row['allele_frequency']);require(math.isfinite(af) and af==markers[key][1],'AF binary64 identity differs')
        total=sum(counts);value=dict(contig=key[0],position=key[1],counts=counts,total=total,panel_AF_literal=markers[key][0],table_AF_literal=row['allele_frequency']);payload.append(value);by_key[key]=value
        if total>10:depths.append(total)
    mean=median=None
    if depths:
        base=sum(depths)/len(depths);correction=0.0
        for value in depths:correction+=value-base
        mean=base+correction/len(depths);ordered=sorted(depths);n=len(ordered);median=ordered[n//2] if n%2 else (ordered[n//2-1]+ordered[n//2])/2
    membership=[];counts=Counter();per_chrom={chrom:Counter() for chrom,_ in markers}
    for key in markers:
        value=by_key.get(key)
        category='no_emission' if value is None else 'emitted_zero' if value['total']==0 else 'emitted_1_to_10' if value['total']<=10 else 'coverage_retained' if value['total']>.5*median and value['total']<3*mean else 'coverage_excluded'
        counts[category]+=1;per_chrom[key[0]][category]+=1;membership.append(dict(contig=key[0],position=key[1],category=category,total_read_observations=None if value is None else value['total']))
    denominators=dict(panel_sites=len(markers),emitted=len(rows),no_emission=counts['no_emission'],emitted_zero=counts['emitted_zero'],emitted_1_to_10=counts['emitted_1_to_10'],preliminary_gt10=len(depths),coverage_retained=counts['coverage_retained'],coverage_excluded=counts['coverage_excluded'],covered_mean=mean,covered_median=median)
    require(sum(counts.values())==len(markers),'marker denominator lost')
    return denominators,payload,membership,{k:dict(v) for k,v in per_chrom.items()}
def interpret(rows,sample):
    require(len(rows)==1 and rows[0]['sample']==sample,'one exact native sample output required');row=rows[0];value=float(row['contamination']);error=float(row['error'])
    require(math.isfinite(value) and math.isfinite(error) and 0<=value<=1 and error>=0,'invalid numerical model output')
    fallback=value==0 and error==1
    return dict(classification='insufficient_information' if fallback else 'limited_model_output_available',raw_contamination=value,raw_error=error,interpreted_estimate=None if fallback else value,fallback_0_1=fallback,clinical_pass=None,zero_means_absence=False)
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',type=Path,required=True);p.add_argument('--expected-config-sha256',required=True);a=p.parse_args();pins=[]
    def pin(path,expected=None):
        result=fingerprint(path,expected);pins.append(result);return result
    config_pin=pin(a.config,a.expected_config_sha256);cfg=parse(a.config.read_text());require(cfg['ready'] is True and type(cfg['fixture_only']) is bool,'reviewed ready configuration required')
    require(cfg['patient_assay_authorized'] is (not cfg['fixture_only']),'authorization scope differs')
    root=Path(cfg['output']);run=root/'producer';out=root/'verification';require(not out.exists(),'new audit output required')
    status=parse((run/'status.json').read_text());summary=parse((run/'summary.json').read_text());require(status['status']==summary['status']=='complete','producer incomplete');pin(run/'summary.json',status['summary_sha256']);pin(run/'status.json')
    require(summary['config']['sha256']==config_pin['sha256'],'producer configuration binding')
    require(set(summary['inputs'])==set(cfg['inputs']) and set(summary['runtime'])==set(cfg['runtime']),'input/runtime inventory')
    for kind in ['inputs','runtime']:
        for name,expected in cfg[kind].items():
            observed=pin(expected['path'],expected['sha256']);declared=summary[kind][name];require(all(observed[k]==declared[k] for k in ['path','bytes','sha256']),'producer input identity differs')
    require(set(summary['sources'])==set(cfg['source_sha256']),'source inventory omitted')
    for name,digest in cfg['source_sha256'].items():
        observed=pin(Path(__file__).parent/name,digest);require(all(observed[k]==summary['sources'][name][k] for k in ['path','bytes','sha256']),'executed source linkage differs')
    if not cfg['fixture_only']:
        q=cfg['qualification'];pin(q['path'],q['sha256']);controls=parse(Path(q['path']).read_text());require(controls['status']=='complete' and type(controls['tests_run']) is int and controls['tests_run']>=10 and controls['skipped']==[] and all(type(controls[k]) is int and controls[k]==0 for k in ['errors','failures','exit_code']),'qualification not complete')
        same(controls['executed_source_sha256'],cfg['source_sha256'],'qualification/source mismatch')
        same(controls['runtime_sha256'],{k:v['sha256'] for k,v in cfg['runtime'].items()},'qualification/runtime mismatch')
        import sys
        pin(sys.executable,controls['python']['sha256'])
    actual={str(p.relative_to(run)) for p in run.rglob('*') if p.is_file()};require(actual==set(summary['outputs'])|{'summary.json','status.json'},'raw output inventory differs')
    for name,value in summary['outputs'].items():
        observed=pin(run/name,value['sha256']);require(observed['bytes']==value['bytes'],'raw output size differs')
    import pysam
    with Path(cfg['inputs']['alignment']['path']).open('rb') as stream:require(stream.read(6)==b'CRAM\x03\x00','CRAM version changed')
    with pysam.AlignmentFile(cfg['inputs']['alignment']['path'],'rc',reference_filename=cfg['inputs']['reference']['path'],index_filename=cfg['inputs']['alignment_index']['path'],format_options=[b'decode_md=0']) as alignment:header=alignment.header.to_dict()
    same(header,parse((run/'header.json').read_text()),'original header payload differs');groups=header.get('RG',[]);ids=[x.get('ID') for x in groups]
    require(len(groups)==cfg['expected_read_groups'] and len(set(ids))==len(ids) and all(type(x) is str and x for x in ids) and {x.get('SM') for x in groups}=={cfg['sample']},'header sample/RG scope')
    reference_rows=[(x.split('\t')[0],int(x.split('\t')[1])) for x in Path(cfg['inputs']['reference_index']['path']).read_text().splitlines()]
    same([(x['SN'],x['LN']) for x in header['SQ']],reference_rows,'reference header geometry')
    dictionary=[]
    for line in Path(cfg['inputs']['reference_dictionary']['path']).read_text().splitlines():
        if line.startswith('@SQ\t'):
            fields=dict(item.split(':',1) for item in line.split('\t')[1:]);dictionary.append((fields['SN'],int(fields['LN'])))
    same(dictionary,reference_rows,'GATK dictionary geometry')
    markers=marker_map(cfg['inputs']['panel']['path'],cfg['expected_panel_markers'],cfg['contigs']);rows=read_rows(run/'pileups.tsv',['contig','position','ref_count','alt_count','other_alt_count','allele_frequency'],cfg['sample']);denominators,payload,membership,per_chrom=reconcile(rows,markers)
    same(denominators,summary['denominators'],'producer count/coverage denominators differ')
    producer_membership=parse((run/'pileup-membership.json').read_text());same(producer_membership['denominators'],denominators,'count types differ');same(producer_membership['emitted_rows'],payload,'full emitted row evidence differs')
    same(producer_membership['no_emission'],[{'contig':x['contig'],'position':x['position']} for x in membership if x['category']=='no_emission'],'missing-row evidence differs')
    require(type(summary['estimator_executed']) is bool and summary['estimator_executed']==bool(denominators['coverage_retained']),'estimator eligibility gate differs')
    if not denominators['coverage_retained']:
        require(not (run/'contamination.tsv').exists() and not (run/'minor-allele-fraction-segments.tsv').exists() and summary['estimator_not_run_reason']=='zero_native_coverage_retained_sites','ineligible estimator ran')
        result=dict(classification='insufficient_information',raw_contamination=None,raw_error=None,interpreted_estimate=None,fallback_0_1=False,clinical_pass=None,zero_means_absence=False);segments=[]
    else:
        require(summary['estimator_not_run_reason'] is None,'unexpected skip reason');result=interpret(read_rows(run/'contamination.tsv',['sample','contamination','error']),cfg['sample'])
        segments=read_rows(run/'minor-allele-fraction-segments.tsv',['contig','start','end','minor_allele_fraction'],cfg['sample']);lengths={line.split('\t')[0]:int(line.split('\t')[1]) for line in Path(cfg['inputs']['reference_index']['path']).read_text().splitlines()}
        for row in segments:
            value=float(row['minor_allele_fraction']);require(row['contig'] in cfg['contigs'] and 0<int(row['start'])<=int(row['end'])<=lengths[row['contig']] and math.isfinite(value) and 0<=value<=.5,'malformed raw segmentation')
    labels=['pileup']+(['estimator'] if denominators['coverage_retained'] else []);require([x['label'] for x in summary['commands']]==labels,'native command scope differs')
    prefix=[cfg['runtime']['java']['path'],'-Xmx2g','-XX:ActiveProcessorCount=2','-jar',cfg['runtime']['gatk']['path']];inputs=cfg['inputs']
    expected_commands={'pileup':prefix+['GetPileupSummaries','-R',inputs['reference']['path'],'-I',inputs['alignment']['path'],'-V',inputs['panel']['path'],'-L',inputs['panel']['path'],'--minimum-mapping-quality','50','--max-depth-per-sample','0','-O',str(run/'pileups.tsv')],'estimator':prefix+['CalculateContamination','-I',str(run/'pileups.tsv'),'-O',str(run/'contamination.tsv'),'--tumor-segmentation',str(run/'minor-allele-fraction-segments.tsv')]}
    for command in summary['commands']:
        require(type(command['exit_code']) is int and command['exit_code']==0,'native command failed');same(command,parse((run/'commands'/(command['label']+'.json')).read_text()),'native command receipt mismatch');same(command['argv'],expected_commands[command['label']],'native command/input/filter changed')
    for observed in pins:
        now=Path(observed['path']).stat();require((now.st_size,now.st_mtime_ns)==(observed['bytes'],observed['mtime_ns']),'verified input/output changed')
    out.mkdir();document=dict(status='complete',completed_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),source_sha256=fingerprint(__file__)['sha256'],verified_files=pins,denominators=denominators,by_chromosome=per_chrom,marker_membership=membership,model_result=result,segment_rows=len(segments),AF_semantics='Both original panel and serialized native table literals retained; finite IEEE754 binary64 values must match exactly.',scope='Independent saved-output/source/denominator/model-result reconciliation, no producer imports. No independent read counter or biological mixture accuracy claim.',limitations=['Read observations have no base-quality threshold or overlap collapse.','No pure/clinical-negative/cell-fraction/mosaicism/diagnosis inference; boundary zero does not imply absence.','Native internal homozygous-site selection or strategy is not inferred from the count ledger.'])
    with (out/'verification.json').open('x') as f:f.write(json.dumps(document,indent=2,sort_keys=True,allow_nan=False)+'\n')
if __name__=='__main__':main()
