"""Strict bounded patient-assay contract; no execution at import."""
import csv,gzip,hashlib,io,json,math,statistics
from pathlib import Path
ROOT=Path('/Users/jeweled/Documents/mva-hackathon')
REAL_INPUTS={'alignment':'c3ef27cc9bde0010da9d6e872f88d52b542cbb435ee1d6f85f7f88509ce6f26b','alignment_index':'4c04806b8ffe3a5dcf48599385b221e4688430cb6843d09cd357cddf5fd6a692','reference':'187917c00fa7373a054f9768009d19e4f6e5fc87d3a6c125ce131164ca69fcd6','reference_index':'d81f672f40cfabf2510f37e3f790f2c79665023f4ac950bb9cc6496a40cebf32','reference_dictionary':'c7ca757292292168446b572776853e809c0547bd00bee4db048d735c455d7e40','panel':'4837d35d3b4c2d793aec31070755b90ad203b9a066add16d3c40ade7c7c01d40','panel_index':'4a425b9fe8926f802a5efa114f76009ca4982e6077d35c1565ce015cf054387a'}
RUNTIME={'java':'ec64ac719290ad080c7f8a4f5d83402e4152f98f486fe14150cf734a46ebcbb1','gatk':'b2824999fafe4dc83036e334b9765d4a7de4f4ff4bec1e929d77bea04a1560c7'}
LIMITS=dict(total_wall_seconds=5400,maximum_group_RSS_bytes=4*1024**3,maximum_fixture_total_bytes=2*1024**3,minimum_free_bytes=40*1024**3,per_command_wall_seconds=3600)
FIELDS=['contig','position','ref_count','alt_count','other_alt_count','allele_frequency']
SOURCE_NAMES={'patient_common.py','run-patient-assay.py','audit-patient-output.py','guard-patient-assay.py','guard-fictional-pilot.py'}
def require(v,m):
    if not v:raise ValueError(m)
def parse(text):
    def pairs(items):
        d={}
        for k,v in items:require(k not in d,'duplicate JSON key');d[k]=v
        return d
    def bad(v):raise ValueError('nonfinite JSON '+v)
    return json.loads(text,object_pairs_hook=pairs,parse_constant=bad)
def fingerprint(path,expected=None):
    p=Path(path).resolve();before=p.stat();h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    after=p.stat();require((before.st_size,before.st_mtime_ns)==(after.st_size,after.st_mtime_ns),'unstable input')
    if expected:require(h.hexdigest()==expected,'SHA differs '+str(p))
    return dict(path=str(p),bytes=after.st_size,mtime_ns=after.st_mtime_ns,sha256=h.hexdigest())
def save(path,data):
    with Path(path).open('x') as f:f.write(json.dumps(data,indent=2,sort_keys=True,allow_nan=False)+'\n')
def config(path,expected):
    pin=fingerprint(path,expected);c=parse(Path(path).read_text());require(type(c['fixture_only']) is bool and c['ready'] is True,'ready typed configuration required')
    require(json.dumps(c['limits'],sort_keys=True)==json.dumps(LIMITS,sort_keys=True) and json.dumps(c['timeouts'],sort_keys=True)==json.dumps({'pileup':3600,'estimator':600},sort_keys=True),'resource contract changed')
    require(set(c['inputs'])==set(REAL_INPUTS) and set(c['runtime'])==set(RUNTIME),'input inventory changed')
    require(type(c['expected_read_groups']) is int and c['expected_read_groups']>0 and type(c['expected_panel_markers']) is int and c['expected_panel_markers']>0,'invalid scope counts')
    if c['fixture_only']:
        require(c['patient_assay_authorized'] is False and c['sample'].startswith('FICTIONAL_'),'fictional scope')
        base=Path(c['fixture_root']).resolve();require(base!=ROOT and not any(base.is_relative_to(ROOT/x) for x in ['align','calls','ref','results']),'invalid fictional root')
        for v in c['inputs'].values():require(Path(v['path']).resolve().is_relative_to(base),'fixture input escaped')
        require(Path(c['output']).resolve().is_relative_to(base),'fixture output escaped')
    else:
        require(c['patient_assay_authorized'] is True and c['sample']=='WGS_EX2312012' and c['expected_read_groups']==4 and c['expected_panel_markers']==37267 and c['contigs']==list(map(str,range(1,23))),'real scope changed')
        require({k:v['sha256'] for k,v in c['inputs'].items()}==REAL_INPUTS,'real input pins changed')
        require(Path(c['output'])==ROOT/'results/nuclear-mixture-assay-20260920/assay-v1' and Path(c['launch'])==ROOT/'methods/nuclear-mixture-assay-20260920/launch-v1','real output changed')
    require({k:v['sha256'] for k,v in c['runtime'].items()}==RUNTIME,'runtime pins changed')
    require(c['inputs']['alignment_index']['path']==c['inputs']['alignment']['path']+'.crai' and c['inputs']['reference_index']['path']==c['inputs']['reference']['path']+'.fai' and c['inputs']['panel_index']['path']==c['inputs']['panel']['path']+'.tbi','index path binding')
    require(c['inputs']['reference_dictionary']['path']==str(Path(c['inputs']['reference']['path']).with_suffix('.dict')),'reference dictionary binding')
    require(set(c['source_sha256'])==SOURCE_NAMES,'executed source inventory differs');sources={}
    for name,digest in c['source_sha256'].items():sources[name]=fingerprint(Path(__file__).parent/name,digest)
    if not c['fixture_only']:
        q=c['qualification'];fingerprint(q['path'],q['sha256']);controls=parse(Path(q['path']).read_text())
        require(controls['status']=='complete' and type(controls['tests_run']) is int and controls['tests_run']>=10 and controls['skipped']==[] and all(type(controls[k]) is int and controls[k]==0 for k in ['errors','failures','exit_code']),'patient preparation not qualified')
        require(set(controls['executed_source_sha256'])==SOURCE_NAMES and controls['executed_source_sha256']==c['source_sha256'],'qualified source differs')
        require(controls['runtime_sha256']==RUNTIME,'qualified native runtime differs')
        import sys
        require(fingerprint(sys.executable)['sha256']==controls['python']['sha256'],'qualified Python differs')
    return c,pin,sources
def panel(path,expected,contigs):
    result={}
    with gzip.open(path,'rt') as f:
        for line in f:
            if line.startswith('#'):continue
            cols=line.rstrip('\n').split('\t');require(len(cols)==8,'panel must be sites only');chrom,pos,_,ref,alt,_,_,info=cols;pos=int(pos);key=(chrom,pos)
            attrs=dict(x.split('=',1) for x in info.split(';') if '=' in x);literal=attrs['AF'];af=float(literal)
            require(key not in result and chrom in contigs and pos>0 and len(ref)==len(alt)==1 and ref in 'ACGT' and alt in 'ACGT' and ref!=alt and math.isfinite(af) and .01<af<.2,'invalid panel marker')
            result[key]=dict(ref=ref,alt=alt,AF_literal=literal,AF_binary64=af)
    require(len(result)==expected and {k[0] for k in result}==set(contigs),'panel scope differs');return result
def table(path,sample):
    lines=Path(path).read_text().splitlines();metadata=[x for x in lines if x.startswith('#<METADATA>SAMPLE=')];require(metadata==['#<METADATA>SAMPLE='+sample],'table sample metadata')
    reader=csv.DictReader(io.StringIO('\n'.join(x for x in lines if not x.startswith('#'))),delimiter='\t');require(reader.fieldnames==FIELDS,'table columns');rows=list(reader)
    require(all(set(row)==set(FIELDS) and all(type(v) is str and v!='' for v in row.values()) for row in rows),'malformed table row width/value');return rows
def inspect_table(path,markers,sample):
    rows=table(path,sample);depths=[];seen=set();payload=[]
    for row in rows:
        key=(row['contig'],int(row['position']));require(key in markers and key not in seen,'unknown/duplicate table marker');seen.add(key)
        counts=[int(row[k]) for k in FIELDS[2:5]];require(all(x>=0 for x in counts),'negative counts')
        af=float(row['allele_frequency']);require(math.isfinite(af) and af==markers[key]['AF_binary64'],'native AF numeric identity changed')
        total=sum(counts);payload.append(dict(contig=key[0],position=key[1],counts=counts,total=total,panel_AF_literal=markers[key]['AF_literal'],table_AF_literal=row['allele_frequency']))
        if total>10:depths.append(total)
    # Apache Commons Math Mean uses initial sum/n followed by a correction pass.
    mean=None;median=None;retained=[]
    if depths:
        mean=sum(depths)/len(depths);correction=0.0
        for depth in depths:correction+=depth-mean
        mean+=correction/len(depths);median=statistics.median(depths)
        retained=[r for r in payload if r['total']>10 and r['total']>.5*median and r['total']<3*mean]
    return dict(panel_sites=len(markers),emitted=len(rows),no_emission=len(markers)-len(rows),emitted_zero=sum(r['total']==0 for r in payload),emitted_1_to_10=sum(1<=r['total']<=10 for r in payload),preliminary_gt10=len(depths),coverage_retained=len(retained),coverage_excluded=len(depths)-len(retained),covered_mean=mean,covered_median=median),payload
