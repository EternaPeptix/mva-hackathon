#!/usr/bin/env python3
"""Independent finite saved-depth audit; no producer imports or alignment access."""
from collections import Counter, defaultdict
from pathlib import Path
import bisect
import datetime
import gzip
import hashlib
import json
import platform
import re
import sys
from urllib.parse import unquote

ROOT=Path('/Users/jeweled/Documents/mva-hackathon')
METHOD=ROOT/'methods/secondary-depth-gap-map-20260920'
RUN=ROOT/'results/secondary-depth-gap-map-20260920/map-v1'
MODES=('MQ20_BQ20','MQ60_BQ20')
CONFIG_SHA='471b97cfc09c6d3cecc9d5c4550d6942ae67a1b5e5cf7034a16342e81d2eaafc'
SUMMARY_SHA='db6d50cc8e115cbdca697f581d7ce5a754a419b50ff282276b8f9a52ba828eb4'
STATUS_SHA='c4560b81a9597715aad4bf532e0e6fcd4d588e994249ecc8d7882783adc81e19'
GUARD_SHA='fd0fc7eaf1c07e6632e267e421c9135ef17ed437a1520d8a77bc848eebb74760'


def need(value,message):
    if not value:raise ValueError(message)


def canon(value):return json.dumps(value,sort_keys=True,allow_nan=False,separators=(',',':'))
def equal(a,b,message):need(canon(a)==canon(b),message)
def read(path):return json.loads(Path(path).read_text())


def fingerprint(path):
    path=Path(path);before=path.stat();h=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    after=path.stat();need((before.st_size,before.st_mtime_ns)==(after.st_size,after.st_mtime_ns),'file changed while hashing')
    return dict(path=str(path),bytes=after.st_size,mtime_ns=after.st_mtime_ns,sha256=h.hexdigest())


def verify_pin(pin,expected,label):
    for key in ('path','bytes','sha256'):equal(pin[key],expected[key],label+' '+key)


def merged(intervals):
    result=[]
    for start,end in sorted({tuple(row) for row in intervals}):
        need(type(start) is int and type(end) is int and 0<=start<end,'bad interval')
        if result and start<=result[-1][1]:result[-1][1]=max(result[-1][1],end)
        else:result.append([start,end])
    return result


def category(d):return 'zero' if d==0 else '1-19' if d<20 else 'at-least-20'


def metrics(depths):
    h=Counter(depths)
    return dict(CDS_union_bases=sum(h.values()),zero_depth_bases=h[0],below20_bases=sum(n for d,n in h.items() if d<20),
                sum_depth=sum(d*n for d,n in h.items()),histogram={str(d):n for d,n in sorted(h.items()) if n})


def run(source_sha):
    need(platform.system()=='Darwin' and platform.machine()=='arm64','Mac native host required')
    import subprocess
    need(subprocess.check_output(['sysctl','-n','hw.model'],text=True).strip()=='Mac15,14','512S2 identity')
    pins={}
    def add(path,expected=None):
        pin=fingerprint(path)
        if expected:
            if isinstance(expected,str):equal(pin['sha256'],expected,'expected SHA')
            else:verify_pin(pin,expected,'expected pin')
        if str(path) in pins:equal(pin,pins[str(path)],'duplicate input changed')
        pins[str(path)]=pin;return pin
    config_path=METHOD/'gap-map-approved-v1.json';add(config_path,CONFIG_SHA);config=read(config_path)
    need(config['ready'] is True and len(config['inputs'])==13,'exact approved inputs')
    expected_roles={'gff','fai','public_design','public_status','public_summary','selection','summary','status','verification'}|{m+'_'+k for m in MODES for k in ('depth','summary')}
    need(set(config['inputs'])==expected_roles,'input roles')
    data={}
    for role,info in config['inputs'].items():
        add(info['path'],info)
        if role not in ('gff','fai') and not role.endswith('_depth'):data[role]=read(info['path'])
    add(RUN/'summary.json',SUMMARY_SHA);add(RUN/'status.json',STATUS_SHA)
    summary=read(RUN/'summary.json');status=read(RUN/'status.json')
    launch=RUN.with_name('map-v1-launch');add(launch/'status.json',GUARD_SHA);guard=read(launch/'status.json')
    for name,info in guard['outputs'].items():add(launch/name,info)
    binding=read(launch/'worker-binding.json');result=read(launch/'worker-result.json')
    need(summary['status']==status['status']==guard['status']==result['status']=='complete','terminal completion')
    need(status['reads_opened'] is False and guard['exit_code']==0 and guard['limit_reason'] is None,'terminal guard')
    verify_pin(pins[str(config_path)],binding['config'],'executed config')
    verify_pin(pins[str(RUN/'status.json')],result['mapper_status'],'worker mapper status')
    for info in binding['sources'].values():add(info['path'],info)
    for rel,digest in config['source_sha256'].items():equal(binding['sources']['scripts/'+rel]['sha256'],digest,'qualified producer source')
    add(binding['qualification']['path'],binding['qualification']);q=read(binding['qualification']['path'])
    need(q['status']=='complete' and q['tests_run']==39 and q['errors']==q['failures']==q['exit_code']==0 and q['skipped']==[],'qualification')
    need(set(q['source_files'])==set(binding['sources']),'qualification source names')
    for label,info in binding['sources'].items():verify_pin(pins[info['path']],q['source_files'][label],'qualified source')
    add(binding['runtime']['path'],binding['runtime']);verify_pin(pins[binding['runtime']['path']],config['python'],'runtime')
    for label,info in summary['inputs'].items():
        add(info['path'],info)
        verify_pin(pins[info['path']],status['inputs'][label],'summary/status input')
    for name,info in status['outputs'].items():add(RUN/name,info)
    original=data['verification']
    need(original['status']=='complete' and original['all_verified_inputs_stable'] and original['public_target_reconciled'],'original independent depth verification')
    for role in ('status','summary','selection','public_design','public_status','gff','fai')+tuple(m+'_'+k for m in MODES for k in ('depth','summary')):
        info=config['inputs'][role];verify_pin(pins[info['path']],original['verified_files'][info['path']],'original audit '+role)
    for role in ('selection','summary')+tuple(m+'_'+k for m in MODES for k in ('depth','summary')):
        info=config['inputs'][role];verify_pin(pins[info['path']],data['status']['outputs'][Path(info['path']).name],'original producer output')
    equal(data['status']['summary_sha256'],config['inputs']['summary']['sha256'],'original terminal summary')
    for role in ('public_design','public_summary'):
        info=config['inputs'][role];verify_pin(pins[info['path']],data['public_status']['outputs'][Path(info['path']).name],'public output')

    # Rebuild membership directly from original annotation lines, without target helper.
    selection=data['selection'];genes=selection['genes'];public=data['public_design']
    need(len(genes)==84 and len({g['id'] for g in genes})==84,'gene identity count')
    equal(genes,[{k:g[k] for k in genes[0]} for g in public['genes']],'public full gene context')
    lengths={}
    for line in Path(config['inputs']['fai']['path']).read_text().splitlines():
        f=line.split('\t');need(f[0] not in lengths,'duplicate contig');lengths[f[0]]=int(f[1])
    wanted_gene={g['ordinal']:g for g in genes};wanted_cds=defaultdict(list);tx_expected={}
    for g in genes:
        for t in g['transcript_ids']:
            need(t not in tx_expected,'transcript in multiple genes');tx_expected[t]=g['id']
        for row in g['cds_rows']:wanted_cds[row['ordinal']].append(row)
    found_gene=set();found_cds=set();found_tx=set();all_selected_cds=[]
    primary={str(n) for n in range(1,23)}|{'X','Y'}
    with gzip.open(config['inputs']['gff']['path'],'rt') as stream:
        for ordinal,line in enumerate(stream,1):
            if line.startswith('#') or not line.strip():continue
            f=line.rstrip('\n').split('\t')
            need(len(f)==9,'bad GFF columns')
            if f[0] not in primary:continue
            attrs={k:unquote(v) for x in f[8].split(';') if '=' in x for k,v in [x.split('=',1)]}
            ident=attrs.get('ID','');parents=attrs.get('Parent','').split(',')
            if ordinal in wanted_gene:
                g=wanted_gene[ordinal];equal((ident,attrs.get('Name'),attrs.get('biotype'),f[0],int(f[3])-1,int(f[4]),f[6]),(g['id'],g['name'],g['biotype'],g['chrom'],g['start'],g['end'],g['strand']),'original gene line');found_gene.add(ordinal)
            if ident in tx_expected:
                equal(parents,[tx_expected[ident]],'original transcript parent');found_tx.add(ident)
            if f[2]=='CDS':
                for t in parents:
                    if t in tx_expected:
                        row=dict(chrom=f[0],start=int(f[3])-1,end=int(f[4]),strand=f[6],ordinal=ordinal,parents=parents,phase=int(f[7]),transcript_id=t)
                        need(0<=row['start']<row['end']<=lengths[row['chrom']],'CDS reference span')
                        need(row in wanted_cds.get(ordinal,[]),'unpreserved selected CDS row');all_selected_cds.append(row)
                if ordinal in wanted_cds:
                    for row in wanted_cds[ordinal]:need(row['transcript_id'] in parents,'CDS parent missing')
                    found_cds.add(ordinal)
    need(found_gene==set(wanted_gene) and found_cds==set(wanted_cds) and found_tx==set(tx_expected),'annotation completeness')
    equal(sorted(all_selected_cds,key=lambda x:(x['ordinal'],x['transcript_id'])),sorted([r for g in genes for r in g['cds_rows']],key=lambda x:(x['ordinal'],x['transcript_id'])),'all CDS payloads')

    # Independent set/merged-union accounting, rather than producer event sweep.
    groups={};metadata={};panel_by_chrom=defaultdict(list)
    for g in genes:
        intervals=merged((r['start'],r['end']) for r in g['cds_rows']);groups[g['id']]=(g['chrom'],intervals)
        equal(sum(e-s for s,e in intervals),g['CDS_union_bases'],'gene union denominator')
        panel_by_chrom[g['chrom']].extend(intervals)
        for tx in g['transcript_ids']:
            rows=[r for r in g['cds_rows'] if r['transcript_id']==tx]
            groups[tx]=(g['chrom'],merged((r['start'],r['end']) for r in rows))
            metadata[tx]=dict(gene_id=g['id'],gene_name=g['name'],chrom=g['chrom'],strand=g['strand'],cds_rows=rows)
    panel={chrom:merged(v) for chrom,v in panel_by_chrom.items()}
    panel_positions=[(chrom,pos) for chrom in lengths if chrom in panel for s,e in panel[chrom] for pos in range(s,e)]
    need(len(panel_positions)==393880 and len(tx_expected)==2578,'scope counts')
    saved_positions=[(w['chrom'],pos) for w in selection['windows'] for pos in range(w['start'],w['end'])]
    equal(saved_positions,panel_positions,'saved windows versus independent panel union')
    need(len(selection['windows'])==2112,'window count')
    depths=[]
    for m in MODES:
        d={}
        with gzip.open(config['inputs'][m+'_depth']['path'],'rt') as stream:
            for expected in panel_positions:
                f=next(stream,'').rstrip('\n').split('\t');need(len(f)==3 and f[0]==expected[0] and f[1]==str(expected[1]+1) and re.fullmatch(r'0|[1-9][0-9]*',f[2]),'depth position/order/type');d[expected]=int(f[2])
            need(next(stream,None) is None,'extra saved depth')
        depths.append(d)
    expected_genes=[];expected_tx=[]
    for g in genes:
        chrom,intervals=groups[g['id']]
        modes={m:metrics(d[(chrom,pos)] for s,e in intervals for pos in range(s,e)) for m,d in zip(MODES,depths)}
        expected_genes.append(dict(gene_id=g['id'],name=g['name'],chrom=chrom,strand=g['strand'],modes=modes))
    for t in sorted(metadata):
        chrom,intervals=groups[t]
        modes={m:metrics(d[(chrom,pos)] for s,e in intervals for pos in range(s,e)) for m,d in zip(MODES,depths)}
        expected_tx.append(dict(transcript_id=t,**metadata[t],modes=modes))
    equal(expected_genes,summary['genes'],'all gene summary payloads')
    equal(expected_tx,summary['transcripts'],'all transcript summary payloads')
    expected_panel={m:metrics(d.values()) for m,d in zip(MODES,depths)}
    equal(expected_panel,summary['panel_union'],'complete panel summary')
    pairs=Counter((category(depths[0][p]),category(depths[1][p])) for p in panel_positions)
    cells=[dict(MQ20_BQ20=a,MQ60_BQ20=b,bases=pairs[a,b]) for a in ('zero','1-19','at-least-20') for b in ('zero','1-19','at-least-20')]
    equal(cells,summary['paired_categories'],'paired categories')
    increases=[p for p in panel_positions if depths[1][p]>depths[0][p]]
    equal(len(increases),summary['MQ60_depth_greater_than_MQ20_bases'],'numeric increases')
    # Resolve only gap-position CDS membership using interval bisection; no sweep helper.
    gap_positions=[p for p in panel_positions if depths[0][p]<20 or depths[1][p]<20]
    by_chrom=defaultdict(list)
    for chrom,pos in gap_positions:by_chrom[chrom].append(pos)
    membership=defaultdict(set)
    for g in genes:
        positions=by_chrom[g['chrom']]
        for row in g['cds_rows']:
            lo=bisect.bisect_left(positions,row['start']);hi=bisect.bisect_left(positions,row['end'])
            for pos in positions[lo:hi]:membership[g['chrom'],pos].add((g['id'],row['transcript_id'],row['ordinal']))
    intervals=[];last=None;last_signature=None
    for key in gap_positions:
        chrom,pos=key;members=sorted(membership[key]);need(members,'unassigned gap')
        values=[d[key] for d in depths];cats=[category(v) for v in values];signature=(chrom,tuple(members),tuple(cats))
        if last is None or last['end']!=pos or signature!=last_signature:
            if last is not None:intervals.append(last)
            last=dict(chrom=chrom,start=pos,end=pos,coordinates='0-based half-open',gene_ids=sorted({m[0] for m in members}),transcript_ids=sorted({m[1] for m in members}),source_CDS_members=[dict(gene_id=g,transcript_id=t,gff_line_ordinal=o) for g,t,o in members],categories=cats,depth_modes=list(MODES),depth_statistics=[dict(minimum=v,maximum=v,sum_depth=0) for v in values])
            last_signature=signature
        last['end']=pos+1
        for st,v in zip(last['depth_statistics'],values):st['minimum']=min(st['minimum'],v);st['maximum']=max(st['maximum'],v);st['sum_depth']+=v
    if last is not None:intervals.append(last)
    with gzip.open(RUN/'gap-intervals.jsonl.gz','rt') as stream:observed=[json.loads(line) for line in stream]
    equal(intervals,observed,'all gap interval fields/order/boundaries')
    equal(len(intervals),summary['gap_intervals'],'interval count')
    need(sum(r['end']-r['start'] for r in intervals)==len(gap_positions),'gap union')
    # Bind original full summaries and independently recompute every original window histogram.
    for m,d in zip(MODES,depths):
        old=data[m+'_summary'];need(old['position_count']==393880,'original position count')
        equal(old['mode'],next(x for x in config['depth_modes'] if x['name']==m),'original mode')
        equal(data['summary']['depth_modes'][m],{k:old[k] for k in ('mode','position_count','groups')},'saved aggregate summary')
        for w,ow in zip(selection['windows'],old['windows']):
            equal({k:ow[k] for k in w},w,'original window identity')
            got=metrics(d[w['chrom'],pos] for pos in range(w['start'],w['end']))
            for key in ('span_depth','ACGT_depth'):
                equal({k:v for k,v in ow[key]['histogram'].items() if v},got['histogram'],'original window histogram')
        for label,oldgroup in old['groups'].items():
            got=expected_panel[m] if label=='panel_union' else next(g for g in expected_genes if g['gene_id']==label)['modes'][m]
            for key in ('span','ACGT'):
                equal({k:v for k,v in oldgroup[key]['histogram'].items() if v},got['histogram'],'original group histogram')
    for path,oldpin in list(pins.items()):equal(fingerprint(path),oldpin,'final stability')
    return dict(status='complete',verified_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),verifier_source_sha256=source_sha,
        host=platform.node(),source_imports='stdlib only; no producer/helper import',verified_files=list(pins.values()),
        all_verified_files_stable=True,input_data_roles=sorted(expected_roles),genes=84,transcripts=2578,windows=2112,panel_union_bases=393880,
        gap_intervals=len(intervals),gap_union_bases=len(gap_positions),paired_categories=cells,numerically_increased_MQ60_bases=len(increases),
        upward_category_crossings=sum(n for (a,b),n in pairs.items() if ('zero','1-19','at-least-20').index(b)>('zero','1-19','at-least-20').index(a)),
        panel_metrics=expected_panel,affected_genes={m:sum(g['modes'][m]['below20_bases']>0 for g in expected_genes) for m in MODES},
        affected_transcripts={m:sum(t['modes'][m]['below20_bases']>0 for t in expected_tx) for m in MODES},
        strict_zero_genes=[dict(name=g['name'],bases=g['modes'][MODES[1]]['zero_depth_bases']) for g in expected_genes if g['modes'][MODES[1]]['zero_depth_bases']],
        independently_checked=['Fresh original input/source/terminal hashes and qualification chains','Original GFF gene/CDS/transcript membership','Merged-set gene/transcript/panel denominators and all histograms','Every saved depth coordinate/count, explicit zero versus missing','All paired cells and numeric increases',f'All {len(intervals)} gap row fields, source memberships, exact order and category boundaries','Original saved depth summaries and output status bindings'],
        not_repeated=['No alignment/FASTA/caller/depth execution','Original aligned-read extraction and same-engine depth recount remain supported by prior independent audit only','No causal mechanism assigned to numeric increases','No clinical callability or locus-origin conclusion'])


if __name__=='__main__':
    print(json.dumps(run(hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2,allow_nan=False))
