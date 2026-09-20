#!/usr/bin/env python3
"""Own one bounded fictional pilot process group, including Python and Java."""
import argparse,datetime,hashlib,json,os,platform,re,shlex,shutil,signal,subprocess,time
from pathlib import Path
ROOT=Path('/Users/jeweled/Documents/mva-hackathon')
EXPECTED_CALLER_ARGV=[str(ROOT/'tools/bio/bin/python'),str(ROOT/'repo/workflow/scripts/call_genome.py'),'--root',str(ROOT),'--reference',str(ROOT/'ref/realignment/GRCh38_nochr.fa'),'--alignment',str(ROOT/'align/gatk-compatible-cram30-20260920/WGS_EX2312012.cram'),'--sample','WGS_EX2312012','--output',str(ROOT/'calls/from-fastq-cram30-20260920'),'--workers','4']
DEFINITION_SHA='6bb6aca4c6e281ccbe345b6f7a4e7a38f340da8a9138b3cd8faec87ccc532267'
TERMINAL_STATUS_SHA='d58b058470ab989ac27257699705401a58a6a81c587b09c4f3018ee64a6788e0'
CALLER_RUN=ROOT/'calls/from-fastq-cram30-20260920'
LIMITS=dict(maximum_fixture_total_bytes=2*1024**3,minimum_free_bytes=40*1024**3,per_command_wall_seconds=300,maximum_group_RSS_bytes=4*1024**3,total_wall_seconds=1800)

def require(v,m):
    if not v:raise ValueError(m)
def stamp():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def read_json(path):
    def pairs(items):
        d={}
        for k,v in items:require(k not in d,'duplicate JSON key');d[k]=v
        return d
    def bad(v):raise ValueError('nonfinite JSON '+v)
    return json.loads(Path(path).read_text(),object_pairs_hook=pairs,parse_constant=bad)
def sha(path):
    p=Path(path);before=p.stat();h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    after=p.stat();require((before.st_size,before.st_mtime_ns)==(after.st_size,after.st_mtime_ns),'input changed while hashing');return h.hexdigest()
def save(path,value):
    tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n');tmp.replace(path)
def parse_caller(line):
    fields=line.strip().split(None,7);require(len(fields)==8,'original caller not live or malformed identity')
    result=dict(pid=int(fields[0]),pgid=int(fields[1]),started=' '.join(fields[2:7]),argv=shlex.split(fields[7]))
    require(result['pid']==79898 and result['argv']==EXPECTED_CALLER_ARGV,'original caller exact command/run differs')
    return result
def validate_terminal_metadata(state,definition,shards,methods,run):
    require(state['status']=='complete' and state['pid']==79898 and state['definition_sha256']==DEFINITION_SHA,'terminal caller run identity')
    require(type(state['expected_shards']) is int and state['expected_shards']==64 and len(definition['shards'])==64,'terminal expected shard scope')
    require(state['completed_shards']==list(range(64)) and all(type(v) is int for v in state['completed_shards']),'terminal shard membership')
    require(bool(state.get('completed_at')) and bool(state.get('started_at')) and set(state['outputs'])=={'gatk','bcftools'},'terminal completion fields')
    require(len(shards)==64 and definition['sample']=='WGS_EX2312012','terminal metadata source scope')
    totals={'gatk':0,'bcftools':0}
    for receipt in shards:
        require(bool(receipt.get('completed_at')) and re.fullmatch('[0-9a-f]{64}',receipt.get('commands_sha256','')),'shard completion identity')
        for caller in totals:
            count=receipt[caller]['counts']['records'];require(type(count) is int and count>=0,'invalid shard record count');totals[caller]+=count
    for caller in totals:
        result=state['outputs'][caller];require(result==methods[caller],'final method receipt differs from status')
        require(result['path']==str(run/('WGS_EX2312012.'+caller+'.normalized.raw.vcf.gz')),'final published path differs')
        require(bool(re.fullmatch('[0-9a-f]{64}',result['sha256'])),'invalid final output digest declaration')
        require(type(result['counts']['records']) is int and result['counts']['records']==totals[caller],'final count/shard sum differs')
    return totals

def caller_identity():
    # Read only small controller metadata. Never open the patient alignment or VCF outputs.
    pins={}
    def metadata(path):
        raw=path.read_bytes();value=json.loads(raw);pins[str(path)]=dict(path=str(path),bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest());return value
    definition=metadata(CALLER_RUN/'run-definition.json');require(pins[str(CALLER_RUN/'run-definition.json')]['sha256']==DEFINITION_SHA,'original definition changed')
    state=metadata(CALLER_RUN/'status.json')
    require(state.get('pid')==79898 and state.get('definition_sha256')==DEFINITION_SHA and type(state.get('expected_shards')) is int and state['expected_shards']==64,'original state identity differs')
    identity=dict(pid=79898,definition_sha256=DEFINITION_SHA,run=str(CALLER_RUN),started_at=state['started_at'],expected_shards=64)
    proc=subprocess.run(['ps','-p','79898','-o','pid=,pgid=,lstart=,command='],text=True,capture_output=True)
    if state['status']=='running':
        require(proc.returncode==0,'original caller absent without terminal complete status')
        return dict(mode='running',identity=identity,process=parse_caller(proc.stdout),metadata_pins=pins)
    require(state['status']=='complete','original caller failed or unknown state')
    require(pins[str(CALLER_RUN/'status.json')]['sha256']==TERMINAL_STATUS_SHA,'reviewed terminal caller status SHA differs')
    shards=[metadata(CALLER_RUN/('shard-%04d'%i)/'complete.json') for i in range(64)]
    methods={c:metadata(CALLER_RUN/(c+'-complete.json')) for c in ['gatk','bcftools']}
    totals=validate_terminal_metadata(state,definition,shards,methods,CALLER_RUN)
    return dict(mode='complete_metadata_only',identity=identity,process=parse_caller(proc.stdout) if proc.returncode==0 else None,metadata_pins=pins,record_count_sums=totals,scope='Controller completion metadata only; output VCF bytes not accessed or independently verified.')

def exact_identity(before,after):
    # Process identity may disappear only through the approved natural-completion metadata gate.
    if 'mode' not in before:
        require(before==after,'original caller PID/start/group/command identity changed');return
    require(before['identity']==after['identity'],'original caller run/definition/PID/start identity changed')
    if before['mode']=='running' and after['mode']=='running':require(before['process']==after['process'],'original live process identity changed')
    elif before['mode']=='running' and after['mode']=='complete_metadata_only':
        if after['process'] is not None:require(before['process']==after['process'],'terminal live process differs')
    elif before['mode']=='complete_metadata_only' and after['mode']=='complete_metadata_only':
        require(before['metadata_pins']==after['metadata_pins'],'completed caller metadata changed')
        if before['process'] is None:require(after['process'] is None,'completed caller PID restarted')
        elif after['process'] is not None:require(before['process']==after['process'],'completed caller process identity changed')
    else:raise ValueError('completed caller restarted or invalid transition')

def group(pgid):
    rows=[]
    for line in subprocess.check_output(['ps','-axo','pid=,ppid=,pgid=,rss='],text=True).splitlines():
        p=line.split()
        if len(p)==4 and int(p[2])==pgid:rows.append(dict(pid=int(p[0]),ppid=int(p[1]),pgid=int(p[2]),RSS_bytes=int(p[3])*1024))
    return rows
def size(root):
    if not root.exists():return 0
    total=0
    for p in root.rglob('*'):
        require(not p.is_symlink(),'unexpected output symlink')
        if p.is_file():total+=p.stat().st_size
    return total
def limit_reason(elapsed,rss,total,free,limits):
    for key,value in [('total_wall_seconds',elapsed),('maximum_group_RSS_bytes',rss),('maximum_fixture_total_bytes',total)]:
        if value>=limits[key]:return key
    return 'minimum_free_bytes' if free<limits['minimum_free_bytes'] else None
def stop_group(process):
    require(process.pid>1 and process.pid!=os.getpgrp(),'unsafe owned process group')
    try:os.killpg(process.pid,signal.SIGTERM)
    except ProcessLookupError:pass
    deadline=time.monotonic()+5
    while time.monotonic()<deadline:
        process.poll()
        if not group(process.pid):break
        time.sleep(.05)
    process.poll()
    if group(process.pid):
        try:os.killpg(process.pid,signal.SIGKILL)
        except ProcessLookupError:pass
    process.wait(timeout=5)
def supervise(command,out,launch,limits,disk_root):
    require(not out.exists() and not launch.exists(),'no retry or overwrite')
    require(shutil.disk_usage(disk_root).free>=limits['minimum_free_bytes'],'minimum free before launch')
    launch.mkdir(parents=True);state=dict(status='launching',started_at=stamp(),command=command,limits=limits,sampling_seconds=.25);save(launch/'status.json',state);process=None;start=time.monotonic()
    try:
        with (launch/'stdout.log').open('xb') as stdout,(launch/'stderr.log').open('xb') as stderr:
            process=subprocess.Popen(command,stdout=stdout,stderr=stderr,start_new_session=True)
            state.update(status='running',pid=process.pid,pgid=process.pid);save(launch/'status.json',state)
            while True:
                rows=group(process.pid);rss=sum(r['RSS_bytes'] for r in rows);total=size(out)+size(launch);free=shutil.disk_usage(disk_root).free;elapsed=time.monotonic()-start
                state.update(checked_at=stamp(),elapsed_seconds=elapsed,RSS_bytes=rss,output_bytes=total,free_bytes=free,peak_sampled_RSS_bytes=max(state.get('peak_sampled_RSS_bytes',0),rss),peak_output_bytes=max(state.get('peak_output_bytes',0),total))
                reason=limit_reason(elapsed,rss,total,free,limits)
                if reason:state['limit_triggered']=reason;raise ValueError('whole-pilot guard '+reason)
                code=process.poll()
                if code is not None:
                    require(not group(process.pid),'owned descendants remain after runner');require(code==0,'pilot runner exit '+str(code));break
                save(launch/'status.json',state);time.sleep(.25)
        state.update(status='complete',completed_at=stamp(),exit_code=0)
        return state
    except BaseException as e:
        if process is not None:stop_group(process)
        state.update(status='failed',failed_at=stamp(),error_type=type(e).__name__,error=str(e),exit_code=process.returncode if process else None);raise
    finally:save(launch/'status.json',state)
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);p.add_argument('--expected-plan-sha256',required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--launch',type=Path,required=True);args=p.parse_args();os.nice(10)
    require(platform.node().split('.')[0]=='512S2','Mac Studio512S2 required');require(sha(args.plan)==args.expected_plan_sha256,'plan SHA mismatch');plan=read_json(args.plan)
    require(plan['ready_for_native_pilot'] is True and plan['case_execution_authorized'] is False and plan['scope']=='public_marker_derivation_and_fictional_controls_only','reviewed fictional-only ready configuration required')
    require(plan['resource_limits']==LIMITS and plan['root']==str(ROOT),'resource/root contract differs')
    require(plan['expected_original_caller_argv']==EXPECTED_CALLER_ARGV,'expected caller command differs')
    require(plan['original_caller_definition_sha256']==DEFINITION_SHA and plan['original_caller_terminal_sha256']==TERMINAL_STATUS_SHA,'reviewed caller metadata pins differ')
    source=Path(__file__).resolve().parent
    require(set(plan['source_sha256'])=={'guard-fictional-pilot.py','run-fictional-pilot.py','fictional_controls.py','test_pilot_controls.py','test_pilot_guard.py'},'exact source inventory')
    for name,digest in plan['source_sha256'].items():require(sha(source/name)==digest,'source hash differs '+name)
    before=caller_identity();out=args.output.resolve();launch=args.launch.resolve()
    require(out.is_relative_to(ROOT/'results/nuclear-mixture-pilot-20260920') and launch.is_relative_to(ROOT/'methods/nuclear-mixture-pilot-20260920'),'isolated output/launch scope')
    require(not out.is_relative_to(launch) and not launch.is_relative_to(out),'separate output/guard scope')
    command=[str(ROOT/'tools/bio/bin/python'),'-B',str(source/'run-fictional-pilot.py'),'--plan',str(args.plan.resolve()),'--expected-plan-sha256',args.expected_plan_sha256,'--output',str(out),'--outer-guard-pid',str(os.getpid())]
    state=None
    try:
        state=supervise(command,out,launch,LIMITS,ROOT)
        terminal=read_json(out/'status.json');require(terminal['status']=='complete','pilot terminal incomplete');require(sha(out/'summary.json')==terminal['summary_sha256'],'pilot summary link')
        after=caller_identity();exact_identity(before,after)
        state.update(original_caller_before=before,original_caller_after=after,source_sha256=sha(__file__),plan_sha256=args.expected_plan_sha256,pilot_status_sha256=sha(out/'status.json'),pilot_summary_sha256=sha(out/'summary.json'))
    except BaseException as e:
        if launch.exists():
            state=read_json(launch/'status.json');state.update(status='failed',failed_at=stamp(),error_type=type(e).__name__,error=str(e),original_caller_before=before)
            try:state['original_caller_after']=caller_identity()
            except ValueError as identity_error:state['original_caller_after_error']=str(identity_error)
        raise
    finally:
        if state is not None:save(launch/'status.json',state)
if __name__=='__main__':main()
