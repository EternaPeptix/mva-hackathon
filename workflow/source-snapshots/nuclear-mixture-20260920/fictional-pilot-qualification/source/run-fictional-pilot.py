#!/usr/bin/env python3
"""Finite public-method fictional qualification; no project sample input access."""
import argparse,csv,datetime,importlib.util,hashlib,json,math,os,platform,shlex,shutil,signal,subprocess,time
from pathlib import Path
from fictional_controls import count_fixture,write_tables,save,require
ROOT=Path('/Users/jeweled/Documents/mva-hackathon')

def stamp():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def fingerprint(p,expected=None):
    p=Path(p).resolve();before=p.stat();h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    after=p.stat();require((before.st_size,before.st_mtime_ns)==(after.st_size,after.st_mtime_ns),'input changed')
    v=dict(path=str(p),bytes=after.st_size,mtime_ns=after.st_mtime_ns,sha256=h.hexdigest())
    if expected:require(v['sha256']==expected,'input hash mismatch '+p.name)
    return v

def read_plan(path,digest):
    pin=fingerprint(path,digest);plan=json.loads(Path(path).read_text())
    require(plan['scope']=='public_marker_derivation_and_fictional_controls_only' and plan['case_execution_authorized'] is False,'nonfictional scope forbidden')
    require(plan['ready_for_native_pilot'] is True,'preparation only; reviewed ready configuration required')
    require(plan['root']==str(ROOT) and plan['host']=='512S2' and platform.node().split('.')[0]=='512S2','wrong execution host')
    require(plan['native_count_policy']==dict(maximum_depth_per_sample=0,minimum_mapping_quality=50,overlap_collapse=False,per_base_quality_threshold=None,unit='read observation'),'count policy drift')
    sources={n:fingerprint(Path(__file__).resolve().parent/n,h) for n,h in plan['source_sha256'].items()}
    require(set(sources)=={'guard-fictional-pilot.py','run-fictional-pilot.py','fictional_controls.py','test_pilot_controls.py','test_pilot_guard.py'},'exact executed qualification source scope')
    runtime=plan['runtime'];require(runtime['expected_version']=='4.6.1.0' and runtime['maximum_heap_gb']==2 and runtime['active_processors']==2,'runtime policy drift')
    design=plan['table_design'];require(design['markers']==10000 and design['population_ALT_frequency']=='0.1' and design['seed']==20260920 and design['genotypes_per_100_individuals']=={'0/0':81,'0/1':18,'1/1':1} and design['fallback_contamination']==0 and design['fallback_error']==1,'table design drift')
    require(plan['table_design']['diploid_clean_acceptance_absolute_error']==.01 and plan['table_design']['diploid_mixture_acceptance_max_absolute_or_relative']==[.01,.2],'fixed tolerance drift')
    limits=plan['resource_limits'];require(limits==dict(maximum_fixture_total_bytes=2*1024**3,minimum_free_bytes=40*1024**3,per_command_wall_seconds=300,maximum_group_RSS_bytes=4*1024**3,total_wall_seconds=1800),'resource limit drift')
    tools={k:fingerprint(ROOT/runtime[k+'_path'],runtime[k+'_sha256']) for k in ['java','gatk']}
    return plan,pin,sources,tools

def group(pgid):
    return [dict(pid=int(parts[0]),rss_bytes=int(parts[2])*1024) for line in subprocess.check_output(['ps','-axo','pid=,pgid=,rss='],text=True).splitlines() if len(parts:=line.split())==3 and int(parts[1])==pgid]
def rss(pgid):return sum(r['rss_bytes'] for r in group(pgid))
def size(path):return sum(p.stat().st_size for p in path.rglob('*') if p.is_file())
def stop(p):
    # The outer guard owns the whole inherited process group. Stop only this Java child here.
    if p.poll() is None:p.terminate()
    try:p.wait(timeout=5)
    except subprocess.TimeoutExpired:p.kill();p.wait(timeout=5)
def caller_module():
    spec=importlib.util.spec_from_file_location('qualified_pilot_guard',Path(__file__).with_name('guard-fictional-pilot.py'))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module

def command_run(command,label,out,limits,deadline,expected_failure=False):
    logs=out/'commands';logs.mkdir(exist_ok=True);started=time.monotonic();record=dict(label=label,command=command,started_at=stamp(),expected_failure=expected_failure);process=None
    try:
        with (logs/(label+'.stdout.log')).open('xb') as stdout,(logs/(label+'.stderr.log')).open('xb') as stderr:
            process=subprocess.Popen(command,stdout=stdout,stderr=stderr,start_new_session=False);record['pid']=process.pid;peak=0
            while True:
                current=time.monotonic();mem=rss(os.getpgrp());peak=max(peak,mem)
                require(current-started<limits['per_command_wall_seconds'],'per-command timeout');require(current<deadline,'pilot timeout')
                require(mem<limits['maximum_group_RSS_bytes'],'own group RSS limit');require(size(out)<limits['maximum_fixture_total_bytes'],'output size limit');require(shutil.disk_usage(ROOT).free>=limits['minimum_free_bytes'],'minimum free disk')
                code=process.poll()
                if code is not None:
                    break
                time.sleep(.25)
        record.update(exit_code=code,elapsed_seconds=time.monotonic()-started,peak_sampled_group_RSS_bytes=peak,completed_at=stamp())
        require((code!=0) if expected_failure else(code==0),'unexpected command result '+label)
        record['status']='complete';return record
    except BaseException as e:
        if process is not None:stop(process)
        record.update(status='failed',error_type=type(e).__name__,error=str(e),elapsed_seconds=time.monotonic()-started);raise
    finally:save(logs/(label+'.json'),record)

def read_table(path):
    with Path(path).open() as f:return list(csv.DictReader((line for line in f if not line.startswith('#')),delimiter='\t'))
def compare_counts(path,expected):
    rows=read_table(path);got=[]
    for r in rows:
        require(set(r)=={'contig','position','ref_count','alt_count','other_alt_count','allele_frequency'},'count columns')
        got.append(dict(contig=r['contig'],position=int(r['position']),ref_count=int(r['ref_count']),alt_count=int(r['alt_count']),other_alt_count=int(r['other_alt_count']),allele_frequency=str(float(r['allele_frequency']))))
    require(got==expected,'count/panel membership mismatch')
    return dict(status='pass',sites=len(got),total_read_observations=sum(r['ref_count']+r['alt_count']+r['other_alt_count'] for r in got))
def assess_estimator(path,spec):
    rows=read_table(path);require(len(rows)==1,'estimator output row count');row=rows[0]
    require('contamination' in row and 'error' in row,'estimator columns');c,e=float(row['contamination']),float(row['error'])
    require(math.isfinite(c) and math.isfinite(e) and 0<=c<=1 and e>=0,'estimator nonfinite/range')
    fallback=c==0 and e==1;kind=spec['kind'];truth=0 if kind=='same_donor' else float(spec['alpha'])
    if kind in ['empty','uninformative']:
        return dict(classification='ineligible_low_information',contamination=c,error=e,fallback_0_1=fallback,numerical_pass=None)
    if kind in ['diploid','same_donor']:
        bound=.01 if truth==0 else max(.01,.2*truth);require(not fallback,'uninformative fallback on informative diploid');require(abs(c-truth)<=bound,'fixed diploid numerical tolerance failed')
        return dict(classification='component_numeric_pass',contamination=c,error=e,declared_fraction=truth,tolerance=bound,fallback_0_1=fallback)
    return dict(classification='descriptive_confounding_control',contamination=c,error=e,fallback_0_1=fallback,numerical_pass=None)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);p.add_argument('--expected-plan-sha256',required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--outer-guard-pid',type=int,required=True);args=p.parse_args();require(args.outer_guard_pid==os.getppid() and args.outer_guard_pid>1 and os.getpgrp()==os.getpid(),'dedicated outer guard parent/group required')
    plan,planpin,sources,tools=read_plan(args.plan,args.expected_plan_sha256);out=args.output.resolve();require(out.is_relative_to(ROOT/'results/nuclear-mixture-pilot-20260920'),'isolated output scope');require(not out.exists(),'no retry/overwrite')
    limits=plan['resource_limits'];require(shutil.disk_usage(ROOT).free>=limits['minimum_free_bytes'],'free space');out.mkdir(parents=True);started=time.monotonic();deadline=started+limits['total_wall_seconds'];commands=[]
    guard=caller_module();require(plan['expected_original_caller_argv']==guard.EXPECTED_CALLER_ARGV,'caller command contract drift');caller=guard.caller_identity()
    state=dict(status='running',started_at=stamp(),scope=plan['scope'],plan=planpin,sources=sources,runtime=tools,original_caller_before=caller)
    try:
        prefix=[tools['java']['path'],'-Xmx2g','-XX:ActiveProcessorCount=2','-jar',tools['gatk']['path']]
        commands.append(command_run(prefix+['--version'],'version',out,limits,deadline))
        require('4.6.1.0' in (out/'commands/version.stdout.log').read_text()+(out/'commands/version.stderr.log').read_text(),'installed version differs')
        fixture=count_fixture(out/'count-fixture');require(fixture['decoded_payload_equal'],'fictional BAM/CRAM full payload changed')
        counts={}
        for mode in ['bam','cram']:
            table=out/(mode+'-pileups.tsv');cmd=prefix+['GetPileupSummaries','-R',fixture['reference'],'-I',fixture[mode],'-V',fixture['sites'],'-L',fixture['sites'],'--minimum-mapping-quality','50','--max-depth-per-sample','0','-O',str(table)]
            commands.append(command_run(cmd,'pileups-'+mode,out,limits,deadline));counts[mode]=compare_counts(table,fixture['expected_counts'])
        badcmd=prefix+['GetPileupSummaries','-R',fixture['reference'],'-I',fixture['bam'],'-V',fixture['missing_AF_sites'],'-L',fixture['missing_AF_sites'],'-O',str(out/'missing-AF-pileups.tsv')]
        commands.append(command_run(badcmd,'missing-AF-refusal',out,limits,deadline,expected_failure=True))
        require('allele frequency' in (out/'commands/missing-AF-refusal.stderr.log').read_text(),'missing AF failed for unrelated reason')
        scenarios=write_tables(out/'table-fixtures');estimator=[]
        for spec in scenarios:
            table=out/(spec['name']+'-contamination.tsv');segments=out/(spec['name']+'-segments.tsv');cmd=prefix+['CalculateContamination','-I',spec['path'],'-O',str(table),'--tumor-segmentation',str(segments)]
            if spec['kind']=='empty':
                # Tool may refuse an empty table; either result is retained as ineligible, never clean.
                try:commands.append(command_run(cmd,spec['name'],out,limits,deadline))
                except ValueError:
                    receipt=json.loads((out/'commands'/(spec['name']+'.json')).read_text())
                    require(receipt.get('exit_code',0)!=0 and receipt.get('error','').startswith('unexpected command result'),'empty control failed outside expected tool boundary')
                    estimator.append(spec|dict(classification='ineligible_low_information',tool_refused=True));continue
            else:commands.append(command_run(cmd,spec['name'],out,limits,deadline))
            estimator.append(spec|assess_estimator(table,spec))
        for pin in [planpin]+list(sources.values())+list(tools.values()):fingerprint(pin['path'],pin['sha256'])
        caller_after=guard.caller_identity();guard.exact_identity(caller,caller_after)
        save(out/'execution-end.json',dict(checked_at=stamp(),original_caller_before=caller,original_caller_after=caller_after,free_bytes=shutil.disk_usage(ROOT).free))
        state.update(original_caller_after=caller_after,status='complete',completed_at=stamp(),elapsed_seconds=time.monotonic()-started,count_controls=counts,table_controls=estimator,commands=commands,limits=limits,limitations=['Fictional counter controls and estimator tables are separate components, not end-to-end donor read mixtures.','No project sample alignment/genotype or case-derived resource accessed; no actual nuclear estimate.','Same-donor, ancestry, CN/LOH and allelic-imbalance limits remain. No clinical sensitivity, causal or cell-fraction claim.'])
        state['outputs']={str(f.relative_to(out)):fingerprint(f) for f in sorted(out.rglob('*')) if f.is_file()};save(out/'summary.json',state);save(out/'status.json',dict(status='complete',completed_at=stamp(),summary_sha256=fingerprint(out/'summary.json')['sha256']))
    except BaseException as e:
        save(out/'status.json',dict(status='failed',failed_at=stamp(),error_type=type(e).__name__,error=str(e),scope=plan['scope'],retained_without_retry=True));raise
    finally:
        if not (out/'execution-end.json').exists():
            record=dict(checked_at=stamp(),free_bytes=shutil.disk_usage(ROOT).free)
            try:record['original_caller_after']=guard.caller_identity()
            except ValueError as error:record['original_caller_identity_error']=str(error)
            save(out/'execution-end.json',record)
if __name__=='__main__':main()
