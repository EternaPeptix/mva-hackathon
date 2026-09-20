#!/usr/bin/env python3
"""Bind a frozen public CDS target to the unchanged qualified depth producer.

Only a separately approved configuration can start a new owned process group.
This wrapper does not change extraction, quality gates or depth calculation.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time

SOURCE = Path(__file__).resolve()
VENDOR = SOURCE.parent / 'vendor/guard_helpers.py'
VENDOR_SHA = 'e30e5ad4819d488e2509feffff1027bab6102e473bb44ed33bcf4dbfc734ffb6'
if hashlib.sha256(VENDOR.read_bytes()).hexdigest() != VENDOR_SHA:
    raise ValueError('guard helper source differs')
spec = importlib.util.spec_from_file_location('secondary_depth_guard_helpers', VENDOR)
guard = importlib.util.module_from_spec(spec)
spec.loader.exec_module(guard)
require = guard.require
LIMITS = dict(wall_seconds=5400, group_RSS_bytes=16*1024**3,
              output_bytes=20*1024**3, minimum_free_bytes=40*1024**3)
ROOT = Path('/Users/jeweled/Documents/mva-hackathon')
REAL_OUTPUT = ROOT/'results/secondary-panel-cds-depth-20260920/review-v2'
QUALIFIED = {
 'measure_coding_depth.py':'6ae51dfa66571f81cfd3d5b55696fd844f12bab45cdd2ac4c7f5052b92389749',
 'coding_depth_targets.py':'1b7df219b0552479fd45836ce30fa1a2e9b95c2b315384d71da83bc108418278',
 'survey_chromosome_depth.py':'403718a03c506f70f84d0164d4bb062a461f6a62f6d1c143e85e4829482d2c29',
 'prepare_gatk_compatible_cram.py':'6e1c5f8a830bb48fba7258d810a0f54166f40d471b21c1b3b4f185e9bed650c9'}
MODES = [dict(name='MQ0_BQ1', minimum_MQ=0, minimum_BQ=1),
         dict(name='MQ20_BQ20', minimum_MQ=20, minimum_BQ=20),
         dict(name='MQ60_BQ20', minimum_MQ=60, minimum_BQ=20)]
PUBLIC_NAMES = {'config','status','launch','reference-design.json','target-summary.json','regions.bed'}
CONTROL_PATHS = {'scripts/launch_secondary_depth_v2.py','tests/test_secondary_depth_guard_v2.py',
                 'scripts/qualify_secondary_depth_guard_v2.py','scripts/vendor/guard_helpers.py'}
CONTROLS = 21


def encoded(value):
    return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n').encode()


def same(a, b, message):
    require(encoded(a) == encoded(b), message)


def core_design(design):
    extras = {'concatenated_CDS_union_sha256','original_gene_feature','reference_ACGT_bases','reference_non_ACGT_bases'}
    return dict(genes=[{k:v for k,v in g.items() if k not in extras} for g in design['genes']],
                windows=[{k:v for k,v in w.items() if not k.startswith('reference_')} for w in design['windows']],
                panel_union_bases=design['panel_union_bases'], gene_union_bases_sum=design['gene_union_bases_sum'])


def shape(cfg, output, authorized=True):
    keys = {'schema_version','ready','preparation_note','fixture_only','output','limits',
            'producer_configuration','producer_config_sha256','qualified_sources','producer_qualification',
            'public_preflight','runtime','completed_caller_status','caller_definition','expected_target'}
    same(sorted(cfg), sorted(keys), 'configuration membership differs')
    same(cfg['schema_version'], 1, 'schema differs')
    require(type(cfg['ready']) is bool and (cfg['ready'] or not authorized), 'case depth not authorized')
    require(type(cfg['fixture_only']) is bool, 'fixture flag differs')
    require(type(cfg['preparation_note']) is str and cfg['preparation_note'].strip(), 'authorization note missing')
    same(cfg['limits'], LIMITS, 'resource limits differ')
    same(str(Path(cfg['output']).resolve()), str(output.resolve()), 'output differs')
    require(not output.exists() and not output.with_name(output.name+'-launch').exists(), 'existing attempt preserved')
    c = cfg['producer_configuration']
    same(hashlib.sha256(encoded(c)).hexdigest(), cfg['producer_config_sha256'], 'inner config bytes differ')
    same(c['schema_version'],1,'producer schema differs');same(c['assembly'],'GRCh38','assembly differs')
    same(c['depth_modes'],MODES,'quality modes differ')
    require(type(c['expected_read_groups']) is int and c['expected_read_groups'] > 0, 'read groups differ')
    require(type(c['genes']) is list and c['genes'] == sorted(set(c['genes'])) and c['genes'], 'gene list differs')
    same(sorted(c['inputs']), sorted(['alignment','alignment_index','reference','reference_index','gff','samtools']), 'producer input membership differs')
    same(sorted(cfg['public_preflight']), sorted(PUBLIC_NAMES), 'public target receipt set differs')
    same(sorted(cfg['qualified_sources']), sorted(QUALIFIED), 'source membership differs')
    same(c['source_sha256'],{k:v['sha256'] for k,v in cfg['qualified_sources'].items()},'source config linkage differs')
    for p in [cfg['runtime'], cfg['completed_caller_status'],cfg['caller_definition'],cfg['producer_qualification'],
              *cfg['public_preflight'].values(),*cfg['qualified_sources'].values(),*c['inputs'].values(),c['gene_list_source']]:
        require(type(p) is dict and {'path','sha256'} <= set(p) <= {'path','sha256','bytes'}, 'pin shape differs')
        require(type(p['path']) is str and Path(p['path']).is_absolute(), 'absolute pin path required')
        require(type(p['sha256']) is str and len(p['sha256'])==64 and all(x in '0123456789abcdef' for x in p['sha256']), 'invalid hash')
        if 'bytes' in p:require(type(p['bytes']) is int and p['bytes'] >= 0,'invalid byte count')
    if not cfg['fixture_only']:
        same(c['source_sha256'],QUALIFIED,'unchanged qualified source required')
        same(str(output.resolve()),str(REAL_OUTPUT),'real output differs')
        same(platform.node().split('.')[0],'512S2','512S2 only')
        same(cfg['expected_target'],dict(genes=84,windows=2112,panel_union_bases=393880,gene_union_bases_sum=393880,panel_ACGT_bases=393880,panel_non_ACGT_bases=0),'measured scope differs')
        same(cfg['public_preflight']['status']['sha256'],'3673bafb827a6b1d6f6ac90b08b837c7053373fa5c3c21e4645f113b2ecf5fc2','frozen public status differs')
        same(cfg['public_preflight']['reference-design.json']['sha256'],'f9b6ab0835a979150ac2d1e32125ed83cec2ac5c430102b5fa7eba46784e5775','frozen design differs')
        same(cfg['producer_qualification']['sha256'],'74398ab873453551e1b2b3aafdf42f032949fcdbb8cc9ca5a9357a1f63f43fa6','historical controls differ')
        same(cfg['completed_caller_status']['sha256'],'d58b058470ab989ac27257699705401a58a6a81c587b09c4f3018ee64a6788e0','completion identity differs')
        same(cfg['caller_definition']['sha256'],'6bb6aca4c6e281ccbe345b6f7a4e7a38f340da8a9138b3cd8faec87ccc532267','caller definition differs')


def bind(cfg):
    """Read only source/runtime and public/completion metadata before producer opens case data."""
    pins = {}
    def read(name, pin):
        pins[name] = guard.fingerprint(pin['path'],pin)
        return guard.parse(Path(pin['path']).read_text())
    c = cfg['producer_configuration']
    pub = {k:read('public/'+k,p) for k,p in cfg['public_preflight'].items() if k != 'regions.bed'}
    bedpin = cfg['public_preflight']['regions.bed'];pins['public/regions.bed']=guard.fingerprint(bedpin['path'],bedpin)
    state, design, summary, publiccfg = (pub[k] for k in ('status','reference-design.json','target-summary.json','config'))
    require(state['status']=='complete' and state['reads_inspected'] is False, 'public preparation not complete')
    require(summary['status']=='complete' and summary['reads_inspected'] is False and design['reads_inspected'] is False, 'public target is not reference only')
    require(pub['launch']['status']=='complete' and type(pub['launch']['exit_code']) is int and pub['launch']['exit_code']==0,'public guard incomplete')
    guard.fingerprint(cfg['public_preflight']['status']['path'],pub['launch']['result_status'])
    guard.fingerprint(cfg['public_preflight']['config']['path'],state['inputs']['config'])
    for name in ('reference-design.json','target-summary.json','regions.bed'):
        guard.fingerprint(cfg['public_preflight'][name]['path'],state['outputs'][name])
    same(sorted(state['outputs']),sorted(['reference-design.json','target-summary.json','regions.bed']),'public output set differs')
    same(c['genes'],publiccfg['genes'],'public gene selection differs')
    same(c['genes'],[g['name'] for g in design['genes']],'target gene membership differs')
    same(c['gene_list_source']['sha256'],state['inputs']['panel']['sha256'],'panel linkage differs')
    same(str(Path(c['gene_list_source']['path']).resolve()),state['inputs']['panel']['path'],'panel path differs')
    for name in ('gff','reference','reference_index'):
        for field in ('path','sha256'):same(c['inputs'][name][field],state['inputs'][name][field],'public input differs '+name)
    same({k:summary[k] for k in cfg['expected_target']},cfg['expected_target'],'public denominator differs')
    same(len(design['genes']),cfg['expected_target']['genes'],'design gene count differs')
    same(len(design['windows']),cfg['expected_target']['windows'],'design window count differs')
    for name in ('panel_union_bases','gene_union_bases_sum'):same(design[name],cfg['expected_target'][name],'design denominator differs')
    same(Path(bedpin['path']).read_text(),''.join(f"{w['chrom']}\t{w['start']}\t{w['end']}\n" for w in design['windows']),'public BED differs')
    q=read('historical_qualification',cfg['producer_qualification'])
    require(q['status']=='complete' and type(q['exit_code']) is int and q['exit_code']==0,'producer qualification failed')
    for name,pin in cfg['qualified_sources'].items():
        same(q['source_sha256']['workflow/scripts/'+name],pin['sha256'],'producer qualification link differs')
        pins['source/'+name]=guard.fingerprint(pin['path'],pin)
    dirs={str(Path(v['path']).resolve().parent) for v in cfg['qualified_sources'].values()}
    require(len(dirs)==1 and all(Path(v['path']).name==k for k,v in cfg['qualified_sources'].items()),'qualified source layout differs')
    caller=read('completed_caller',cfg['completed_caller_status']);definition=read('caller_definition',cfg['caller_definition'])
    require(caller['status']=='complete','original caller not complete')
    same(caller['expected_shards'],64,'expected caller scope differs');same(caller['completed_shards'],list(range(64)),'completed shard IDs differ')
    same(len(definition['shards']),64,'definition scope differs');same(caller['definition_sha256'],cfg['caller_definition']['sha256'],'definition link differs')
    for key,target in [('alignment_sha256','alignment'),('reference_sha256','reference')]:same(definition[key],c['inputs'][target]['sha256'],'original source differs')
    same(definition['sample'],c['sample'],'original sample differs')
    pins['runtime']=guard.fingerprint(cfg['runtime']['path'],cfg['runtime'])
    same(str(Path(sys.executable).resolve()),pins['runtime']['path'],'runtime path differs')
    return design,pins


def qualification(path, expected):
    pin=guard.fingerprint(path,{'sha256':expected});q=guard.parse(Path(path).read_text())
    same(q['tests_run'],CONTROLS,'wrapper control count differs')
    require(q['status']=='complete' and q['skipped']==[] and all(type(q[k]) is int and q[k]==0 for k in ('errors','failures','exit_code')),'wrapper controls not wholly passing')
    same(sorted(q['source_files']),sorted(CONTROL_PATHS),'wrapper qualification membership differs')
    workflow=SOURCE.parents[1]
    for name,value in q['source_files'].items():guard.fingerprint(workflow/name,value)
    return pin


def terminal(cfg, output, inner_pin, design):
    s=guard.parse((output/'status.json').read_text())
    require(s['status']=='complete','depth producer incomplete')
    guard.fingerprint(output/'summary.json',{'sha256':s['summary_sha256']})
    summary=guard.parse((output/'summary.json').read_text())
    require(summary['status']=='complete','depth summary incomplete')
    same(s['inputs'],summary['inputs'],'producer input receipts differ')
    guard.fingerprint(inner_pin['path'],s['inputs']['config'])
    same(s['inputs']['config']['sha256'],inner_pin['sha256'],'producer configuration differs')
    names={'selection.json','regions.bed','selected.bam','selected.bam.bai','commands.json','extraction-counts.json','extract-preserving-records.log','quickcheck.log','index.log','summary.json'}
    names.update(m['name']+suffix for m in MODES for suffix in ('.json','.depth.tsv.gz','.log'))
    same(sorted(s['outputs']),sorted(names),'depth output set differs')
    for name,pin in s['outputs'].items():guard.fingerprint(output/name,pin)
    selected=guard.parse((output/'selection.json').read_text())
    same(selected['config'],cfg['producer_configuration'],'selection config differs')
    same(core_design(selected),core_design(design),'producer target drift')
    require((output/'regions.bed').read_bytes()==Path(cfg['public_preflight']['regions.bed']['path']).read_bytes(),'producer BED drift')
    for key in ('panel_union_bases','gene_union_bases_sum'):same(summary[key],cfg['expected_target'][key],'depth summary denominator differs')
    same(sorted(summary['depth_modes']),sorted(m['name'] for m in MODES),'summary modes differ')
    for m in MODES:
        v=summary['depth_modes'][m['name']];same(v['mode'],m,'mode metadata differs');same(v['position_count'],cfg['expected_target']['panel_union_bases'],'position count differs')
    return dict(status=guard.fingerprint(output/'status.json'),summary=guard.fingerprint(output/'summary.json'))


def worker(args, cfg, config_pin, launch):
    os.nice(10)  # The unchanged producer adds its own increment before extraction.
    design,pins=bind(cfg)
    qualification(args.qualification,args.expected_qualification_sha256)
    inner=launch/'producer-config.json';inner.write_bytes(encoded(cfg['producer_configuration']))
    inner_pin=guard.fingerprint(inner,{'sha256':cfg['producer_config_sha256']})
    pins['outer_config']=config_pin;pins['wrapper_source']=guard.fingerprint(SOURCE)
    binding=dict(status='bound',before_case_access=True,inputs=pins,producer_config=inner_pin,expected_target=cfg['expected_target'])
    guard.save(launch/'worker-binding.json',binding)
    cmd=[cfg['runtime']['path'],'-B',cfg['qualified_sources']['measure_coding_depth.py']['path'],'--config',str(inner),'--expected-config-sha256',inner_pin['sha256'],'--output',str(args.output.resolve())]
    completed=subprocess.run(cmd,check=False)
    require(completed.returncode==0,'unchanged producer failed '+str(completed.returncode))
    results=terminal(cfg,args.output.resolve(),inner_pin,design)
    # Fresh hashes after work; preserve producer's own input and output receipts.
    for pin in pins.values():guard.fingerprint(pin['path'],pin)
    binding.update(status='complete',completed_at_utc=guard.stamp(),producer_results=results,command=cmd)
    guard.save(launch/'worker-binding.json',binding)


def run(args):
    guard.fingerprint(SOURCE,{'sha256':args.expected_source_sha256})
    config_pin=guard.fingerprint(args.config,{'sha256':args.expected_config_sha256})
    cfg=guard.parse(args.config.read_text());output=args.output.resolve();shape(cfg,output)
    draft=guard.parse(args.proposal.read_text());guard.fingerprint(args.proposal,{'sha256':args.expected_proposal_sha256})
    guard.approved_delta(draft,cfg)
    qpin=qualification(args.qualification,args.expected_qualification_sha256)
    # Parent performs only small source/config reads. All public/input hashing is in the guarded worker.
    disk_root=output.parent
    while not disk_root.exists():disk_root=disk_root.parent
    require(shutil.disk_usage(disk_root).free>=LIMITS['minimum_free_bytes'],'insufficient free disk')
    output.parent.mkdir(parents=True,exist_ok=True)
    launch=output.with_name(output.name+'-launch');launch.mkdir()
    state=dict(status='launching',started_at_utc=guard.stamp(),config=config_pin,source=guard.fingerprint(SOURCE),qualification=qpin,limits=LIMITS,sampled_every_seconds=2,before_processes=guard.snapshot())
    guard.save(launch/'status.json',state)
    command=[sys.executable,'-B',str(SOURCE),*sys.argv[1:],'--internal-worker']
    env=os.environ.copy();env['MVA_SECONDARY_DEPTH_PARENT']=str(os.getpid())
    process=None;started=time.monotonic();peaks=dict(group_RSS_bytes=0,output_bytes=0)
    try:
        with (launch/'stdout.log').open('xb') as stdout,(launch/'stderr.log').open('xb') as stderr:
            process=subprocess.Popen(command,stdout=stdout,stderr=stderr,env=env,start_new_session=True)
            state.update(status='running',pid=process.pid,pgid=process.pid,command=command)
            while True:
                rows=guard.snapshot();elapsed=time.monotonic()-started;rss=guard.owned_rss(rows,process.pid)
                size=guard.tree_bytes(output)+guard.tree_bytes(launch);free=shutil.disk_usage(output.parent).free
                peaks={k:max(peaks[k],v) for k,v in [('group_RSS_bytes',rss),('output_bytes',size)]}
                state.update(elapsed_seconds=elapsed,group_RSS_bytes=rss,output_bytes=size,free_disk_bytes=free,peaks=peaks)
                reason=guard.limit_reason(elapsed,rss,size,free,LIMITS);require(reason is None,'guard limit '+str(reason))
                code=process.poll()
                if code is not None:
                    require(not any(r['pgid']==process.pid and r['pid']!=process.pid for r in rows),'owned children remain')
                    require(code==0,'worker exit '+str(code));break
                guard.save(launch/'status.json',state);time.sleep(2)
        done=guard.parse((launch/'worker-binding.json').read_text());require(done['status']=='complete','worker binding incomplete')
        state.update(status='complete',completed_at_utc=guard.stamp(),exit_code=0,worker_binding=guard.fingerprint(launch/'worker-binding.json'),producer_status=guard.fingerprint(output/'status.json'))
    except BaseException as error:
        if process is not None:guard.stop_owned(process)
        state.update(status='failed',failed_at_utc=guard.stamp(),error=str(error),exit_code=process.returncode if process else None)
        raise
    finally:
        state['after_processes']=guard.snapshot();guard.save(launch/'status.json',state)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('config','proposal','output','qualification'):p.add_argument('--'+name,type=Path,required=True)
    for name in ('config','proposal','source','qualification'):p.add_argument('--expected-'+name+'-sha256',required=True)
    p.add_argument('--internal-worker',action='store_true',help=argparse.SUPPRESS);args=p.parse_args()
    if args.internal_worker:
        require(os.environ.get('MVA_SECONDARY_DEPTH_PARENT')==str(os.getppid()) and os.getpgrp()==os.getpid(),'worker requires owned parent/session')
        guard.fingerprint(SOURCE,{'sha256':args.expected_source_sha256});pin=guard.fingerprint(args.config,{'sha256':args.expected_config_sha256})
        cfg=guard.parse(args.config.read_text())
        # Existing launch belongs to this parent; still refuse an existing producer output.
        require(not args.output.exists(),'existing producer output preserved')
        worker(args,cfg,pin,args.output.resolve().with_name(args.output.name+'-launch'))
    else:run(args)


if __name__=='__main__':main()
