#!/usr/bin/env python3
"""Guarded public-reference target preparation only; never opens case inputs."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
import time

SOURCE = Path(__file__).resolve()
VENDOR = SOURCE.parent / 'vendor'
PINS = {'guard_helpers.py': 'e30e5ad4819d488e2509feffff1027bab6102e473bb44ed33bcf4dbfc734ffb6',
        'coding_depth_targets.py': '1b7df219b0552479fd45836ce30fa1a2e9b95c2b315384d71da83bc108418278'}
LIMITS = dict(wall_seconds=1800, group_RSS_bytes=16*1024**3,
              output_bytes=2*1024**3, minimum_free_bytes=40*1024**3)
REAL_INPUTS = {'panel': '8089748ab8645336eedab56e4f548f7c902e2f7d163fa3a559b49c9fbd93151d',
              'gff': '0857407b4a8be8792089a69fd961737256fb7129b6b799f8a8920ebea4fb134f',
              'reference': '187917c00fa7373a054f9768009d19e4f6e5fc87d3a6c125ce131164ca69fcd6',
              'reference_index': 'd81f672f40cfabf2510f37e3f790f2c79665023f4ac950bb9cc6496a40cebf32'}
QUALIFIED = {'coding_depth_targets.py': PINS['coding_depth_targets.py'],
             'measure_coding_depth.py': '6ae51dfa66571f81cfd3d5b55696fd844f12bab45cdd2ac4c7f5052b92389749',
             'survey_chromosome_depth.py': '403718a03c506f70f84d0164d4bb062a461f6a62f6d1c143e85e4829482d2c29',
             'prepare_gatk_compatible_cram.py': '6e1c5f8a830bb48fba7258d810a0f54166f40d471b21c1b3b4f185e9bed650c9'}
REAL_OUTPUT = Path('/Users/jeweled/Documents/mva-hackathon/methods/secondary-panel-depth-preparation-20260920/public-target-v1')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def load_vendor(name):
    path = VENDOR / name
    require(hashlib.sha256(path.read_bytes()).hexdigest() == PINS[name], 'vendor source differs: '+name)
    spec = importlib.util.spec_from_file_location('public_targets_'+name.replace('.', '_'), path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


guard = load_vendor('guard_helpers.py')


def same(a, b, message):
    require(json.dumps(a, sort_keys=True, allow_nan=False) == json.dumps(b, sort_keys=True, allow_nan=False), message)


def configuration(path, expected, output):
    config_pin = guard.fingerprint(path, {'sha256': expected})
    cfg = guard.parse(Path(path).read_text())
    keys = {'schema_version', 'public_preflight_authorized', 'fixture_only', 'inputs', 'genes',
            'expected_genes', 'producer_qualification', 'qualified_sources',
            'expected_gene_identity', 'output', 'approval_note'}
    require(set(cfg) == keys, 'configuration field set differs')
    require(type(cfg['schema_version']) is int and cfg['schema_version'] == 1, 'schema differs')
    require(cfg['public_preflight_authorized'] is True, 'public target generation is not authorized')
    require(type(cfg['fixture_only']) is bool, 'fixture flag must be boolean')
    require(type(cfg['expected_genes']) is int and cfg['expected_genes'] > 0, 'invalid gene count')
    genes = cfg['genes']
    require(type(genes) is list and all(type(g) is str and g for g in genes), 'invalid gene list')
    require(genes == sorted(set(genes)) and len(genes) == cfg['expected_genes'], 'panel list differs')
    require(set(cfg['expected_gene_identity']) == set(genes), 'gene identity membership differs')
    require(set(cfg['inputs']) == set(REAL_INPUTS), 'public-only input membership differs')
    require(set(cfg['qualified_sources']) == set(QUALIFIED), 'qualified source membership differs')
    require(type(cfg['approval_note']) is str and cfg['approval_note'].strip(), 'approval note required')
    require(Path(cfg['output']).resolve() == output.resolve(), 'output path differs')
    require(not output.exists(), 'existing output preserved')
    for pin in list(cfg['inputs'].values()) + list(cfg['qualified_sources'].values()) + [cfg['producer_qualification']]:
        require(set(pin) <= {'path', 'sha256', 'bytes'} and {'path', 'sha256'} <= set(pin), 'invalid pin shape')
        require(type(pin['path']) is str and Path(pin['path']).is_absolute(), 'absolute input path required')
        require(type(pin['sha256']) is str and re.fullmatch('[0-9a-f]{64}', pin['sha256']), 'invalid SHA256')
        if 'bytes' in pin:
            require(type(pin['bytes']) is int and pin['bytes'] >= 0, 'invalid pin size')
    if not cfg['fixture_only']:
        same(cfg['expected_genes'], 84, 'frozen gene count differs')
        same({k:v['sha256'] for k,v in cfg['inputs'].items()}, REAL_INPUTS, 'frozen public inputs differ')
        same({k:v['sha256'] for k,v in cfg['qualified_sources'].items()}, QUALIFIED, 'qualified source differs')
        require(cfg['producer_qualification']['sha256'] == '74398ab873453551e1b2b3aafdf42f032949fcdbb8cc9ca5a9357a1f63f43fa6', 'historical qualification differs')
        require(output.resolve() == REAL_OUTPUT, 'real output root differs')
        require(platform.node().split('.')[0] == '512S2', '512S2 required')
    return cfg, config_pin


def bind_inputs(cfg, config_pin):
    # The only data-input loop is the exact four-public-file whitelist.
    pins = {'config': config_pin}
    for key, expected in cfg['inputs'].items():
        pins[key] = guard.fingerprint(expected['path'], expected)
    for key, expected in cfg['qualified_sources'].items():
        pins['qualified_source/'+key] = guard.fingerprint(expected['path'], expected)
    q = cfg['producer_qualification']
    pins['producer_qualification'] = guard.fingerprint(q['path'], q)
    controls = guard.parse(Path(q['path']).read_text())
    require(controls['status'] == 'complete' and type(controls['exit_code']) is int and controls['exit_code'] == 0, 'qualification incomplete')
    for key, expected in cfg['qualified_sources'].items():
        same(controls['source_sha256']['workflow/scripts/'+key], expected['sha256'], 'qualification source linkage differs')
    same(cfg['qualified_sources']['coding_depth_targets.py']['sha256'], PINS['coding_depth_targets.py'], 'unchanged target helper required')
    pins['executed_source'] = guard.fingerprint(SOURCE)
    for name, digest in PINS.items():
        pins['executed_vendor/'+name] = guard.fingerprint(VENDOR/name, {'sha256': digest})
    require(pins['reference_index']['path'] == pins['reference']['path']+'.fai', 'reference index binding differs')
    panel = guard.parse(Path(pins['panel']['path']).read_text())
    rows = panel['rows']
    require(type(panel['total']) is int and type(rows) is list, 'panel count type differs')
    symbols = [row['symbol'] for row in rows]
    same(sorted(symbols), cfg['genes'], 'frozen panel membership differs')
    require(len(set(symbols)) == len(symbols) == panel['total'] == cfg['expected_genes'], 'panel count or duplicates differ')
    for row in rows:
        same(row['hgnc_id'], cfg['expected_gene_identity'][row['symbol']]['hgnc_id'], 'panel HGNC identity differs')
    return pins, rows


def verify_design(design, cfg, fasta):
    same([g['name'] for g in design['genes']], cfg['genes'], 'selected gene membership differs')
    rank = {name: i for i, name in enumerate(fasta.references)}
    previous = None
    for window in design['windows']:
        identity = (rank[window['chrom']], window['start'])
        require(previous is None or identity > previous[0], 'window order differs')
        if previous and window['chrom'] == previous[1]['chrom']:
            require(window['start'] >= previous[1]['end'], 'overlapping target windows')
        previous = (identity, window)
        sequence = fasta.fetch(window['chrom'], window['start'], window['end']).upper()
        require(len(sequence) == window['end']-window['start'], 'reference window length differs')
        window['reference_ACGT_bases'] = sum(sequence.count(base) for base in 'ACGT')
        window['reference_non_ACGT_bases'] = len(sequence)-window['reference_ACGT_bases']
        window['reference_sequence_sha256'] = hashlib.sha256(sequence.encode()).hexdigest()
    for gene in design['genes']:
        same(gene['id'], cfg['expected_gene_identity'][gene['name']]['gene_id'], 'gene identity differs')
        windows = [w for w in design['windows'] if gene['id'] in w['gene_ids']]
        same(sum(w['end']-w['start'] for w in windows), gene['CDS_union_bases'], 'gene denominator differs')
        sequence = ''.join(fasta.fetch(gene['chrom'], a, b).upper() for a,b in gene['cds_union'])
        gene['reference_ACGT_bases'] = sum(sequence.count(base) for base in 'ACGT')
        gene['reference_non_ACGT_bases'] = len(sequence)-gene['reference_ACGT_bases']
        gene['concatenated_CDS_union_sha256'] = hashlib.sha256(sequence.encode()).hexdigest()
    same(sum(w['end']-w['start'] for w in design['windows']), design['panel_union_bases'], 'panel union denominator differs')
    same(sum(g['CDS_union_bases'] for g in design['genes']), design['gene_union_bases_sum'], 'gene sum differs')


def prepare(cfg, config_pin, output):
    pins, rows = bind_inputs(cfg, config_pin)
    targets = load_vendor('coding_depth_targets.py')
    import pysam
    with pysam.FastaFile(pins['reference']['path']) as fasta:
        design = targets.prepare_targets(pins['gff']['path'], cfg['genes'], dict(zip(fasta.references, fasta.lengths)))
        verify_design(design, cfg, fasta)
    # Capture the literal gene feature, including version and HGNC description.
    import gzip
    from urllib.parse import unquote
    ordinals = {g['ordinal']:g for g in design['genes']}
    with gzip.open(pins['gff']['path'], 'rt') as stream:
        for ordinal, line in enumerate(stream, 1):
            if ordinal not in ordinals:
                continue
            gene = ordinals[ordinal]
            fields = line.rstrip('\n').split('\t')
            attrs = dict((k,unquote(v)) for k,v in (item.split('=',1) for item in fields[8].split(';') if '=' in item))
            hgnc = re.findall(r'Acc:(HGNC:\d+)', attrs.get('description',''))
            same(hgnc, [cfg['expected_gene_identity'][gene['name']]['hgnc_id']], 'GFF HGNC identity differs')
            gene['original_gene_feature'] = line.rstrip('\n')
    require(all('original_gene_feature' in g for g in design['genes']), 'gene source row missing')
    design.update(status='prepared', prepared_at_utc=guard.stamp(), coordinates='0-based half-open',
                  reads_inspected=False, panel_rows=rows, inputs=pins)
    output.mkdir(parents=True, exist_ok=False)
    guard.save(output/'reference-design.json', design)
    (output/'regions.bed').write_text(''.join(f"{w['chrom']}\t{w['start']}\t{w['end']}\n" for w in design['windows']))
    summary = dict(status='complete', completed_at_utc=guard.stamp(), reads_inspected=False,
        genes=len(design['genes']), windows=len(design['windows']), panel_union_bases=design['panel_union_bases'],
        gene_union_bases_sum=design['gene_union_bases_sum'],
        panel_ACGT_bases=sum(w['reference_ACGT_bases'] for w in design['windows']),
        panel_non_ACGT_bases=sum(w['reference_non_ACGT_bases'] for w in design['windows']),
        gene_summary=[{k:g[k] for k in ['name','id','chrom','CDS_union_bases','reference_ACGT_bases','reference_non_ACGT_bases']} for g in design['genes']],
        scope='Public all-transcript annotated CDS union only; no alignment/index/callset, depth, genotype, functional consequence or clinical callability.')
    guard.save(output/'target-summary.json', summary)
    for pin in pins.values():
        guard.fingerprint(pin['path'], pin)
    guard.save(output/'status.json', dict(status='complete', completed_at_utc=guard.stamp(), reads_inspected=False,
        inputs=pins, outputs={p.name:guard.fingerprint(p) for p in sorted(output.iterdir())},
        runtime=dict(host=platform.node(),python=sys.version,executable=sys.executable,pysam=pysam.__version__)))


def run(args):
    guard.fingerprint(SOURCE, {'sha256': args.expected_source_sha256})
    output = args.output.resolve()
    cfg, config_pin = configuration(args.config, args.expected_config_sha256, output)
    launch = output.with_name(output.name+'-launch')
    require(not launch.exists(), 'existing launch attempt preserved')
    require(shutil.disk_usage(output.parent).free >= LIMITS['minimum_free_bytes'], 'insufficient disk')
    launch.mkdir()
    state = dict(status='launching', started_at_utc=guard.stamp(), limits=LIMITS, sampled_every_seconds=2,
                 config=config_pin, source=guard.fingerprint(SOURCE), before_processes=guard.snapshot(), reads_inspected=False)
    guard.save(launch/'status.json', state)
    command = [sys.executable,'-B',str(SOURCE),'--config',str(args.config),'--expected-config-sha256',args.expected_config_sha256,
               '--expected-source-sha256',args.expected_source_sha256,'--output',str(output),'--internal-worker']
    environment = os.environ.copy(); environment['MVA_PUBLIC_PREFLIGHT_PARENT'] = str(os.getpid())
    process = None; started = time.monotonic()
    try:
        with (launch/'stdout.log').open('xb') as stdout, (launch/'stderr.log').open('xb') as stderr:
            process = subprocess.Popen(command, stdout=stdout, stderr=stderr, env=environment, start_new_session=True)
            state.update(status='running', pid=process.pid, pgid=process.pid, command=command)
            while True:
                rows=guard.snapshot(); elapsed=time.monotonic()-started; rss=guard.owned_rss(rows,process.pid)
                size=guard.tree_bytes(output)+guard.tree_bytes(launch); free=shutil.disk_usage(output.parent).free
                state.update(elapsed_seconds=elapsed,group_RSS_bytes=rss,output_bytes=size,free_disk_bytes=free)
                reason=guard.limit_reason(elapsed,rss,size,free,LIMITS)
                require(reason is None, 'guard limit '+str(reason))
                code=process.poll()
                if code is not None:
                    require(not any(r['pgid']==process.pid and r['pid']!=process.pid for r in rows), 'live group remains')
                    require(code==0, 'worker exit '+str(code)); break
                guard.save(launch/'status.json',state); time.sleep(2)
        terminal=guard.parse((output/'status.json').read_text())
        require(terminal['status']=='complete' and terminal['reads_inspected'] is False, 'target receipt incomplete')
        same(terminal['inputs']['config'],config_pin,'worker configuration binding differs')
        for name,pin in terminal['outputs'].items():guard.fingerprint(output/name,pin)
        state.update(status='complete',completed_at_utc=guard.stamp(),exit_code=0,result_status=guard.fingerprint(output/'status.json'))
    except BaseException as error:
        if process is not None:guard.stop_owned(process)
        state.update(status='failed',failed_at_utc=guard.stamp(),error=str(error),exit_code=process.returncode if process else None)
        raise
    finally:
        state['after_processes']=guard.snapshot();guard.save(launch/'status.json',state)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',required=True,type=Path);parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--expected-config-sha256',required=True);parser.add_argument('--expected-source-sha256',required=True)
    parser.add_argument('--internal-worker',action='store_true',help=argparse.SUPPRESS)
    args=parser.parse_args()
    if args.internal_worker:
        require(os.environ.get('MVA_PUBLIC_PREFLIGHT_PARENT')==str(os.getppid()) and os.getpgrp()==os.getpid(), 'worker requires owned parent/session')
        guard.fingerprint(SOURCE,{'sha256':args.expected_source_sha256})
        cfg,pin=configuration(args.config,args.expected_config_sha256,args.output)
        os.nice(10);prepare(cfg,pin,args.output)
    else:run(args)


if __name__=='__main__':main()
