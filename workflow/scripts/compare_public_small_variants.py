#!/usr/bin/env python3
"""Run the frozen, genotype-aware RTG public comparison after all evidence gates.

This compares complete raw chr20 callsets from the unchanged GATK/bcftools
workflow. No call/GT/QUAL edits or benchmark-specific parameter tuning occur.
"""
import argparse
from collections import Counter
import json
import os
from pathlib import Path
import shutil
import time

import pysam

from rtg_small_variant_helpers import audit_vcf, comparison_command, count_outputs, execute, fingerprint, runtime, save, stamp


def verified(entry):
    current = fingerprint(Path(entry['path']))
    if current['sha256'] != entry['sha256']: raise ValueError('completed artifact changed: ' + entry['path'])
    return current


def identities(path):
    rows = Counter()
    with pysam.VariantFile(str(path)) as vcf:
        for record in vcf:
            sample = next(iter(record.samples))
            rows[(record.contig,record.pos,record.ref,record.alts,record.samples[sample].get('GT'))] += 1
    return rows


def reconcile_native(result):
    """Unthresholded native row must equal independent split-VCF record counts."""
    expected = [result['counts'][k] for k in ('tp_baseline','tp_calls','fp','fn')]
    if expected == [0,0,0,0] and result['native_summary'].strip() == '0 total baseline variants, no summary statistics available':
        return {'native_empty_scope_notice_matches_empty_split_VCFs': True,
                'requested_QUAL_threshold': 0, 'native_reported_observed_score': None}
    rows = [line.split() for line in result['native_summary'].splitlines() if line.strip().startswith('None')]
    if len(rows) != 1: raise ValueError('one native unthresholded summary row required')
    if [int(value) for value in rows[0][1:5]] != expected:
        raise ValueError('native summary and independent VCF counts disagree')
    # A fixed-score row should also be unthresholded because QUAL >= 0 was checked.
    numeric_rows = []
    for line in result['native_summary'].splitlines():
        cells = line.split()
        if not cells: continue
        try: float(cells[0])
        except ValueError: continue
        numeric_rows.append(cells)
    if not numeric_rows and expected[1] + expected[2] == 0:
        return {'native_unthresholded_matches_split_VCF_counts': True,
                'requested_QUAL_threshold': 0, 'native_reported_observed_score': None,
                'no_scored_calls': True}
    if len(numeric_rows) != 1 or [int(v) for v in numeric_rows[0][1:5]] != expected:
        raise ValueError('fixed QUAL threshold does not include exactly all scored calls')
    return {'native_unthresholded_and_fixed_score_match_split_VCF_counts': True,
            'requested_QUAL_threshold': 0, 'native_reported_observed_score': numeric_rows[0][0]}


def run(root, output, config_path, prepared, controls, pilot):
    if output.resolve().is_relative_to(Path(__file__).resolve().parents[2]):
        raise ValueError('large comparison files must remain outside the repository')
    # Missing dependencies do not create a misleading running/failed experiment directory.
    for folder in (prepared,controls,pilot):
        if not (folder/'summary.json').is_file(): raise ValueError('dependency not complete: '+str(folder))
    config = json.loads(config_path.read_text()); config_hash = fingerprint(config_path)['sha256']
    reference_run = json.loads((prepared/'summary.json').read_text())
    controls_run = json.loads((controls/'summary.json').read_text())
    pilot_run = json.loads((pilot/'summary.json').read_text())
    for entry in (reference_run, controls_run, pilot_run):
        if entry['status'] != 'complete': raise ValueError('dependency is not complete')
    if any(entry['config']['sha256'] != config_hash for entry in (reference_run,controls_run)):
        raise ValueError('comparison configuration differs from preparation/controls')
    if controls_run['cases_passed'] != 12 or controls_run['cases_total'] != 12 or controls_run['failed_cases']:
        raise ValueError('all 12 native controls must pass')
    if not reference_run['SDF_roundtrip_chr20_all_64444167_bases_identical']:
        raise ValueError('complete chromosome reference identity not established')
    output.mkdir(parents=True,exist_ok=False); started=time.monotonic()
    state={'status':'running','pid':os.getpid(),'started_at':stamp(),'stage':'verify completed inputs'}
    save(output/'status.json',state)
    try:
        rtg, runtime_files = runtime(root,config)
        reference=Path(reference_run['reference']['path']); baseline=Path(reference_run['truth']['path'])
        bed=Path(reference_run['bed']['path']); sdf=prepared/'reference.sdf'
        inputs={name:verified(reference_run[name]) for name in ('reference','truth','truth_index','bed')}
        for entry in reference_run['sdf_files'].values(): verified(entry)
        for name in ('gatk','bcftools'):
            inputs[name]=verified(pilot_run['outputs'][name])
        for entry in controls_run['runtime'].values(): verified(entry)
        sources={name:fingerprint(folder/'summary.json') for name,folder in [('preparation',prepared),('controls',controls),('native_calling',pilot)]}
        shutil.copyfile(config_path,output/'frozen-comparison-config.json')
        shutil.copyfile(Path(__file__),output/'executed-compare_public_small_variants.py')
        baseline_audit=audit_vcf(baseline,reference,'HG002','chr20')
        if baseline_audit!=reference_run['truth_audit']: raise ValueError('truth audit changed')
        baseline_ids=identities(baseline); results={}
        for method in ('gatk','bcftools'):
            state['stage']='genotype-aware '+method+' comparison';save(output/'status.json',state)
            calls=Path(inputs[method]['path']); call_audit=audit_vcf(calls,reference,'HG002','chr20',require_raw_qual=True)
            target=output/method
            command=execute(comparison_command(rtg,config,baseline,calls,sdf,bed,target,'chr20'),output,method+'-compare')
            result=count_outputs(target);native_reconciliation=reconcile_native(result)
            tpbase=identities(target/'tp-baseline.vcf.gz');fn=identities(target/'fn.vcf.gz')
            tpcall=identities(target/'tp.vcf.gz');fp=identities(target/'fp.vcf.gz');call_ids=identities(calls)
            if (tpbase+fn)-baseline_ids or (tpcall+fp)-call_ids:
                raise ValueError('comparison changed or added original variant/genotype identities')
            if tpbase & fn or tpcall & fp: raise ValueError('record appears in both success and failure categories')
            result.update(command=command,input_audit=call_audit,native_reconciliation=native_reconciliation,
                          original_variant_and_GT_identity_preserved=True,
                          input_calls_not_in_scored_outputs=sum(call_ids.values())-sum((tpcall+fp).values()),
                          input_truth_not_in_scored_outputs=sum(baseline_ids.values())-sum((tpbase+fn).values()),
                          ignored_scope_note='Records outside evaluation regions or ignored by native FILTER/GT/allele semantics are not accuracy denominators. Exact native logs and input class/GT/FILTER counts are retained.')
            results[method]=result;save(output/(method+'-summary.json'),result)
        # Both methods must have used the same truth record universe, not just equal counts.
        universes={method:identities(output/method/'tp-baseline.vcf.gz')+identities(output/method/'fn.vcf.gz') for method in results}
        if universes['gatk']!=universes['bcftools']:
            raise ValueError('the two comparator truth universes differ; investigate boundary handling before cross-caller comparison')
        summary={'status':'complete','completed_at':stamp(),'elapsed_seconds':time.monotonic()-started,
                 'scope':config['scope'],'config':fingerprint(config_path),'inputs':inputs,
                 'source_receipts':sources,'runtime':runtime_files,'truth_audit':baseline_audit,
                 'benchmark_bed_intervals':reference_run['bed_intervals'],'benchmark_bed_bases':reference_run['bed_bases'],
                 'callers':results,'identical_scored_truth_universe':True,
                 'limits':config['limits'],'code':fingerprint(Path(__file__)),
                 'helpers':fingerprint(Path(__file__).with_name('rtg_small_variant_helpers.py'))}
        save(output/'summary.json',summary);state.update(status='complete',completed_at=summary['completed_at']);save(output/'status.json',state)
        return summary
    except BaseException as error:
        state.update(status='failed',failed_at=stamp(),error_type=type(error).__name__,error=str(error));save(output/'status.json',state);raise


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('root','output','config','prepared','controls','pilot'):parser.add_argument('--'+key,required=True,type=Path)
    args=parser.parse_args();result=run(args.root,args.output,args.config,args.prepared,args.controls,args.pilot)
    print(json.dumps({'status':result['status'],'counts':{k:v['counts'] for k,v in result['callers'].items()}}))
