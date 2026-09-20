#!/usr/bin/env python3
"""Exercise the guarded structural runner on retained synthetic reads and known masks."""
import argparse
import json
import os
from pathlib import Path
import time

import pysam

from build_structural_exclusions import fingerprint, parse_exclusions, query_rows, save, stamp, write_bed
from run_structural_variants import NATIVE_OPTIONS, run_analysis


def fixture_config(experiment, mask_folder, scenario, delly, samtools):
    reference = experiment / 'fixture/reference.fa'
    truth = json.loads((experiment / 'fixture/truth.json').read_text())
    original = json.loads((experiment / 'summary.json').read_text())
    if truth['scope'] != 'Synthetic diploid simple-SV reference; no patient sequence or observed biological effect.':
        raise ValueError('retained synthetic fixture required')
    if (original['truth_sha256'] != fingerprint(experiment / 'fixture/truth.json')['sha256']
            or truth['files']['reference.fa']['sha256'] != fingerprint(reference)['sha256']):
        raise ValueError('retained synthetic truth/reference changed')
    with pysam.FastaFile(str(reference)) as fasta: dictionary = dict(zip(fasta.references, fasta.lengths))
    choices = {
        'control': ('3\t0\t1\n', ['1', '2', '3']),
        'deletion_window': ('1\t140000\t165000\n', ['1', '2', '3']),
        'exclude_partner': ('2\n', ['1', '3']),
        'exclude_changed_chromosome': ('1\n', ['2', '3']),
    }
    text, chroms = choices[scenario]; mask_folder.mkdir(parents=True, exist_ok=False)
    native = mask_folder / 'native-exclusions.tsv'; native.write_text(text)
    intervals, whole = parse_exclusions(native, dictionary)
    write_bed(mask_folder / 'query-intervals.bed', query_rows(dictionary, intervals, whole))
    write_bed(mask_folder / 'primary-template-intervals.bed', [[chrom,a,b] for chrom in dictionary
                                                            if chrom in chroms for a,b in intervals[chrom]])
    summary = {'scope': 'synthetic_query_exclusions', 'chromosomes': chroms, 'reference_dictionary': dictionary,
               'inputs': {'reference': fingerprint(reference), 'reference_index': fingerprint(Path(str(reference) + '.fai'))},
               'outputs': {name:fingerprint(mask_folder/name) for name in
                           ('native-exclusions.tsv','query-intervals.bed','primary-template-intervals.bed')}}
    save(mask_folder / 'summary.json', summary)
    config = {'scope':'synthetic','sample':'SYNTH_SV','chromosomes':chroms,
              'reference_sha256':fingerprint(reference)['sha256'],
              'reference_index_sha256':fingerprint(Path(str(reference)+'.fai'))['sha256'],
              'exclusion_summary_sha256':fingerprint(mask_folder/'summary.json')['sha256'],
              'delly_sha256':fingerprint(delly)['sha256'],'samtools_sha256':fingerprint(samtools)['sha256'],
              'threads':2,'native_options':list(NATIVE_OPTIONS),
              'read_policy':'native_delly_sr_unmodified_alignment_reject_unknown_mapped_MQ',
              'ploidy_policy':'native_diploid_discovery_sex_genotypes_uninterpreted'}
    return config


def probe(experiment, output, delly, samtools, alignment_format='bam'):
    output.mkdir(parents=True, exist_ok=False); start = time.monotonic()
    original = json.loads((experiment/'summary.json').read_text())
    alignment = experiment/'aligned.bam' if alignment_format=='bam' else experiment/'cram-check/aligned.cram'
    expected_hash = original['output_sha256']['aligned.bam'] if alignment_format=='bam' else json.loads((experiment/'cram-check/summary.json').read_text())['outputs']['aligned.cram']
    if fingerprint(alignment)['sha256'] != expected_hash: raise ValueError('original synthetic alignment changed')
    scenarios = ['control','deletion_window','exclude_partner','exclude_changed_chromosome'] if alignment_format=='bam' else ['control']
    records = []
    expected_types = {'control':{'DEL':1,'DUP':1,'INV':2,'BND':2},
                      'deletion_window':{'DUP':1,'INV':2,'BND':2},
                      'exclude_partner':{'DEL':1,'DUP':1,'INV':2}, 'exclude_changed_chromosome':{}}
    for scenario in scenarios:
        folder = output/scenario; folder.mkdir()
        config=fixture_config(experiment, folder/'mask', scenario, delly, samtools);save(folder/'config.json',config)
        result=run_analysis(alignment,experiment/'fixture/reference.fa',folder/'mask',config,delly,samtools,folder/'analysis')
        observed=result['native_discovery']['types'];passed=observed==expected_types[scenario]
        records.append({'scenario':scenario,'expected_types':expected_types[scenario],'observed_types':observed,
                        'passed':passed,'summary':fingerprint(folder/'analysis/summary.json'),
                        'elapsed_seconds':result['elapsed_seconds'],'analysis_status':result['analysis_status'],
                        'alignment_census':result['alignment_census']})
        if not passed:
            save(output/'unexpected-result.json',records);raise ValueError('whole-caller exclusion result differs from prespecified expectation')
        if scenario=='control':
            rows=[json.loads(line) for line in (folder/'analysis/calls-review.jsonl').read_text().splitlines()]
            for before,after in zip(original['native_calls'],rows):
                for key in ('chrom','pos_1based','end_chrom','end_position','type','connection_type','filter','qual','evidence'):
                    if before[key]!=after[key]:raise ValueError('control parsed call field differs: '+key)
                for key,value in before['genotype'].items():
                    if value!=after['genotype'].get(key):raise ValueError('control genotype evidence differs')
    result={'completed_at':stamp(),'elapsed_seconds':time.monotonic()-start,'alignment_format':alignment_format,
            'probe_source_sha256':fingerprint(Path(__file__))['sha256'],
            'runner_source_sha256':fingerprint(Path(__file__).with_name('run_structural_variants.py'))['sha256'],
            'original_summary_sha256':fingerprint(experiment/'summary.json')['sha256'],
            'scenarios':records,'scope':'Known synthetic events and original aligned reads through the guarded runner. Exclusion effects are assessed for wide windows/whole contigs, not exact breakpoint edges or clinical sensitivity.'}
    save(output/'summary.json',result);return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('experiment','output','delly','samtools'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--alignment-format',choices=('bam','cram'),default='bam')
    a=p.parse_args();os.nice(10);r=probe(a.experiment,a.output,a.delly,a.samtools,a.alignment_format)
    print(json.dumps({'completed_at':r['completed_at'],'elapsed_seconds':r['elapsed_seconds'],
                      'scenarios':[{k:row[k] for k in ('scenario','observed_types','passed','analysis_status')} for row in r['scenarios']]},indent=2))


if __name__=='__main__':main()
