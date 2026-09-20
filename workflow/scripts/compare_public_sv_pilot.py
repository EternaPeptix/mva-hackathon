#!/usr/bin/env python3
"""Compare the pinned public pilot with full-context GIAB SV truth and MAFFT.

Keep raw calls and all derivations. Pre-resolve symbolic DEL/DUP using the pinned
Truvari implementation so full-context refinement does not silently omit them.
"""
import argparse
from collections import Counter
import importlib.metadata
import json
from pathlib import Path
import time

import pysam
import truvari

from build_structural_exclusions import fingerprint, save, stamp
from prepare_public_sv_benchmark import REFERENCE
from probe_sv_comparator import compare, environment, run_command, write_vcf


def resolve_comparison(source, output, reference):
    params = truvari.VariantParams(dup_to_ins=True, max_resolve=50000, write_resolved=True)
    counts = Counter(); rows = []
    with pysam.FastaFile(str(reference)) as ref, truvari.VariantFile(str(source), params=params) as reader:
        with truvari.VariantFile(str(output), 'wz', header=reader.header, params=params) as writer:
            for entry in reader:
                was_resolved=entry.is_resolved(); original_type=str(entry.info.get('SVTYPE','absent'))
                size=entry.var_size(); counts['input_records']+=1
                if not was_resolved: entry.resolve(ref)
                resolved=entry.is_resolved()
                if not resolved and 0 <= size <= 50000:
                    raise ValueError('in-scope symbolic call cannot be sequence resolved: ' + str(entry.id))
                if resolved:
                    if ref.fetch(entry.chrom, entry.start, entry.start+len(entry.get_ref())).upper()!=entry.get_ref().upper():
                        raise ValueError('resolved REF disagrees with reference')
                    counts['resolved_records']+=1
                else:counts['unresolved_above_size_ceiling']+=1
                if resolved and not was_resolved:counts['resolved_'+original_type]+=1
                rows.append({'id':entry.id,'original_type':original_type,'native_size':size,
                             'was_sequence_resolved':was_resolved,'is_sequence_resolved':resolved,
                             'resolved_sequence_length_difference':len(entry.get_alt())-len(entry.get_ref()) if resolved else None})
                writer.write(entry)
    pysam.tabix_index(str(output),preset='vcf')
    return {'counts':dict(counts),'records':rows,'input':fingerprint(source),'output':fingerprint(output),
            'policy':'Truvari 5.5.0 reference resolution, max 50000 bp; DUP assumes perfect reference copy at start. Native INFO fields retained.'}


def conversion_controls(root, output, config, previous):
    output.mkdir(exist_ok=False)
    ref=previous/'control.fa';bed=previous/'regions.bed'
    with pysam.FastaFile(str(ref)) as f:seq=f.fetch('chr1')
    base=output/'base.vcf.gz';source=output/'symbolic.vcf.gz';resolved=output/'resolved.vcf.gz'
    write_vcf(base,seq,[('COPY_AS_INS',8000,8000,seq[7999:8120])])
    with pysam.VariantFile(str(base)) as reader, pysam.VariantFile(str(source),'wz',header=reader.header) as writer:
        entry=writer.new_record(contig='chr1',start=7999,stop=8120,id='DUP_CONTROL',alleles=(seq[7999],'<DUP>'))
        entry.filter.add('PASS');entry.info['SVTYPE']='DUP';entry.info['SVLEN']=120
        entry.samples['CONTROL']['GT']=(1,0);entry.samples['CONTROL'].phased=True;writer.write(entry)
    pysam.tabix_index(str(source),preset='vcf')
    resolved_receipt=resolve_comparison(source,resolved,ref)
    got=compare(root,base,resolved,ref,bed,output/'duplication',config)
    passed=all(got['bench'][k]==v for k,v in {'TP-base':1,'TP-comp':1,'FP':0,'FN':0}.items())
    if not passed or got['refined']['FP'] or got['refined']['FN']:raise ValueError('DUP conversion native control failed')
    unknown=output/'unresolved.vcf.gz'
    with pysam.VariantFile(str(source)) as reader, pysam.VariantFile(str(unknown),'wz',header=reader.header) as writer:
        for entry in reader:entry.alts=('<INS>',);entry.info['SVTYPE']='INS';writer.write(entry)
    pysam.tabix_index(str(unknown),preset='vcf');rejected=False
    try:resolve_comparison(unknown,output/'must-reject.vcf.gz',ref)
    except ValueError as error:rejected='cannot be sequence resolved' in str(error)
    if not rejected:raise ValueError('unresolved symbolic insertion was not rejected')
    result={'perfect_duplication_conversion_passed':True,'unresolved_in_scope_insertion_rejected':True,
            'conversion':resolved_receipt,'native_comparison':got}
    save(output/'summary.json',result);return result


def run(root, pilot, output, config_path):
    if output.resolve().is_relative_to(Path(__file__).resolve().parents[2]):raise ValueError('raw comparison stays outside repository')
    output.mkdir(parents=True,exist_ok=False);start=time.monotonic();state={'status':'running','started_at':stamp()}
    save(output/'status.json',state)
    try:
        config=json.loads(config_path.read_text())
        if importlib.metadata.version('truvari')!=config['truvari_version']:raise ValueError('Truvari version mismatch')
        previous=root/'results/sv-comparator-controls-20260919-v2';controls=json.loads((previous/'summary.json').read_text())
        if not controls['all_passed'] or controls['configuration']['sha256']!=fingerprint(config_path)['sha256']:
            raise ValueError('matching comparator controls must pass before public comparison')
        if json.loads((pilot/'status.json').read_text())['status']!='complete':raise ValueError('public pilot incomplete')
        ref=root/'ref/giab-sv-20260919/inputs'/REFERENCE
        extra=conversion_controls(root,output/'conversion-controls',config,previous)
        baseline=output/'truth-chr20.no-star.full-context.vcf.gz'
        prep=run_command([root/'tools/bio/bin/bcftools','view','-e','ALT="*"','-Oz','-o',baseline,pilot/'truth-chr20.full-context.vcf.gz'],
                         output/'filter-stars.log',environment(root))
        pysam.tabix_index(str(baseline),preset='vcf')
        comp=output/'comp.sequence-resolved.full-context.vcf.gz'
        converted=resolve_comparison(pilot/'del-ins-dup.vcf.gz',comp,ref);save(output/'conversion.json',converted)
        counts=Counter()
        with pysam.VariantFile(str(baseline)) as reader:
            for entry in reader:
                counts['records_retained']+=1
                if '*' in (entry.alts or ()):raise ValueError('star allele remains after publisher filter')
        state['stage']='native full-context public comparison';save(output/'status.json',state)
        result=compare(root,baseline,comp,ref,pilot/'official-chr20-benchmark.bed',output/'benchmark',config)
        receipt={'scope':'public_HG002_chr20_prealigned_caller_component','completed_at':stamp(),
                 'elapsed_seconds':time.monotonic()-start,'comparison':result,'conversion_counts':converted['counts'],
                 'truth_context_counts':dict(counts),'prior_controls':fingerprint(previous/'summary.json'),
                 'additional_conversion_controls':extra,'preparation_command':prep,
                 'inputs':{name:fingerprint(path) for name,path in {'baseline':baseline,'comp':comp,'reference':ref,
                           'benchmark_bed':pilot/'official-chr20-benchmark.bed','pilot_summary':pilot/'summary.json'}.items()},
                 'code':{name:fingerprint(Path(__file__).with_name(name)) for name in ('compare_public_sv_pilot.py','probe_sv_comparator.py')},
                 'config':fingerprint(config_path),
                 'interpretation_status':'algorithmic comparison only; class/size counts and discordance review required',
                 'limitations':['Single chromosome and prealigned real reads; does not validate patient raw-read alignment.',
                                'Truth is v5.0q draft and includes complex alleles; apparent FP/FN can include representation or benchmark errors.',
                                'Only DEL and insertion-like INS/DUP, 50-50000 bp baseline and >=30 bp comparison; balanced events outside scope.',
                                'Genotype agreement is a separate result and is not a clinical probability.']}
        save(output/'summary.json',receipt);state.update(status='complete',completed_at=receipt['completed_at']);save(output/'status.json',state)
        return receipt
    except BaseException as error:
        state.update(status='failed',error_type=type(error).__name__,error=str(error),failed_at=stamp());save(output/'status.json',state);raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('root','pilot','output','config'):p.add_argument('--'+name,required=True,type=Path)
    a=p.parse_args();r=run(a.root,a.pilot,a.output,a.config)
    print(json.dumps({'elapsed_seconds':r['elapsed_seconds'],'bench':r['comparison']['bench'],'refined':r['comparison']['refined']}))
