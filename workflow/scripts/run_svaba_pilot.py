#!/usr/bin/env python3
"""Run the pinned native SvABA build on two fixed nonpatient pilot inputs."""
import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

import pysam

from build_structural_exclusions import fingerprint, save, stamp
from prepare_public_sv_benchmark import ALIGNMENT
from convert_svaba_pilot import reference_status

SYNTH_BAM_SHA = 'cb3c815adbe16131d3e76cae7b25b947c44d1c96df89176e31ef43efde6d43f9'
SYNTH_REF_SHA = '21ae55bcdc0a516f852b6fa99f5bcd63d7cb0c4479cf76382295ea41c87dcd60'
PUBLIC_BAM_SHA = '34ac157739e1feeb590f6eb7e11046ccc2aa3277fd55a3ce0942e774d931ed81'
EXE_SHA = '1e408ade23375b1b1e1281ec1d26b2527f32f126f13924afdaee0c9c12f1a3de'
OPTIONS = ['--germline', '-p', '2', '-C', '100', '--num-assembly-rounds', '3', '--lod', '8']


def inspect_outputs(output, reference, alignment):
    summaries = {}
    with pysam.FastaFile(str(reference)) as ref:
        for kind in ('indel', 'sv'):
            unfiltered = output / ('svaba.svaba.unfiltered.' + kind + '.vcf')
            filtered = output / ('svaba.svaba.' + kind + '.vcf')
            all_records = {}; pass_ids = set(); counts = Counter()
            with pysam.VariantFile(str(unfiltered)) as reader:
                if list(reader.header.samples) != [str(alignment)]:
                    raise ValueError('native sample column must be the exact input BAM path')
                for row in reader:
                    if row.id in all_records: raise ValueError('duplicate native record ID')
                    if len(row.alts or ()) != 1: raise ValueError('unexpected multiallelic native record')
                    verification=reference_status(row,ref);counts['reference:'+verification]+=1
                    counts['records'] += 1
                    counts['filter:' + ';'.join(row.filter)] += 1
                    counts['GT:' + str(row.samples[str(alignment)].get('GT'))] += 1
                    all_records[row.id] = str(row)
                    if list(row.filter) == ['PASS']: pass_ids.add(row.id)
            filtered_ids = set()
            with pysam.VariantFile(str(filtered)) as reader:
                for row in reader:
                    if row.id in filtered_ids or all_records.get(row.id) != str(row):
                        raise ValueError('filtered records differ from unfiltered source')
                    filtered_ids.add(row.id)
            if filtered_ids != pass_ids: raise ValueError('filtered VCF is not the exact PASS subset')
            summaries[kind] = {'counts': dict(counts), 'filtered_records': len(filtered_ids),
                               'filtered_exact_PASS_subset': True, 'all_interior_REF_match': True,
                               'all_REF_match':not counts['reference:virtual_boundary_no_reference_base']}
    return summaries


def run(root, mode, output, config_path):
    if output.resolve().is_relative_to(Path(__file__).resolve().parents[2]):
        raise ValueError('raw pilot artifacts must stay outside Git')
    output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    state = {'status': 'running', 'pid': os.getpid(), 'started_at': stamp(), 'mode': mode}
    save(output/'status.json', state)
    try:
        if shutil.disk_usage(output).free < 350 * 2**30: raise ValueError('insufficient reserve for active patient workflow')
        config = json.loads(config_path.read_text())
        if config['common_options'] != OPTIONS or config['executable_sha256'] != EXE_SHA:
            raise ValueError('pilot configuration is not the frozen configuration')
        native = root/'tools/svaba-1.2.0/bin/svaba'
        if fingerprint(native)['sha256'] != EXE_SHA: raise ValueError('native executable changed')
        if mode == 'synthetic':
            folder = root/'results/structural-synthetic-20260919-v2'
            alignment = folder/'aligned.bam'; reference = folder/'fixture/reference.fa'; sample = 'SYNTH_SV'
            expected_bam = SYNTH_BAM_SHA; expected_ref = SYNTH_REF_SHA; region_args = []
        elif mode == 'public':
            folder = root/'ref/giab-sv-20260919'
            alignment = folder/'inputs'/ALIGNMENT; reference = folder/'assembly-reference/GRCh38_no_alt.fa'; sample = 'HG002'
            prep = json.loads((reference.parent/'summary.json').read_text())
            if prep['status'] != 'complete' or not prep['all_reference_sequences_identical']:
                raise ValueError('public assembly reference preparation incomplete')
            expected_ref = prep['outputs']['']['sha256']; expected_bam = PUBLIC_BAM_SHA; region_args = ['-k', 'chr20']
            if prep['inputs']['compressed_reference']['sha256'] != '3c8def6d325c5d1e934b2dd530c4d1709f027677a15859467071fddf3fff2026':
                raise ValueError('public reference preparation used a different reference')
        else: raise ValueError('only fixed nonpatient pilots supported')
        inputs = {name: fingerprint(path) for name, path in {'alignment': alignment, 'reference': reference,
                  'alignment_index': Path(str(alignment)+'.bai'), 'native': native, 'config': config_path,
                  'code': Path(__file__)}.items()}
        if inputs['alignment']['sha256'] != expected_bam or inputs['reference']['sha256'] != expected_ref:
            raise ValueError('pilot alignment or reference identity mismatch')
        index = {suffix: fingerprint(Path(str(reference)+suffix)) for suffix in ('.fai','.amb','.ann','.bwt','.pac','.sa')}
        if mode == 'public' and any(v['sha256'] != prep['outputs'][suffix]['sha256'] for suffix, v in index.items()):
            raise ValueError('public reference index changed after preparation')
        with pysam.AlignmentFile(str(alignment)) as reader, pysam.FastaFile(str(reference)) as ref:
            header = reader.header.to_dict()
            if not reader.has_index() or header.get('HD',{}).get('SO') != 'coordinate': raise ValueError('indexed sorted BAM required')
            if {x.get('SM') for x in header.get('RG',[])} != {sample}: raise ValueError('input sample mismatch')
            if list(zip(reader.references,reader.lengths)) != list(zip(ref.references,ref.lengths)):
                raise ValueError('complete BAM/reference dictionary mismatch')
            if mode == 'public':
                if any(x.mapped and x.contig != 'chr20' for x in reader.get_index_statistics()): raise ValueError('public input is not chr20 only')
                if hashlib.sha256(ref.fetch('chr20').upper().encode()).hexdigest() != '2442701218d3b2f7e2d7cacd78774e1628fd48a7a13b58bee53ab8b7cb3d1af4':
                    raise ValueError('public chr20 sequence mismatch')
        save(output/'inputs.json', {'files': inputs, 'reference_indexes': index, 'sample_from_read_groups': sample,
                                    'sample_output_policy': 'Native VCF uses BAM path; this mapping is preserved explicitly'})
        args = [str(native), 'run', '-t', str(alignment), '-G', str(reference), '-a', 'svaba', *OPTIONS, *region_args]
        env = dict(os.environ); env['OMP_NUM_THREADS'] = '2'; env['OPENBLAS_NUM_THREADS'] = '2'
        state.update(stage='native assembly and calling', command=args); save(output/'status.json', state)
        tick = time.monotonic()
        with (output/'native.log').open('w') as log:
            proc = subprocess.run(args, cwd=output, env=env, stdout=log, stderr=subprocess.STDOUT)
        command = {'command': args, 'exit_code': proc.returncode, 'elapsed_seconds': time.monotonic()-tick}
        save(output/'command.json', command)
        if proc.returncode: raise RuntimeError('native SvABA failed; original log retained')
        summaries = inspect_outputs(output, reference, alignment)
        if fingerprint(native)['sha256'] != EXE_SHA: raise ValueError('executable changed during pilot')
        artifacts = {p.name: fingerprint(p) for p in sorted(output.iterdir()) if p.is_file() and p.name != 'status.json'}
        result = {'status': 'complete', 'mode': mode, 'started_at': state['started_at'], 'completed_at': stamp(),
                  'elapsed_seconds': time.monotonic()-started, 'native_command': command, 'inputs': inputs,
                  'configuration': config, 'native_outputs': summaries, 'artifacts': artifacts,
                  'comparison_status': 'not yet compared; paired-breakend conversion and controls required',
                  'scope': 'Nonpatient caller-component validation only; misses and extra calls must be retained'}
        save(output/'summary.json', result); state.update(status='complete',completed_at=result['completed_at']); save(output/'status.json', state)
        return result
    except BaseException as error:
        state.update(status='failed',error_type=type(error).__name__,error=str(error),failed_at=stamp())
        save(output/'status.json', state); raise


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('root','output','config'): p.add_argument('--'+name, required=True, type=Path)
    p.add_argument('--mode', required=True, choices=('synthetic','public'))
    a = p.parse_args(); os.nice(10); r = run(a.root,a.mode,a.output,a.config)
    print(json.dumps({'mode':r['mode'],'elapsed_seconds':r['elapsed_seconds'],'native_outputs':r['native_outputs']}))
