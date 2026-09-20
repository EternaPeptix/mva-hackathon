#!/usr/bin/env python3
"""Run synthetic haplotypes through paired FASTQ alignment, duplicate marking and DELLY.

The experiment reports missing events and extra calls rather than redefining truth.
Patient structural analysis is a separate workflow and is not run by this script.
"""
import argparse
from collections import Counter
import datetime
import json
import os
from pathlib import Path
import subprocess
import time

import pysam

from simulate_sv_reads import simulate, sha


def stamp(): return datetime.datetime.now(datetime.timezone.utc).isoformat()


def save(path, value): path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def read_calls(path, sample):
    calls = []
    with pysam.VariantFile(str(path)) as source:
        if list(source.header.samples) != [sample]: raise ValueError('synthetic caller sample identity mismatch')
        for record in source:
            svtype = record.info.get('SVTYPE')
            end_chrom = record.info.get('CHR2', record.contig)
            end = record.info.get('POS2') if svtype == 'BND' else record.stop
            calls.append({'id': record.id, 'chrom': record.contig, 'pos_1based': record.pos,
                'end_chrom': end_chrom, 'end_position': end, 'type': svtype,
                'connection_type': record.info.get('CT'), 'length': record.info.get('SVLEN'),
                'precise': record.info.get('PRECISE', False), 'filter': list(record.filter),
                'qual': record.qual,
                'evidence': {key: record.info.get(key) for key in ('PE', 'SR', 'MAPQ', 'SRMAPQ', 'CIPOS', 'CIEND')},
                'genotype': {key: record.samples[sample].get(key) for key in ('GT', 'GQ', 'GL', 'DR', 'DV', 'RR', 'RV', 'RDCN')}})
    return calls


def compare_truth(truth, calls, tolerance=150):
    comparisons = []; matched_ids = set()
    for event in truth['events_on_one_haplotype']:
        matches = []; normalized = []
        for call in calls:
            if call['type'] != event['type']: continue
            chrom, pos = call['chrom'], call['pos_1based']
            partner, partner_pos = call['end_chrom'], call['end_position']
            orientation = call['connection_type']
            if call['type'] == 'BND' and chrom == event['partner_chrom'] and partner == event['chrom']:
                chrom, partner = partner, chrom; pos, partner_pos = partner_pos, pos
                ends = orientation.split('to') if orientation else []
                if len(ends) != 2 or any(end not in ('3','5') for end in ends):
                    raise ValueError('unrecognized breakend orientation')
                orientation = ends[1] + 'to' + ends[0]
            if chrom != event['chrom'] or pos is None or abs(pos - event['left_cut_0based']) > tolerance: continue
            end = event.get('partner_cut_0based', event.get('right_cut_0based'))
            if partner != event.get('partner_chrom', event['chrom']): continue
            if call['type'] != 'INS' and (partner_pos is None or abs(partner_pos - end) > tolerance): continue
            if call['type'] == 'INS' and (call['length'] is None or abs(abs(call['length']) - event['length_bp']) > 10): continue
            matches.append(call); matched_ids.add(call['id'])
            normalized.append({'id':call['id'],'chrom':chrom,'position':pos,'end_chrom':partner,
                               'end_position':partner_pos,'connection_type':orientation})
        required_ct = set(event.get('expected_connection_types', []))
        observed_ct = {row['connection_type'] for row in normalized}
        comparisons.append({'event_id': event['id'], 'type': event['type'], 'truth': event,
            'site_detected_within_tolerance': bool(matches),
            'required_connection_types_detected': required_ct <= observed_ct,
            'matched_ids': [row['id'] for row in matches],
            'matched_endpoint_views': normalized,
            'matched_genotypes': [row['genotype']['GT'] for row in matches],
            'all_matched_genotypes_heterozygous': bool(matches) and all(tuple(row['genotype'].get('GT') or ()) in ((0,1),(1,0)) for row in matches),
            'all_matched_PASS': bool(matches) and all(row['filter'] == ['PASS'] for row in matches)})
    return {'coordinate_tolerance_bp': tolerance, 'insertion_length_tolerance_bp': 10,
        'coordinate_note': 'Native VCF break positions compared with known zero-based interbase cuts within the declared tolerance; no exact-normalization claim.',
        'events': comparisons, 'unmatched_calls': [row for row in calls if row['id'] not in matched_ids],
        'control_chromosome_calls': [row for row in calls if truth['unaltered_control_chromosome'] in (row['chrom'], row['end_chrom'])],
        'scope': 'Five simple engineered events, required inversion/translocation connections and one unchanged chromosome; not a general sensitivity/specificity estimate.'}


def run(output, bwa, samtools, delly, coverage=30, seed=20260919):
    output.mkdir(parents=True, exist_ok=False); started = time.monotonic()
    state = {'status': 'running', 'pid': os.getpid(), 'started_at': stamp(), 'stage': 'generate independent synthetic paired reads'}
    save(output/'status.json', state); commands = []
    tools = {name: {'path': str(path), 'sha256': sha(path)} for name, path in [('bwa',bwa),('samtools',samtools),('delly',delly)]}
    def command(label, args):
        state.update(stage=label); save(output/'status.json', state)
        log = output / ('step-%02d.log' % len(commands)); tick = time.monotonic()
        with log.open('w') as stream: result = subprocess.run(list(map(str,args)), stdout=stream, stderr=subprocess.STDOUT)
        commands.append({'stage': label, 'command': list(map(str,args)), 'exit_code': result.returncode,
                         'elapsed_seconds': time.monotonic()-tick, 'log_sha256': sha(log)})
        save(output/'commands.json', commands)
        if result.returncode: raise RuntimeError('native stage failed: '+label)
    try:
        truth = simulate(output/'fixture', coverage, seed); reference = output/'fixture/reference.fa'
        command('reference index for BWA', [bwa,'index',reference])
        command('reference FASTA index', [samtools,'faidx',reference])
        state.update(stage='paired read alignment and name sort'); save(output/'status.json',state)
        name_bam = output/'names.bam'
        align_args = [str(bwa),'mem','-Y','-K','100000000','-t','2','-R',
                      r'@RG\tID:SYNTH.1\tSM:SYNTH_SV\tLB:SYNTH\tPL:ILLUMINA',str(reference),
                      str(output/'fixture/SYNTH_SV_R1.fastq.gz'),str(output/'fixture/SYNTH_SV_R2.fastq.gz')]
        sort_args = [str(samtools),'sort','-n','-@','1','-m','128M','-o',str(name_bam),'-']
        tick = time.monotonic()
        with (output/'bwa.log').open('w') as bwa_log, (output/'name-sort.log').open('w') as sort_log:
            process = subprocess.Popen(align_args, stdout=subprocess.PIPE, stderr=bwa_log)
            try:
                sorted_result = subprocess.run(sort_args, stdin=process.stdout, stdout=sort_log, stderr=subprocess.STDOUT)
                process.stdout.close()
                if sorted_result.returncode and process.poll() is None: process.terminate()
                code = process.wait()
            finally:
                if process.poll() is None: process.terminate(); process.wait()
        commands.append({'stage': 'paired read alignment and name sort','pipeline': [align_args,sort_args],
            'exit_codes': [code,sorted_result.returncode], 'elapsed_seconds': time.monotonic()-tick,
            'log_sha256': {'bwa': sha(output/'bwa.log'),'name_sort': sha(output/'name-sort.log')}})
        save(output/'commands.json',commands)
        if code or sorted_result.returncode: raise RuntimeError('synthetic alignment pipeline failed')
        command('mate metadata', [samtools,'fixmate','-m','-@','1',name_bam,output/'mates.bam'])
        command('coordinate sort', [samtools,'sort','-@','1','-m','128M','-o',output/'coordinates.bam',output/'mates.bam'])
        alignment=output/'aligned.bam'
        command('duplicate marking', [samtools,'markdup','-@','1','-f',output/'markdup.stats',output/'coordinates.bam',alignment])
        command('alignment index', [samtools,'index',alignment])
        command('alignment format check', [samtools,'quickcheck','-v',alignment])
        counts=Counter()
        with pysam.AlignmentFile(str(alignment)) as reader:
            if reader.header.to_dict()['RG'][0]['SM'] != 'SYNTH_SV': raise ValueError('alignment sample identity changed')
            for read in reader:
                counts['alignment_records']+=1
                if not read.flag & 0x900: counts['primary_read_records']+=1
                if read.is_duplicate: counts['duplicate_flag_records']+=1
                if read.is_supplementary: counts['supplementary_records']+=1
                if read.is_unmapped: counts['unmapped_records']+=1
                if read.mapping_quality==255: counts['unknown_MQ_255']+=1
        if counts['primary_read_records'] != truth['read_records']: raise ValueError('FASTQ to primary-read accounting mismatch')
        command('native structural discovery', [delly,'sr','-h','2','-t','ALL','-q','20','-u','20','-r','20',
            '-z','3','-s','9','-c','25','-m','25','-n','40','-p','20','-a','250','-g',reference,
            '-o',output/'calls.bcf',alignment])
        if not (output/'calls.bcf.csi').exists(): raise ValueError('native variant index missing')
        calls=read_calls(output/'calls.bcf','SYNTH_SV');comparison=compare_truth(truth,calls)
        for identity in tools.values():
            if sha(Path(identity['path']))!=identity['sha256']: raise ValueError('native executable changed during experiment')
        result={'started_at':state['started_at'],'completed_at':stamp(),'elapsed_seconds':time.monotonic()-started,
            'tools':tools,'probe_source_sha256':sha(Path(__file__)),'simulator_source_sha256':sha(Path(__file__).with_name('simulate_sv_reads.py')),
            'truth_sha256':sha(output/'fixture/truth.json'),'commands':commands,'alignment_counts':dict(counts),
            'native_calls':calls,'comparison':comparison,
            'output_sha256':{p.name:sha(p) for p in (alignment,Path(str(alignment)+'.bai'),output/'calls.bcf',output/'calls.bcf.csi',output/'markdup.stats')},
            'scope':'Synthetic diploid haplotypes to paired FASTQs, real BWA alignment/duplicate marking and native SV discovery. Misses, unmatched calls and genotype discrepancies are retained. No patient analysis or clinical validation.'}
        save(output/'summary.json',result);state.update(status='complete',stage='complete',completed_at=result['completed_at']);save(output/'status.json',state)
        return result
    except BaseException as error:
        state.update(status='failed',failed_at=stamp(),error_type=type(error).__name__);save(output/'status.json',state);raise


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('output','bwa','samtools','delly'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--coverage',type=float,default=30);p.add_argument('--seed',type=int,default=20260919)
    a=p.parse_args();os.nice(10);r=run(a.output,a.bwa,a.samtools,a.delly,a.coverage,a.seed)
    print(json.dumps({'completed_at':r['completed_at'],'elapsed_seconds':r['elapsed_seconds'],
        'alignment_counts':r['alignment_counts'],'native_call_count':len(r['native_calls']),
        'event_comparisons':r['comparison']['events'],'unmatched_calls':len(r['comparison']['unmatched_calls']),
        'control_chromosome_calls':len(r['comparison']['control_chromosome_calls'])},indent=2))


if __name__=='__main__':main()
