#!/usr/bin/env python3
"""Resume a failed native SV call using an already record-verified CRAM 3.0 copy.

Preserves the failed run. Reuses its complete census only after pinning that
artifact, both alignments, the conversion proof, original runner and configuration.
The native caller/options/mask and decoded alignment evidence remain unchanged.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import time

import pysam
from build_structural_exclusions import fingerprint, save, stamp
from run_structural_variants import validate_config, verify_completion, load_mask, review_calls, PATIENT_SCOPE


def require(value, message):
    if not value:
        raise ValueError(message)


def pinned(path, digest):
    record=fingerprint(path);require(record['sha256']==digest,'identity differs: '+str(path));return record


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('failed-run','conversion','completion-state','reference','exclusions','config','control','output'):
        p.add_argument('--'+n,type=Path,required=True)
    for n in ('census','conversion','config','control','original-runner'):
        p.add_argument('--expected-'+n+'-sha256',required=True)
    a=p.parse_args();os.nice(10);start=time.monotonic()
    out=a.output.resolve();require(not out.is_relative_to(Path(__file__).resolve().parents[2]),'case outputs must remain outside Git')
    out.mkdir(parents=True,exist_ok=False)
    state={'status':'running','pid':os.getpid(),'started_at':stamp(),'stage':'verify preserved census and equivalent copy'}
    save(out/'status.json',state);commands=[]
    try:
        cfgid=pinned(a.config,a.expected_config_sha256);cfg=json.loads(a.config.read_text());validate_config(cfg)
        require(cfg['scope']==PATIENT_SCOPE,'full primary-nuclear scope required')
        sourceid=pinned(Path(__file__).with_name('run_structural_variants.py'),a.expected_original_runner_sha256)
        proofid=pinned(a.conversion/'summary.json',a.expected_conversion_sha256)
        proof=json.loads((a.conversion/'summary.json').read_text());conversion_state=json.loads((a.conversion/'status.json').read_text())
        require(conversion_state['status']=='complete' and conversion_state['summary_sha256']==proofid['sha256'],'format conversion incomplete')
        require(proof['header_equal'] and proof['index_readable'] and proof['cram_version']==[3,0],'format proof incomplete')
        verified=proof['decoded_record_verification']
        require(all(verified.get(k) is True for k in ('full_flags_equal','all_auxiliary_types_values_equal','absolute_record_order_equal','all_compared_fields_equal')),'decoded evidence equality is incomplete')
        original,alignment=Path(proof['source']),Path(proof['cram_path'])
        origid=pinned(original,proof['source_sha256']);alignid=pinned(alignment,proof['cram_sha256'])
        indexid=pinned(Path(str(alignment)+'.crai'),proof['index_sha256'])
        completion=verify_completion(cfg,original,origid,a.completion_state)
        censusid=pinned(a.failed_run/'alignment-census.json',a.expected_census_sha256)
        census=json.loads((a.failed_run/'alignment-census.json').read_text());counts=census['counts']
        require(counts['alignment_records']==verified['records_compared'],'census/record comparison totals differ')
        require(counts['alignment_records']==counts['mapped_records']+counts['unmapped_records']==sum(census['records_by_read_group'].values()),'read census does not reconcile')
        require(counts['mapped_records']==sum(census['mapped_MQ_histogram'].values()),'MQ census differs')
        require(not counts.get('mapped_unknown_MQ_255') and not census['mapped_MQ_histogram'].get('255'),'unknown mapped MQ requires review')
        for k in ('duplicate_records','secondary_records','supplementary_records','unmapped_records'):
            require(counts.get(k,0)==verified[k],'census differs from decoded-record proof: '+k)
        oldstate=json.loads((a.failed_run/'status.json').read_text())
        oldcmds=json.loads((a.failed_run/'commands.json').read_text())
        require(oldstate['status']=='failed' and oldstate['stage']=='native structural discovery','not a native-discovery failure')
        require(len(oldcmds)==2 and oldcmds[0]['exit_code']==0 and oldcmds[1]['exit_code']!=0,'unexpected prior stage outcomes')
        oldcmd=oldcmds[1]['command'];delly=Path(oldcmd[0]);samtools=Path(oldcmds[0]['command'][0])
        expected=[str(delly),'sr','-h',str(cfg['threads']),*cfg['native_options'],'-x',str(a.exclusions/'native-exclusions.tsv'),
                  '-g',str(a.reference),'-o',str(a.failed_run/'calls.bcf'),str(original)]
        require(oldcmd==expected,'failed native command differs from frozen options/inputs')
        inputs={'original_alignment':origid,'alignment':alignid,'alignment_index':indexid,
            'reference':pinned(a.reference,cfg['reference_sha256']),
            'reference_index':pinned(Path(str(a.reference)+'.fai'),cfg['reference_index_sha256']),
            'delly':pinned(delly,cfg['delly_sha256']),'samtools':pinned(samtools,cfg['samtools_sha256'])}
        require(proof['reference_sha256']==cfg['reference_sha256'] and proof['samtools_sha256']==cfg['samtools_sha256'],'conversion used different reference/tool')
        controlid=pinned(a.control,a.expected_control_sha256);control=json.loads(a.control.read_text())
        require(control['cram_version']==[3,0] and control['all_decoded_fields_and_typed_tags_equal'] and control['all_compared_call_fields_equal']
                and control['native_call_count']==6 and control['delly']['sha256']==cfg['delly_sha256'],'native format control incomplete')
        with pysam.AlignmentFile(str(original),reference_filename=str(a.reference)) as old,pysam.AlignmentFile(str(alignment),reference_filename=str(a.reference)) as new,pysam.FastaFile(str(a.reference)) as fasta:
            dictionary=dict(zip(fasta.references,fasta.lengths))
            require(list(dictionary.items())==list(zip(new.references,new.lengths)),'reference dictionary differs')
            for section in ('HD','SQ','RG'):
                require(old.header.to_dict()[section]==new.header.to_dict()[section],'alignment header section differs: '+section)
            groups=new.header.to_dict()['RG']
            require(new.has_index() and new.header.to_dict()['HD']['SO']=='coordinate','copy is not indexed/sorted')
            require(len(groups)==4 and all(g['SM']==cfg['sample']==proof['sample'] for g in groups),'copy sample/read groups differ')
            require(set(census['records_by_read_group'])=={g['ID'] for g in groups},'census groups differ')
            maskid,intervals,whole=load_mask(a.exclusions,dictionary,cfg,inputs['reference'],inputs['reference_index'])
        provenance={'created_at':stamp(),'scope':'Same-record CRAM compatibility recovery; no realignment, filtering or new census.',
            'config':cfgid,'original_runner':sourceid,'recovery_script':fingerprint(Path(__file__)),
            'conversion':proofid,'alignment_completion':completion,'preserved_census':censusid,'census':census,
            'native_control':controlid,'inputs':inputs,'exclusion_summary':maskid,'failed_run':str(a.failed_run)}
        save(out/'provenance.json',provenance)
        def run(label,cmd):
            state.update(stage=label);save(out/'status.json',state);log=out/('step-%02d.log'%len(commands));t=time.monotonic()
            with log.open('x') as stream:r=subprocess.run(list(map(str,cmd)),stdout=stream,stderr=subprocess.STDOUT)
            commands.append({'stage':label,'command':list(map(str,cmd)),'exit_code':r.returncode,'elapsed_seconds':time.monotonic()-t,'log':fingerprint(log)})
            save(out/'commands.json',commands);require(r.returncode==0,'native stage failed: '+label)
        run('alignment format check',[samtools,'quickcheck','-v',alignment])
        cmd=list(oldcmd);cmd[-1]=str(alignment);cmd[cmd.index('-o')+1]=str(out/'calls.bcf')
        save(out/'native-command.json',cmd)
        run('native structural discovery on record-verified CRAM 3.0',cmd)
        require((out/'calls.bcf.csi').is_file(),'native index missing')
        reviewed=review_calls(out/'calls.bcf',out/'calls-review.jsonl',cfg['sample'],dictionary,cfg['chromosomes'],intervals,whole)
        for rec in inputs.values():
            s=Path(rec['path']).stat();require((s.st_size,s.st_mtime_ns)==(rec['bytes'],rec['mtime_ns']),'input changed during resumed discovery')
        load_mask(a.exclusions,dictionary,cfg,inputs['reference'],inputs['reference_index'])
        result={'completed_at':stamp(),'elapsed_seconds':time.monotonic()-start,'scope':cfg['scope'],'config':cfg,
            'provenance':fingerprint(out/'provenance.json'),'commands':commands,'native_discovery':reviewed,
            'analysis_status':'complete_with_calls' if reviewed['counts'].get('records') else 'complete_without_calls',
            'outputs':{n:fingerprint(out/n) for n in ('calls.bcf','calls.bcf.csi','calls-review.jsonl')},
            'limits':['Reuses completed metadata census through direct equality of every decoded record; no independent biological replication.',
                      'Native diploid discovery remains exploratory; no sex-genotype, phase, mosaicism, clinical or EPCR conclusion.',
                      'Primary nuclear query mask and previous public sensitivity limitations remain unchanged.']}
        save(out/'summary.json',result);state.update(status='complete',stage='complete',completed_at=result['completed_at'],summary_sha256=fingerprint(out/'summary.json')['sha256']);save(out/'status.json',state)
    except BaseException as e:
        state.update(status='failed',failed_at=stamp(),error_type=type(e).__name__,error=str(e));save(out/'status.json',state);raise


if __name__=='__main__':main()
