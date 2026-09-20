#!/usr/bin/env python3
"""One guarded CRAM pileup and conditional GATK estimator; immutable inputs."""
import argparse,datetime,json,os,subprocess,sys,time
from pathlib import Path
import patient_common as c
def stamp():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def header_check(cfg):
    import pysam
    inputs=cfg['inputs'];alignment=Path(inputs['alignment']['path'])
    with alignment.open('rb') as f:c.require(f.read(6)==b'CRAM\x03\x00','CRAM3.0 required')
    with pysam.FastaFile(inputs['reference']['path']) as reference,pysam.AlignmentFile(str(alignment),'rc',reference_filename=inputs['reference']['path'],index_filename=inputs['alignment_index']['path'],format_options=[b'decode_md=0']) as reads:
        h=reads.header.to_dict();groups=h.get('RG',[]);ids=[x.get('ID') for x in groups]
        c.require(len(groups)==cfg['expected_read_groups'] and all(isinstance(x,str) and x for x in ids) and len(set(ids))==len(ids),'read group identity/count')
        c.require({x.get('SM') for x in groups}=={cfg['sample']},'one exact nonempty sample required')
        lengths=dict(zip(reference.references,reference.lengths));sq=[(x['SN'],x['LN']) for x in h['SQ']]
        c.require(sq==list(lengths.items()),'CRAM/reference dictionary order/length differs')
        dictionary=[]
        for line in Path(inputs['reference_dictionary']['path']).read_text().splitlines():
            if line.startswith('@SQ\t'):
                row=dict(x.split(':',1) for x in line.split('\t')[1:]);dictionary.append((row['SN'],int(row['LN'])))
        c.require(dictionary==sq,'GATK reference dictionary differs')
        return h
def run(cfg,config_pin,sources,out):
    c.require(not out.exists(),'new output required');out.mkdir(parents=True);started=time.monotonic();state=dict(status='running',started_at=stamp(),fixture_only=cfg['fixture_only'],config=config_pin,sources=sources);commands=[]
    try:
        inputs={k:c.fingerprint(v['path'],v['sha256']) for k,v in cfg['inputs'].items()};runtime={k:c.fingerprint(v['path'],v['sha256']) for k,v in cfg['runtime'].items()}
        header=header_check(cfg);c.save(out/'header.json',header)
        markers=c.panel(inputs['panel']['path'],cfg['expected_panel_markers'],cfg['contigs']);logs=out/'commands';logs.mkdir()
        prefix=[runtime['java']['path'],'-Xmx2g','-XX:ActiveProcessorCount=2','-jar',runtime['gatk']['path']]
        def command(label,args):
            argv=prefix+args;begin=time.monotonic();record=dict(label=label,argv=argv,started_at=stamp())
            try:
                with (logs/(label+'.stdout.log')).open('xb') as stdout,(logs/(label+'.stderr.log')).open('xb') as stderr:result=subprocess.run(argv,stdout=stdout,stderr=stderr,timeout=cfg['timeouts'][label])
                record.update(exit_code=result.returncode,elapsed_seconds=time.monotonic()-begin,completed_at=stamp());commands.append(record);c.save(logs/(label+'.json'),record);c.require(result.returncode==0,'native command failed '+label)
            except subprocess.TimeoutExpired:
                record.update(status='timeout',elapsed_seconds=time.monotonic()-begin);c.save(logs/(label+'.json'),record);raise
        command('pileup',['GetPileupSummaries','-R',inputs['reference']['path'],'-I',inputs['alignment']['path'],'-V',inputs['panel']['path'],'-L',inputs['panel']['path'],'--minimum-mapping-quality','50','--max-depth-per-sample','0','-O',str(out/'pileups.tsv')])
        denominators,rows=c.inspect_table(out/'pileups.tsv',markers,cfg['sample']);observed={(r['contig'],r['position']) for r in rows};c.save(out/'pileup-membership.json',dict(denominators=denominators,emitted_rows=rows,no_emission=[dict(contig=k[0],position=k[1]) for k in markers if k not in observed]))
        if denominators['coverage_retained']:
            command('estimator',['CalculateContamination','-I',str(out/'pileups.tsv'),'-O',str(out/'contamination.tsv'),'--tumor-segmentation',str(out/'minor-allele-fraction-segments.tsv')])
            reason=None
        else:reason='zero_native_coverage_retained_sites'
        for pin in list(inputs.values())+list(runtime.values()):
            now=Path(pin['path']).stat();c.require((now.st_size,now.st_mtime_ns)==(pin['bytes'],pin['mtime_ns']),'input changed during producer')
        state.update(status='complete',completed_at=stamp(),elapsed_seconds=time.monotonic()-started,inputs=inputs,runtime=runtime,denominators=denominators,commands=commands,estimator_executed=bool(denominators['coverage_retained']),estimator_not_run_reason=reason,scope='Raw native model inputs/outputs only; independent checker determines eligibility. No clinical, genotype, purity or cell-fraction claim.')
        state['outputs']={str(p.relative_to(out)):c.fingerprint(p) for p in sorted(out.rglob('*')) if p.is_file()}
        c.save(out/'summary.json',state);c.save(out/'status.json',dict(status='complete',summary_sha256=c.fingerprint(out/'summary.json')['sha256']))
    except BaseException as e:
        c.save(out/'failure.json',dict(status='failed',failed_at=stamp(),error_type=type(e).__name__,error=str(e),commands=commands));raise
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',type=Path,required=True);p.add_argument('--expected-config-sha256',required=True);p.add_argument('--guard-pid',type=int,required=True);a=p.parse_args()
    c.require(a.guard_pid==os.getppid() and os.getpgrp()==os.getpid(),'owned guard parent/group required');cfg,pin,sources=c.config(a.config,a.expected_config_sha256);root=Path(cfg['output']);run(cfg,pin,sources,root/'producer')
    command=[sys.executable,'-B',str(Path(__file__).with_name('audit-patient-output.py')),'--config',str(a.config.resolve()),'--expected-config-sha256',a.expected_config_sha256]
    with (root/'audit.stdout.log').open('xb') as stdout,(root/'audit.stderr.log').open('xb') as stderr:result=subprocess.run(command,stdout=stdout,stderr=stderr,timeout=900)
    c.require(result.returncode==0,'independent output audit failed')
    c.save(root/'status.json',dict(status='complete',completed_at=stamp(),producer_status_sha256=c.fingerprint(root/'producer/status.json')['sha256'],verification_sha256=c.fingerprint(root/'verification/verification.json')['sha256']))
if __name__=='__main__':main()
