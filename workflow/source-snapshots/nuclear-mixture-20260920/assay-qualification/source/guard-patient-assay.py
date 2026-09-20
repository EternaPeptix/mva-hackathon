#!/usr/bin/env python3
"""Reviewed patient or fictional assay in one owned bounded process group."""
import argparse,hashlib,importlib.util,os,platform,sys
from pathlib import Path
import patient_common as c
G=Path(__file__).with_name('guard-fictional-pilot.py')
if hashlib.sha256(G.read_bytes()).hexdigest()!='61eb3dffa51c23f30e203036d43d9fbdee76d314b5e869daf9fefbbd7aaacb18':raise ValueError('qualified guard helper differs')
spec=importlib.util.spec_from_file_location('assay_owned_guard',G);g=importlib.util.module_from_spec(spec);spec.loader.exec_module(g)
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',type=Path,required=True);p.add_argument('--expected-config-sha256',required=True);a=p.parse_args();cfg,pin,sources=c.config(a.config,a.expected_config_sha256)
    c.require(platform.node().split('.')[0]=='512S2','Mac Studio512S2 required');out=Path(cfg['output']).resolve();launch=Path(cfg['launch']).resolve();c.require(not out.is_relative_to(launch) and not launch.is_relative_to(out),'output/launch overlap')
    if cfg['fixture_only']:c.require(launch.is_relative_to(Path(cfg['fixture_root']).resolve()),'fixture launch escaped');before=None
    else:before=g.caller_identity()
    os.nice(10);state=None
    try:
        command=[sys.executable,'-B',str(Path(__file__).with_name('run-patient-assay.py')),'--config',str(a.config.resolve()),'--expected-config-sha256',a.expected_config_sha256,'--guard-pid',str(os.getpid())]
        state=g.supervise(command,out,launch,c.LIMITS,c.ROOT)
        terminal=c.parse((out/'status.json').read_text());c.require(terminal['status']=='complete','assay/checker terminal incomplete')
        c.require(c.fingerprint(out/'producer/status.json')['sha256']==terminal['producer_status_sha256'] and c.fingerprint(out/'verification/verification.json')['sha256']==terminal['verification_sha256'],'terminal links changed')
        after=None if cfg['fixture_only'] else g.caller_identity()
        if before is not None:g.exact_identity(before,after)
        state.update(config=pin,sources=sources,fixture_only=cfg['fixture_only'],result_status=c.fingerprint(out/'status.json'),original_caller_before=before,original_caller_after=after)
    except BaseException as error:
        if launch.exists():state=g.read_json(launch/'status.json');state.update(status='failed',error=str(error),original_caller_before=before)
        raise
    finally:
        if state is not None:g.save(launch/'status.json',state)
if __name__=='__main__':main()
