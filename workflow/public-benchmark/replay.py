#!/usr/bin/env python3
"""Download and replay the fixed PUBLIC HG002 chr20 benchmark in an exclusive tree.

The default is the full chromosome, not a toy. --controls-only tests installation
and the unchanged twelve RTG semantic controls; it is not a benchmark result.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import datetime
import gzip
import hashlib
import importlib
import json
import os
from pathlib import Path
import platform
import shutil
import shlex
import signal
import subprocess
import sys
import time
import urllib.request
import zipfile

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
GIB = 2**30
LIMITS = dict(seconds=21600, rss_bytes=32*GIB, output_bytes=30*GIB, minimum_free_bytes=40*GIB)


def stamp(): return datetime.datetime.now(datetime.timezone.utc).isoformat()

def save(path, obj):
    p = Path(path)
    tmp = p.with_name(p.name + '.tmp')
    tmp.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')
    os.replace(tmp, p)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(8*1024*1024), b''): h.update(b)
    return h.hexdigest()


def pin(path):
    p = Path(path).resolve()
    return dict(path=str(p), bytes=p.stat().st_size, sha256=digest(p))


def require(ok, message):
    if not ok: raise ValueError(message)


def system_platform():
    key = (platform.system(), platform.machine())
    result = {('Darwin','arm64'):'osx-arm64', ('Linux','x86_64'):'linux-64',
              ('Linux','aarch64'):'linux-aarch64'}.get(key)
    require(result is not None, 'Supported platforms: macOS arm64, Linux x86_64/aarch64 (glibc).')
    return result


def fetch(url, path, sha256=None, size=None):
    """One HTTPS transfer; no silent overwrite, fallback data, or automatic retry."""
    path = Path(path)
    require(url.startswith('https://'), 'HTTPS input required')
    require(not path.exists(), 'Existing download preserved: ' + str(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + '.partial')
    with urllib.request.urlopen(url, timeout=120) as response, partial.open('xb') as dest:
        count = 0
        while True:
            block = response.read(8*1024*1024)
            if not block: break
            count += len(block)
            require(count <= (size if size is not None else 2*GIB), 'Download larger than declared limit')
            dest.write(block)
    require(size is None or count == size, 'Download byte count differs: ' + path.name)
    require(sha256 is None or digest(partial) == sha256, 'Download checksum differs: ' + path.name)
    os.replace(partial, path)
    return dict(url=url, **pin(path))


def command(args, folder, label, env=None):
    args = list(map(str,args)); start = time.monotonic()
    print(stamp(), label, flush=True)
    with (folder/(label+'.stdout')).open('w') as out, (folder/(label+'.stderr')).open('w') as err:
        result = subprocess.run(args, stdout=out, stderr=err, env=env)
    receipt = dict(command=args, exit_code=result.returncode, elapsed_seconds=time.monotonic()-start)
    save(folder/(label+'.command.json'),receipt)
    require(result.returncode == 0, 'Command failed; inspect retained log: ' + label)
    return receipt


def source_check(profile):
    for rel, expected in profile['source_sha256'].items():
        require(digest(REPO/rel) == expected, 'Qualified source differs: ' + rel)


def bootstrap(root, mode, workers):
    os.nice(10)
    profile = json.loads((HERE/'profile.json').read_text())
    source_check(profile)
    subdir = system_platform()
    lock = HERE/(subdir+'.explicit.txt')
    require(digest(lock) == profile['locks'][subdir]['sha256'], 'Runtime lock differs')
    save(root/'profile.json',profile)
    save(root/'source-provenance.json',{'repository_sources':[pin(REPO/x) for x in profile['source_sha256']],
         'replay_source':pin(Path(__file__)), 'launcher':pin(HERE/'run.sh'), 'environment_lock':pin(lock),
         'platform':subdir, 'mode':mode, 'workers':workers})
    logs = root/'logs'; logs.mkdir()
    bootstrap_dir = root/'bootstrap'; bootstrap_dir.mkdir()
    mamba = bootstrap_dir/'micromamba'
    transfer = fetch(profile['micromamba']['url_prefix']+subdir,mamba,profile['micromamba']['sha256'][subdir])
    mamba.chmod(0o755)
    prefix = root/'tools/bio'
    env = dict(os.environ)
    # Avoid activating or modifying any user's environment or shell files.
    env.pop('CONDA_PREFIX',None); env.pop('MAMBA_ROOT_PREFIX',None)
    env.pop('PYTHONPATH',None); env.pop('PYTHONHOME',None)
    command([mamba,'--no-rc','create','--yes','--root-prefix',root/'package-cache',
             '--prefix',prefix,'--file',lock],logs,'install-exact-environment',env)
    save(root/'bootstrap.json',{'status':'complete','bootstrap':transfer,'environment_lock':pin(lock),
         'conda_packages':[pin(p) for p in sorted((prefix/'conda-meta').glob('*.json'))]})
    env.update(PATH=str(prefix/'bin')+os.pathsep+env.get('PATH',''),
               JAVA_HOME=str(prefix/'lib/jvm'), PYTHONNOUSERSITE='1',
               OMP_NUM_THREADS='2', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    # The pinned environment Python executes the same source in the owned process group.
    os.execve(str(prefix/'bin/python'),[str(prefix/'bin/python'),str(Path(__file__)),str(root),
              '--pipeline-worker','--mode',mode,'--workers',str(workers)],env)


def pipeline(root, mode, workers):
    import pysam
    profile = json.loads((root/'profile.json').read_text()); source_check(profile)
    logs = root/'logs'; tools = root/'tools'; downloads = root/'downloads'; downloads.mkdir()
    java = tools/'bio/lib/jvm/bin/java'
    require(java.is_file(), 'Pinned OpenJDK layout not found')
    for name in ('samtools','bcftools'):
        command([tools/'bio/bin'/name,'--version'],logs,name+'-version')
    command([java,'-version'],logs,'java-version')
    require(pysam.__version__ == '0.24.1','pysam version differs')
    rtg = profile['rtg']
    fetch(rtg['url'],downloads/'rtg.zip',rtg['sha256'],rtg['bytes'])
    with zipfile.ZipFile(downloads/'rtg.zip') as archive:
        for info in archive.infolist():
            target = tools/info.filename
            require(target.resolve().is_relative_to(tools.resolve()), 'Unsafe runtime archive member')
        archive.extractall(tools)
    rtgdir = tools/'rtg-tools-3.13'
    config = json.loads((REPO/'workflow/configs/public-small-variant-comparison.json').read_text())
    for name in ('rtg','RTG.jar'):
        require(digest(rtgdir/name) == config['runtime_sha256'][name], 'RTG release bytes differ')
    (rtgdir/'rtg').chmod(0o755)
    # Runtime-only relocation. All scientific comparator arguments remain frozen.
    (rtgdir/'rtg.cfg').write_text('RTG_JAVA='+shlex.quote(str(java))+'\nRTG_MEM="4G"\nRTG_DEFAULT_THREADS=2\n'
                                'RTG_JAVA_OPTS="-XX:ActiveProcessorCount=2"\nRTG_TALKBACK=false\nRTG_USAGE=false\n')
    config['runtime_sha256']['rtg.cfg'] = digest(rtgdir/'rtg.cfg')
    config_path = root/'comparison-config.json'; save(config_path,config)
    runtime_pins = {name:pin(rtgdir/name) for name in ('rtg','RTG.jar','rtg.cfg')}
    runtime_pins['java'] = pin(java)
    save(root/'portable-runtime.json',{'files':runtime_pins,'platform':system_platform(),
         'note':'Same RTG wrapper/jar and comparison arguments. Java/build pins are from this platform lock; rtg.cfg only relocates Java.'})
    sys.path.insert(0,str(REPO/'workflow/scripts'))
    helpers = importlib.import_module('rtg_small_variant_helpers')
    def portable_runtime(requested_root, requested_config):
        require(requested_root.resolve()==root.resolve(),'Unexpected public root')
        require(requested_config == config,'Comparator configuration changed')
        host = subprocess.check_output(['hostname','-s'],text=True).strip()
        require(not Path('/etc/rtg.cfg').exists() and not (rtgdir/('rtg.'+host+'.cfg')).exists(),
                'Unreviewed RTG configuration override')
        for name,entry in runtime_pins.items():
            require(digest(entry['path'])==entry['sha256'],'Runtime changed: '+name)
        return rtgdir/'rtg', runtime_pins
    # Override only the historical Mac-specific runtime locator in memory. The
    # comparator, reference checks and twelve semantic controls remain unchanged.
    helpers.runtime = portable_runtime
    controls = root/'controls'
    result = importlib.import_module('probe_rtg_comparator').run(root,controls,config_path)
    require(result['cases_passed']==12 and not result['failed_cases'],'RTG controls failed')
    if mode == 'controls-only':
        save(root/'summary.json',{'status':'complete','mode':mode,'platform':system_platform(),
             'RTG_controls_passed':12,'benchmark_executed':False,'completed_at':stamp()})
        return
    gatk = profile['gatk']
    fetch(gatk['url'],downloads/'gatk.zip',size=gatk['archive_bytes'])
    jar = root/gatk['jar_member'];jar.parent.mkdir()
    with zipfile.ZipFile(downloads/'gatk.zip') as archive, archive.open(gatk['jar_member']) as src, jar.open('xb') as out:
        shutil.copyfileobj(src,out,8*1024*1024)
    require(digest(jar)==gatk['jar_sha256'],'GATK jar checksum differs')
    def get(item):
        folder = root/('ref/giab-small-variants-20260919/inputs' if item['kind']=='truth' else 'ref/giab-sv-20260919/inputs')
        return fetch(item['url'],folder/item['name'],item['sha256'],item['bytes'])
    with ThreadPoolExecutor(2) as pool: receipts = list(pool.map(get,profile['inputs']))
    save(root/'input-receipt.json',{'status':'complete','files':receipts,'public_only':True})
    truth_dir = root/'ref/giab-small-variants-20260919/inputs'
    save(truth_dir/'input-receipt.json',{'status':'complete','files':[r for r in receipts if str(truth_dir) in r['path']]})
    assembly = root/'ref/giab-sv-20260919/assembly-reference';assembly.mkdir()
    reference = assembly/'GRCh38_no_alt.fa'
    reference_gz = next(Path(r['path']) for r in receipts if r['path'].endswith('.fasta.gz'))
    with gzip.open(reference_gz,'rb') as src, reference.open('xb') as dst: shutil.copyfileobj(src,dst,8*1024*1024)
    require(digest(reference)=='9cce8b926416dd96b152deea85188495b75f7ac8d634cc723a017067be8702b7','Reference differs')
    pysam.faidx(str(reference))
    alignment = next(Path(r['path']) for r in receipts if r['path'].endswith('.bam'))
    with pysam.AlignmentFile(str(alignment)) as bam, pysam.FastaFile(str(reference)) as fasta:
        require(list(zip(bam.references,bam.lengths))==list(zip(fasta.references,fasta.lengths)), 'Public dictionaries differ')
        require(bam.header.to_dict().get('HD',{}).get('SO')=='coordinate' and bam.check_index(),'Indexed sorted BAM required')
        groups=bam.header.to_dict().get('RG',[])
        require(len(groups)==1 and groups[0].get('SM')=='HG002' and groups[0].get('ID')=='HG002.pcr-free','Unexpected public sample')
        require({s.contig:s.mapped for s in bam.get_index_statistics() if s.mapped}=={'chr20':19430746},'Unexpected public read scope')
    prepared = root/'prepared'
    importlib.import_module('prepare_public_small_variant_benchmark').run(root,prepared,config_path)
    pilot=root/'pilot';pilot.mkdir()
    command([sys.executable,REPO/'workflow/scripts/call_genome.py','--root',root,'--reference',reference,
             '--alignment',alignment,'--sample','HG002','--output',pilot/'calls','--workers',workers],logs,'call-complete-public-chromosome')
    native=json.loads((pilot/'calls/status.json').read_text());require(native['status']=='complete','Calling did not finish')
    save(pilot/'summary.json',{'status':'complete','outputs':native['outputs'],'native_completion':pin(pilot/'calls/status.json'),
         'scope':'Public prealigned chr20; unchanged all-reference-contig sharding, no non-chr20 mapped reads.'})
    comparison=root/'comparison'
    report=importlib.import_module('compare_public_small_variants').run(root,comparison,config_path,prepared,controls,pilot)
    command([sys.executable,REPO/'workflow/scripts/audit_public_small_variants.py','--comparison',comparison,
             '--output',root/'audit'],logs,'audit-public-record-accounting')
    for entry in receipts:require(digest(entry['path'])==entry['sha256'],'Public input changed')
    source_check(profile)
    summary={'status':'complete','completed_at':stamp(),'platform':system_platform(),'mode':mode,
             'public_only':True,'RTG_controls_passed':12,'identical_scored_truth_universe':report['identical_scored_truth_universe'],
             'counts':{k:v['counts'] for k,v in report['callers'].items()},'comparison':pin(comparison/'summary.json'),
             'native_completion':pin(pilot/'calls/status.json'),'audit':pin(root/'audit/summary.json'),
             'inputs':pin(root/'input-receipt.json'),'runtime':pin(root/'portable-runtime.json'),
             'scope':profile['scope'],'limits':config['limits'],
             'normalization':'Frozen primary bcftools norm -m -any representation; post-observation unsplit sensitivity is separate.'}
    save(root/'summary.json',summary)
    print(json.dumps(summary['counts']),flush=True)


def group_rss(pgid):
    total=0
    text=subprocess.check_output(['ps','-axo','pid=,rss='],text=True)
    for line in text.splitlines():
        parts=line.split()
        if len(parts)!=2:continue
        try:
            if os.getpgid(int(parts[0]))==pgid:total+=int(parts[1])*1024
        except ProcessLookupError:pass
    return total


def output_size(root):
    total = 0
    for path in root.rglob('*'):
        try:
            if path.is_file() and not path.is_symlink(): total += path.stat().st_size
        except FileNotFoundError:
            pass  # An atomic receipt replacement can remove a temporary file.
    return total


def guard(root, mode, workers):
    require(not root.exists(),'Output already exists; it is preserved. Choose a new directory.')
    require(not root.is_relative_to(REPO),'Output must be outside the clone')
    require(shutil.disk_usage(root.parent).free>=LIMITS['minimum_free_bytes'],'At least40GiB free required')
    root.mkdir()
    start=time.monotonic();peak=0;maxbytes=0;reason=None
    with (root/'worker.stdout').open('w') as out, (root/'worker.stderr').open('w') as err:
        child=subprocess.Popen([sys.executable,str(Path(__file__)),str(root),'--bootstrap-worker','--mode',mode,
                                '--workers',str(workers)],stdout=out,stderr=err,start_new_session=True)
        state={'status':'running','started_at':stamp(),'parent_pid':os.getpid(),'owned_pgid':child.pid,
               'platform':system_platform(),'mode':mode,'limits':LIMITS,'scope':'Public downloads and benchmark only'}
        save(root/'status.json',state)
        print(json.dumps({'status':'running','output':str(root),'owned_pgid':child.pid}),flush=True)
        try:
            while child.poll() is None:
                rss=group_rss(child.pid);peak=max(peak,rss);size=output_size(root);maxbytes=max(maxbytes,size)
                if time.monotonic()-start>LIMITS['seconds']:reason='timeout'
                elif rss>LIMITS['rss_bytes']:reason='group RSS limit'
                elif size>LIMITS['output_bytes']:reason='output limit'
                elif shutil.disk_usage(root).free<LIMITS['minimum_free_bytes']:reason='free disk limit'
                if reason:raise RuntimeError(reason)
                time.sleep(5)
            require(child.returncode==0,'Public replay worker failed; logs retained')
            summary=json.loads((root/'summary.json').read_text());require(summary['status']=='complete','Missing completion')
            state.update(status='complete',summary=pin(root/'summary.json'),exit_code=child.returncode)
        except BaseException as error:
            if child.poll() is None:
                os.killpg(child.pid,signal.SIGTERM)
                try:child.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    os.killpg(child.pid,signal.SIGKILL);child.wait()
            state.update(status='failed',error=str(error),exit_code=child.returncode)
            raise
        finally:
            state.update(completed_at=stamp(),elapsed_seconds=time.monotonic()-start,
                         peak_sampled_group_rss_bytes=peak,peak_sampled_output_bytes=maxbytes,
                         limits_are_sampled_not_instantaneous=True)
            save(root/'status.json',state)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('output',type=Path)
    p.add_argument('--controls-only',action='store_true',help='Install pinned runtime and run12RTG controls; no public benchmark claim')
    p.add_argument('--workers',type=int,default=2,choices=range(1,5))
    p.add_argument('--bootstrap-worker',action='store_true',help=argparse.SUPPRESS)
    p.add_argument('--pipeline-worker',action='store_true',help=argparse.SUPPRESS)
    p.add_argument('--mode',choices=['full','controls-only'],help=argparse.SUPPRESS)
    a=p.parse_args();root=a.output.expanduser().resolve();mode=a.mode or ('controls-only' if a.controls_only else 'full')
    require(sys.version_info>=(3,9),'System Python3.9+ required')
    if a.bootstrap_worker:bootstrap(root,mode,a.workers)
    elif a.pipeline_worker:pipeline(root,mode,a.workers)
    else:guard(root,mode,a.workers)

if __name__=='__main__':main()
