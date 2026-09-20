#!/usr/bin/env python3
"""Native, resumable diploid small-variant discovery on every reference contig.

This produces raw discovery callsets, not clinical classifications. Sex-chromosome,
mitochondrial, mosaic, structural-variant and callability analyses are separate gates.
Whole contigs are assigned to shards so no assembly interval is cut mid-contig.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import fcntl
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

import pysam


def stamp():
    return time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as source:
        for chunk in iter(lambda: source.read(8*1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def save(path, data):
    temporary = Path(str(path)+'.tmp')
    temporary.write_text(json.dumps(data, indent=2)+'\n')
    os.replace(temporary, path)


def shard_contigs(contigs, maximum_bases=20_000_000, maximum_contigs=64):
    """Preserve reference order and assign each complete contig exactly once."""
    shards, current, bases = [], [], 0
    for name, length in contigs:
        if length <= 0:
            raise ValueError('invalid contig length')
        if current and (bases+length > maximum_bases or len(current) >= maximum_contigs):
            shards.append(current)
            current, bases = [], 0
        current.append((name, length))
        bases += length
        if bases >= maximum_bases:
            shards.append(current)
            current, bases = [], 0
    if current:
        shards.append(current)
    return shards


def inspect_vcf(path, sample):
    counts = {'records': 0, 'snv_records': 0, 'indel_records': 0, 'other_records': 0}
    if not any(Path(str(path)+extension).is_file() for extension in ('.tbi','.csi')):
        raise ValueError('VCF index missing')
    with pysam.VariantFile(str(path)) as vcf:
        if list(vcf.header.samples) != [sample]:
            raise ValueError('unexpected VCF sample header')
        if vcf.index is None:
            raise ValueError('VCF index could not be read')
        # The indexed iterator requires a working index, including for an empty valid VCF.
        for record in vcf.fetch():
            counts['records'] += 1
            alts = record.alts or ()
            if alts and len(record.ref) == 1 and all(len(a) == 1 and a != '*' for a in alts):
                counts['snv_records'] += 1
            elif alts and all(a != '*' and not a.startswith('<') for a in alts) and any(len(a) != len(record.ref) for a in alts):
                counts['indel_records'] += 1
            else:
                counts['other_records'] += 1
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--alignment', type=Path, required=True)
    parser.add_argument('--sample', required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--shard-max-bases', type=int, default=20_000_000)
    args = parser.parse_args()
    if not 1 <= args.workers <= 4:
        parser.error('use one to four workers, at most eight advertised JVM CPUs')
    if args.shard_max_bases < 1:
        parser.error('shard maximum must be positive')
    root, reference, alignment, output = [p.resolve() for p in (args.root,args.reference,args.alignment,args.output)]
    output.mkdir(parents=True, exist_ok=True)
    runtime = root/'tools/bio'
    java = runtime/'lib/jvm/bin/java'
    jar = root/'gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar'
    env = dict(os.environ, PATH=str(runtime/'bin')+os.pathsep+os.environ.get('PATH',''))
    os.nice(10)
    with (output/'.calling.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        def run(command, log):
            with log.open('a') as stream:
                subprocess.run([str(x) for x in command], env=env, stdout=stream,
                               stderr=subprocess.STDOUT, check=True, pass_fds=(lock.fileno(),))
        run(['samtools','quickcheck','-v',alignment], output/'preflight.log')
        with pysam.FastaFile(str(reference)) as fasta:
            contigs = list(zip(fasta.references, fasta.lengths))
        with pysam.AlignmentFile(str(alignment),reference_filename=str(reference)) as bam:
            if list(zip(bam.references,bam.lengths)) != contigs:
                raise ValueError('alignment and reference dictionaries differ')
            bam.check_index()
            if {g.get('SM') for g in bam.header.to_dict().get('RG',[])} != {args.sample}:
                raise ValueError('alignment read groups do not name exactly the expected sample')
        shards = shard_contigs(contigs,args.shard_max_bases)
        definition = {'reference_sha256':digest(reference),'alignment_sha256':digest(alignment),
                      'gatk_jar_sha256':digest(jar),'script_sha256':digest(Path(__file__)),
                      'sample':args.sample,'contigs':len(contigs),'reference_bases':sum(n for _,n in contigs),
                      'shards':shards,'ploidy':2,'gatk_pair_hmm':'LOGLESS_CACHING',
                      'bcftools_minimum_mapping_quality':20,'bcftools_minimum_base_quality':20,
                      'bcftools_maximum_depth':10000,'bcftools_excluded_flags':'0xF04',
                      'scope':'All reference contigs, diploid small-variant discovery; raw calls, no BQSR or clinical interpretation.'}
        definition_path = output/'run-definition.json'
        if definition_path.exists():
            if json.loads(definition_path.read_text()) != json.loads(json.dumps(definition)):
                raise ValueError('inputs or method changed; existing run preserved, choose a new output')
        elif any(p.name not in ('.calling.lock','preflight.log') for p in output.iterdir()):
            raise ValueError('unlabelled existing output preserved')
        else:
            save(definition_path, definition)
        dictionary = reference.with_suffix('.dict')
        if not dictionary.exists():
            run(['samtools','dict','-o',str(dictionary)+'.tmp',reference], output/'preflight.log')
            os.replace(str(dictionary)+'.tmp',dictionary)
        state_path = output/'status.json'
        if state_path.exists():
            prior = json.loads(state_path.read_text())
            if prior.get('status') == 'complete':
                for result in prior['outputs'].values():
                    if digest(result['path']) != result['sha256']:
                        raise ValueError('completed final callset changed')
                    inspect_vcf(Path(result['path']),args.sample)
                print('[done] existing verified callsets retained',flush=True)
                return
        state = {'status':'running','pid':os.getpid(),'started_at':stamp(),'workers':args.workers,
                 'scope':definition['scope'],'completed_shards':[],'expected_shards':len(shards),
                 'definition_sha256':digest(definition_path)}
        if state_path.exists():
            save(output/('previous-status-%d.json'%time.time_ns()),json.loads(state_path.read_text()))
        save(state_path,state)
        base_java = [java,'-XX:ActiveProcessorCount=2','-Xmx6g','-jar',jar]

        def process_shard(index, intervals):
            folder = output/('shard-%04d'%index)
            folder.mkdir(exist_ok=True)
            receipt = folder/'complete.json'
            if receipt.exists():
                prior = json.loads(receipt.read_text())
                for method in ('gatk','bcftools'):
                    p = folder/prior[method]['path']
                    if digest(p) != prior[method]['sha256']:
                        raise ValueError('completed shard changed')
                    inspect_vcf(p,args.sample)
                return index,prior
            attempt = folder/('attempt-%d'%time.time_ns())
            attempt.mkdir()
            bed = attempt/'regions.bed'
            bed.write_text(''.join('%s\t0\t%d\n'%(c,n) for c,n in intervals))
            commands = []
            def execute(command,name):
                commands.append([str(x) for x in command])
                save(attempt/'commands.json',commands)
                run(command,attempt/(name+'.log'))
            started = time.monotonic()
            try:
                raw = attempt/'gatk.raw.vcf.gz'
                execute(base_java+['HaplotypeCaller','-R',reference,'-I',alignment,'-O',raw,'-L',bed,
                        '--sample-ploidy','2','--pair-hmm-implementation','LOGLESS_CACHING'], 'gatk')
                execute(['bcftools','mpileup','-Ob','-o',attempt/'pileup.bcf','-f',reference,'-R',bed,
                        '-a','AD,DP,SP','-q','20','-Q','20','-d','10000','--ff','0xF04',alignment], 'mpileup')
                execute(['bcftools','call','-mv','--ploidy','2','-Ob','-o',attempt/'bcftools.raw.bcf',attempt/'pileup.bcf'], 'bcftools-call')
                result = {'completed_at':stamp(),'commands_sha256':digest(attempt/'commands.json')}
                for method,source in [('gatk',raw),('bcftools',attempt/'bcftools.raw.bcf')]:
                    p = attempt/(method+'.normalized.raw.vcf.gz')
                    execute(['bcftools','norm','--check-ref','e','-f',reference,'-m','-any','-Oz','-o',p,source], method+'-normalize')
                    execute(['bcftools','index','-t',p],method+'-index')
                    result[method] = {'path':str(p.relative_to(folder)),'sha256':digest(p),'counts':inspect_vcf(p,args.sample)}
                result['commands_sha256'] = digest(attempt/'commands.json')
                result['elapsed_seconds'] = round(time.monotonic()-started,2)
                save(receipt,result)
                return index,result
            except BaseException as error:
                save(attempt/'failure.json',{'failed_at':stamp(),'error_type':type(error).__name__})
                raise

        try:
            completed = {}
            with ThreadPoolExecutor(max_workers=args.workers) as executor:
                futures = [executor.submit(process_shard,i,intervals) for i,intervals in enumerate(shards)]
                for future in as_completed(futures):
                    try:
                        index,result = future.result()
                    except BaseException:
                        for queued in futures:
                            queued.cancel()
                        raise
                    completed[index] = result
                    state['completed_shards'] = sorted(completed)
                    save(state_path,state)
                    print('[shard] %d/%d complete'%(len(completed),len(shards)),flush=True)
            final = {}
            for method in ('gatk','bcftools'):
                method_receipt = output/(method+'-complete.json')
                if method_receipt.exists():
                    prior = json.loads(method_receipt.read_text())
                    if digest(prior['path']) != prior['sha256']:
                        raise ValueError('published callset changed')
                    inspect_vcf(Path(prior['path']),args.sample)
                    final[method] = prior
                    continue
                paths = [output/('shard-%04d'%i)/completed[i][method]['path'] for i in range(len(shards))]
                list_file = output/(method+'-inputs.list')
                list_file.write_text(''.join(str(p)+'\n' for p in paths))
                published = output/(args.sample+'.'+method+'.normalized.raw.vcf.gz')
                if published.exists():
                    raise ValueError('existing final callset preserved; inspect before republishing')
                gather = output/('gather-%s-%d'%(method,time.time_ns()))
                gather.mkdir()
                temporary = gather/(args.sample+'.'+method+'.normalized.raw.partial.vcf.gz')
                run(['bcftools','concat','-f',list_file,'-Oz','-o',temporary],output/(method+'-gather.log'))
                run(['bcftools','index','-t',temporary],output/(method+'-gather.log'))
                counts = inspect_vcf(temporary,args.sample)
                expected = sum(completed[i][method]['counts']['records'] for i in completed)
                if counts['records'] != expected:
                    raise ValueError('gather changed record count')
                os.replace(temporary,published)
                os.replace(str(temporary)+'.tbi',str(published)+'.tbi')
                final[method] = {'path':str(published),'sha256':digest(published),'counts':counts}
                save(method_receipt,final[method])
            state.update(status='complete',completed_at=stamp(),outputs=final,
                         next_gates=['Quality filtering and callability','Ploidy-aware sex/mitochondrial review',
                                     'Caller concordance, independent prioritisation, classification and phenotype review',
                                     'Structural/mosaic variation and BQSR sensitivity where justified'])
            save(state_path,state)
            print('[done] both raw discovery callsets validated; downstream gates remain',flush=True)
        except BaseException as error:
            state.update(status='failed',failed_at=stamp(),error_type=type(error).__name__)
            save(state_path,state)
            raise


if __name__ == '__main__':
    main()
