#!/usr/bin/env python3
"""Bounded mtDNA mixture screen using pinned mutserve and haplocheck runtimes.

Preserves original whole-genome mappings and reports both supported SNPs and a
hotspot-excluded sensitivity arm. No numeric result is transferred to Mutect2 or
interpreted as nuclear contamination, absence of contamination, or diagnosis.
"""
import argparse
from collections import Counter
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess

import pysam
from prepare_gatk_compatible_cram import digest,save,stamp
from prepare_haplocheck_input import qualify,verify_import

MUTSERVE_SHA='bb741acf57ed23388082260358aa44794b00006506297b9b115fcc3d58096387'
HAPLOCHECK_SHA='2b30f573f98395593076d32d9e377818c3feb4505a08a53151a39656cf075782'


def select_read(read):
    if read.flag & 0xF04:
        return 'flag_excluded_0xF04'
    if read.mapping_quality == 255:
        return 'mapping_quality_unknown_255'
    if not read.has_tag('AS') or not isinstance(read.get_tag('AS'), int):
        return 'alignment_score_missing_or_noninteger'
    return None


def check_reports(out, label, expected_sample):
    reports={}
    for suffix in ['.txt','.raw.txt','.html']:
        p=out/(label+suffix)
        if not p.is_file() or p.stat().st_size==0:
            raise ValueError('missing/empty fresh haplocheck report: '+p.name)
    for suffix in ['.txt','.raw.txt']:
        with (out/(label+suffix)).open() as handle:
            rows=list(csv.DictReader(handle,delimiter='\t'))
        if len(rows)!=1 or rows[0].get('Sample')!=expected_sample:
            raise ValueError('haplocheck report sample accounting differs')
        row=rows[0]
        for key in ['Contamination Status','Contamination Level','Distance','Sample Coverage']:
            if key not in row:
                raise ValueError('haplocheck report schema differs')
        if row['Contamination Status'] not in ['YES','NO','ND']:
            raise ValueError('unknown haplocheck status')
        if row['Contamination Status']=='YES':
            value=float(row['Contamination Level'])
            if not math.isfinite(value) or not 0<=value<=1:
                raise ValueError('invalid haplocheck level')
        elif row['Contamination Level']!='ND':
            raise ValueError('unexpected numeric level on NO/ND report')
        reports[suffix]=row
    for key in ['Contamination Status','Contamination Level','Distance','Sample Coverage']:
        if reports['.txt'][key]!=reports['.raw.txt'][key]:
            raise ValueError('haplocheck report views disagree')
    return reports


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['root','reference','alignment','config','qualification','parser-probe','output']:
        parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args()
    os.nice(10)
    root,out=args.root.resolve(),args.output.resolve()
    if out.exists():
        raise ValueError('existing run preserved')
    if '.' in str(out):
        raise ValueError('mutserve derives output prefixes at the first dot; use a directory path without dots')
    config=json.loads(args.config.read_text())
    control=json.loads(args.qualification.read_text())
    if control['status']!='qualified_with_explicit_exclusions':
        raise ValueError('native method qualification unavailable')
    for rel,info in control['files'].items():
        if digest(root/rel)!=info['sha256']:
            raise ValueError('native method control changed: '+rel)
    mutserve=root/'tools/mutserve-2.0.0-rc12/mutserve.jar'
    haplocheck=root/'tools/haplocheck-1.3.3/haplocheck.jar'
    bundled_ref=root/'tools/mutserve-2.0.0-rc12/rCRS.fasta'
    for path,expected in [(mutserve,MUTSERVE_SHA),(haplocheck,HAPLOCHECK_SHA),
                          (args.reference,config['reference_sha256']),
                          (args.alignment,config['alignment_sha256'])]:
        if digest(path)!=expected:
            raise ValueError('frozen input identity differs: '+str(path))
    seq=''.join(x for x in bundled_ref.read_text().splitlines() if not x.startswith('>')).upper()
    if len(seq)!=16569 or hashlib.sha256(seq.encode()).hexdigest()!=config['mitochondrial_sequence_sha256']:
        raise ValueError('bundled mitochondrial reference differs')
    with pysam.FastaFile(str(args.reference)) as ref,pysam.AlignmentFile(str(args.alignment),reference_filename=str(args.reference)) as bam:
        if ref.fetch(config['input_contig']).upper()!=seq:
            raise ValueError('mitochondrial sequences differ')
        if list(zip(ref.references,ref.lengths))!=list(zip(bam.references,bam.lengths)):
            raise ValueError('alignment/reference dictionaries differ')
        bam.check_index()
        rg=bam.header.to_dict().get('RG',[])
        if len(rg)!=config['expected_read_groups'] or {x['SM'] for x in rg}!={config['sample']}:
            raise ValueError('read group/sample mismatch')
    out.mkdir(parents=True,exist_ok=False)
    java=root/'tools/bio/lib/jvm/bin/java'
    source_paths=[Path(__file__),Path(__file__).with_name('prepare_haplocheck_input.py'),
                  Path(__file__).with_name('prepare_gatk_compatible_cram.py'),args.parser_probe,
                  args.config,args.qualification,mutserve,haplocheck,bundled_ref,java]
    state=dict(status='running',pid=os.getpid(),started_at=stamp(),config=config,
               source_and_tool_sha256={str(p):digest(p) for p in source_paths},
               scope='Exploratory mtDNA mixture screen on original whole-genome mappings, with native-format and hotspot sensitivity controls.')
    commands=[]
    def run(command,label):
        state['current_stage']=label;save(out/'status.json',state)
        row=dict(label=label,command=list(map(str,command)),started_at=stamp())
        commands.append(row);save(out/'commands.json',commands)
        with (out/(label+'.log')).open('w') as stream:
            code=subprocess.run(row['command'],cwd=out,stdout=stream,stderr=subprocess.STDOUT).returncode
        row.update(exit_code=code,completed_at=stamp());save(out/'commands.json',commands)
        log=(out/(label+'.log')).read_text()
        if code or re.search(r'\bERROR\b|Exception|failed|terminated',log,re.I):
            raise RuntimeError('native tool failure or error log: '+label)
    try:
        state['current_stage']='select-original-mitochondrial-reads';save(out/'status.json',state)
        selected=out/'input.bam';counts=Counter();groups=Counter();mapping=Counter();mates=Counter()
        with pysam.AlignmentFile(str(args.alignment),reference_filename=str(args.reference)) as src,pysam.AlignmentFile(str(selected),'wb',header=src.header) as dest:
            tid=src.get_tid(config['input_contig'])
            for read in src.fetch(config['input_contig']):
                counts['all_input_contig_records']+=1
                reason=select_read(read)
                if reason:
                    counts[reason]+=1;continue
                dest.write(read);counts['selected_records']+=1
                groups[read.get_tag('RG')]+=1
                mapping[str(read.mapping_quality)]+=1
                mates['unmapped' if read.mate_is_unmapped else 'same_contig' if read.next_reference_id==tid else 'different_contig']+=1
        if not counts['selected_records']:
            raise ValueError('no selected reads')
        if set(groups)-{x['ID'] for x in rg}:
            raise ValueError('unexpected selected read group')
        pysam.index(str(selected))
        selection=dict(counts=dict(counts),read_groups=dict(groups),mapping_quality=dict(mapping),mate_location=dict(mates),
                       selected_bam_sha256=digest(selected),
                       scope='Original alignments, sequence, qualities, flags and scores retained. No mitochondrial-only remapping or same-contig mate requirement. Flags 0xF04, MQ255 and absent/noninteger AS excluded; native MQ30/BQ20/AS30/read-length filters run next. Overlapping mates remain separate observations.')
        save(out/'read-selection.json',selection)
        ref=out/'reference.fasta';shutil.copy2(bundled_ref,ref);pysam.faidx(str(ref))
        run([java,'-XX:ActiveProcessorCount=2','-Xmx2g','-jar',mutserve,'call','--level','0.01','--reference',ref,
             '--mapQ','30','--baseQ','20','--alignQ','30','--deletions','--write-raw','--contig-name',config['input_contig'],
             '--output',out/'calls.vcf.gz','--threads','2','--no-ansi',selected],'mutserve')
        for name in ['calls.vcf.gz','calls.txt','calls_raw.txt']:
            if not (out/name).is_file() or (out/name).stat().st_size==0:
                raise ValueError('missing native caller output: '+name)
        pysam.tabix_index(str(out/'calls.vcf.gz'),preset='vcf',force=False)
        shutil.copy2(args.parser_probe,out/'HaplocheckParserProbe.java')
        run([java.with_name('javac'),'-cp',haplocheck,out/'HaplocheckParserProbe.java'],'compile-parser-probe')
        arms={};sample='input.bam'
        for label,masked in [('supported',False),('hotspots-excluded',True)]:
            vcf=out/(label+'.vcf')
            qualification=qualify(out/'calls.vcf.gz',out/'calls.txt',seq,sample,vcf,masked)
            save(out/(label+'-qualification.json'),qualification)
            if qualification['accepted_records']==0:
                arms[label]=dict(qualification=qualification,status='uninformative_no_qualified_variants')
                continue
            run([java,'-XX:ActiveProcessorCount=2','-Xmx2g','-cp',str(haplocheck)+':'+str(out),'HaplocheckParserProbe',vcf,out/(label+'-parser.json')],label+'-parser')
            imported=verify_import(json.loads((out/(label+'-parser.json')).read_text()),qualification)
            run([java,'-XX:ActiveProcessorCount=2','-Xmx2g','-jar',haplocheck,'--raw','--out',out/(label+'.txt'),vcf],label+'-haplocheck')
            arms[label]=dict(qualification=qualification,import_verification=imported,reports=check_reports(out,label,sample))
        raw={}
        with (out/'calls_raw.txt').open() as handle:
            for row in csv.DictReader(handle,delimiter='\t'):
                pos=int(row['POS'])
                if row['SAMPLE']!=sample or pos in raw or not 1<=pos<=len(seq) or row['REF']!=seq[pos-1]:
                    raise ValueError('raw depth sample/position/reference mismatch')
                raw[pos]=[int(row[k]) for k in ['COV-FWD','COV-REV','COV-TOTAL']]
                if raw[pos][0]+raw[pos][1]!=raw[pos][2]:
                    raise ValueError('raw depth accounting differs')
        concrete=[p for p,b in enumerate(seq,1) if b in 'ACGT']
        depth=[raw.get(p,[0,0,0])[2] for p in concrete]
        coverage=dict(concrete_positions=len(concrete),positions_with_raw_entry=len(raw),
                      mean_read_observation_depth=sum(depth)/len(depth),min_depth=min(depth),max_depth=max(depth),
                      no_qualified_depth=[p for p in concrete if raw.get(p,[0,0,0])[2]==0],
                      below_10_per_strand=[p for p in concrete if min(raw.get(p,[0,0,0])[:2])<10],
                      scope='Native MQ30/BQ20/AS30 pileup read observations, including its deletion counting. Overlapping mates are not collapsed; depth is not molecule count or a detection limit.')
        save(out/'summary.json',dict(completed_at=stamp(),host='512S2',sample_mapping={sample:config['sample']},
             read_selection=selection,coverage=coverage,arms=arms,
             limitations=['Native homoplasmic is an operational caller class, not demonstrated molecular homogeneity.',
             'Native SNP classification at a 1 percent setting has no independently measured clinical detection limit for this library.',
             'The linear original whole-genome mapping can reduce origin-spanning coverage; missing/unsupported sites are not reference evidence.',
             'Mitochondrial haplotype mixing does not measure nuclear-genome contamination or validate a low-fraction pathogenic variant.',
             'NO or ND does not mean zero contamination; same/indistinguishable haplotypes and mixtures below informative-marker sensitivity remain possible.',
             'No filter contamination estimate or causal ranking is changed by this screen.']))
        state.update(status='complete',completed_at=stamp())
        state['output_sha256']={p.name:digest(p) for p in sorted(out.iterdir()) if p.is_file() and p.name!='status.json'}
        save(out/'status.json',state)
        print(json.dumps(dict(status='complete',output=str(out),reports={k:v.get('reports') for k,v in arms.items()})))
    except Exception as exc:
        state.update(status='failed',failed_at=stamp(),error_type=type(exc).__name__,error=str(exc));save(out/'status.json',state)
        raise


if __name__=='__main__':
    main()
