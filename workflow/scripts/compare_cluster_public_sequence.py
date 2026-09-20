#!/usr/bin/env python3
"""Compare retained cluster reads with an independently archived public clone.

The small comparison reference does not supply genome-wide mapping confidence.
Full-reference alignment of the public clone provides coordinate context only.
"""
import argparse
from collections import Counter, defaultdict
import datetime
import json
from pathlib import Path
import re
import subprocess

import pysam

from review_cnv_context import fingerprint, require, save, sha
from review_paralog_candidate_pair import remap_details


def parse_embl(text):
    lines = text.splitlines()
    match = re.fullmatch(r'ID\s+(\w+); SV (\d+);.*; (\d+) BP\.',lines[0])
    require(match is not None and sum(line.startswith('ID ') for line in lines)==1,'Expected one EMBL DNA record')
    require(lines[-1]=='//','Incomplete EMBL record')
    sq = [i for i,line in enumerate(lines) if line.startswith('SQ ')]
    require(len(sq)==1,'Missing or duplicate sequence section')
    sequence = ''.join(''.join(c for c in line if c.isalpha()) for line in lines[sq[0]+1:-1]).upper()
    require(len(sequence)==int(match.group(3)) and set(sequence)<=set('ACGTN'),'Sequence length/alphabet differs')
    return match.group(1)+'.'+match.group(2),sequence


def complete_exact(row):
    return not row['flag'] & 4 and row['query_alignment_length']==row['query_length'] and row['tags'].get('NM')==0


def score_difference(a,b):
    if a['flag'] & 4 or b['flag'] & 4 or 'AS' not in a['tags'] or 'AS' not in b['tags']:
        return None
    return a['tags']['AS']-b['tags']['AS']


def plausible_pair(a,b,maximum_span):
    if a['flag'] & 4 or b['flag'] & 4 or a['chrom']!=b['chrom']:
        return False
    if any(r['query_alignment_length']<0.9*r['query_length'] for r in (a,b)):
        return False
    if bool(a['flag'] & 16)==bool(b['flag'] & 16):
        return False
    forward,reverse = (b,a) if a['flag'] & 16 else (a,b)
    span=max(a['end_exclusive'],b['end_exclusive'])-min(a['start_0based'],b['start_0based'])
    return forward['start_0based']<=reverse['start_0based'] and max(a['query_length'],b['query_length'])<=span<=maximum_span


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['root','config','output']:
        parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args();root=args.root.resolve();output=args.output.resolve()
    require(not output.exists() and not output.is_relative_to(Path(__file__).resolve().parents[2]),'New external output required')
    config=json.loads(args.config.read_text());parent=root/config['cluster_directory'];public=root/config['public_embl']
    require(sha(parent/'summary.json')==config['cluster_summary_sha256'],'Cluster summary differs')
    prior=json.loads((parent/'summary.json').read_text())
    for name in ['queries.fa','read-review.json','remapping.json']:
        require(sha(parent/name)==prior['outputs'][name]['sha256'],'Cluster input differs')
    require(sha(public)==config['public_embl_sha256'],'Public record differs')
    name,sequence=parse_embl(public.read_text())
    require(name==config['accession_version'],'Unexpected accession version')
    bwa=root/config['bwa_binary'];reference=root/config['genome_reference']
    require(sha(bwa)==prior['inputs']['bwa']['sha256'],'BWA differs')
    require(sha(reference)==prior['inputs']['reference']['sha256'],'Genome reference differs')
    output.mkdir(parents=True)
    fasta=output/'public-clone.fa';fasta.write_text('>'+name+'\n'+sequence+'\n')
    commands=[]
    def run(command,label,sam=None):
        commands.append(list(map(str,command)));save(output/'commands.json',commands)
        with (output/(label+'.log')).open('w') as log:
            if sam is None:subprocess.run(commands[-1],stdout=log,stderr=subprocess.STDOUT,check=True)
            else:
                with sam.open('w') as stream:subprocess.run(commands[-1],stdout=stream,stderr=log,check=True)
    run([bwa,'index',fasta],'index-public-clone')
    run([bwa,'mem','-a','-Y','-t','2',fasta,parent/'queries.fa'],'map-reads-to-clone',output/'reads-to-clone.sam')
    run([bwa,'mem','-a','-Y','-t','2',reference,fasta],'map-clone-to-genome',output/'clone-to-genome.sam')
    remaps=defaultdict(list)
    with pysam.AlignmentFile(str(output/'reads-to-clone.sam'),'r') as sam:
        for read in sam:
            row=remap_details(read);row['query_length']=read.infer_read_length()
            remaps[read.query_name].append(row)
    primary={}
    for identifier,rows in remaps.items():
        eligible=[r for r in rows if not r['flag'] & (256|2048)]
        require(len(eligible)==1,'Unexpected primary census')
        primary[identifier]=eligible[0]
    old=json.loads((parent/'remapping.json').read_text())
    require(set(primary)==set(old),'Query identities differ')
    reads=json.loads((parent/'read-review.json').read_text());contexts=[]
    for site,configuration in zip(prior['configuration']['sites'],reads['configurations']['20']):
        groups={}
        for label in configuration['counts']:
            fragments=[f for f in configuration['fragments'] if f['label']==label]
            ids={i for f in fragments for i in f['observation_ids']}
            values=[primary[i] for i in ids]
            concordant=[];exact_mates=[]
            for fragment in fragments:
                pair=False;exact=False
                for identifier in fragment['observation_ids']:
                    mid=reads['observations'][identifier].get('mate_id')
                    if mid:
                        m=reads['mates'][mid]
                        if m['flag'] & 0xF04:continue
                        pair |= plausible_pair(primary[identifier],primary[mid],config['maximum_fragment_span'])
                        exact |= complete_exact(primary[identifier]) and complete_exact(primary[mid])
                concordant.append(pair);exact_mates.append(exact)
            groups[label]={'fragments':len(fragments),'read_observations':len(ids),
                           'full_length_exact_read_observations':sum(complete_exact(v) for v in values),
                           'fragments_with_any_full_length_exact_read':sum(any(complete_exact(primary[i]) for i in f['observation_ids']) for f in fragments),
                           'fragments_with_both_ends_full_length_exact':sum(exact_mates),
                           'fragments_with_plausible_separately_mapped_pair':sum(concordant),
                           'mapped_read_NM':dict(Counter(str(v['tags'].get('NM')) for v in values)),
                           'mapped_read_CIGAR':dict(Counter(str(v['cigar']) for v in values)),
                           'public_clone_minus_full_reference_primary_AS':dict(Counter(str(score_difference(primary[i],old[i]['primary']))
                                if score_difference(primary[i],old[i]['primary']) is not None else 'unreported_comparison' for i in ids)),
                           'public_clone_primary_intervals_0based':sorted({(v['start_0based'],v['end_exclusive']) for v in values if not v['flag'] & 4})}
        contexts.append({'site':site,'groups':groups})
    with pysam.AlignmentFile(str(output/'clone-to-genome.sam'),'r') as sam:
        genome=[remap_details(read) for read in sam]
    save(output/'read-mappings.json',dict(remaps))
    summary={'completed_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'configuration':config,
             'inputs':{n:fingerprint(p) for n,p in {'public_embl':public,'cluster_summary':parent/'summary.json','queries':parent/'queries.fa',
                       'source':Path(__file__),'helper':Path(__file__).with_name('review_paralog_candidate_pair.py'),
                       'io_helper':Path(__file__).with_name('review_cnv_context.py'),'config':args.config,'bwa':bwa,'genome_reference':reference}.items()},
             'public_sequence_length':len(sequence),'public_sequence_sha256':sha(fasta),'query_count':len(primary),
             'per_site':contexts,'public_clone_to_genome':genome,
             'outputs':{p.name:fingerprint(p) for p in output.iterdir() if p.is_file()},
             'limits':['Small-reference MQ is not genome-wide confidence; score differences are aligner heuristics, not likelihood ratios.',
                       'Archived clone sequence is independent public material; it is not a patient assembly, parental haplotype or proof of locus origin.',
                       'Plausible pair means both ends >=90% aligned, inward orientation and configured span bound after separate mapping, not an independent paired-end aligner call.',
                       'Per-site groups overlap. A sequence match does not establish its chromosomal location, copy number, a gene-conversion mechanism or clinical significance.']}
    save(output/'summary.json',summary)
    print(json.dumps({'status':'complete','queries':len(primary),'summary_sha256':sha(output/'summary.json')}))


if __name__=='__main__':main()
