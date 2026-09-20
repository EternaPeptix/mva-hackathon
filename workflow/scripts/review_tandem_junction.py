#!/usr/bin/env python3
"""Review one precise duplication's literal junction and conditional transcript geometry.

All read-scale/sequence outputs must stay outside the source repository. Signatures
from the same input are not biological replication, calibrated VAF or a diagnosis.
"""
import argparse
from collections import Counter, defaultdict
import datetime
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

from review_cnv_context import fingerprint, overlap, require, save, sha, verify


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGTNacgtn', 'TGCANtgcan'))[::-1]


def occurrences(sequence, needle):
    start = 0
    while True:
        i = sequence.find(needle, start)
        if i < 0:
            return
        yield i
        start = i + 1


def decompose_consensus(consensus, windows, minimum_anchor):
    """Windows are (0-based start, forward reference sequence); retain every exact solution."""
    out = []
    for direction, sequence in [('native', consensus), ('reverse_complement', reverse_complement(consensus))]:
        for split in range(minimum_anchor, len(sequence) - minimum_anchor + 1):
            for first, second in [('right', 'left'), ('left', 'right')]:
                a, left = windows[first]; b, right = windows[second]
                for i in occurrences(left, sequence[:split]):
                    for j in occurrences(right, sequence[split:]):
                        out.append({'orientation': direction, 'order': [first, second], 'split_offset': split,
                                    'prefix_last_reference_base_1based': a + i + split,
                                    'suffix_first_reference_base_1based': b + j + 1})
    return out


def exact_high_quality(sequence, qualities, needle, minimum_quality):
    if not sequence or qualities is None:
        return False
    require(len(sequence) == len(qualities), 'Sequence/quality lengths differ')
    return any(min(qualities[i:i+len(needle)]) >= minimum_quality
               for query in (needle, reverse_complement(needle)) for i in occurrences(sequence, query))


def query_span(cigar):
    """Oriented original-query coordinates, retaining both soft and hard clips."""
    require(bool(cigar), 'Missing CIGAR')
    length = sum(n for op, n in cigar if op in (0, 1, 4, 5, 7, 8))
    front = back = 0
    for op, n in cigar:
        if op not in (4, 5): break
        front += n
    for op, n in reversed(cigar):
        if op not in (4, 5): break
        back += n
    require(front + back < length, 'No aligned query segment')
    return front, length-back, length


def transcript_geometry(transcript, exons, coding, left, right):
    """left/right are 1-based literal suffix-start/prefix-end, not normalized VCF alleles."""
    ranked = sorted(exons, key=lambda e: int(e['attributes']['rank']))
    require([int(e['attributes']['rank']) for e in ranked] == list(range(1, len(ranked)+1)), 'Exon ranks incomplete')
    ordered = sorted(ranked, key=lambda e: e['start'])
    require(all(a['end'] < b['start'] for a,b in zip(ordered, ordered[1:])), 'Overlapping transcript exons')
    require(all(any(e['start'] <= c['start'] <= c['end'] <= e['end'] for e in exons) for c in coding), 'CDS outside exon')
    introns = []
    for a,b in zip(ordered, ordered[1:]):
        if a['end'] < left <= b['start']:
            introns.append({'start_1based': a['end']+1, 'end_inclusive': b['start']-1,
                            'adjacent_transcription_exon_ranks': sorted([int(a['attributes']['rank']), int(b['attributes']['rank'])])})
    return {'exon_count': len(exons), 'CDS_length_bp': sum(c['end']-c['start']+1 for c in coding),
            'copied_CDS_bp_under_literal_tandem_model': sum(overlap(left,right,c['start'],c['end']) for c in coding),
            'complete_exon_ranks_inside_literal_copy_interval': [int(e['attributes']['rank']) for e in ranked if left <= e['start'] and e['end'] <= right],
            'partial_exon_ranks_inside_literal_copy_interval': [int(e['attributes']['rank']) for e in ranked if overlap(left,right,e['start'],e['end']) and not (left <= e['start'] and e['end'] <= right)],
            'whole_transcript_span_copied_under_literal_tandem_model': left <= transcript['start'] and transcript['end'] <= right,
            'left_boundary_introns': introns,
            'original_transcript_span_uninterrupted_in_simple_tandem_model': not (transcript['start'] <= right < transcript['end']),
            'limit': 'Geometry conditional on a simple direct tandem copy inserted after the right endpoint. No RNA, splice, expression, complete allele or loss/gain-of-function inference.'}


def transcript_context(gff, chrom, left, right, fasta):
    features = []
    with gzip.open(gff, 'rt') as stream:
        for line in stream:
            if line.startswith('#'): continue
            a = line.rstrip().split('\t')
            if len(a) != 9 or a[0] != chrom: continue
            features.append({'type':a[2], 'start':int(a[3]), 'end':int(a[4]), 'strand':a[6], 'phase':a[7],
                             'attributes':dict(p.split('=',1) for p in a[8].split(';') if '=' in p)})
    genes = [f for f in features if f['type']=='gene' and overlap(left,right,f['start'],f['end'])]
    ids = {g['attributes']['ID'] for g in genes}
    transcripts = [f for f in features if f['attributes'].get('Parent') in ids]
    children = defaultdict(list)
    for f in features:
        for parent in f['attributes'].get('Parent','').split(','): children[parent].append(f)
    rows=[]
    for t in transcripts:
        kids=children[t['attributes']['ID']]; exons=[k for k in kids if k['type']=='exon']; cds=[k for k in kids if k['type']=='CDS']
        if not exons: continue
        s=''.join(fasta.fetch(chrom,c['start']-1,c['end']).upper() for c in sorted(cds,key=lambda c:c['start']))
        if t['strand']=='-': s=reverse_complement(s)
        rows.append({'transcript':t, 'geometry':transcript_geometry(t,exons,cds,left,right),
                     'reference_CDS':{'length':len(s), 'sha256':hashlib.sha256(s.encode()).hexdigest(), 'starts_ATG':s.startswith('ATG'),
                                      'ends_stop':s[-3:] in ('TAA','TAG','TGA'), 'length_modulo_three':len(s)%3}, 'exons':exons, 'CDS':cds})
    return {'genes':genes, 'transcripts':rows, 'MANE_Select':[r for r in rows if 'MANE_Select' in r['transcript']['attributes'].get('tag','').split(',')],
            'limit':'Pinned GFF tag/version and reference sequence; no public transcript-sequence API comparison or patient expression assay.'}


def inspect_reads(bam, chrom, left, right, needles, config):
    import pysam
    records=[]; by_fragment=defaultdict(list); by_read=defaultdict(list); tags=defaultdict(set)
    with pysam.AlignmentFile(str(bam)) as reader:
        groups=reader.header.to_dict().get('RG',[])
        declared={g['ID'] for g in groups}
        for i,r in enumerate(reader):
            require(r.has_tag('RG') and r.get_tag('RG') in declared, 'Missing or undeclared read group')
            rg=r.get_tag('RG'); fragment=(rg,r.query_name); readkey=fragment+(1 if r.is_read1 else 2 if r.is_read2 else 0,)
            reasons=[]
            if r.flag & config['excluded_flags']: reasons.append('excluded_flag')
            if r.mapping_quality==255 or r.mapping_quality<config['minimum_mapping_quality']:reasons.append('mapping_quality')
            matches={name:exact_high_quality(r.query_sequence,r.query_qualities,needle,config['minimum_base_quality']) for name,needle in needles.items()}
            row={'index':i,'RG':rg,'query_name':r.query_name,'read_number':readkey[-1],'flag':r.flag,'chrom':r.reference_name,
                 'start_0based':r.reference_start,'end_exclusive':r.reference_end,'mapping_quality':r.mapping_quality,
                 'cigar':r.cigarstring,'query_span':query_span(r.cigartuples) if r.cigartuples else None,
                 'reverse':r.is_reverse,'supplementary':r.is_supplementary,'secondary':r.is_secondary,
                 'sequence_quality_matches':matches,'exclusions':reasons,'SA':r.get_tag('SA') if r.has_tag('SA') else None}
            records.append(row);by_fragment[fragment].append(row);by_read[readkey].append(row)
            if not reasons:
                for name,match in matches.items():
                    if match:tags[name].add(fragment)
    tolerance=config['split_coordinate_tolerance']; anchor=config['primary_read_anchor_length']
    split_details=[];pair_details=[]
    for readkey,rows in by_read.items():
        good=[r for r in rows if not r['exclusions'] and r['chrom']==chrom and r['query_span']]
        for a in good:
            for b in good:
                if a['index']==b['index'] or a['reverse']!=b['reverse']:continue
                qa,qb=a['query_span'],b['query_span']
                if (abs(a['end_exclusive']-right)<=tolerance and abs(b['start_0based']-(left-1))<=tolerance
                    and qa[2]==qb[2] and abs(qa[1]-qb[0])<=tolerance and qa[1]-qa[0]>=anchor and qb[1]-qb[0]>=anchor):
                    tags['two_observed_alignment_segments'].add(readkey[:2])
                    split_details.append({'read_key':readkey,'right_segment_index':a['index'],'left_segment_index':b['index'],'query_gap':qb[0]-qa[1]})
    for fragment,rows in by_fragment.items():
        primary=[r for r in rows if not r['supplementary'] and not r['secondary']]
        if len(primary)!=2 or sorted(r['read_number'] for r in primary)!=[1,2] or any(r['exclusions'] for r in primary):continue
        a,b=sorted(primary,key=lambda r:r['start_0based']);distance=config['paired_endpoint_distance']
        if (a['chrom']==b['chrom']==chrom and a['reverse'] and not b['reverse']
            and -tolerance<=a['start_0based']-(left-1)<=distance and -tolerance<=right-b['end_exclusive']<=distance):
            tags['outward_primary_pair'].add(fragment);pair_details.append({'fragment_key':fragment,'left_index':a['index'],'right_index':b['index']})
    qualified=set().union(*tags.values()) if tags else set()
    tag_summary={name:{'fragments':len(keys),'by_read_group':dict(Counter(k[0] for k in keys))} for name,keys in tags.items()}
    primary_alt=tags.get('junction_'+str(anchor),set()); split=tags.get('two_observed_alignment_segments',set()); pairs=tags.get('outward_primary_pair',set())
    return {'records':records,'split_details':split_details,'pair_details':pair_details,
            'fragment_evidence_tags':[{'RG':rg,'query_name':q,'tags':[name for name,keys in tags.items() if (rg,q) in keys]} for rg,q in sorted(qualified)],
            'summary':{'extracted_alignment_records':len(records),'fragments_in_windows':len(by_fragment),'read_groups':groups,
                       'records_by_read_group':dict(Counter(r['RG'] for r in records)), 'excluded_record_reasons':dict(Counter(reason for r in records for reason in r['exclusions'])),
                       'evidence':tag_summary,'primary_alt_fragment_union':len(primary_alt|split|pairs),
                       'primary_exact_junction_and_split':len(primary_alt&split),'primary_exact_junction_and_outward_pair':len(primary_alt&pairs),
                       'all_three':len(primary_alt&split&pairs),
                       'limit':'Unique RG/QNAME fragments, not additive evidence counts or calibrated allele fractions. Exact k-mer support can miss mismatching alleles; bounded extraction can omit distant/low-quality alternatives.'}}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('analysis','gff','config','samtools','output'):parser.add_argument('--'+name,type=Path,required=True)
    for name in ('summary','gff'):parser.add_argument('--expected-'+name+'-sha256',required=True)
    parser.add_argument('--record-id',required=True);args=parser.parse_args();tick=time.monotonic();os.nice(10)
    import pysam
    require(not args.output.exists() and not args.output.resolve().is_relative_to(Path(__file__).resolve().parents[2]), 'Use a new output outside source tree')
    args.output.mkdir(parents=True);state={'status':'running','pid':os.getpid(),'started_at':datetime.datetime.now(datetime.timezone.utc).isoformat()};save(args.output/'status.json',state)
    try:
        require(sha(args.analysis/'summary.json')==args.expected_summary_sha256,'Native summary differs')
        native=json.loads((args.analysis/'summary.json').read_text());verify(args.analysis/'provenance.json',native['provenance'])
        proof=json.loads((args.analysis/'provenance.json').read_text());require(native['analysis_status']=='complete_with_calls','Native run incomplete')
        for name in ('calls.bcf','calls.bcf.csi'):verify(args.analysis/name,native['outputs'][name])
        source=Path(proof['inputs']['alignment']['path']);reference=Path(proof['inputs']['reference']['path'])
        index_paths=[Path(proof['inputs'][name]['path']) for name in ('alignment_index','reference_index')]
        monitored=[source,reference]+index_paths
        before={str(p):(p.stat().st_size,p.stat().st_mtime_ns) for p in monitored}
        verify(source,proof['inputs']['alignment']);verify(reference,proof['inputs']['reference'])
        for name in ('alignment_index','reference_index'):verify(Path(proof['inputs'][name]['path']),proof['inputs'][name])
        require(sha(args.gff)==args.expected_gff_sha256,'GFF differs');config=json.loads(args.config.read_text())
        with pysam.VariantFile(str(args.analysis/'calls.bcf')) as reader:
            selected=[r.copy() for r in reader if r.id==args.record_id]
        require(len(selected)==1,'Record ID not unique');record=selected[0];info=dict(record.info)
        require(info.get('SVTYPE')=='DUP' and info.get('CT')=='5to3' and info.get('PRECISE') and info.get('INSLEN')==0,'Focused review requires precise 5to3 DUP without native insertion')
        with pysam.FastaFile(str(reference)) as fa:
            radius=config['reference_window_radius'];chrom=record.contig
            windows={name:(max(0,p-radius-1),fa.fetch(chrom,max(0,p-radius-1),p+radius).upper()) for name,p in [('left',record.pos),('right',record.stop)]}
            solutions=decompose_consensus(info['CONSENSUS'],windows,config['consensus_minimum_anchor']);save(args.output/'consensus-solutions.json',solutions)
            require(len(solutions)==1 and solutions[0]['order']==['right','left'] and solutions[0]['orientation']=='native','Literal junction ambiguous or different; preserve for manual review')
            sol=solutions[0];left=sol['suffix_first_reference_base_1based'];right=sol['prefix_last_reference_base_1based']
            needles={}
            for n in config['read_anchor_lengths']:
                needles['junction_'+str(n)]=fa.fetch(chrom,right-n,right).upper()+fa.fetch(chrom,left-1,left-1+n).upper()
            n=config['primary_read_anchor_length']
            for name,cut in [('reference_left',left-1),('reference_right',right)]:needles[name]=fa.fetch(chrom,cut-n,cut+n).upper()
            require(all('N' not in s and len(s)>=2*n for name,s in needles.items() if name.startswith('reference_')),'Reference k-mer not resolved')
            contexts=transcript_context(args.gff,chrom,left,right,fa)
        save(args.output/'transcript-context.json',contexts);save(args.output/'junction-needles.json',needles)
        save(args.output/'native-record.json',{'id':record.id,'chrom':record.contig,'POS':record.pos,'END':record.stop,'REF':record.ref,'ALT':record.alts,'FILTER':list(record.filter),'INFO':info,'samples':{s:dict(record.samples[s]) for s in record.samples}})
        radius=config['read_window_radius'];regions=[[chrom,max(0,p-radius-1),p+radius] for p in (left,right)]
        require(regions[0][2]<regions[1][1],'Focused windows overlap')
        bed=args.output/'regions.bed';bed.write_text(''.join('\t'.join(map(str,row))+'\n' for row in regions));bam=args.output/'junction-windows.bam'
        commands=[[str(args.samtools),'view','-@','1','--no-PG','-T',str(reference),'-b','-M','-L',str(bed),'-o',str(bam),str(source)],
                  [str(args.samtools),'quickcheck','-v',str(bam)],[str(args.samtools),'index',str(bam)]]
        for i,cmd in enumerate(commands):
            with (args.output/f'command-{i}.log').open('w') as log:subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
        save(args.output/'commands.json',commands)
        reads=inspect_reads(bam,chrom,left,right,needles,config);save(args.output/'read-evidence.json',reads)
        require(before=={str(p):(p.stat().st_size,p.stat().st_mtime_ns) for p in monitored},'Input changed during review')
        summary={'completed_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'elapsed_seconds':time.monotonic()-tick,
                 'inputs':{'native_summary':fingerprint(args.analysis/'summary.json'),'native_provenance':native['provenance'],'gff':fingerprint(args.gff),'config':fingerprint(args.config),'script':fingerprint(Path(__file__)),'samtools':fingerprint(args.samtools)},
                 'whole_alignment_and_reference_SHA256_rechecked':True,'pysam_version':pysam.__version__,'configuration':config,
                 'native_record_id':record.id,'native_POS':record.pos,'native_END':record.stop,'native_SVLEN':info.get('SVLEN'),'native_HOMLEN':info.get('HOMLEN'),
                 'exact_consensus_solution':sol,'consensus_length':len(info['CONSENSUS']),
                 'literal_simple_tandem_copy_interval_length':right-left+1,'length_note':'Native SVLEN and literal suffix-start through prefix-end inclusive differ in convention; native record is not rewritten into an asserted normalized allele.',
                 'reads':reads['summary'],'transcript_context':{'genes':contexts['genes'],'transcript_count':len(contexts['transcripts']),'MANE_Select':contexts['MANE_Select']},
                 'outputs':{p.name:fingerprint(p) for p in args.output.iterdir() if p.is_file() and p.name not in ('status.json','summary.json')},
                 'limitations':['Read signatures reuse the same specimen, reads and alignment; consensus is caller-derived. No independent biological validation.',
                                'One literal junction does not establish full haplotype, parent of origin, insertion site completeness, copy count, mosaic fraction or transcript expression.',
                                'Conditional simple tandem geometry is not a gene loss-of-function or clinical classification. Genome-wide junction uniqueness was not established.']}
        save(args.output/'summary.json',summary);state.update(status='complete',completed_at=summary['completed_at'],summary_sha256=sha(args.output/'summary.json'));save(args.output/'status.json',state)
        print(json.dumps({'elapsed_seconds':summary['elapsed_seconds'],'literal_junction':sol,'reads':reads['summary'],'MANE_transcripts':[r['transcript']['attributes']['transcript_id'] for r in contexts['MANE_Select']]}))
    except BaseException as error:
        state.update(status='failed',error_type=type(error).__name__,error=str(error));save(args.output/'status.json',state);raise


if __name__=='__main__':main()
