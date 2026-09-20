#!/usr/bin/env python3
"""Check synthetic novel junctions and full paired FASTQ structure independently."""
import argparse
from collections import Counter
import gzip
import json
from pathlib import Path
import time

import pysam

from simulate_sv_reads import reverse_complement, sha


def audit(fixture, output):
    started = time.monotonic(); truth = json.loads((fixture/'truth.json').read_text())
    if not truth['scope'].startswith('Synthetic diploid simple-SV reference;'):
        raise ValueError('this audit is restricted to the synthetic fixture')
    for name, identity in truth['files'].items():
        p = fixture/name
        if p.stat().st_size != identity['bytes'] or sha(p) != identity['sha256']:
            raise ValueError('fixture file differs from saved identity')
    def fasta(path):
        with pysam.FastxFile(str(path)) as source: return {r.name:r.sequence for r in source}
    reference = fasta(fixture/'reference.fa'); alternate = fasta(fixture/'alternate-haplotype.fa')
    flank = 31; queries = []
    for event in truth['events_on_one_haplotype']:
        seq = reference[event['chrom']]; left = event['left_cut_0based']; kind = event['type']
        if kind == 'BND':
            partner = reference[event['partner_chrom']]; right = event['partner_cut_0based']
            parts = [('first_connection',seq[left-flank:left]+partner[right:right+flank]),
                     ('reciprocal_connection',partner[right-flank:right]+seq[left:left+flank])]
        else:
            right = event['right_cut_0based']
            if kind == 'DEL': parts = [('deletion_join',seq[left-flank:left]+seq[right:right+flank])]
            elif kind == 'DUP': parts = [('tandem_join',seq[right-flank:right]+seq[left:left+flank])]
            elif kind == 'INV': parts = [('left_join',seq[left-flank:left]+reverse_complement(seq[right-flank:right])),
                                        ('right_join',reverse_complement(seq[left:left+flank])+seq[right:right+flank])]
            elif kind == 'INS': parts = [('left_join',seq[left-flank:left]+event['inserted_sequence'][:flank]),
                                        ('right_join',event['inserted_sequence'][-flank:]+seq[left:left+flank])]
            else: raise ValueError('unexpected synthetic event class')
        for label, query in parts:
            rc = reverse_complement(query)
            locations = [{'chrom':c,'start_0based':sequence.find(query)} for c,sequence in alternate.items() if query in sequence]
            occurrences = sum(sequence.count(query)+sequence.count(rc) for sequence in alternate.values())
            ref_count = sum(sequence.count(query)+sequence.count(rc) for sequence in reference.values())
            if occurrences != 1 or ref_count != 0: raise AssertionError('novel junction construction is not uniquely supported')
            queries.append({'event_id':event['id'],'junction':label,'sequence_sha256':__import__('hashlib').sha256(query.encode()).hexdigest(),
                            'alternate_occurrences_both_strands':occurrences,'reference_occurrences_both_strands':ref_count,
                            'alternate_forward_location':locations,'query':query,'reverse_query':rc})
    if alternate[truth['unaltered_control_chromosome']] != reference[truth['unaltered_control_chromosome']]:
        raise AssertionError('control chromosome changed')
    counts=Counter(); junction_reads=Counter()
    with gzip.open(fixture/'SYNTH_SV_R1.fastq.gz','rt') as r1, gzip.open(fixture/'SYNTH_SV_R2.fastq.gz','rt') as r2:
        while True:
            h1,h2=r1.readline(),r2.readline()
            if not h1 and not h2:break
            if not h1 or not h2:raise ValueError('FASTQ mates have different lengths')
            one=[r1.readline().rstrip('\n') for _ in range(3)];two=[r2.readline().rstrip('\n') for _ in range(3)]
            if h1.rstrip()[:-2] != h2.rstrip()[:-2] or not h1.rstrip().endswith('/1') or not h2.rstrip().endswith('/2'):
                raise ValueError('paired FASTQ names do not match')
            for seq,plus,quality in (one,two):
                if len(seq)!=truth['read_length'] or plus!='+' or quality!='?'*truth['read_length'] or set(seq)-set('ACGT'):
                    raise ValueError('malformed synthetic read or quality data')
                for query in queries:
                    if query['query'] in seq or query['reverse_query'] in seq:
                        junction_reads[query['event_id']+'/'+query['junction']]+=1
                counts['read_records']+=1;counts['bases']+=len(seq)
            counts['read_pairs']+=1
    if counts['read_records']!=truth['read_records'] or counts['read_pairs']!=truth['read_pairs'] or counts['bases']!=truth['read_bases']:
        raise AssertionError('full FASTQ census differs from claimed counts')
    for query in queries:
        query['exact_junction_read_records']=junction_reads[query['event_id']+'/'+query['junction']]
        del query['query'];del query['reverse_query']
    result={'elapsed_seconds':time.monotonic()-started,'truth_sha256':sha(fixture/'truth.json'),
        'audit_source_sha256':sha(Path(__file__)),'junction_flank_bp':flank,'junctions':queries,'full_FASTQ_census':dict(counts),
        'control_chromosome_sequence_unchanged':True,
        'scope':'Independent junction-sequence and full paired-FASTQ checks of a synthetic fixture. Exact junction read counts are not independent molecules, mapped support or a clinical detection threshold.'}
    output.write_text(json.dumps(result,indent=2)+'\n');return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--fixture',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();print(json.dumps(audit(a.fixture,a.output),indent=2))
