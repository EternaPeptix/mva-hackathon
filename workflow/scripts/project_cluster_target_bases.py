#!/usr/bin/env python3
"""Check whether each original target base participates in a public-sequence match."""
import argparse
from collections import Counter
import datetime
import json
from pathlib import Path

import pysam

from compare_cluster_public_sequence import parse_embl
from review_cnv_context import fingerprint, require, save, sha


def reoriented_query_position(position, length, old_reverse, new_reverse):
    require(0 <= position < length,'Query position outside read')
    return length-1-position if bool(old_reverse) != bool(new_reverse) else position


def complement(base):
    return base.translate(str.maketrans('ACGTN','TGCAN'))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['root','comparison','output']:
        parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args();root=args.root.resolve();run=args.comparison.resolve();output=args.output.resolve()
    require(not output.exists(),'New output file required')
    comparison=json.loads((run/'summary.json').read_text());config=comparison['configuration']
    parent=root/config['cluster_directory'];cluster=json.loads((parent/'summary.json').read_text())
    require(sha(parent/'summary.json')==config['cluster_summary_sha256'],'Cluster differs')
    for name in ['read-review.json','queries.fa']:
        require(sha(parent/name)==cluster['outputs'][name]['sha256'],'Original evidence differs')
    require(sha(run/'reads-to-clone.sam')==comparison['outputs']['reads-to-clone.sam']['sha256'],'Comparison differs')
    public=root/config['public_embl'];require(sha(public)==config['public_embl_sha256'],'Public sequence differs')
    accession,sequence=parse_embl(public.read_text())
    review=json.loads((parent/'read-review.json').read_text())
    original={}
    target=root/cluster['configuration']['target_directory']/'targets.bam'
    target_state=json.loads((target.parent/'status.json').read_text())
    require(sha(target)==target_state['output_sha256']['targets.bam'],'Original target BAM differs')
    with pysam.AlignmentFile(str(target),'rb') as bam:
        for i,read in enumerate(bam.fetch(cluster['configuration']['chrom'],*cluster['configuration']['read_window_0based'])):
            identifier='read_'+str(i);record=review['observations'][identifier]
            require((read.query_name,read.get_tag('RG'),read.reference_start,read.flag)==
                    (record['query_name'],record['RG'],record['start_0based'],record['flag']),'Original record order differs')
            original[identifier]=read
    remapped={}
    with pysam.AlignmentFile(str(run/'reads-to-clone.sam'),'r') as sam:
        for read in sam:
            if not read.is_secondary and not read.is_supplementary:remapped[read.query_name]=read
    projections=[];summaries=[]
    for index,site in enumerate(cluster['configuration']['sites']):
        groups={}
        configuration=review['configurations']['20'][index]
        for label in ['0','1']:
            fragments=[f for f in configuration['fragments'] if f['label']==label];fragment_matches=0;states=Counter()
            for fragment in fragments:
                matched=False
                for identifier in fragment['observation_ids']:
                    old=original[identifier];new=remapped[identifier]
                    q=next(q for q,r in old.get_aligned_pairs() if r==site[0]-1 and q is not None)
                    n=reoriented_query_position(q,old.query_length,old.is_reverse,new.is_reverse)
                    pos=None if new.is_unmapped else dict(new.get_aligned_pairs()).get(n)
                    base=old.query_sequence[q].upper()
                    oriented=complement(base) if old.is_reverse != new.is_reverse else base
                    exact=pos is not None and new.reference_name==accession and sequence[pos]==oriented
                    status='mapped_target_matches_public_base' if exact else 'target_unmapped_or_clipped' if pos is None else 'mapped_target_differs'
                    matched |= exact;states[status]+=1
                    projections.append({'site':site,'group':label,'observation_id':identifier,'original_query_position_0based':q,
                                        'comparison_query_position_0based':n,'public_position_1based':pos+1 if pos is not None else None,
                                        'original_base':base,'public_oriented_base':oriented,'public_reference_base':sequence[pos] if pos is not None else None,
                                        'status':status,'original_reverse':old.is_reverse,'comparison_reverse':new.is_reverse})
                fragment_matches+=matched
            groups[label]={'fragments':len(fragments),'fragments_with_any_target_matched':fragment_matches,'read_target_statuses':dict(states)}
        summaries.append({'site':site,'groups':groups})
    save(output,{'completed_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'public_accession':accession,
                 'inputs':{n:fingerprint(p) for n,p in {'comparison_summary':run/'summary.json','cluster_summary':parent/'summary.json',
                           'original_bam':target,'remapped_sam':run/'reads-to-clone.sam','public_record':public,'script':Path(__file__),
                           'public_sequence_helper':Path(__file__).with_name('compare_cluster_public_sequence.py'),
                           'io_helper':Path(__file__).with_name('review_cnv_context.py')}.items()},
                 'summary':summaries,'projections':projections,
                 'limits':['Checks the original selected base itself, not only another matching segment of its read.',
                           'A public transcript match supports sequence similarity, not genomic locus, gene conversion, expression, physical phase or a patient assembly.']})
    print(json.dumps({'output':str(output),'summary':summaries,'sha256':sha(output)}))


if __name__=='__main__':main()
