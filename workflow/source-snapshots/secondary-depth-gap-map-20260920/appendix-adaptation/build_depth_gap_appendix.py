#!/usr/bin/env python3
"""Readable derivative of completed gap-map outputs; no original depths or reads."""
from collections import defaultdict
import argparse
import os
import gzip
import hashlib
import json
from pathlib import Path

MODES=('MQ20_BQ20','MQ60_BQ20')


def pin(path):
    data=path.read_bytes()
    return dict(path=str(path),bytes=len(data),sha256=hashlib.sha256(data).hexdigest())


def merge(values):
    out=[]
    for row,ordinal in values:
        if out and out[-1]['chrom']==row['chrom'] and out[-1]['end']==row['start']:
            out[-1]['end']=row['end']
            out[-1]['atomic_row_numbers'].append(ordinal)
            out[-1]['transcript_ids_overlapping_any_part']=sorted(set(out[-1]['transcript_ids_overlapping_any_part'])|set(row['transcript_ids']))
        else:
            assert not out or out[-1]['chrom']!=row['chrom'] or row['start']>=out[-1]['end']
            out.append(dict(chrom=row['chrom'],start=row['start'],end=row['end'],atomic_row_numbers=[ordinal],
                            transcript_ids_overlapping_any_part=list(row['transcript_ids'])))
    return out


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--expected-summary-sha256',required=True)
    parser.add_argument('--result-link',default='RESULT.md')
    args=parser.parse_args()
    source=args.input_dir.resolve();output_dir=args.output_dir.resolve()
    assert not output_dir.exists(),'existing output directory preserved'
    summary_path=source/'summary.json';interval_path=source/'gap-intervals.jsonl.gz'
    assert pin(summary_path)['sha256']==args.expected_summary_sha256
    summary=json.loads(summary_path.read_text())
    status=json.loads((source/'status.json').read_text())
    assert pin(interval_path)['sha256']==status['outputs']['gap-intervals.jsonl.gz']['sha256']
    with gzip.open(interval_path,'rt') as f:rows=[json.loads(line) for line in f]
    assert len(rows)==summary['gap_intervals']
    selected=defaultdict(list)
    for ordinal,row in enumerate(rows,1):
        for mode,category in zip(row['depth_modes'],row['categories']):
            for criterion in ('below20','zero'):
                if (criterion=='zero' and category=='zero') or (criterion=='below20' and category!='at-least-20'):
                    for key in ['panel_union']+row['gene_ids']+row['transcript_ids']:
                        selected[(key,mode,criterion)].append((row,ordinal))
    def with_locations(key,metadata):
        out=dict(metadata)
        out['intervals']={mode:{criterion:merge(selected[(key,mode,criterion)]) for criterion in ('below20','zero')} for mode in MODES}
        for mode in MODES:
            for criterion in ('below20','zero'):
                assert sum(r['end']-r['start'] for r in out['intervals'][mode][criterion])==metadata['modes'][mode][criterion+'_depth_bases' if criterion=='zero' else 'below20_bases']
        return out
    genes=[with_locations(g['gene_id'],g) for g in summary['genes']]
    transcripts=[with_locations(t['transcript_id'],{k:v for k,v in t.items() if k!='cds_rows'}) for t in summary['transcripts']]
    panel=with_locations('panel_union',{'modes':summary['panel_union']})
    value=dict(status='complete',scope='Readable derivative of saved map outputs only; not an independent scientific audit.',
               coordinates='0-based half-open',inputs=[pin(summary_path),pin(interval_path)],
               interval_membership_note='Merged intervals preserve the union of transcript IDs overlapping any part; consult atomic_row_numbers in the original JSONL for exact changing membership. Zero intervals are a subset of below20, not an additive category.',
               panel=panel,genes=genes,transcripts=transcripts,paired_categories=summary['paired_categories'],
               strict_depth_numeric_increase_bases=summary['MQ60_depth_greater_than_MQ20_bases'])
    output_dir.mkdir(parents=True,exist_ok=False)
    output=output_dir/'gap-locations.json';assert not output.exists();output.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
    md=['# Complete low-depth location appendix','',
        'All coordinates below are **GRCh38, 0-based half-open**. Rows merge adjacent gene-union positions under one criterion; they do not bridge an untargeted interval. Depth <20 includes zero. The two modes reuse the same reads. These are saved aligned-evidence locations, not clinical callability or variant findings.','',
        f'[All {len(summary["transcripts"]):,} transcript summaries and interval memberships](gap-locations.json) · [Original {len(rows):,} membership/category rows]({os.path.relpath(interval_path,output_dir)}) · [Main result]({args.result_link})','',
        'Transcript sets in merged intervals mean overlap with **any part**, not coverage of the complete interval by every listed transcript. Atomic row numbers preserve exact membership changes. The original summary retains all source CDS rows and ordinals.','',
        '| Gene | Mode | Criterion | Chr | Start (included) | End (excluded) | Bases | Atomic row numbers |',
        '|---|---|---|---|---:|---:|---:|---|']
    for g in genes:
        for mode in MODES:
            for criterion in ('below20','zero'):
                for r in g['intervals'][mode][criterion]:
                    md.append(f"| {g['name']} | {mode} | {criterion} | {r['chrom']} | {r['start']} | {r['end']} | {r['end']-r['start']} | {', '.join(map(str,r['atomic_row_numbers']))} |")
    out=output_dir/'GAP-LOCATIONS.md';assert not out.exists();out.write_text('\n'.join(md)+'\n')
    compact=dict(inputs=value['inputs'],derived_json=pin(output),appendix=pin(out),
                 affected_genes={m:sum(g['modes'][m]['below20_bases']>0 for g in genes) for m in MODES},
                 affected_transcripts={m:sum(t['modes'][m]['below20_bases']>0 for t in transcripts) for m in MODES},
                 gene_interval_counts={m:{c:sum(len(g['intervals'][m][c]) for g in genes) for c in ('below20','zero')} for m in MODES})
    out=output_dir/'location-derivation.json';assert not out.exists();out.write_text(json.dumps(compact,indent=2,sort_keys=True)+'\n')
    print(json.dumps(compact,indent=2))


if __name__=='__main__':main()
