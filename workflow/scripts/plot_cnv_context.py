#!/usr/bin/env python3
"""Plot selected internal CNV depth/read-balance contexts; no clinical inference.

Patient-derived figures and point tables must remain outside the source export.
"""
import argparse
import datetime
import gzip
import json
import os
from pathlib import Path

os.environ.setdefault('MPLBACKEND', 'Agg')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')

from review_cnv_context import allele_fraction, fingerprint, require, sha, verify, save


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('analysis', 'context', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--expected-context-sha256', required=True)
    parser.add_argument('--call', action='append', required=True, help='Call ID=display label')
    args = parser.parse_args(); os.nice(10)
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from cnvpytor.io import IO, FLAG_GC_CORR
    require(not args.output.exists(), 'Existing figures preserved')
    require(sha(args.context / 'summary.json') == args.expected_context_sha256, 'Context summary differs')
    context = json.loads((args.context / 'summary.json').read_text())
    for name in ('all-call-context.json.gz', 'all-marker-context.jsonl.gz'):
        verify(args.context / name, context['outputs'][name])
    with gzip.open(args.context / 'all-call-context.json.gz', 'rt') as f:
        calls = {r['call_id']: r for r in json.load(f)['calls']}
    review = json.loads(Path(context['inputs']['cnv_review_summary']['path']).read_text())
    require(sha(args.analysis / 'summary.json') == review['inputs']['native_summary']['sha256'], 'Native CNV summary differs')
    native = json.loads((args.analysis / 'summary.json').read_text())
    verify(args.analysis / 'depth.pytor', native['outputs']['depth.pytor'])
    markers = {}
    with gzip.open(args.context / 'all-marker-context.jsonl.gz', 'rt') as f:
        for line in f:
            row = json.loads(line); markers.setdefault(row['chrom'], []).append(row)
    chosen = []
    for value in args.call:
        ident, label = value.split('=', 1)
        require(ident in calls and label, 'Unknown figure call')
        chosen.append((calls[ident], label))
    require(len({r['call_id'] for r, _ in chosen}) == len(chosen), 'Repeated figure call')
    args.output.mkdir(parents=True)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.titlesize': 10,
                         'svg.hashsalt': 'mva-cnv-context-20260920'})
    fig, axes = plt.subplots(len(chosen), 2, figsize=(13, 2.45 * len(chosen) + 1.1), squeeze=False)
    points = []
    io = IO(str(args.analysis / 'depth.pytor'), ro=True, create=False)
    try:
        for row_index, (call, label) in enumerate(chosen):
            chrom = call['chrom']; start = call['flanks']['left']['start_1based']; end = call['flanks']['right']['end_inclusive']
            a, b = call['start_1based'], call['end_inclusive']
            size = 1000
            mean = native['normalization'][str(size)]['fitted_mean_read_midpoints_per_bin']
            values = np.asarray(io.get_signal(chrom, size, 'RD', FLAG_GC_CORR), dtype=float) / mean
            lo, hi = (start-1)//size, min(len(values), (end-1)//size+1)
            x = (np.arange(lo, hi) + .5) * size
            y = values[lo:hi]
            require(np.all(np.isfinite(y)), 'Nonfinite plotted depth')
            left, right = axes[row_index]
            left.scatter(x/1e6, y, s=7, color='#225E87', alpha=.8, linewidths=0)
            left.axhline(1, color='#555555', lw=.8, ls='--')
            left.set_ylim(0, max(2.15, float(y.max())*1.06))
            left.set_ylabel('Relative read density\n(1 kb GC-corrected bins)')
            left.set_title(f'{label} | {call["type"]} | {call["call_id"]}', loc='left', fontweight='bold')
            selected = [r for r in markers[chrom] if start <= r['pos'] <= end]
            inside = [r for r in selected if a <= r['pos'] <= b]
            plotted = []
            for r in selected:
                fraction = allele_fraction(r['AD'])
                enough = fraction is not None and sum(r['AD']) >= context['config']['marker_minimum_allelic_depth']
                if enough:
                    right.scatter(r['pos']/1e6, fraction, s=25 if r['usable'] else 42,
                                  marker='o' if r['usable'] else 'x', color='#225E87' if r['usable'] else '#B23A34', linewidths=.9)
                    plotted.append(dict(r, ALT_fraction=fraction, figure_status='usable' if r['usable'] else 'previously_excluded'))
                else:
                    right.scatter(r['pos']/1e6, -.075, marker='v', s=24, color='#888888')
                    plotted.append(dict(r, ALT_fraction=fraction, figure_status='missing_or_below_AD_minimum'))
            for fraction in (1/3, .5, 2/3):
                right.axhline(fraction, color='#AAAAAA', lw=.6, ls='--' if fraction == .5 else ':', zorder=0)
            right.set_ylim(-.12, 1.06)
            right.set_yticks([0, 1/3, .5, 2/3, 1], ['0', '1/3', '1/2', '2/3', '1'])
            right.set_ylabel('ALT / (REF + ALT reads)')
            hets = sum(r['category'] == 'heterozygous' for r in inside)
            right.set_title(f'{len(inside)} panel markers inside; {hets} native diploid heterozygotes', loc='left')
            if not inside:
                right.text((a+b)/2/1e6, .8, 'No panel marker\ninside call', ha='center', va='center', fontsize=9,
                           bbox={'facecolor':'white', 'alpha':.8, 'edgecolor':'none'})
            for axis in (left, right):
                axis.axvspan(a/1e6, b/1e6, color='#DAA245', alpha=.18, zorder=-1)
                axis.set_xlim(start/1e6, end/1e6)
                axis.set_xlabel(f'Chromosome {chrom}, GRCh38 position (Mb)')
                axis.ticklabel_format(axis='x', style='plain', useOffset=False)
                axis.spines[['top', 'right']].set_visible(False)
            points.append({'call_id': call['call_id'], 'label': label, 'chrom': chrom,
                'plot_bounds_1based_inclusive': [start,end], 'call_bounds_1based_inclusive':[a,b],
                'bin_midpoints':x.tolist(), 'relative_read_density':y.tolist(), 'marker_points':plotted})
    finally:
        io.file.close()
    fig.suptitle('Selected copy-number candidates: depth and common-marker read balance', fontsize=15, y=.994, va='top')
    fig.text(.5, .965, 'Internal research figure. Same specimen and reads; neither panel establishes allele structure or clinical significance.',
             ha='center', va='top', fontsize=9)
    handles = [Patch(facecolor='#DAA245',alpha=.25,label='Nominal CNV interval'),
        Line2D([],[],ls='',marker='o',color='#225E87',label='Marker retained for prior ROH analysis'),
        Line2D([],[],ls='',marker='x',color='#B23A34',label='Marker excluded from prior ROH analysis'),
        Line2D([],[],ls='',marker='v',color='#888888',label='No usable AD or total AD < 10')]
    fig.legend(handles=handles,loc='lower center',ncol=2,frameon=False,bbox_to_anchor=(.5,.014),fontsize=9)
    fig.tight_layout(rect=(0,.065,1,.955), h_pad=1.5, w_pad=2.5)
    for suffix in ('png','svg'):
        fig.savefig(args.output / ('selected-context.' + suffix), dpi=180)
    plt.close(fig)
    save(args.output / 'plotted-points.json', {'points':points,
        'selection':'Five descriptive examples including sparse/zero-marker regions, not a causal ranking.',
        'limits':['Read-midpoint depth includes low-MQ reads as in the native method.',
                  'Fractions at common ascertained SNPs use native AD; not a copy-state, molecular phase or mosaic-fraction estimate.',
                  'The 1/3, 1/2 and 2/3 guides are visual references, not a fitted hypothesis test.']})
    summary = {'completed_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'context_summary':fingerprint(args.context / 'summary.json'), 'native_summary':fingerprint(args.analysis / 'summary.json'),
        'script':fingerprint(Path(__file__)), 'calls':[r['call_id'] for r, _ in chosen],
        'matplotlib':matplotlib.__version__, 'numpy':np.__version__,
        'outputs':{name:fingerprint(args.output/name) for name in ('selected-context.png','selected-context.svg','plotted-points.json')},
        'publication_scope':'Patient-derived internal artifacts. Exclude all figure outputs and point tables from source repository export.'}
    save(args.output / 'summary.json', summary)
    print(json.dumps({'status':'complete','calls':summary['calls'],'outputs':list(summary['outputs'])}))


if __name__ == '__main__':
    main()
