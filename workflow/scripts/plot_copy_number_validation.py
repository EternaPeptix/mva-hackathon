#!/usr/bin/env python3
"""Plot retained synthetic CNVpytor profiles and truth; never patient profiles."""
import argparse
import hashlib
import json
import os
from pathlib import Path

os.environ.setdefault('MPLBACKEND', 'Agg')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pysam
from cnvpytor.io import IO, FLAG_GC_CORR


def sha256(path):
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            value.update(chunk)
    return value.hexdigest()


def load_fixture(root):
    summary_path = root / 'analysis/summary.json'
    validation_path = root / 'synthetic-validation.json'
    truth_path = root / 'fixture/truth.json'
    summary = json.loads(summary_path.read_text())
    validation = json.loads(validation_path.read_text())
    truth = json.loads(truth_path.read_text())
    if summary['config']['scope'] != 'synthetic' or summary['config']['sample'] != 'SYNTH':
        raise ValueError('only explicitly synthetic profiles may be plotted')
    if sha256(summary_path) != validation['result_summary_sha256']:
        raise ValueError('analysis summary differs from validation receipt')
    depth = root / 'analysis/depth.pytor'
    if sha256(depth) != summary['outputs']['depth.pytor']['sha256']:
        raise ValueError('profile differs from completed analysis')
    reference = root / 'fixture/synthetic.fa'
    if sha256(reference) != summary['inputs']['reference']['sha256']:
        raise ValueError('synthetic reference differs from analysis')
    with pysam.FastaFile(str(reference)) as fasta:
        lengths = dict(zip(fasta.references, fasta.lengths))
    inputs = {str(path): sha256(path) for path in
              (summary_path, validation_path, truth_path, depth, reference)}
    return summary, truth, lengths, depth, inputs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fine-root', required=True, type=Path)
    parser.add_argument('--coarse-root', required=True, type=Path)
    parser.add_argument('--output-stem', required=True, type=Path)
    args = parser.parse_args()
    outputs = [args.output_stem.with_suffix(suffix) for suffix in ('.png', '.svg', '.json')]
    if any(path.exists() for path in outputs):
        raise FileExistsError('choose a new output stem')
    outputs[0].parent.mkdir(parents=True, exist_ok=True)
    fine = load_fixture(args.fine_root.resolve())
    coarse = load_fixture(args.coarse_root.resolve())
    plt.rcParams.update({'font.size': 11, 'axes.titlesize': 12, 'axes.labelsize': 11,
                         'font.family': 'DejaVu Sans', 'svg.hashsalt': 'mva-cnv-synthetic-v3'})
    fig, axes = plt.subplots(3, 2, figsize=(13, 9.6))
    rows = [(fine, 1000), (fine, 10000), (coarse, 100000)]
    for row, (fixture, bin_size) in enumerate(rows):
        summary, truth, lengths, depth, _ = fixture
        if bin_size not in summary['config']['bin_sizes']:
            raise ValueError('requested resolution absent from completed analysis')
        reader = IO(str(depth), ro=True, create=False)
        try:
            mean = summary['normalization'][str(bin_size)]['fitted_mean_read_midpoints_per_bin']
            for col, chrom in enumerate(('1', '2')):
                axis = axes[row, col]
                y = np.asarray(reader.get_signal(chrom, bin_size, 'RD', FLAG_GC_CORR)) / mean
                x = (np.arange(len(y)) + 0.5) * bin_size
                visible = x < lengths[chrom]
                if not np.all(np.isfinite(y[visible])):
                    raise ValueError('nonfinite plotted density')
                axis.plot(x[visible] / 1e6, y[visible], color='#266a9d', lw=0.7,
                          marker='.', markersize=2.3, alpha=0.9)
                expected_x = [0.0]
                expected_y = [1.0]
                for change in sorted((r for r in truth['truth'] if r['chrom'] == chrom),
                                     key=lambda r: r['start_0based']):
                    a, b = change['start_0based'] / 1e6, change['end_exclusive'] / 1e6
                    expected_x.extend([a, a, b, b])
                    expected_y.extend([1.0, change['relative_depth'], change['relative_depth'], 1.0])
                    axis.axvspan(a, b, color='#c75560' if change['type'] == 'deletion' else '#3c957c',
                                 alpha=0.13, lw=0)
                expected_x.append(lengths[chrom] / 1e6); expected_y.append(1.0)
                axis.plot(expected_x, expected_y, color='#222222', linestyle='--', lw=1.1)
                title = 'Simulated loss and gain' if chrom == '1' else 'Unchanged control'
                axis.set_title(f'{title} · {bin_size // 1000:,} kb windows', loc='left', pad=9)
                axis.set_xlim(0, lengths[chrom] / 1e6)
                axis.set_ylim(0, max(2.0, float(y[visible].max()) * 1.05))
                axis.set_yticks([0, 0.5, 1, 1.5, 2])
                axis.set_xlabel('Position on artificial chromosome (Mb)')
                if col == 0:
                    axis.set_ylabel('Relative read density')
                axis.grid(axis='y', color='#e4e8eb', linewidth=0.6)
                axis.spines[['top', 'right']].set_visible(False)
        finally:
            reader.file.close()
    fig.suptitle('Copy-number workflow: synthetic signal recovery', x=0.08, y=0.975,
                 ha='left', fontsize=19, weight='bold')
    fig.text(0.08, 0.939,
             'Simulation only · 2 Mb chromosomes at 1/10 kb windows; 10 Mb chromosomes at 100 kb windows',
             fontsize=11, color='#444444')
    handles = [Line2D([], [], color='#266a9d', label='Measured GC-corrected density'),
               Line2D([], [], color='#222222', ls='--', label='Expected density'),
               Patch(facecolor='#c75560', alpha=0.25, label='Simulated loss'),
               Patch(facecolor='#3c957c', alpha=0.25, label='Simulated gain')]
    fig.legend(handles=handles, loc='lower left', bbox_to_anchor=(0.073, 0.04),
               ncol=4, frameon=False, fontsize=10)
    fig.text(0.08, 0.021,
             'Directly mapped read-density fixtures. No realignment or realistic mapping bias; not clinical detection limits.',
             fontsize=10, color='#444444')
    fig.subplots_adjust(left=0.08, right=0.98, top=0.89, bottom=0.135, hspace=0.53, wspace=0.2)
    fig.savefig(outputs[0], dpi=170, facecolor='white')
    fig.savefig(outputs[1], facecolor='white', metadata={'Date': None})
    plt.close(fig)
    receipt = {'scope': 'synthetic only; no patient profile',
               'script_sha256': sha256(Path(__file__)), 'matplotlib': matplotlib.__version__,
               'inputs': dict(fine[4], **coarse[4]),
               'outputs': {str(path.resolve()): sha256(path) for path in outputs[:2]}}
    outputs[2].write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'outputs': [str(p) for p in outputs], 'scope': receipt['scope']}))


if __name__ == '__main__':
    main()
