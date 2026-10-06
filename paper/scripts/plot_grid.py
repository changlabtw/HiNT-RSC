"""Accuracy vs number of CNV calls: unsmoothed lambda sweep vs mean-filter window sweeps at several lambdas.

Figure 3 (K562, WGS BIC-seq2 reference), Figure 6 (PANC-1, DepMap WGS and array-state references) and
Supplementary Figure S5 (Caki2, DepMap WGS reference).

Run from paper/, e.g.:
  python scripts/plot_grid.py --cell k562                     # expected/k562_MboI.tsv -> output/figures/figure3_k562_grid.png
  python scripts/plot_grid.py --cell panc1 --results <PANC-1 grid TSV> -o output/figures/figure6_panc1_grid.png
  python scripts/plot_grid.py --cell caki2 --results <Caki2 grid TSV> -o output/figures/figureS5_caki2_grid.png
A dagger marks segmentations in which BIC-seq2 dropped at least one chromosome; this is read from the segmentation
files ('seg' column, relative to --seg-root), so the daggers are drawn only when those files are present.
"""
import argparse
import os
import sys

import matplotlib
import pandas as pd

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

LAM_COLORS = [(3, 'tab:purple'), (2, 'tab:orange'), (1, 'tab:green'), (0.5, 'tab:red')]
CELLS = {
    'k562': [('K562, WGS BIC-seq2 reference', 'SCC_bicseq2_l4', 'iou50_hits_trg', 'iou50_P', 'iou50_calls')],
    'panc1': [('PANC-1, DepMap WGS reference', 'depmap_SCC', 'depmap_hits', 'depmap_P', 'depmap_calls'),
              ('PANC-1, array-state reference', 'array_SCC', 'arraystate_hits', 'arraystate_P', 'arraystate_calls')],
    'caki2': [('Caki2, DepMap WGS reference', 'depmap_SCC', 'depmap_hits', 'depmap_P', 'depmap_calls')],
}
DEFAULT_OUT = {'k562': 'output/figures/figure3_k562_grid.png', 'panc1': 'output/figures/figure6_panc1_grid.png',
               'caki2': 'output/figures/figureS5_caki2_grid.png'}
N_CHROM = 23  # chr1-22, chrX


def seg_path(p, seg_root):
    return p if os.path.isabs(p) or os.path.exists(p) else os.path.join(seg_root, p)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--cell', choices=list(CELLS), required=True)
    ap.add_argument('--results', help='results TSV; default for k562: expected/k562_MboI.tsv')
    ap.add_argument('--seg-root', default='output', help="directory the 'seg' column is relative to; default: %(default)s")
    ap.add_argument('-o', '--output', help='default: ' + ', '.join(f'{k}: {v}' for k, v in DEFAULT_OUT.items()))
    args = ap.parse_args()
    if args.results is None:
        if args.cell != 'k562':
            ap.error('--results is required for --cell panc1/caki2')
        args.results = 'expected/k562_MboI.tsv'
    args.output = args.output or DEFAULT_OUT[args.cell]
    os.makedirs(os.path.dirname(args.output) or '.', exist_ok=True)
    d = pd.read_csv(args.results, sep='\t')
    d = d[(d.setting == 'complete') & d.method.isin(['none', 'mean'])]
    if 'segmenter' in d.columns:
        d = d[d.segmenter == 'bicseq2']
    d['lambda'] = d['lambda'].astype(float)
    paths = [seg_path(p, args.seg_root) for p in d['seg']]
    if all(os.path.exists(p) for p in paths):
        d['incomplete'] = [len({l.split('\t')[0] for l in open(p).readlines()[1:]}) < N_CHROM for p in paths]
    else:
        print(f'warning: segmentation files not found under {args.seg_root}/; chromosome-absence daggers are not drawn', file=sys.stderr)
        d['incomplete'] = False
    rows = CELLS[args.cell]
    fig, axes = plt.subplots(len(rows), 3, figsize=(15, 4.6 * len(rows) + 1.2), squeeze=False)
    for ri, (title, scc, hits, prec, calls) in enumerate(rows):
        un = d[d.method == 'none'].sort_values(calls)
        for ax, (ycol, yname) in zip(axes[ri], [(scc, 'segment-level SCC'), (hits, '# matched calls'), (prec, 'precision')]):
            ax.plot(un[calls], un[ycol], 'o-', color='tab:blue')
            for _, r in un.iterrows():
                ax.annotate(f"{r['lambda']:g}", (r[calls], r[ycol]), fontsize=7, color='tab:blue', xytext=(2, -8), textcoords='offset points')
            for lam, col in LAM_COLORS:
                s = d[(d.method == 'mean') & (d['lambda'] == lam)].sort_values('param')
                ax.plot(s[calls], s[ycol], 's--', color=col, ms=5)
                for _, r in s.iterrows():
                    dag = '†' if r['incomplete'] else ''
                    ax.annotate(f"{r['param']:g}{dag}", (r[calls], r[ycol]), fontsize=7, xytext=(2, 2), textcoords='offset points')
            base = d[(d.method == 'none') & (d['lambda'] == 3)].iloc[0]
            w9 = d[(d.method == 'mean') & (d.param == 9) & (d['lambda'] == 3)].iloc[0]
            ax.scatter([base[calls]], [base[ycol]], s=170, facecolors='none', edgecolors='k', zorder=5)
            ax.scatter([w9[calls]], [w9[ycol]], s=170, marker='D', facecolors='none', edgecolors='k', zorder=5)
            ax.set_xlabel('# CNV calls')
            ax.set_ylabel(yname)
        axes[ri][0].set_title(title, loc='left', fontsize=10)
    handles = [Line2D([], [], color='tab:blue', marker='o', label='unsmoothed, lambda sweep (blue labels = lambda)')]
    handles += [Line2D([], [], color=c, marker='s', ls='--', label=f'mean filter w = 3..41, lambda = {l:g} (black labels = w)') for l, c in LAM_COLORS]
    handles += [Line2D([], [], color='k', marker='o', ls='', mfc='none', ms=10, label='baseline (unsmoothed, lambda = 3)'),
                Line2D([], [], color='k', marker='D', ls='', mfc='none', ms=9, label='HiNT-RSC (w = 9, lambda = 3, frozen)')]
    if d['incomplete'].any():
        handles.append(Line2D([], [], ls='', label='† = at least one chromosome absent from the segmentation'))
    fig.legend(handles=handles, loc='upper right', ncol=2, fontsize=8, frameon=True, bbox_to_anchor=(0.995, 0.995))
    plt.tight_layout(rect=(0, 0, 1, 1 - 1.1 / (4.6 * len(rows) + 1.2)))
    plt.savefig(args.output, dpi=150)
    print('wrote', args.output)


if __name__ == '__main__':
    main()
