"""Supplementary Figure S6: bin-level (a) and segment-level (b) agreement with WGS as contacts are thinned.

Reads the per-seed depth summary written by depth_summary.py; points are means over seeds, bars min-max.

Run from paper/:
  python scripts/figS6_depth.py --per-seed output/analyses/depth_per_seed.tsv --out output/figures/figS6_depth
  (use --per-seed expected/depth_per_seed.tsv to redraw the figure from the released values)
"""
import argparse
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.ticker import FixedFormatter, FixedLocator, LogLocator, NullFormatter  # noqa: E402

INK = '#52514e'
BLUE = '#2a78d6'
ORANGE = '#eb6834'
AQUA = '#1baf7a'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--per-seed', default='output/analyses/depth_per_seed.tsv')
    ap.add_argument('--out', default='output/figures/figS6_depth', help='output prefix (.png and .pdf are written)')
    args = ap.parse_args()

    d = pd.read_csv(args.per_seed, sep='\t')
    g = d.groupby('p')
    P = np.array(sorted(d.p.unique()))

    plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
                         'font.size': 7.5, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6,
                         'ytick.major.width': 0.6, 'xtick.minor.width': 0.4, 'axes.spines.top': False, 'axes.spines.right': False,
                         'legend.frameon': False, 'pdf.fonttype': 42, 'ps.fonttype': 42})

    def series(col):
        m = g[col].mean().reindex(P).values
        lo = g[col].min().reindex(P).values
        hi = g[col].max().reindex(P).values
        return m, m - lo, hi - m

    fig, axes = plt.subplots(1, 2, figsize=(180 / 25.4, 70 / 25.4))
    ax = axes[0]
    for col, lab, c, mk in (('binSCC_w0', 'Unsmoothed', INK, 'o'), ('binSCC_w3', 'Mean filter, w = 3', BLUE, 's'),
                            ('binSCC_w5', 'Mean filter, w = 5', ORANGE, '^'), ('binSCC_w9', 'Mean filter, w = 9', AQUA, 'D')):
        m, el, eh = series(col)
        ax.errorbar(P, m, yerr=[el, eh], color=c, marker=mk, ms=3.5, lw=1.2, capsize=1.5, elinewidth=0.7, label=lab)
    ax.set_ylabel('Bin-level SCC (vs WGS BIC-seq2 bins)')
    ax.set_ylim(0.35, 1.0)
    ax.legend(loc='upper right', ncol=2, fontsize=6.5, handlelength=1.8, columnspacing=1.0)

    ax = axes[1]
    for col, lab, c, mk, ls in (('best_u_SCC', 'Unsmoothed, best λ', INK, 'o', '-'),
                                ('best_s_SCC', 'Best smoothed setting (w, λ)', BLUE, 's', '-'),
                                ('l3_SCC_w9', 'Mean filter w = 9, λ = 3', AQUA, 'D', '--'),
                                ('l3_SCC_w0', 'Unsmoothed, λ = 3', INK, 'o', '--')):
        m, el, eh = series(col)
        ax.errorbar(P, m, yerr=[el, eh], color=c, marker=mk, ms=3.5, lw=1.2, ls=ls, capsize=1.5, elinewidth=0.7, label=lab,
                    mfc='white' if ls == '--' else c)
    ax.set_ylabel('Segment-level SCC (vs WGS BIC-seq2 λ = 4)')
    ax.set_ylim(0.70, 0.95)
    ax.legend(loc='upper right', ncol=2, fontsize=6.5, handlelength=2.2, columnspacing=1.0)

    ticks = [1e-4, 1e-3, 1e-2, 1e-1, 1]
    for a, lab in zip(axes, 'ab'):
        a.set_xscale('log')
        a.set_xlim(6e-5, 1.6)
        a.xaxis.set_major_locator(FixedLocator(ticks))
        a.xaxis.set_major_formatter(FixedFormatter(['0.0001', '0.001', '0.01', '0.1', '1']))
        a.xaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10) * 0.1, numticks=20))
        a.xaxis.set_minor_formatter(NullFormatter())
        a.set_xlabel('Fraction of contacts retained')
        a.grid(axis='y', color='#e4e3df', lw=0.5)
        a.set_axisbelow(True)
        a.text(-0.16, 1.04, lab, transform=a.transAxes, fontsize=10, fontweight='bold', va='bottom')
    fig.tight_layout(w_pad=2.5)
    os.makedirs(os.path.dirname(args.out) or '.', exist_ok=True)
    fig.savefig(args.out + '.png', dpi=300)
    fig.savefig(args.out + '.pdf')
    print('wrote', args.out + '.png', args.out + '.pdf')


if __name__ == '__main__':
    main()
