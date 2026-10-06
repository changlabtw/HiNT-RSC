"""Table 1 and Supplementary Figure S2: agreement among the K562 WGS-derived CNV references (BIC-seq2, Control-FREEC, CNVpytor).

Run from paper/:  python scripts/table1_wgs_callers.py [--outdir output]
Writes <outdir>/tables/k562/table1.tsv (+ .md) and <outdir>/figures/figureS2{a,b,c}_wgs_callers.png
(chromosomes 1-8, 9-16, 17-X).

Each caller is mapped to genome-wide 50-kb hg19 bins: BIC-seq2 (lambda = 4) segment log2 copy ratio; CNVpytor
50-kb read-depth copy number (bins without a call = 1); Control-FREEC per-window copy number (bins without a window = 0).
Table 1 is the Spearman correlation between these bin profiles. For plotting only, values are capped at the display
range (BIC-seq2 [-2, 2], CNVpytor 3, Control-FREEC 6).
"""
import argparse
import contextlib
import io
import math
import os

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

from load import chr_bin_number, load_seg_to_bin, resolution  # noqa: E402
from metric import calculate_dict_correlation  # noqa: E402

BICSEQ2 = 'data/k562/wgs/bicseq2/K562_WGS_CNV_lambda4.txt'
CNVPYTOR = 'data/k562/wgs/cnvpytor/K562_CNV_50000.tsv'
FREEC = 'data/k562/wgs/freec/G15509.K-562.2.bam_ratio.txt'
PANELS = [('a', 'chr1', 'chr8'), ('b', 'chr9', 'chr16'), ('c', 'chr17', 'chrX')]


def load_cnvpytor(path):
    d = {c: np.ones(n) for c, n in chr_bin_number.items()}
    for line in open(path):
        f = line.strip().split()
        chrom, rng = f[1].split(':')
        start, end = (int(x) for x in rng.split('-'))
        if chrom in ('Y', 'M'):
            continue
        value = float(f[3])
        for b in range(math.floor(start / resolution), math.floor(end / resolution) + 1):
            if b >= chr_bin_number[f'chr{chrom}']:
                break
            d[f'chr{chrom}'][b] = value
    return d


def load_freec(path):
    d = {c: np.zeros(n) for c, n in chr_bin_number.items()}
    for line in open(path).readlines()[1:]:
        f = line.strip().split()
        chrom, start, cn = f[0], int(float(f[1])), int(f[4])
        if chrom in ('Y', 'M'):
            continue
        b = math.floor(start / resolution)
        if b >= chr_bin_number[f'chr{chrom}']:
            break
        d[f'chr{chrom}'][b] = cn
    return d


def spearman(a, b):
    with contextlib.redirect_stdout(io.StringIO()):
        return calculate_dict_correlation(a, b, pearson=False, spearman=True)[1]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--outdir', default='output', help='default: %(default)s')
    args = ap.parse_args()
    tdir, fdir = os.path.join(args.outdir, 'tables', 'k562'), os.path.join(args.outdir, 'figures')
    os.makedirs(tdir, exist_ok=True)
    os.makedirs(fdir, exist_ok=True)

    prof = {'BIC-seq2': load_seg_to_bin(BICSEQ2), 'Control-FREEC': load_freec(FREEC), 'CNVpytor': load_cnvpytor(CNVPYTOR)}
    names = list(prof)
    tab = pd.DataFrame([[1.0 if a == b else spearman(prof[a], prof[b]) for b in names] for a in names], index=names, columns=names)
    tab.index.name = 'Reference profile'
    tab.to_csv(os.path.join(tdir, 'table1.tsv'), sep='\t', float_format='%.3f')
    md = ['**Table 1. Consistency among WGS-derived CNV reference profiles (Spearman correlation, K562, 50-kb bins)**', '',
          '| Reference profile | ' + ' | '.join(names) + ' |', '|' + '---|' * (len(names) + 1)]
    md += [f'| {a} | ' + ' | '.join(f'{tab.loc[a, b]:.3f}' for b in names) + ' |' for a in names]
    with open(os.path.join(tdir, 'table1.md'), 'w') as f:
        f.write('\n'.join(md) + '\n')
    print('\n'.join(md))

    # Supplementary Figure S2 (display caps as in the original notebook)
    bic = {c: np.where(v >= 2, 1.98, np.where(v <= -2, -1.98, v)) for c, v in prof['BIC-seq2'].items()}
    cnvp = {c: np.where(v >= 3, 2.98, v) for c, v in prof['CNVpytor'].items()}
    freec = {c: np.where(v >= 6, 5.98, v) for c, v in prof['Control-FREEC'].items()}
    chroms = list(chr_bin_number)
    for tag, first, last in PANELS:
        sel = chroms[chroms.index(first):chroms.index(last) + 1]
        fig, ax = plt.subplots(len(sel), 3, figsize=(20, 80 * len(sel) / 23), squeeze=False)
        for i, c in enumerate(sel):
            ax[i][0].scatter(range(len(bic[c])), bic[c], marker='s', s=1, linewidths=0)
            ax[i][0].set_ylim([-2, 2])
            ax[i][0].set_ylabel(c, fontsize=20)
            ax[i][0].axhline(y=0, color='grey')
            ax[i][0].axhline(y=0.2, color='grey', linestyle='--')
            ax[i][0].axhline(y=-0.2, color='grey', linestyle='--')
            ax[i][1].scatter(range(len(cnvp[c])), cnvp[c], color='red', marker='s', s=1, linewidths=0)
            ax[i][1].set_ylim([0, 3])
            ax[i][1].axhline(y=1, color='grey')
            ax[i][2].scatter(range(len(freec[c])), freec[c], color='green', marker='s', s=1, linewidths=0)
            ax[i][2].set_ylim([0, 6])
            ax[i][2].axhline(y=2, color='grey')
            if i == 0:
                ax[i][0].set_title('BIC-seq2', fontdict={'fontsize': 24})
                ax[i][1].set_title('CNVpytor', fontdict={'fontsize': 24})
                ax[i][2].set_title('Control-FREEC', fontdict={'fontsize': 24})
        out = os.path.join(fdir, f'figureS2{tag}_wgs_callers.png')
        fig.savefig(out, dpi=100, bbox_inches='tight')
        plt.close(fig)
        print('wrote', out)


if __name__ == '__main__':
    main()
