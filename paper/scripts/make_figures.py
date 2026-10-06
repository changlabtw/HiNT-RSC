#!/usr/bin/env python
"""K562 chromosome-track figures from the k562_rerun.py output (MboI covariate, 50-kb bins, hg19):

  Figure 2                  chr6: HiNT residual and segments, mean filter w = 9 and segments, WGS BIC-seq2 segments
  Figure 4                  chr6: K562 and GM12878 GAM residuals, control-corrected residuals (f = 0.3, 0.5, 0.7)
  Figure 5                  chr6: complete versus intra-chromosomal-only residual
  Supplementary Figure S1   chr5 and chr6: WGS BIC-seq2 log2 ratio and segments, HiNT log2 ratio and segments
  Supplementary Figure S3   chr6: residual before and after mean filtering, w = 3..41

Run from paper/ after workflow/10_run_k562.sh:   python scripts/make_figures.py [--run output/k562] [--outdir output/figures]

Conventions
- bin index = genomic start // 50000 (0-based starts)
- residual of a bin  = obs - exp          (per-chromosome .bin files)
- log2 copy ratio    = log2(obs / exp)    (as load.load_bin_dir_to_ratio);
  bins with obs <= 0 or exp <= 0 have no defined log2 ratio and are not drawn
- bins absent from a .bin file are missing and are not drawn
- segment value = log2.copyRatio, drawn as a horizontal line from
  start//50000 to end//50000; WGS BIC-seq2 segments are 1-based and are
  shifted by -1 first (load.load_seg_pair(has_offset=True)).
"""
import argparse
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from binlib import chr_bin_number  # hg19 bin numbers (HINT_RSC_REF / HINT_RSC_GENOME, default hg19)

RES = 50000
DPI = 200

WGS_SEG = 'data/k562/wgs/bicseq2/K562_WGS_CNV_lambda4.txt'
WGS_BIN = 'data/k562/wgs/bicseq2/bin'
SEG_NAME = 'bicseq_lambda3.seg'
OUT = 'output/figures'
# run-directory layout of k562_rerun.py; set by set_run()
COMPLETE_NONE = COMPLETE_MEAN = INTRA_NONE = CONTROL = GAM_K562 = GAM_GM12878 = None


def set_run(run):
    global COMPLETE_NONE, COMPLETE_MEAN, INTRA_NONE, CONTROL, GAM_K562, GAM_GM12878
    COMPLETE_NONE = f'{run}/complete/none_0'
    COMPLETE_MEAN = {w: f'{run}/complete/mean_{w}' for w in (3, 5, 7, 9, 11, 21, 41)}
    INTRA_NONE = f'{run}/intra/none_0'
    CONTROL = {f: f'{run}/control/f{f}_{f}' for f in ('0.3', '0.5', '0.7')}
    GAM_K562 = f'{run}/gam/control_target'
    GAM_GM12878 = f'{run}/gam/control'


# track colours (as in the manuscript figures)
BLUE = '#00bce5'
PURPLE = '#b658d1'
RED = '#ff4847'
GREY = '#a0a0a0'

USED = {}       # figure -> list of (track label, data file)
NOTES = []      # problems / caveats encountered


def _use(fig, label, path):
    USED.setdefault(fig, []).append((label, path))


# ----------------------------------------------------------------- loaders
def read_bin_file(path):
    df = pd.read_csv(path, sep='\t')
    df.columns = ['start', 'end', 'obs', 'exp', 'var'][:len(df.columns)]
    df['bin'] = df['start'].astype(int) // RES
    return df


def residual_track(bin_dir, chrom):
    path = f'{bin_dir}/_{chrom}.bin'
    df = read_bin_file(path)
    return df['bin'].values, (df['obs'] - df['exp']).values, path


def log2_track(path):
    df = read_bin_file(path)
    ok = (df['obs'] > 0) & (df['exp'] > 0)
    if (~ok).any():
        NOTES.append(f'{path}: {int((~ok).sum())} bin(s) with obs<=0 or exp<=0 '
                     f'have undefined log2(obs/exp) and are not drawn')
    df = df[ok]
    return df['bin'].values, np.log2(df['obs'] / df['exp']).values


def hint_log2_track(bin_dir, chrom):
    path = f'{bin_dir}/_{chrom}.bin'
    x, y = log2_track(path)
    return x, y, path


def wgs_log2_track(chrom):
    # WGS bin starts are 1-based (50001, 100001, ...): subtract 1 before // RES
    path = f'{WGS_BIN}/{chrom}.bin'
    df = pd.read_csv(path, sep='\t')
    df.columns = ['start', 'end', 'obs', 'exp', 'var'][:len(df.columns)]
    df['bin'] = (df['start'].astype(int) - 1) // RES
    ok = (df['obs'] > 0) & (df['exp'] > 0)
    if (~ok).any():
        NOTES.append(f'{path}: {int((~ok).sum())} bin(s) with obs<=0 '
                     f'have undefined log2(obs/exp) and are not drawn')
    df = df[ok]
    return df['bin'].values, np.log2(df['obs'] / df['exp']).values, path


def read_segments(path, chrom, has_offset=False):
    df = pd.read_csv(path, sep='\t')
    df = df[df['chrom'] == chrom]
    start = df['start'].astype(int).values
    end = df['end'].astype(int).values
    if has_offset:
        start = np.maximum(start - 1, 0)
        end = end - 1
    return start // RES, end // RES, df['log2.copyRatio'].astype(float).values


def gam_residual_track(gam_dir, chrom):
    bias = pd.read_csv(f'{gam_dir}/bias.txt', sep='\t')
    res = np.loadtxt(f'{gam_dir}/residual.txt')
    if len(res) != len(bias):
        raise ValueError(f'{gam_dir}: residual.txt has {len(res)} lines, '
                         f'bias.txt has {len(bias)} data rows')
    m = (bias['chrom'] == chrom).values
    return bias.loc[m, 'bin_num'].astype(int).values, res[m], \
        f'{gam_dir}/residual.txt (+ bias.txt for chrom/bin_num)'


# ----------------------------------------------------------------- drawing
def style_axis(ax, chrom, ylabel):
    ax.set_xlim(0, chr_bin_number[chrom])
    ax.set_ylim(-2, 2)
    ax.set_yticks([-2, -1, 0, 1, 2])
    ax.set_ylabel(ylabel, fontsize=12)
    ax.tick_params(labelsize=10)


def scatter_track(ax, chrom, x, y, color, ylabel):
    ax.scatter(x, y, s=4, marker='s', color=color, linewidths=0, zorder=1)
    ax.axhline(0, color=GREY, lw=2, zorder=2)
    style_axis(ax, chrom, ylabel)


def segment_track(ax, chrom, segs, color, ylabel):
    s, e, v = segs
    ax.axhline(0, color=GREY, lw=2, zorder=1)
    ax.hlines(v, s, e, color=color, lw=3, zorder=2)
    style_axis(ax, chrom, ylabel)


def new_stack(n, figsize, chrom, hspace=0.3):
    fig, axes = plt.subplots(n, 1, figsize=figsize)
    fig.subplots_adjust(left=0.07, right=0.985, top=0.95, bottom=0.03, hspace=hspace)
    axes[0].set_title(chrom, fontsize=20)
    return fig, axes


def save(fig, name):
    path = f'{OUT}/{name}'
    fig.savefig(path, dpi=DPI, bbox_inches='tight', pad_inches=0.05, facecolor='white')
    plt.close(fig)
    return path


# ----------------------------------------------------------------- figures
def fig2():
    chrom, name = 'chr6', 'figure2_chr6_smoothing.png'
    fig, ax = new_stack(5, (13.5, 10.7), chrom)
    x, y, p = residual_track(COMPLETE_NONE, chrom)
    scatter_track(ax[0], chrom, x, y, BLUE, 'HiNT residual'); _use(name, 'HiNT residual', p)
    p = f'{COMPLETE_NONE}/{SEG_NAME}'
    segment_track(ax[1], chrom, read_segments(p, chrom), BLUE, 'HiNT CNVs'); _use(name, 'HiNT CNVs', p)
    x, y, p = residual_track(COMPLETE_MEAN[9], chrom)
    scatter_track(ax[2], chrom, x, y, PURPLE, 'residual, win=9'); _use(name, 'residual, win=9', p)
    p = f'{COMPLETE_MEAN[9]}/{SEG_NAME}'
    segment_track(ax[3], chrom, read_segments(p, chrom), PURPLE, 'CNVs, win=9'); _use(name, 'CNVs, win=9', p)
    segment_track(ax[4], chrom, read_segments(WGS_SEG, chrom, has_offset=True), RED, 'BICseq2 CNVs')
    _use(name, 'BICseq2 CNVs', WGS_SEG)
    return save(fig, name)


def fig4():
    chrom, name = 'chr6', 'figure4_chr6_control.png'
    fig, ax = new_stack(5, (13.6, 10.85), chrom)
    x, y, p = gam_residual_track(GAM_K562, chrom)
    scatter_track(ax[0], chrom, x, y, BLUE, 'K562 residual'); _use(name, 'K562 residual', p)
    x, y, p = gam_residual_track(GAM_GM12878, chrom)
    scatter_track(ax[1], chrom, x, y, BLUE, 'GM12878 residual'); _use(name, 'GM12878 residual', p)
    for i, f in enumerate(('0.3', '0.5', '0.7')):
        x, y, p = residual_track(CONTROL[f], chrom)
        scatter_track(ax[2 + i], chrom, x, y, PURPLE, f'f = {f}'); _use(name, f'f = {f}', p)
    return save(fig, name)


def fig5():
    chrom, name = 'chr6', 'figure5_chr6_intra.png'
    fig, ax = new_stack(2, (15.7, 6.7), chrom, hspace=0.2)
    fig.subplots_adjust(top=0.92, bottom=0.05)
    x, y, p = residual_track(COMPLETE_NONE, chrom)
    scatter_track(ax[0], chrom, x, y, BLUE, 'residual'); _use(name, 'residual', p)
    x, y, p = residual_track(INTRA_NONE, chrom)
    scatter_track(ax[1], chrom, x, y, PURPLE, 'residual (intra)'); _use(name, 'residual (intra)', p)
    return save(fig, name)


def figS1():
    name = 'figureS1_chr5_chr6.png'
    fig = plt.figure(figsize=(13.5, 23.5))
    subfigs = fig.subfigures(2, 1, hspace=0.04)
    for sf, chrom, tag in zip(subfigs, ('chr5', 'chr6'), ('(a)', '(b)')):
        ax = sf.subplots(4, 1)
        sf.subplots_adjust(left=0.08, right=0.985, top=0.9, bottom=0.04, hspace=0.3)
        ax[0].set_title(chrom, fontsize=20)
        sf.text(0.01, 0.985, tag, fontsize=28, ha='left', va='top')
        x, y, p = wgs_log2_track(chrom)
        scatter_track(ax[0], chrom, x, y, RED, 'BICseq2 log2_copy_ratio')
        _use(name, f'{chrom} BICseq2 log2_copy_ratio', p)
        segment_track(ax[1], chrom, read_segments(WGS_SEG, chrom, has_offset=True), RED, 'BICseq2 CNVs')
        _use(name, f'{chrom} BICseq2 CNVs', WGS_SEG)
        x, y, p = hint_log2_track(COMPLETE_NONE, chrom)
        scatter_track(ax[2], chrom, x, y, BLUE, 'HiNT log2_copy_ratio')
        _use(name, f'{chrom} HiNT log2_copy_ratio', p)
        p = f'{COMPLETE_NONE}/{SEG_NAME}'
        segment_track(ax[3], chrom, read_segments(p, chrom), BLUE, 'HiNT CNVs')
        _use(name, f'{chrom} HiNT CNVs', p)
    return save(fig, name)


def figS3():
    chrom, name = 'chr6', 'figureS3_chr6_windows.png'
    wins = (3, 5, 7, 9, 11, 21, 41)
    fig, ax = new_stack(1 + len(wins), (11.7, 11.35), chrom, hspace=0.25)
    x, y, p = residual_track(COMPLETE_NONE, chrom)
    scatter_track(ax[0], chrom, x, y, BLUE, 'residual'); _use(name, 'residual', p)
    for i, w in enumerate(wins):
        x, y, p = residual_track(COMPLETE_MEAN[w], chrom)
        scatter_track(ax[1 + i], chrom, x, y, PURPLE, f'win = {w}'); _use(name, f'win = {w}', p)
    return save(fig, name)


def main():
    global OUT
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--run', default='output/k562', help='k562_rerun.py output directory; default: %(default)s')
    ap.add_argument('--outdir', default=OUT, help='default: %(default)s')
    args = ap.parse_args()
    set_run(args.run)
    OUT = args.outdir
    os.makedirs(OUT, exist_ok=True)
    outs = [fig2(), fig4(), fig5(), figS1(), figS3()]
    for o in outs:
        print('wrote', o, os.path.getsize(o), 'bytes')
    print('\n# data files per figure')
    for f, items in USED.items():
        print(f)
        for label, path in items:
            print(f'  {label:32s} {path}')
    if NOTES:
        print('\n# notes')
        for n in dict.fromkeys(NOTES):
            print(' -', n)


if __name__ == '__main__':
    main()
