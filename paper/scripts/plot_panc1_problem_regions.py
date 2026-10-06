"""Supplementary Figure S4: PANC-1 genome-wide Hi-C profile before and after hg38 problem-region exclusion, with DepMap WGS.

Top: DepMap WGS log2 copy ratio. Middle: Hi-C log2 ratio per 50-kb bin and BIC-seq2 segments (lambda = 3) without exclusion;
bins >= 25 % in centromere bands, assembly gaps or ENCODE blacklist v2 regions in orange. Bottom: the same after exclusion.
Inputs are the outputs of workflow/20_run_independent.sh (run from paper/). Genome: hg38 (data/ref_genome/hg38).
"""
import argparse
import os

if os.path.basename(os.path.normpath(os.environ.get('HINT_RSC_REF', ''))) != 'hg38':
    os.environ['HINT_RSC_REF'] = 'data/ref_genome/hg38'  # PANC-1 and Caki2 are analysed on hg38

import matplotlib  # noqa: E402
import numpy as np  # noqa: E402

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

from load import chr_bin_number, load_bin_dir_to_ratio, load_seg_to_bin  # noqa: E402
from panc1_rerun import ANNOTATIONS, exclusion_mask  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--noexcl-bins', default='output/panc1/noexcl/complete/w0', help='default: %(default)s')
    ap.add_argument('--noexcl-seg', default='output/panc1/noexcl/complete/w0/w0_lambda3.seg', help='default: %(default)s')
    ap.add_argument('--excl-bins', default='output/panc1/excl/complete/none_0', help='default: %(default)s')
    ap.add_argument('--excl-seg', default='output/panc1/excl/complete/none_0/bicseq_lambda3.seg', help='default: %(default)s')
    ap.add_argument('--depmap-seg', default='output/panc1/excl/depmap_PR-Tm9gPE.seg', help='default: %(default)s')
    ap.add_argument('--annotations', default='data/annotations/hg38', help='default: %(default)s')
    ap.add_argument('--min-frac', type=float, default=0.25)
    ap.add_argument('-o', '--output', default='output/figures/FigureS4_PANC1_problem_regions.png', help='default: %(default)s')
    args = ap.parse_args()

    mask = exclusion_mask(*[os.path.join(args.annotations, f) for f in ANNOTATIONS], args.min_frac)
    cat = lambda d: np.concatenate([d[c] for c in chr_bin_number])  # noqa: E731
    x = np.arange(sum(chr_bin_number.values()))
    edges = np.cumsum([0] + list(chr_bin_number.values()))
    m = cat({c: mask[c].astype(float) for c in chr_bin_number}).astype(bool)
    dm = cat(load_seg_to_bin(args.depmap_seg))
    b0, s0 = cat(load_bin_dir_to_ratio(args.noexcl_bins)), cat(load_seg_to_bin(args.noexcl_seg))
    b1, s1 = cat(load_bin_dir_to_ratio(args.excl_bins)), cat(load_seg_to_bin(args.excl_seg))
    fig, ax = plt.subplots(3, 1, figsize=(18, 9), sharex=True)
    ax[0].plot(x, dm, color='#ff4847', lw=1)
    ax[0].set_ylabel('DepMap WGS\nlog2 ratio', fontsize=9)
    ax[1].scatter(x[~m], b0[~m], s=0.3, c='#9bbcd6')
    ax[1].scatter(x[m], b0[m], s=0.6, c='#e8a33d', label='bins ≥ 25 % in centromere / gap / blacklist regions')
    ax[1].plot(x, s0, color='#1f77b4', lw=1)
    ax[1].set_ylabel('Hi-C, no exclusion\n(bins and segments)', fontsize=9)
    ax[1].legend(loc='lower left', fontsize=8, markerscale=8)
    ax[2].scatter(x, b1, s=0.3, c='#c9a6df')
    ax[2].plot(x, s1, color='#8e44ad', lw=1)
    ax[2].set_ylabel('Hi-C, problem regions\nexcluded', fontsize=9)
    for a in ax:
        for e in edges:
            a.axvline(e, c='lightgrey', lw=0.5)
        for t in (0.3, -0.3):
            a.axhline(t, c='grey', ls=':', lw=0.5)
        a.set_ylim(-1.5, 1.5)
    ax[-1].set_xticks((edges[:-1] + edges[1:]) / 2)
    ax[-1].set_xticklabels([c.replace('chr', '') for c in chr_bin_number], fontsize=8)
    ax[-1].set_xlabel('chromosome (50-kb bins, hg38)')
    plt.tight_layout()
    os.makedirs(os.path.dirname(args.output) or '.', exist_ok=True)
    plt.savefig(args.output, dpi=150)
    print('wrote', args.output)


if __name__ == '__main__':
    main()
