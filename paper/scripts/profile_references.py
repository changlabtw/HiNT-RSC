"""Bin-level agreement of the K562 Hi-C residual profiles with three WGS references (Supplementary Table S1).

For every residual profile written by the K562 run (one BIC-seq2 bin directory per setting; taken from the `seg`
column of the run's results table; relative paths not found from the working directory are resolved against --seg-root), the genome-wide 50-kb log2 profile (obs/exp of the bin files) is correlated
(Spearman SCC, Pearson PCC) with
  * the K562 WGS BIC-seq2 bin profile (data/k562/wgs/bicseq2/bin),
  * Control-FREEC bin ratios (data/k562/wgs/freec), and
  * CNVpytor 50-kb calls (data/k562/wgs/cnvpytor).

Outputs
  --out   Supplementary Table S1: unsmoothed (HiNT-CNV reproduced) and mean filter w = 3 ... 41 (complete Hi-C),
          SCC / PCC against each reference.
  --long  (optional) the same correlations for every distinct profile in the results table (all filters and settings).
The agreement between the three references themselves is printed.

Run from paper/ after the K562 run (workflow/10_run_k562.sh):
  HINT_RSC_REF=data/ref_genome/hg19 python scripts/profile_references.py \
      --results output/k562/results_MboI.tsv --out output/analyses/profile_references.tsv
"""
import argparse
import os

import numpy as np
import pandas as pd
from scipy import stats

from load import load_bin_dir_to_ratio, load_cnvpytpr_to_bin, load_freec_to_bin
from metric import chr_bin_number

REF_LABELS = {'wgs_bins': 'WGS BIC-seq2 bins', 'freec': 'Control-FREEC', 'cnvpytor': 'CNVpytor'}


def concat(d, chroms):
    return np.concatenate([np.asarray(d[c], float) for c in chroms])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--results', default='output/k562/results_MboI.tsv', help='K562 run results table (seg column -> bin directory)')
    ap.add_argument('--seg-root', default='output', help="relative 'seg' paths not found from the working directory are resolved against this directory")
    ap.add_argument('--wgs-bins', default='data/k562/wgs/bicseq2/bin')
    ap.add_argument('--freec', default='data/k562/wgs/freec/G15509.K-562.2.bam_ratio.txt')
    ap.add_argument('--cnvpytor', default='data/k562/wgs/cnvpytor/K562_CNV_50000.tsv')
    ap.add_argument('--out', default='output/analyses/profile_references.tsv', help='Supplementary Table S1')
    ap.add_argument('--long', help='optional: correlations for every distinct profile in --results')
    args = ap.parse_args()

    chroms = list(chr_bin_number.keys())
    refs = {'wgs_bins': concat(load_bin_dir_to_ratio(args.wgs_bins), chroms),
            'freec': concat(load_freec_to_bin(args.freec), chroms),
            'cnvpytor': concat(load_cnvpytpr_to_bin(args.cnvpytor), chroms)}

    d = pd.read_csv(args.results, sep='\t')
    d['bindir'] = [os.path.dirname(s if os.path.isabs(s) or os.path.exists(s) else os.path.join(args.seg_root, s)) for s in d.seg]
    u = d.drop_duplicates('bindir')[['setting', 'method', 'param', 'bindir']]
    rows = []
    for _, r in u.iterrows():
        x = concat(load_bin_dir_to_ratio(r.bindir), chroms)
        o = dict(setting=r.setting, method=r.method, param=r.param)
        for k, y in refs.items():
            o[f'S_{k}'] = round(stats.spearmanr(x, y)[0], 3)
            o[f'P_{k}'] = round(stats.pearsonr(x, y)[0], 3)
        rows.append(o)
    out = pd.DataFrame(rows)

    ks = list(refs)
    for i in range(len(ks)):
        for j in range(i + 1, len(ks)):
            print(f'reference agreement (SCC) {REF_LABELS[ks[i]]} ~ {REF_LABELS[ks[j]]}: {stats.spearmanr(refs[ks[i]], refs[ks[j]])[0]:.3f}')

    if args.long:
        os.makedirs(os.path.dirname(args.long) or '.', exist_ok=True)
        out.to_csv(args.long, sep='\t', index=False)

    t = out[(out.setting == 'complete') & out.method.isin(['none', 'mean'])].copy()
    t = t.sort_values(['method', 'param'], key=lambda s: s.map({'none': 0, 'mean': 1}) if s.name == 'method' else s)
    t.insert(0, 'Setting', ['HiNT-CNV reproduced (unsmoothed)' if m == 'none' else f'w = {int(p)}' for m, p in zip(t.method, t.param)])
    cols = {'Setting': 'Setting'}
    for k, lab in REF_LABELS.items():
        cols[f'S_{k}'] = f'SCC vs {lab}'
        cols[f'P_{k}'] = f'PCC vs {lab}'
    T = t[list(cols)].rename(columns=cols)
    os.makedirs(os.path.dirname(args.out) or '.', exist_ok=True)
    T.to_csv(args.out, sep='\t', index=False, float_format='%.3f')
    print(T.to_string(index=False))
    print('wrote', args.out)


if __name__ == '__main__':
    main()
