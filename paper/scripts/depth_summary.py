"""Summarise the depth experiment (depth_downsampling.py) into per-seed values and Supplementary Table S21.

Per (p, seed): residual noise (MAD of first differences, unsmoothed and w = 3, 5, 9), median contacts per bin,
bin-level SCC per window, best segment SCC over lambda (unsmoothed; any smoothed w x lambda; per window), and the
lambda = 3 values per window (segment SCC vs WGS BIC-seq2 lambda 4, calls, matched calls, precision, and matched
reference events >= 20 Mb = size strata 20-50 Mb + >= 50 Mb of iou50_recall_by_refsize).
Table S21 reports mean [min-max] over seeds (p = 1 has a single unthinned fit).

Run from paper/:
  python scripts/depth_summary.py --results output/depth/results_depth.tsv --noise output/depth/noise_depth.tsv \
      --per-seed output/analyses/depth_per_seed.tsv --table output/analyses/tableS21_depth.tsv
"""
import argparse
import os

import numpy as np
import pandas as pd

WINDOWS = (0, 3, 5, 9)


def strata(s, idx):
    parts = [tuple(map(int, x.split('/'))) for x in s.split()]
    return sum(parts[i][0] for i in idx), sum(parts[i][1] for i in idx)


def mmm(x, nd=3):
    x = np.asarray(x, float)
    if len(x) == 1:
        return f'{x[0]:.{nd}f}'
    return f'{x.mean():.{nd}f} [{x.min():.{nd}f}–{x.max():.{nd}f}]'


def per_seed(r, nz):
    rows = []
    for (p, seed), g in r.groupby(['p', 'seed']):
        o = dict(p=p, seed=seed)
        n = nz[(nz.p == p) & (nz.seed == seed)].iloc[0]
        for w in WINDOWS:
            o[f'noise_w{w}'] = n[f'mad_diff_w{w}']
        o['rowsum_median'] = n.rowsum_median
        u = g[g.w == 0]
        bu = u.SCC_bicseq2_l4.max()
        o['best_u_SCC'] = bu
        o['best_u_lam'] = ','.join(f'{x:g}' for x in u[u.SCC_bicseq2_l4 == bu]['lambda'])
        sm = g[g.w > 0]
        bs = sm.SCC_bicseq2_l4.max()
        o['best_s_SCC'] = bs
        o['best_s_set'] = ','.join(f'w{int(a)} l{b:g}' for a, b in sm[sm.SCC_bicseq2_l4 == bs][['w', 'lambda']].values)
        for w in WINDOWS[1:]:
            o[f'best_w{w}_SCC'] = g[g.w == w].SCC_bicseq2_l4.max()
        for w in WINDOWS:
            x = g[g.w == w]
            assert x.bin_SCC.nunique() == 1
            o[f'binSCC_w{w}'] = x.bin_SCC.iloc[0]
            l3 = x[x['lambda'] == 3].iloc[0]
            o[f'l3_SCC_w{w}'] = l3.SCC_bicseq2_l4
            o[f'l3_calls_w{w}'] = l3.iou50_calls
            o[f'l3_hits_w{w}'] = l3.iou50_hits_trg
            o[f'l3_P_w{w}'] = l3.iou50_P
            o[f'l3_L20_w{w}'] = l3.L20
        rows.append(o)
    return pd.DataFrame(rows).sort_values(['p', 'seed'], ascending=[False, True])


def table_s21(r, dd):
    t = []
    for p, g in dd.groupby('p', sort=False):
        for w, lab in ((0, 'unsmoothed'), (3, 'w = 3'), (5, 'w = 5'), (9, 'w = 9')):
            o = {'Fraction retained': f'{p:g}', 'Median contacts per bin': f'{g.rowsum_median.mean():.0f}', 'Seeds': len(g), 'Setting': lab,
                 'Residual noise (MAD of first differences)': mmm(g[f'noise_w{w}'], 4),
                 'Bin-level SCC': mmm(g[f'binSCC_w{w}']),
                 'Segment SCC, λ = 3': mmm(g[f'l3_SCC_w{w}']),
                 'Calls, λ = 3': mmm(g[f'l3_calls_w{w}'], 1), 'Matched, λ = 3': mmm(g[f'l3_hits_w{w}'], 1),
                 'Precision, λ = 3': mmm(g[f'l3_P_w{w}']),
                 '≥20 Mb matched (of 24), λ = 3': mmm(g[f'l3_L20_w{w}'], 1)}
            if w == 0:
                o['Best segment SCC over λ (λ)'] = mmm(g.best_u_SCC) + ' (λ ' + '; '.join(g.best_u_lam) + ')'
            else:
                lam = []
                for _, gg in r[(r.p == p) & (r.w == w)].groupby('seed'):
                    mx = gg.SCC_bicseq2_l4.max()
                    lam.append(','.join(f'{v:g}' for v in gg[gg.SCC_bicseq2_l4 == mx]['lambda']))
                o['Best segment SCC over λ (λ)'] = mmm(g[f'best_w{w}_SCC']) + ' (λ ' + '; '.join(lam) + ')'
            t.append(o)
    return pd.DataFrame(t)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--results', default='output/depth/results_depth.tsv')
    ap.add_argument('--noise', default='output/depth/noise_depth.tsv')
    ap.add_argument('--per-seed', default='output/analyses/depth_per_seed.tsv')
    ap.add_argument('--table', default='output/analyses/tableS21_depth.tsv', help='Supplementary Table S21')
    args = ap.parse_args()

    r = pd.read_csv(args.results, sep='\t')
    nz = pd.read_csv(args.noise, sep='\t')
    r['L20'] = r.iou50_recall_by_refsize.map(lambda s: strata(s, [5, 6])[0])
    n20 = r.iou50_recall_by_refsize.map(lambda s: strata(s, [5, 6])[1])
    assert (n20 == 24).all(), 'expected 24 reference events >= 20 Mb (K562 WGS BIC-seq2 lambda 4)'

    dd = per_seed(r, nz)
    t21 = table_s21(r, dd)
    for f in (args.per_seed, args.table):
        os.makedirs(os.path.dirname(f) or '.', exist_ok=True)
    dd.to_csv(args.per_seed, sep='\t', index=False)
    t21.to_csv(args.table, sep='\t', index=False)
    pd.set_option('display.width', 250, 'display.max_columns', 20)
    print(t21.to_string(index=False))
    print('wrote', args.per_seed, args.table)


if __name__ == '__main__':
    main()
