"""Sequencing-depth experiment (Supplementary Table S21, Supplementary Figure S6): binomial thinning of K562 coverage.

For each retained fraction p the per-bin contact totals (rowsum column of the K562 MboI GAM input table written by
k562_rerun.py, <k562 outdir>/gam/complete/bias.txt) are thinned binomially, rowsum ~ Binomial(round(rowsum), p), with
numpy default_rng(seed). The thinned table is refitted with GAM.R; the residual is segmented unsmoothed and after the
mean filter (w = 3, 5, 9), with HiNT's min-shift transform (binlib.output_bin), by BIC-seq2 over a lambda grid, and
scored against the K562 WGS references (k562_rerun.Scorer). p = 1 is the unthinned table, fitted once (seed 0).

Approximation: each contact contributes to two bins' rowsum; thinning bins independently ignores that correlation
(and keeps the intra/inter mix of each bin fixed in expectation).

Outputs (in --outdir): results_depth.tsv (one row per p, seed, w, lambda), noise_depth.tsv (one row per p, seed:
total/median rowsum and the MAD of first differences of the residual, unsmoothed and per window), runs/ (thinned tables,
GAM fits, bin directories and segmentations; cached, so an interrupted run resumes).

Paper arguments (the defaults):
  --fractions 1,0.5,0.2,0.1,0.05,0.02,0.01,0.005,0.002,0.001,0.0005,0.0002,0.0001 --seeds 1,2
  --windows 3,5,9 (plus unsmoothed) --lambdas 0.5,1,1.5,2,2.5,3,4,5
Run from paper/ after workflow/10_run_k562.sh (about 25 GAM fits and 800 BIC-seq2 runs):
  HINT_RSC_REF=data/ref_genome/hg19 python scripts/depth_downsampling.py \
      --bias output/k562/gam/complete/bias.txt --outdir output/depth --bicseq NBICseq-seg.pl --rscript Rscript --workers 2
"""
import argparse
import os
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd

import k562_rerun as kr

FRACTIONS = '1,0.5,0.2,0.1,0.05,0.02,0.01,0.005,0.002,0.001,0.0005,0.0002,0.0001'


def tag(p, seed):
    return f'p{p}_s{seed}'


def make_residual(p, seed, bias, outdir, rscript):
    d = os.path.join(outdir, 'runs', tag(p, seed))
    os.makedirs(d, exist_ok=True)
    bf = os.path.join(d, 'bias.txt')
    if not os.path.exists(bf):
        b = pd.read_csv(bias, sep='\t')
        if p < 1:
            rng = np.random.default_rng(seed)
            b['rowsum'] = rng.binomial(np.round(b.rowsum.values).astype(np.int64), p).astype(float)
        b.to_csv(bf, sep='\t', index=False)
    t = time.time()
    res = kr.fit_gam(bf, os.path.join(d, 'gam'), rscript)
    print(f'GAM {tag(p, seed)} {time.time() - t:.0f}s', flush=True)
    return d, res


def noise(res, bias_file):
    """MAD of first differences of the residual over bins present in the bias table (within chromosomes)."""
    b = pd.read_csv(bias_file, sep='\t')
    diffs = []
    for c, g in b.groupby('chrom', sort=False):
        idx = np.sort(g.bin_num.values.astype(int))
        diffs.append(np.diff(res[c][idx]))
    dd = np.concatenate(diffs)
    mad = np.median(np.abs(dd - np.median(dd)))
    return mad, mad / 0.6745 / np.sqrt(2)


def parse_list(s, typ):
    return [typ(x) for x in s.split(',') if x != '']


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--bias', default='output/k562/gam/complete/bias.txt', help='K562 MboI GAM input table written by k562_rerun.py')
    ap.add_argument('--outdir', default='output/depth')
    ap.add_argument('--fractions', default=FRACTIONS, help='retained fractions p (1 = unthinned, fitted once with seed 0)')
    ap.add_argument('--seeds', default='1,2', help='thinning seeds for p < 1')
    ap.add_argument('--windows', default='3,5,9', help='mean-filter windows (the unsmoothed residual is always included)')
    ap.add_argument('--lambdas', default='0.5,1,1.5,2,2.5,3,4,5')
    ap.add_argument('--bicseq', default='NBICseq-seg.pl', help='BIC-seq2 segmentation script')
    ap.add_argument('--rscript', default='Rscript')
    ap.add_argument('--workers', type=int, default=2)
    args = ap.parse_args()

    ps = parse_list(args.fractions, float)
    seeds = parse_list(args.seeds, int)
    ws = parse_list(args.windows, int)
    lams = parse_list(args.lambdas, float)
    os.makedirs(args.outdir, exist_ok=True)

    cases = [(1.0, 0)] if 1.0 in ps else []
    cases += [(p, s) for p in ps if p < 1 for s in seeds]
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        gams = list(ex.map(lambda c: make_residual(*c, args.bias, args.outdir, args.rscript), cases))

    jobs, noise_rows = [], []
    for (p, seed), (d, res) in zip(cases, gams):
        bf = os.path.join(d, 'bias.txt')
        tot = pd.read_csv(bf, sep='\t').rowsum
        mad, sig = noise(res, bf)
        nr = {'p': p, 'seed': seed, 'rowsum_total': tot.sum(), 'rowsum_median': tot.median(), 'zero_bins': int((tot == 0).sum()),
              'mad_diff_w0': mad, 'sigma_w0': sig}
        for w in ws:
            sm = kr.smooth_dict(res, 'mean', w)
            nr[f'mad_diff_w{w}'], nr[f'sigma_w{w}'] = noise(sm, bf)
            bd = kr.write_bins(sm, os.path.join(d, f'mean_{w}'))
            jobs += [(p, seed, w, lam, bd) for lam in lams]
        bd = kr.write_bins(res, os.path.join(d, 'none_0'))
        jobs += [(p, seed, 0, lam, bd) for lam in lams]
        noise_rows.append(nr)
    pd.DataFrame(noise_rows).to_csv(os.path.join(args.outdir, 'noise_depth.tsv'), sep='\t', index=False)
    print(f'{len(jobs)} BIC-seq2 jobs', flush=True)

    tmp_root = os.path.join(args.outdir, 'tmp')
    t0 = time.time()

    def run(j):
        p, seed, w, lam, bd = j
        seg = kr.bicseq(bd, f'{lam:g}', args.bicseq, tmp_root)
        print(f'  seg {tag(p, seed)} w={w} lambda={lam:g} (elapsed {time.time() - t0:.0f}s)', flush=True)
        return j, seg
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        done = list(ex.map(run, jobs))

    scorer = kr.Scorer()
    rows = []
    for (p, seed, w, lam, bd), seg in done:
        r = {'p': p, 'seed': seed, 'w': w, 'lambda': lam, 'seg': seg}
        r.update(scorer.score(seg, bd))
        rows.append(r)
    out = os.path.join(args.outdir, 'results_depth.tsv')
    pd.DataFrame(rows).to_csv(out, sep='\t', index=False)
    print('wrote', out, flush=True)


if __name__ == '__main__':
    main()
