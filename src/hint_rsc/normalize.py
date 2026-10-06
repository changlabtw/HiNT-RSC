"""Problem-region exclusion and GAM normalization (Poisson GAM on GC, mappability and restriction sites; R mgcv)."""
import gzip
import os
import subprocess

import numpy as np
import pandas as pd

from .genome import resource


def _open(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def read_problem_regions(files):
    """Regions from UCSC cytoBand ('acen' bands only), UCSC gap tables and BED files (chrom, start, end)."""
    regions = []
    for path in files:
        base = os.path.basename(path)
        for line in _open(path):
            if not line.strip() or line.startswith(('#', 'track', 'browser')):
                continue
            p = line.rstrip('\n').split('\t')
            if base.startswith('cytoBand'):
                if p[4] == 'acen':
                    regions.append((p[0], int(p[1]), int(p[2])))
            elif base.startswith('gap'):
                regions.append((p[1], int(p[2]), int(p[3])))
            else:
                regions.append((p[0], int(p[1]), int(p[2])))
    return regions


def exclusion_mask(regions, genome, min_frac=0.25):
    """{chrom: bool array over bins}: True if >= min_frac of the bin lies in a problem region (1-kb resolution)."""
    per = genome.resolution // 1000
    covered = {c: np.zeros(n * per + 60, dtype=bool) for c, n in genome.nbins.items()}
    for c, s, e in regions:
        if c in covered:
            covered[c][s // 1000:e // 1000 + 1] = True
    return {c: np.array([covered[c][b * per:(b + 1) * per].mean() >= min_frac for b in range(n)]) for c, n in genome.nbins.items()}


def load_table(coverage, covariates):
    cov = pd.read_csv(covariates, sep='\t')
    return cov.merge(pd.read_csv(coverage, sep='\t'), on=['chrom', 'bin_num'])


def fit_gam(table, column, genome, outdir, rscript='Rscript', excluded=None):
    """Fit the GAM on bins with any non-zero covariate and not excluded; return ({chrom: residual}, {chrom: fitted mask}).

    The working residuals of the fit are placed in their bins; all other bins get 0.
    """
    os.makedirs(outdir, exist_ok=True)
    order = {c: i for i, c in enumerate(genome.chroms)}
    table = table[table.chrom.isin(order)]
    table = table.iloc[np.lexsort((table.bin_num.values, table.chrom.map(order).values))]    # genome order
    keep = ~((table.cutSitesNumber == 0) & (table.mapscore == 0) & (table.gcper == 0))
    if excluded is not None:
        keep &= ~np.array([excluded[c][b] for c, b in zip(table.chrom, table.bin_num)])
    bias = table.loc[keep, ['chrom', 'bin_num', 'cutSitesNumber', 'gcper', 'mapscore', column]].rename(columns={column: 'rowsum'})
    bias_file = os.path.join(outdir, 'bias_allchroms_nonzeros.txt')
    res_file = os.path.join(outdir, 'residual.txt')
    text = bias.to_csv(sep='\t', index=False)
    if not (os.path.exists(bias_file) and open(bias_file).read() == text):    # inputs changed: refit
        with open(bias_file, 'w') as f:
            f.write(text)
        if os.path.exists(res_file):
            os.remove(res_file)
    if not os.path.exists(res_file):
        p = subprocess.run([rscript, resource('R', 'GAM.R'), '--input', bias_file, '--log_folder', os.path.join(outdir, 'GAM_log'),
                            '--output', res_file], capture_output=True, text=True)
        if p.returncode != 0:
            raise SystemExit(f'GAM.R failed (needs R with mgcv and optparse):\n{p.stderr[-2000:]}')
    res = np.loadtxt(res_file)
    if len(res) != len(bias):
        raise SystemExit(f'{res_file}: {len(res)} residuals for {len(bias)} bins; delete it and rerun')
    resid = {c: np.zeros(n) for c, n in genome.nbins.items()}
    fitted = {c: np.zeros(n, dtype=bool) for c, n in genome.nbins.items()}
    for (c, b), r in zip(bias[['chrom', 'bin_num']].itertuples(index=False), res):
        resid[c][int(b)] = r
        fitted[c][int(b)] = True
    return resid, fitted
