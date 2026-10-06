"""Regenerate all K562 HiNT-RSC results (Tables 2-6, Supplementary Tables S1-S15, S19, S20; Figure 3).

Run from paper/:  python scripts/k562_rerun.py --workers 2      (full grid; --quick = mean-filter grid at lambda 3 only)

Inputs (relative to paper/):
  data/k562/K562_bias_complete.txt        K562 complete Hi-C (intra + inter) coverage + covariates (GSE63525, 50 kb, hg19)
  data/k562/K562_bias_intra.txt           K562 intra-chromosomal-only coverage
  data/k562/K562_bias_control_target.txt  K562 coverage for the control-correction model (same coverage as complete)
  data/k562/GM12878_bias_control.txt      GM12878 control Hi-C coverage
  data/covariates/hg19_MboI_cov.txt       MboI restriction-site covariate (Section 2; see below)
  data/k562/hint_MboI/                    HiNT-CNV's own MboI run: b50000 residual bins and its segmentation (QC)
  data/k562/wgs/                          WGS references: BIC-seq2 (lambda 3/4/5 segments + bins), Control-FREEC, CNVpytor

Covariate: hg19_MboI_cov.txt (chrom, bin_num, cutSitesNumber, gcper, mapscore) was taken from HiNT-CNV's own MboI
regression inputs (K562_MboI_dataForRegression/regressionInfo_<chrom>_K562_MboI.txt, columns 1-4 = bin, MboI
cut sites, GC, mappability). Rebuilding the MboI cut-site counts from hg19 with make_covariates.py gives r = 1.000000
(99.96 % of bins identical). Only the cutSitesNumber column of the bias tables is replaced; gcper/mapscore are
checked to agree with the covariate table and rowsum (coverage) is unchanged.

Pipeline: GAM.R -> working residual -> (smoothing / control correction) -> binlib.output_bin (HiNT's per-chromosome
min-shift: obs = r - min(r), exp = -min(r)) -> BIC-seq2 (and CBS) -> metrics. For complete Hi-C the residual is taken
from HiNT's own MboI bins (--hint-bins), as in the paper; intra-only and control use GAM.R.
All bins and segmentations are cached under --outdir; the 'seg' column of the results is relative to --seg-root.
"""
import argparse
import contextlib
import io
import math
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
import pywt
from scipy.ndimage import gaussian_filter1d, median_filter
from scipy.signal import savgol_filter

from binlib import chr_bin_number, output_bin
from load import load_bin_dir_to_ratio, load_cnvpytpr_to_bin, load_freec_to_bin, load_seg_pair, load_seg_to_bin
from metric import calculate_dict_correlation

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RES = 50000
ENZYME = 'MboI'
SIZE_EDGES = [1e6, 2e6, 5e6, 1e7, 2e7, 5e7]
REFS = {'bicseq2_l4': 'data/k562/wgs/bicseq2/K562_WGS_CNV_lambda4.txt',     # primary WGS reference
        'bicseq2_l3': 'data/k562/wgs/bicseq2/K562_WGS_CNV_lambda3.txt',     # Supplementary Table S20
        'bicseq2_l5': 'data/k562/wgs/bicseq2/K562_WGS_CNV_lambda5.txt'}
WGS_BIN_DIR = 'data/k562/wgs/bicseq2/bin'
FREEC = 'data/k562/wgs/freec/G15509.K-562.2.bam_ratio.txt'
CNVPYTOR = 'data/k562/wgs/cnvpytor/K562_CNV_50000.tsv'
BIAS = {'complete': 'data/k562/K562_bias_complete.txt',
        'intra': 'data/k562/K562_bias_intra.txt',
        'control_target': 'data/k562/K562_bias_control_target.txt',
        'control': 'data/k562/GM12878_bias_control.txt'}
HINT_MBOI_SEG = 'data/k562/hint_MboI/K562_MboI_CNV_segments.txt'
HINT_MBOI_BINS = 'data/k562/hint_MboI/b50000'
COVARIATES = 'data/covariates/hg19_MboI_cov.txt'

MEAN_W = [3, 5, 7, 9, 11, 21, 41]
LAMBDA_SWEEP = [0.5, 1, 1.5, 2, 2.5, 3, 4, 5]
FILTERS = {'gaussian': MEAN_W, 'median': MEAN_W, 'savgol': [5, 7, 9, 11, 21, 41], 'haar': [1, 2, 3, 4, 5], 'sym4': [1, 2, 3, 4, 5]}
FRACTIONS = [0.3, 0.5, 0.7, 1.0]
INTRA_W = [3, 5, 7, 9]
COMBINED_W = [3, 5, 7, 9]


# ----------------------------------------------------------------------------- residuals
def swap_cutsites(bias_file, cov, out_file):
    """Replace cutSitesNumber in a committed bias table; check gcper/mapscore agree with the covariate table."""
    bias = pd.read_csv(bias_file, sep='\t')
    m = bias.merge(cov[['chrom', 'bin_num', 'cutSitesNumber', 'gcper', 'mapscore']], on=['chrom', 'bin_num'], suffixes=('', '_cov'))
    assert len(m) == len(bias), f'{bias_file}: {len(bias) - len(m)} bins missing from the covariate table'
    assert np.allclose(m.gcper, m.gcper_cov, atol=1e-6) and np.allclose(m.mapscore, m.mapscore_cov, atol=1e-6), 'gcper/mapscore mismatch'
    m['cutSitesNumber'] = m['cutSitesNumber_cov']
    m[['chrom', 'bin_num', 'cutSitesNumber', 'gcper', 'mapscore', 'rowsum']].to_csv(out_file, sep='\t', index=False)
    return m


def read_hint_bins(bin_dir):
    """Residual (obs - exp) from HiNT's own segmentation bins (bins absent from the files = 0), as in the paper."""
    import glob
    import re
    d = {c: np.zeros(n) for c, n in chr_bin_number.items()}
    for f in sorted(glob.glob(os.path.join(bin_dir, '*.bin'))):
        c = 'chr' + re.search(r'chrm_0?([0-9X]+)\.b50000', f).group(1)
        for line in open(f).readlines()[1:]:
            s, e, obs, exp, var = line.split()
            d[c][int(s) // RES] = float(obs) - float(exp)
    return d


def fit_gam(bias_file, outdir, rscript):
    os.makedirs(outdir, exist_ok=True)
    res_file = os.path.join(outdir, 'residual.txt')
    if not os.path.exists(res_file):
        subprocess.run([rscript, os.path.join(SCRIPT_DIR, 'GAM.R'), '--input', bias_file, '--log_folder', os.path.join(outdir, 'GAM_log'), '--output', res_file],
                       check=True, capture_output=True)
    bias = pd.read_csv(bias_file, sep='\t')
    res = np.loadtxt(res_file)
    assert len(res) == len(bias)
    d = {c: np.zeros(n) for c, n in chr_bin_number.items()}
    for (c, b), r in zip(bias[['chrom', 'bin_num']].itertuples(index=False), res):
        d[c][int(b)] = r
    return d


# ----------------------------------------------------------------------------- smoothing
def mean_filter(x, w):
    """Mean filter exactly as in the paper (window truncated at chromosome ends)."""
    n = len(x)
    out = np.zeros(n)
    for i in range(n):
        s = max(0, i - math.floor(w / 2))
        e = min(i + math.floor(w / 2), n)
        out[i] = np.mean(x[s:e + 1])
    return out


def smooth(x, method, p):
    if method == 'mean':
        return mean_filter(x, p)
    if method == 'gaussian':                                   # same variance as a width-p box filter
        return gaussian_filter1d(x, sigma=p / math.sqrt(12), mode='nearest', truncate=4.0)
    if method == 'median':
        return median_filter(x, size=p, mode='nearest')
    if method == 'savgol':                                     # local quadratic fit over p bins
        return savgol_filter(x, p, 2, mode='nearest')
    if method in ('haar', 'sym4'):                             # soft universal threshold on detail levels 1..p
        coeffs = pywt.wavedec(x, method, level=p, mode='symmetric')
        sigma = np.median(np.abs(coeffs[-1])) / 0.6745
        thr = sigma * math.sqrt(2 * math.log(len(x)))
        coeffs = [coeffs[0]] + [pywt.threshold(c, thr, 'soft') for c in coeffs[1:]]
        return pywt.waverec(coeffs, method, mode='symmetric')[:len(x)]
    raise ValueError(method)


def smooth_dict(res, method, p):
    return {c: smooth(v, method, p) for c, v in res.items()}


# ----------------------------------------------------------------------------- segmentation
def write_bins(res_dict, bin_dir):
    if not os.path.exists(os.path.join(bin_dir, '_config.txt')):
        os.makedirs(os.path.dirname(bin_dir), exist_ok=True)
        output_bin(res_dict, bin_dir)
    return bin_dir


def bicseq(bin_dir, lam, bicseq_pl, tmp_root):
    seg = os.path.join(bin_dir, f'bicseq_lambda{lam}.seg')
    if not os.path.exists(seg):
        tmp = os.path.join(tmp_root, os.path.basename(os.path.dirname(bin_dir)) + '_' + os.path.basename(bin_dir) + f'_l{lam}')
        os.makedirs(tmp, exist_ok=True)
        subprocess.run([bicseq_pl, os.path.join(bin_dir, '_config.txt'), seg, '--lambda', str(lam), '--bootstrap', '--tmp', tmp],
                       check=True, capture_output=True)
    return seg


def cbs(bin_dir, rscript):
    seg = os.path.join(bin_dir, 'cbs_alpha0.01.seg')
    if not os.path.exists(seg):
        subprocess.run([rscript, os.path.join(SCRIPT_DIR, 'cbs.R'), bin_dir, seg, '0.01'], check=True, capture_output=True)
    return seg


# ----------------------------------------------------------------------------- metrics
def size_bin(length):
    return sum(length >= e for e in SIZE_EDGES)


def match_calls(seg_file, ref_file, rule):
    """rule: ('iou', t) or ('reciprocal', 0.5). Returns matched target/reference call sets and sizes."""
    trg = load_seg_pair(seg_file, threshold=0.3)
    ref = load_seg_pair(ref_file, threshold=0.3, has_offset=True)
    trg_hit, ref_hit, trg_sizes, ref_sizes = set(), set(), [], []
    for c in chr_bin_number:
        for i, (s, e, d) in enumerate(trg[c]):
            trg_sizes.append(e - s)
            for j, (rs, re_, rd) in enumerate(ref[c]):
                if rd != d:
                    continue
                ov = max(0, min(e, re_) - max(s, rs))
                if rule[0] == 'iou':
                    ok = ov / (max(e, re_) - min(s, rs)) > rule[1]
                else:
                    ok = ov > rule[1] * (e - s) and ov > rule[1] * (re_ - rs)
                if ok:
                    trg_hit.add((c, i))
                    ref_hit.add((c, j))
        for j, (rs, re_, rd) in enumerate(ref[c]):
            ref_sizes.append(re_ - rs)
    return trg, ref, trg_hit, ref_hit


def pr_metrics(seg_file, ref_file, rule, prefix):
    trg, ref, trg_hit, ref_hit = match_calls(seg_file, ref_file, rule)
    n_trg = sum(len(v) for v in trg.values())
    n_ref = sum(len(v) for v in ref.values())
    out = {f'{prefix}_calls': n_trg, f'{prefix}_ref': n_ref, f'{prefix}_hits_trg': len(trg_hit), f'{prefix}_hits_ref': len(ref_hit),
           f'{prefix}_P': round(len(trg_hit) / n_trg, 3) if n_trg else None, f'{prefix}_R': round(len(ref_hit) / n_ref, 3) if n_ref else None}
    rec = [[0, 0] for _ in range(7)]
    prec = [[0, 0] for _ in range(7)]
    for c in chr_bin_number:
        for j, (rs, re_, _) in enumerate(ref[c]):
            b = size_bin(re_ - rs)
            rec[b][1] += 1
            rec[b][0] += (c, j) in ref_hit
        for i, (s, e, _) in enumerate(trg[c]):
            b = size_bin(e - s)
            prec[b][1] += 1
            prec[b][0] += (c, i) in trg_hit
    out[f'{prefix}_recall_by_refsize'] = ' '.join(f'{a}/{b}' for a, b in rec)
    out[f'{prefix}_precision_by_callsize'] = ' '.join(f'{a}/{b}' for a, b in prec)
    return out


class Scorer:
    def __init__(self):
        self.ref_bins = {k: load_seg_to_bin(v) for k, v in REFS.items()}
        self.freec = load_freec_to_bin(FREEC)
        self.cnvpytor = load_cnvpytpr_to_bin(CNVPYTOR)
        self.wgs_bin_ratio = load_bin_dir_to_ratio(WGS_BIN_DIR)

    def corr(self, a, b, pearson=True):
        with contextlib.redirect_stdout(io.StringIO()):
            return calculate_dict_correlation(a, b, pearson=pearson, spearman=True)

    def score(self, seg_file, bin_dir):
        row = {}
        segbins = load_seg_to_bin(seg_file)
        for k, rb in self.ref_bins.items():
            p, s = self.corr(segbins, rb)
            row[f'PCC_{k}'], row[f'SCC_{k}'] = p, s
        row['SCC_freec'] = self.corr(segbins, self.freec, pearson=False)[1]
        row['SCC_cnvpytor'] = self.corr(segbins, self.cnvpytor, pearson=False)[1]
        bp, bs = self.corr(load_bin_dir_to_ratio(bin_dir), self.wgs_bin_ratio)
        row['bin_PCC'], row['bin_SCC'] = bp, bs
        row.update(pr_metrics(seg_file, REFS['bicseq2_l4'], ('iou', 0.5), 'iou50'))
        row.update(pr_metrics(seg_file, REFS['bicseq2_l4'], ('iou', 0.3), 'iou30'))
        row.update(pr_metrics(seg_file, REFS['bicseq2_l4'], ('iou', 0.7), 'iou70'))
        row.update(pr_metrics(seg_file, REFS['bicseq2_l4'], ('reciprocal', 0.5), 'recip50'))
        for k in ('bicseq2_l3', 'bicseq2_l5'):
            r = pr_metrics(seg_file, REFS[k], ('iou', 0.5), k)
            row.update({kk: v for kk, v in r.items() if kk.endswith(('_P', '_R', '_ref', '_hits_trg'))})
        return row


# ----------------------------------------------------------------------------- QC
def compare_seg(a, b, tol=1e-6):
    ra = [l.split()[:7] for l in open(a).readlines()[1:]]
    rb = [l.split()[:7] for l in open(b).readlines()[1:]]
    if len(ra) != len(rb) or any(x[:3] != y[:3] for x, y in zip(ra, rb)):
        return False, f'{len(ra)} vs {len(rb)} segments / boundaries differ'
    d = max(abs(float(x[6]) - float(y[6])) for x, y in zip(ra, rb))
    return d < tol, f'max |log2 diff| = {d:.2e}'


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--covariates', default=COVARIATES, help='restriction-site covariate table (hg19, MboI); default: %(default)s')
    ap.add_argument('--outdir', default='output/k562', help='default: %(default)s')
    ap.add_argument('--seg-root', default='output', help="'seg' paths in the results are written relative to this directory; default: %(default)s")
    ap.add_argument('--bicseq', default='NBICseq-seg.pl', help='BIC-seq2 NBICseq-seg.pl (v0.7.2); default: %(default)s (on PATH)')
    ap.add_argument('--rscript', default='Rscript', help='default: %(default)s (needs mgcv, optparse, DNAcopy)')
    ap.add_argument('--workers', type=int, default=2, help='parallel segmentation jobs; default: %(default)s')
    ap.add_argument('--quick', action='store_true', help='baseline + mean-filter grid at lambda 3 only (quick check)')
    ap.add_argument('--hint-bins', default=HINT_MBOI_BINS,
                    help="HiNT-CNV MboI b50000 bin directory; the complete-Hi-C residual is taken from HiNT's own GAM output "
                         "(as in the paper). Pass '' to use the GAM.R fit instead. Default: %(default)s")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    tmp_root = os.path.join(args.outdir, 'tmp')
    cov = pd.read_csv(args.covariates, sep='\t')
    resid = {}
    for key, bias in BIAS.items():
        d = os.path.join(args.outdir, 'gam', key)
        os.makedirs(d, exist_ok=True)
        m = swap_cutsites(bias, cov, os.path.join(d, 'bias.txt'))
        if key == 'control_target':
            base = pd.read_csv(os.path.join(args.outdir, 'gam', 'complete', 'bias.txt'), sep='\t')
            assert np.allclose(m.rowsum.values, base.rowsum.values), 'control-target coverage != baseline coverage'
        resid[key] = fit_gam(os.path.join(d, 'bias.txt'), d, args.rscript)
        print(f'GAM {key}: {len(m)} bins', flush=True)
    if args.hint_bins:
        hint = read_hint_bins(args.hint_bins)
        diff = max(np.abs(hint[c] - resid['complete'][c]).max() for c in chr_bin_number)
        print(f'complete residual: using HiNT bins from {args.hint_bins} (max |diff| vs GAM.R = {diff:.3g})', flush=True)
        resid['complete'] = hint

    # ---- configurations: (setting, method, param, residual dict)
    configs = [('complete', 'none', 0, resid['complete'])]
    configs += [('complete', 'mean', w, smooth_dict(resid['complete'], 'mean', w)) for w in MEAN_W]
    if not args.quick:
        for meth, grid in FILTERS.items():
            configs += [('complete', meth, p, smooth_dict(resid['complete'], meth, p)) for p in grid]
        configs += [('intra', 'none', 0, resid['intra'])]
        configs += [('intra', 'mean', w, smooth_dict(resid['intra'], 'mean', w)) for w in INTRA_W]
        corrected = {}
        for f in FRACTIONS:
            corrected[f] = {c: resid['control_target'][c] - f * resid['control'][c] for c in chr_bin_number}
            configs.append(('control', f'f{f}', f, corrected[f]))
        configs += [('combined', 'f0.3+mean', w, smooth_dict(corrected[0.3], 'mean', w)) for w in COMBINED_W]

    # ---- segmentation jobs
    jobs = []
    for setting, method, p, res in configs:
        bin_dir = write_bins(res, os.path.join(args.outdir, setting, f'{method}_{p}'))
        lambdas = [3]
        if setting == 'complete' and method == 'none' and not args.quick:
            lambdas = LAMBDA_SWEEP
        elif setting == 'complete' and method == 'mean' and not args.quick:
            lambdas = [1, 2, 3]
        for lam in lambdas:
            jobs.append((setting, method, p, bin_dir, 'bicseq2', lam))
        if setting in ('complete',) and method in ('none', 'mean') and not args.quick:
            jobs.append((setting, method, p, bin_dir, 'cbs', None))
    print(f'{len(jobs)} segmentation jobs', flush=True)

    def run(job):
        setting, method, p, bin_dir, segmenter, lam = job
        seg = bicseq(bin_dir, lam, args.bicseq, tmp_root) if segmenter == 'bicseq2' else cbs(bin_dir, args.rscript)
        return job, seg
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        done = list(ex.map(run, jobs))

    # ---- QC: baseline must equal HiNT-CNV's own MboI segmentation
    base_seg = os.path.join(args.outdir, 'complete', 'none_0', 'bicseq_lambda3.seg')
    ok, msg = compare_seg(base_seg, HINT_MBOI_SEG)
    print(f'QC baseline (complete Hi-C, unsmoothed, lambda 3) vs HiNT-CNV MboI segmentation: {"PASS" if ok else "FAIL"} ({msg})', flush=True)

    # ---- metrics
    scorer = Scorer()
    rows = []
    for (setting, method, p, bin_dir, segmenter, lam), seg in done:
        row = {'enzyme': ENZYME, 'setting': setting, 'method': method, 'param': p, 'segmenter': segmenter, 'lambda': lam,
               'seg': os.path.relpath(seg, args.seg_root)}
        row.update(scorer.score(seg, bin_dir))
        rows.append(row)
        print(f"{setting:9s} {method:10s} {str(p):4s} {segmenter:7s} l={lam}: SCC {row['SCC_bicseq2_l4']:.3f} bin {row['bin_SCC']:.3f} "
              f"calls {row['iou50_calls']} hits {row['iou50_hits_trg']} P {row['iou50_P']} R {row['iou50_R']}", flush=True)
    out = os.path.join(args.outdir, f'results_{ENZYME}.tsv')
    pd.DataFrame(rows).to_csv(out, sep='\t', index=False)
    print('wrote', out)


if __name__ == '__main__':
    main()
