"""Independent cell lines (PANC-1, Caki2) with hg38 problem-region exclusion: primary comparison plus the post hoc grid.

Problem-region exclusion (Supplementary Note S1): a 50-kb bin is dropped from the GAM fit and from segmentation when
>= 25 % of it lies in UCSC hg38 cytoBand 'acen' bands, UCSC hg38 assembly gaps or ENCODE blacklist v2 regions.
Primary comparison (fixed on K562): unsmoothed baseline vs mean filter w = 9, both at lambda = 3 (main-text Table 7).
Post hoc grid (Supplementary Tables S17, S18; Figure 6, Supplementary Figure S5): unsmoothed lambda sweep 0.5-5;
mean filter w = 3..41 at lambda 0.5, 1, 2, 3; Gaussian / median / Savitzky-Golay / Haar / sym4 at lambda 3;
intra-only Hi-C, unsmoothed and mean filter, at lambda 3.
References: DepMap 24Q4 WGS (|log2| >= 0.3) for both lines; for PANC-1 also the ENCODE array-state calls
(ENCFF392GUL; replicate ENCFF155BFQ) and the array logR with |logR| >= 0.3.
Genome: hg38 (data/ref_genome/hg38). Run from paper/, e.g.
  python scripts/panc1_rerun.py --line panc1        # -> output/panc1/excl/results.tsv
  python scripts/panc1_rerun.py --line caki2        # -> output/caki2/excl/results.tsv
"""
import argparse
import gzip
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor

if os.path.basename(os.path.normpath(os.environ.get('HINT_RSC_REF', ''))) != 'hg38':
    os.environ['HINT_RSC_REF'] = 'data/ref_genome/hg38'  # PANC-1 and Caki2 are analysed on hg38

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from binlib import chr_bin_number, output_bin  # noqa: E402
from k562_rerun import FILTERS, LAMBDA_SWEEP, MEAN_W, bicseq, smooth_dict  # noqa: E402
from run_validation import score  # noqa: E402
from score_reference import depmap_to_seg, hint_supported  # noqa: E402

RES = 50000
LINES = {
    'panc1': {'coverage': 'data/panc1/PANC1_coverage_50kb.txt',
              'depmap_csv': 'data/panc1/depmap/PANC1_PR-Tm9gPE_CNsegments.csv', 'profile': 'PR-Tm9gPE',
              'array_prefix': 'data/panc1/array/ENCFF392GUL', 'array2_prefix': 'data/panc1/array/ENCFF155BFQ'},
    'caki2': {'coverage': 'data/caki2/Caki2_coverage_50kb.txt',
              'depmap_csv': 'data/caki2/depmap/Caki2_PR-g3gtC7_CNsegments.csv', 'profile': 'PR-g3gtC7',
              'array_prefix': None, 'array2_prefix': None},
}
ANNOTATIONS = ('cytoBand.txt.gz', 'gap.txt.gz', 'blacklist_v2.bed.gz')


def exclusion_mask(cyto, gap, blacklist, min_frac):
    """Per chromosome, a boolean array over 50-kb bins: True if >= min_frac of the bin is in a problem region (1-kb resolution)."""
    regions = []
    for l in gzip.open(cyto, 'rt'):
        c, s, e, _, st = l.rstrip('\n').split('\t')[:5]
        if st == 'acen':
            regions.append((c, int(s), int(e)))
    for l in gzip.open(gap, 'rt'):
        p = l.rstrip('\n').split('\t')
        regions.append((p[1], int(p[2]), int(p[3])))
    for l in gzip.open(blacklist, 'rt'):
        c, s, e = l.split('\t')[:3]
        regions.append((c, int(s), int(e)))
    covered = {c: np.zeros(n * RES // 1000 + 60, dtype=bool) for c, n in chr_bin_number.items()}
    for c, s, e in regions:
        if c in covered:
            covered[c][s // 1000:e // 1000 + 1] = True
    return {c: np.array([covered[c][b * 50:(b + 1) * 50].mean() >= min_frac for b in range(n)]) for c, n in chr_bin_number.items()}


def absent_chromosomes(seg_file):
    present = {l.split('\t')[0] for l in open(seg_file).readlines()[1:]}
    return [c for c in chr_bin_number if c not in present]


def fit_gam(table, column, outdir, rscript, gam_r, excluded):
    os.makedirs(outdir, exist_ok=True)
    keep = ~((table.cutSitesNumber == 0) & (table.mapscore == 0) & (table.gcper == 0))
    if excluded is not None:
        keep &= ~np.array([excluded[c][b] for c, b in zip(table.chrom, table.bin_num)])
    bias = table.loc[keep, ['chrom', 'bin_num', 'cutSitesNumber', 'gcper', 'mapscore', column]].rename(columns={column: 'rowsum'})
    bias_file = os.path.join(outdir, 'bias_allchroms_nonzeros.txt')
    bias.to_csv(bias_file, sep='\t', index=False)
    res_file = os.path.join(outdir, 'residual.txt')
    if not os.path.exists(res_file):
        subprocess.run([rscript, gam_r, '--input', bias_file, '--log_folder', os.path.join(outdir, 'GAM_log'), '--output', res_file],
                       check=True, capture_output=True)
    res = np.loadtxt(res_file)
    d = {c: np.zeros(n) for c, n in chr_bin_number.items()}
    for (c, b), r in zip(bias[['chrom', 'bin_num']].itertuples(index=False), res):
        d[c][int(b)] = r
    return d, int(keep.sum())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--line', choices=sorted(LINES), required=True, help='cell line; sets the defaults below')
    ap.add_argument('--coverage', help='make_coverage.py output (default: data/<line>/..._coverage_50kb.txt)')
    ap.add_argument('--covariates', default='data/covariates/hg38_HindIII_cov.txt', help='default: %(default)s')
    ap.add_argument('--depmap-csv', help='DepMap segments (OmicsCNSegmentsProfile.csv or the shipped per-profile extract)')
    ap.add_argument('--profile', help='DepMap ProfileID (PANC-1 PR-Tm9gPE, Caki2 PR-g3gtC7)')
    ap.add_argument('--array-prefix', help='array reference prefix (PANC-1 default data/panc1/array/ENCFF392GUL; none for Caki2)')
    ap.add_argument('--array2-prefix', help='array replicate prefix (PANC-1 default data/panc1/array/ENCFF155BFQ)')
    ap.add_argument('--annotations', default='data/annotations/hg38', help='folder with ' + ', '.join(ANNOTATIONS) + ' (default: %(default)s)')
    ap.add_argument('--no-exclude', action='store_true', help='skip problem-region exclusion')
    ap.add_argument('--min-frac', type=float, default=0.25)
    ap.add_argument('--outdir', help='default: output/<line>/excl (output/<line>/noexcl_grid with --no-exclude)')
    ap.add_argument('--bicseq', default='NBICseq-seg.pl', help='BIC-seq2 segmentation script (default: %(default)s on PATH)')
    ap.add_argument('--rscript', default='Rscript')
    ap.add_argument('--gam', default='scripts/GAM.R', help='GAM regression script (default: %(default)s)')
    ap.add_argument('--workers', type=int, default=2, help='parallel BIC-seq2 jobs (default: %(default)s)')
    args = ap.parse_args()
    for k, v in LINES[args.line].items():
        if getattr(args, k) is None:
            setattr(args, k, v)
    if args.outdir is None:
        args.outdir = os.path.join('output', args.line, 'noexcl_grid' if args.no_exclude else 'excl')
    os.makedirs(args.outdir, exist_ok=True)
    depmap_seg = depmap_to_seg(args.depmap_csv, args.profile, os.path.join(args.outdir, f'depmap_{args.profile}.seg'))

    table = pd.read_csv(args.covariates, sep='\t').merge(pd.read_csv(args.coverage, sep='\t'), on=['chrom', 'bin_num'])
    excluded = None if args.no_exclude else exclusion_mask(*[os.path.join(args.annotations, f) for f in ANNOTATIONS], args.min_frac)
    if excluded is not None:
        print(f'excluded bins: {sum(int(v.sum()) for v in excluded.values())}', flush=True)
    resid = {}
    for setting, column in [('complete', 'rowsum_complete'), ('intra', 'rowsum_intra')]:
        resid[setting], n = fit_gam(table, column, os.path.join(args.outdir, 'gam', setting), args.rscript, args.gam, excluded)
        print(f'GAM {setting}: {n} bins', flush=True)

    configs = [('complete', 'none', 0, resid['complete'], LAMBDA_SWEEP)]
    configs += [('complete', 'mean', w, smooth_dict(resid['complete'], 'mean', w), [0.5, 1, 2, 3]) for w in MEAN_W]
    for meth, grid in FILTERS.items():
        configs += [('complete', meth, p, smooth_dict(resid['complete'], meth, p), [3]) for p in grid]
    configs += [('intra', 'none', 0, resid['intra'], [3])]
    configs += [('intra', 'mean', w, smooth_dict(resid['intra'], 'mean', w), [3]) for w in MEAN_W]

    jobs = []
    for setting, method, p, res, lambdas in configs:
        bin_dir = os.path.join(args.outdir, setting, f'{method}_{p}')
        if not os.path.exists(os.path.join(bin_dir, '_config.txt')):
            os.makedirs(os.path.dirname(bin_dir), exist_ok=True)
            output_bin(res, bin_dir)
        jobs += [(setting, method, p, bin_dir, lam) for lam in lambdas]
    print(f'{len(jobs)} segmentation jobs', flush=True)
    tmp_root = os.path.join(args.outdir, 'tmp')
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        segs = list(ex.map(lambda j: bicseq(j[3], j[4], args.bicseq, tmp_root), jobs))

    refs_array = None if not args.array_prefix else {
        'corr': args.array_prefix + '.logR.seg',
        'calls': {'arraystate': args.array_prefix + '.state_calls.seg', 'arraystate_rep2': args.array2_prefix + '.state_calls.seg',
                  'logR0.3': args.array_prefix + '.logR.seg'}}
    refs_depmap = {'corr': depmap_seg, 'calls': {'depmap': depmap_seg}}
    rows = []
    for (setting, method, p, bin_dir, lam), seg in zip(jobs, segs):
        row = {'excluded': excluded is not None, 'setting': setting, 'method': method, 'param': p, 'lambda': lam, 'seg': seg}
        if refs_array:
            a = score(seg, bin_dir, refs_array)
            row.update({f'array_{k}' if k in ('PCC', 'SCC', 'bin_PCC', 'bin_SCC') else k: v for k, v in a.items()})
        dm = score(seg, bin_dir, refs_depmap)
        row.update({f'depmap_{k}' if k in ('PCC', 'SCC', 'bin_PCC', 'bin_SCC') else k: v for k, v in dm.items()})
        sup, n, sup_l, n_l = hint_supported(seg, depmap_seg)
        row.update({'hint_supported': f'{sup}/{n}', 'hint_supported_frac': round(sup / n, 3) if n else None,
                    'hint_large_supported': f'{sup_l}/{n_l}', 'hint_large_frac': round(sup_l / n_l, 3) if n_l else None})
        row['absent_chromosomes'] = ','.join(absent_chromosomes(seg)) or 'none'
        rows.append(row)
        print(f"{setting:8s} {method:8s} {str(p):3s} l={lam}: DepMap SCC {row['depmap_SCC']:.3f} bin {row['depmap_bin_SCC']:.3f} "
              f"calls {row['depmap_calls']} hits {row['depmap_hits']} P {row['depmap_P']} | array hits {row.get('arraystate_hits')} P {row.get('arraystate_P')} "
              f"| HiNT-style {row['hint_supported']} large {row['hint_large_supported']}", flush=True)
    out = os.path.join(args.outdir, 'results.tsv')
    pd.DataFrame(rows).to_csv(out, sep='\t', index=False)
    print('wrote', out)


if __name__ == '__main__':
    main()
