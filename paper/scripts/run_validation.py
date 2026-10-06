"""PANC-1 pre-registered run without problem-region exclusion: GAM baseline -> mean-filter smoothing -> BIC-seq2 -> metrics.

Settings fixed on K562 before the run (Supplementary Note S1): baseline (unsmoothed) at lambda 1, 2, 3 and mean filter
w = 3..41 at lambda 3, for complete and intra-only Hi-C. References: ENCODE array-state calls (primary), array replicate 2,
and the array logR with |logR| >= 0.3. Writes <outdir>/results.tsv (array-state rows of Supplementary Table S16);
score_reference.py adds the DepMap WGS reference. Genome: hg38 (data/ref_genome/hg38).
"""
import argparse
import contextlib
import io
import math
import os
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor

if os.path.basename(os.path.normpath(os.environ.get('HINT_RSC_REF', ''))) != 'hg38':
    os.environ['HINT_RSC_REF'] = 'data/ref_genome/hg38'  # PANC-1 and Caki2 are analysed on hg38

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from binlib import chr_bin_number, output_bin  # noqa: E402
from load import load_bin_dir_to_ratio, load_seg_to_bin  # noqa: E402
from metric import calculate_correlation_with_groundtruth, calculate_dict_correlation, calculate_precision_racall  # noqa: E402


def mean_filter(signal, w):
    """Mean filter with the window truncated at chromosome ends (same rule as k562_rerun.mean_filter)."""
    n = len(signal)
    out = np.zeros(n)
    for i in range(n):
        start = max(0, i - math.floor(w / 2))
        end = min(i + math.floor(w / 2), n)
        out[i] = np.mean(signal[start:end + 1])
    return out


def fit_gam(table, column, outdir, rscript, gam_r):
    os.makedirs(outdir, exist_ok=True)
    keep = ~((table.cutSitesNumber == 0) & (table.mapscore == 0) & (table.gcper == 0))
    bias = table.loc[keep, ['chrom', 'bin_num', 'cutSitesNumber', 'gcper', 'mapscore', column]].rename(columns={column: 'rowsum'})
    bias_file = os.path.join(outdir, 'bias_allchroms_nonzeros.txt')
    bias.to_csv(bias_file, sep='\t', index=False)
    res_file = os.path.join(outdir, 'residual.txt')
    subprocess.run([rscript, gam_r, '--input', bias_file, '--log_folder', os.path.join(outdir, 'GAM_log'), '--output', res_file],
                   check=True, capture_output=True)
    res = np.loadtxt(res_file)
    res_dict = {c: np.zeros(n) for c, n in chr_bin_number.items()}
    for (c, b), r in zip(bias[['chrom', 'bin_num']].itertuples(index=False), res):
        res_dict[c][int(b)] = r
    return res_dict


def segment(bin_dir, lam, seg_file, bicseq, tmp):
    if not os.path.exists(seg_file):
        os.makedirs(tmp, exist_ok=True)
        subprocess.run([bicseq, os.path.join(bin_dir, '_config.txt'), seg_file, '--lambda', str(lam), '--bootstrap', '--tmp', tmp],
                       check=True, capture_output=True)
    return seg_file


def score(seg_file, bin_dir, refs):
    row = {}
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        row['PCC'], row['SCC'] = calculate_correlation_with_groundtruth(seg_file, groundtruth_file=refs['corr'])
        hic_bins = load_bin_dir_to_ratio(bin_dir)
        ref_bins = load_seg_to_bin(refs['corr'])
        row['bin_PCC'], row['bin_SCC'] = calculate_dict_correlation(hic_bins, ref_bins)
    for name, gt in refs['calls'].items():
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            calculate_precision_racall(seg_file, groundtruth_file=gt)
        out = buf.getvalue()
        g = lambda k: float(re.search(k + r'\s*:\s*([\d.]+)', out).group(1))  # noqa: E731
        row[f'{name}_calls'] = int(g('Number of target'))
        row[f'{name}_ref'] = int(g('Number of ground truth'))
        row[f'{name}_hits'] = int(g('Number of hits'))
        row[f'{name}_P'] = g('precision')
        row[f'{name}_R'] = g('recall')
        row[f'{name}_size'] = ' '.join(f'{a}/{b}' for a, b in re.findall(r':\s+(\d+)\s+/\s+(\d+)', out))
    return row


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--coverage', default='data/panc1/PANC1_coverage_50kb.txt', help='make_coverage.py output (default: %(default)s)')
    ap.add_argument('--covariates', default='data/covariates/hg38_HindIII_cov.txt', help='make_covariates.py output (default: %(default)s)')
    ap.add_argument('--ref-prefix', default='data/panc1/array/ENCFF392GUL',
                    help='primary reference prefix (<p>.state_calls.seg, <p>.logR.seg; default: %(default)s)')
    ap.add_argument('--ref2-prefix', default='data/panc1/array/ENCFF155BFQ', help='sensitivity reference prefix (default: %(default)s)')
    ap.add_argument('--outdir', default='output/panc1/noexcl')
    ap.add_argument('--bicseq', default='NBICseq-seg.pl', help='BIC-seq2 segmentation script (default: %(default)s on PATH)')
    ap.add_argument('--rscript', default='Rscript')
    ap.add_argument('--gam', default='scripts/GAM.R', help='GAM regression script (default: %(default)s)')
    ap.add_argument('--windows', default='3,5,7,9,11,21,41')
    ap.add_argument('--baseline-lambdas', default='1,2,3')
    ap.add_argument('--workers', type=int, default=2, help='parallel BIC-seq2 jobs (default: %(default)s)')
    args = ap.parse_args()

    cov = pd.read_csv(args.covariates, sep='\t')
    rows = pd.read_csv(args.coverage, sep='\t')
    table = cov.merge(rows, on=['chrom', 'bin_num'])
    refs = {'corr': args.ref_prefix + '.logR.seg',
            'calls': {'arraystate': args.ref_prefix + '.state_calls.seg',
                      'arraystate_rep2': args.ref2_prefix + '.state_calls.seg',
                      'logR0.3': args.ref_prefix + '.logR.seg'}}
    jobs = []
    for setting, column in [('complete', 'rowsum_complete'), ('intra', 'rowsum_intra')]:
        res = fit_gam(table, column, os.path.join(args.outdir, setting), args.rscript, args.gam)
        runs = [(0, int(l) if float(l).is_integer() else float(l)) for l in args.baseline_lambdas.split(',')]
        runs += [(int(w), 3) for w in args.windows.split(',')]
        for w, lam in runs:
            bin_dir = os.path.join(args.outdir, setting, f'w{w}')
            if not os.path.exists(os.path.join(bin_dir, '_config.txt')):
                smoothed = res if w == 0 else {c: mean_filter(v, w) for c, v in res.items()}
                output_bin(smoothed, bin_dir)
            jobs.append((setting, w, lam, bin_dir, os.path.join(bin_dir, f'w{w}_lambda{lam}.seg')))
    tmp_root = os.path.join(args.outdir, 'tmp')
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        list(ex.map(lambda j: segment(j[3], j[2], j[4], args.bicseq, os.path.join(tmp_root, f'{j[0]}_w{j[1]}_l{j[2]}')), jobs))
    results = []
    for setting, w, lam, bin_dir, seg_file in jobs:
        row = {'setting': setting, 'w': w, 'lambda': lam, **score(seg_file, bin_dir, refs)}
        results.append(row)
        print(f"{setting} w={w} lambda={lam}: SCC {row['SCC']} binSCC {row['bin_SCC']} calls {row['arraystate_calls']} "
              f"hits {row['arraystate_hits']} P {row['arraystate_P']} R {row['arraystate_R']}", flush=True)
    out = os.path.join(args.outdir, 'results.tsv')
    pd.DataFrame(results).to_csv(out, sep='\t', index=False)
    print('wrote', out)


if __name__ == '__main__':
    main()
