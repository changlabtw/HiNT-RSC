"""Post-hoc CBS variant (smooth.CNA + undo.splits = "sdundo", undo.SD = 3; Supplementary Table S19) on the cached K562
bin directories of k562_rerun.py, scored with the same Scorer. Appends rows with segmenter = 'cbs_sdundo' to the results TSV.

Run from paper/ after k562_rerun.py:  python scripts/cbs_variant.py --workers 2
"""
import argparse
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor

import pandas as pd

from k562_rerun import ENZYME, MEAN_W, SCRIPT_DIR, Scorer


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--outdir', default='output/k562', help='k562_rerun.py output directory; default: %(default)s')
    ap.add_argument('--seg-root', default='output', help="'seg' paths are written relative to this directory; default: %(default)s")
    ap.add_argument('--rscript', default='Rscript', help='default: %(default)s (needs DNAcopy)')
    ap.add_argument('--workers', type=int, default=2, help='default: %(default)s')
    args = ap.parse_args()
    configs = [('none', 0)] + [('mean', w) for w in MEAN_W]

    def run(cfg):
        bin_dir = os.path.join(args.outdir, 'complete', f'{cfg[0]}_{cfg[1]}')
        seg = os.path.join(bin_dir, 'cbs_alpha0.01_sdundo3.seg')
        if not os.path.exists(seg):
            subprocess.run([args.rscript, os.path.join(SCRIPT_DIR, 'cbs.R'), bin_dir, seg, '0.01', 'sdundo'], check=True, capture_output=True)
        return cfg, bin_dir, seg
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        done = list(ex.map(run, configs))
    scorer = Scorer()
    rows = []
    for (method, p), bin_dir, seg in done:
        row = {'enzyme': ENZYME, 'setting': 'complete', 'method': method, 'param': p, 'segmenter': 'cbs_sdundo', 'lambda': None,
               'seg': os.path.relpath(seg, args.seg_root)}
        row.update(scorer.score(seg, bin_dir))
        n_seg = sum(1 for _ in open(seg)) - 1
        rows.append(row)
        print(f"{method:5s} {p:3d}: segments {n_seg:4d} | SCC vs BIC-seq2 {row['SCC_bicseq2_l4']:.3f} FREEC {row['SCC_freec']:.3f} CNVpytor {row['SCC_cnvpytor']:.3f} | "
              f"calls {row['iou50_calls']} matched {row['iou50_hits_trg']} P {row['iou50_P']} R {row['iou50_R']}", flush=True)
    res = os.path.join(args.outdir, f'results_{ENZYME}.tsv')
    d = pd.read_csv(res, sep='\t')
    d = d[d.segmenter != 'cbs_sdundo']
    pd.concat([d, pd.DataFrame(rows)]).to_csv(res, sep='\t', index=False)
    print('updated', res)


if __name__ == '__main__':
    main()
