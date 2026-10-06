"""Post-hoc grid extension: mean filter w = 3..41 at extra lambda values (default 0.5) on cached bins, for K562
(Supplementary Table S4, Figure 3) or PANC-1 / Caki2 (Supplementary Tables S17/S18). Replaces the matching rows of --results.

Run from paper/, e.g. after k562_rerun.py:  python scripts/grid_extend.py --cell k562 --workers 2
"""
import argparse
import os
from concurrent.futures import ThreadPoolExecutor

import pandas as pd

from k562_rerun import ENZYME, MEAN_W, bicseq


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--cell', choices=['k562', 'panc1'], required=True, help="'panc1' = an independent line run by panc1_rerun.py (PANC-1 or Caki2)")
    ap.add_argument('--outdir', help='run directory with cached bins; default: output/k562 for --cell k562')
    ap.add_argument('--results', help=f'results TSV to update; default: <outdir>/results_{ENZYME}.tsv for --cell k562')
    ap.add_argument('--seg-root', default='output', help="'seg' paths are written relative to this directory; default: %(default)s")
    ap.add_argument('--bicseq', default='NBICseq-seg.pl', help='default: %(default)s (on PATH)')
    ap.add_argument('--workers', type=int, default=2, help='default: %(default)s')
    ap.add_argument('--lambdas', default='0.5')
    ap.add_argument('--array-prefix')
    ap.add_argument('--array2-prefix')
    ap.add_argument('--depmap-seg')
    args = ap.parse_args()
    if args.cell == 'k562':
        args.outdir = args.outdir or 'output/k562'
        args.results = args.results or os.path.join(args.outdir, f'results_{ENZYME}.tsv')
    elif not (args.outdir and args.results):
        ap.error('--outdir and --results are required for --cell panc1')
    lams = [float(x) for x in args.lambdas.split(',')]
    jobs = [(w, lam, os.path.join(args.outdir, 'complete', f'mean_{w}')) for w in MEAN_W for lam in lams]
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        segs = list(ex.map(lambda j: bicseq(j[2], j[1], args.bicseq, os.path.join(args.outdir, 'tmp')), jobs))
    rows = []
    if args.cell == 'k562':
        from k562_rerun import Scorer
        sc = Scorer()
        for (w, lam, b), seg in zip(jobs, segs):
            row = {'enzyme': ENZYME, 'setting': 'complete', 'method': 'mean', 'param': w, 'segmenter': 'bicseq2', 'lambda': lam,
                   'seg': os.path.relpath(seg, args.seg_root)}
            row.update(sc.score(seg, b))
            rows.append(row)
            print(f"K562 w={w} l={lam}: SCC {row['SCC_bicseq2_l4']:.3f} calls {row['iou50_calls']} matched {row['iou50_hits_trg']} P {row['iou50_P']} R {row['iou50_R']}")
    else:
        from run_validation import score
        from score_reference import hint_supported
        ra = {'corr': args.array_prefix + '.logR.seg', 'calls': {'arraystate': args.array_prefix + '.state_calls.seg',
              'arraystate_rep2': args.array2_prefix + '.state_calls.seg', 'logR0.3': args.array_prefix + '.logR.seg'}} if args.array_prefix else None
        rd = {'corr': args.depmap_seg, 'calls': {'depmap': args.depmap_seg}}
        for (w, lam, b), seg in zip(jobs, segs):
            row = {'excluded': True, 'setting': 'complete', 'method': 'mean', 'param': w, 'lambda': lam, 'seg': os.path.relpath(seg, args.seg_root)}
            if ra:
                a = score(seg, b, ra)
                row.update({f'array_{k}' if k in ('PCC', 'SCC', 'bin_PCC', 'bin_SCC') else k: v for k, v in a.items()})
            dm = score(seg, b, rd)
            row.update({f'depmap_{k}' if k in ('PCC', 'SCC', 'bin_PCC', 'bin_SCC') else k: v for k, v in dm.items()})
            sup, n, sl, nl = hint_supported(seg, args.depmap_seg)
            row.update({'hint_supported': f'{sup}/{n}', 'hint_supported_frac': round(sup / n, 3) if n else None,
                        'hint_large_supported': f'{sl}/{nl}', 'hint_large_frac': round(sl / nl, 3) if nl else None})
            rows.append(row)
            msg = f" | array matched {row['arraystate_hits']} P {row['arraystate_P']} SCC {row['array_SCC']:.3f}" if ra else ''
            print(f"w={w} l={lam}: DepMap SCC {row['depmap_SCC']:.3f} calls {row['depmap_calls']} matched {row['depmap_hits']} P {row['depmap_P']}{msg}")
    d = pd.read_csv(args.results, sep='\t')
    key = d.method.eq('mean') & d.setting.eq('complete') & d['lambda'].astype(float).isin(lams)
    if 'segmenter' in d.columns:
        key &= d.segmenter.eq('bicseq2')
    pd.concat([d[~key], pd.DataFrame(rows)]).to_csv(args.results, sep='\t', index=False)
    print('updated', args.results)


if __name__ == '__main__':
    main()
