"""PANC-1 baseline diagnosis: overlap of unmatched calls and extreme bins with hg38 problem regions (Supplementary Note S1).

Problem regions: UCSC hg38 cytoBand 'acen' bands, UCSC hg38 assembly gaps and ENCODE blacklist v2.
Calls are segments with |log2| >= 0.3; a call is matched if a reference call of the same direction has IoU > 0.5.
Run from paper/ on the unexcluded baseline (complete Hi-C, lambda = 3) against DepMap:
  python scripts/diagnose_problem_regions.py --bins output/panc1/noexcl/complete/w0 --ref output/panc1/noexcl/depmap_PR-Tm9gPE.seg
Genome: hg38 (data/ref_genome/hg38).
"""
import argparse
import collections
import gzip
import os

if os.path.basename(os.path.normpath(os.environ.get('HINT_RSC_REF', ''))) != 'hg38':
    os.environ['HINT_RSC_REF'] = 'data/ref_genome/hg38'  # PANC-1 and Caki2 are analysed on hg38

import numpy as np  # noqa: E402

from load import chr_bin_number, load_bin_dir_to_ratio, load_seg_pair  # noqa: E402

RES = 50000


def problem_regions(folder):
    reg = collections.defaultdict(list)
    for l in gzip.open(os.path.join(folder, 'cytoBand.txt.gz'), 'rt'):
        c, s, e, _, st = l.rstrip('\n').split('\t')[:5]
        if st == 'acen':
            reg[c].append((int(s), int(e), 'acen'))
    for l in gzip.open(os.path.join(folder, 'gap.txt.gz'), 'rt'):
        p = l.rstrip('\n').split('\t')
        reg[p[1]].append((int(p[2]), int(p[3]), 'gap'))
    for l in gzip.open(os.path.join(folder, 'blacklist_v2.bed.gz'), 'rt'):
        c, s, e = l.split('\t')[:3]
        reg[c].append((int(s), int(e), 'blacklist'))
    return reg


def bin_masks(reg):
    """Fraction of each 50-kb bin in problem regions (1-kb resolution) and the region types touching it."""
    frac, kind = {}, {}
    for c, n in chr_bin_number.items():
        cov = np.zeros(n * 50 + 60, dtype=bool)
        kind[c] = [set() for _ in range(n)]
        for s, e, k in reg[c]:
            cov[s // 1000:e // 1000 + 1] = True
            for b in range(s // RES, min(e // RES, n - 1) + 1):
                kind[c][b].add(k)
        frac[c] = np.array([cov[b * 50:(b + 1) * 50].mean() for b in range(n)])
    return frac, kind


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--bins', default='output/panc1/noexcl/complete/w0', help='BIC-seq2 bin folder of the segmentation (default: %(default)s)')
    ap.add_argument('--seg', help='segmentation (default: <bins>/w0_lambda3.seg, else <bins>/bicseq_lambda3.seg)')
    ap.add_argument('--ref', default='output/panc1/noexcl/depmap_PR-Tm9gPE.seg', help='reference seg (default: %(default)s)')
    ap.add_argument('--annotations', default='data/annotations/hg38', help='default: %(default)s')
    ap.add_argument('--min-frac', type=float, default=0.25, help='problem-region fraction for a bin or call (default: %(default)s)')
    ap.add_argument('--extreme', type=float, default=1.5, help='|log2 ratio| of an extreme bin (default: %(default)s)')
    args = ap.parse_args()
    seg = args.seg or next((p for p in (os.path.join(args.bins, 'w0_lambda3.seg'), os.path.join(args.bins, 'bicseq_lambda3.seg'))
                            if os.path.exists(p)), None)
    if seg is None:
        raise SystemExit(f'no lambda-3 segmentation in {args.bins}; pass --seg')

    frac, kind = bin_masks(problem_regions(args.annotations))
    allf = np.concatenate([frac[c] for c in frac])
    print(f'segmentation: {seg}\nreference: {args.ref}')
    print(f'bins with >= {args.min_frac:.0%} in problem regions: {(allf >= args.min_frac).sum()} of {len(allf)} '
          f'({100 * (allf >= args.min_frac).mean():.1f}%)')
    trg = load_seg_pair(seg, 0.3)
    rf = load_seg_pair(args.ref, 0.3, has_offset=True)
    rows = []
    for c, calls in trg.items():
        for s, e, d in calls:
            sup = any(rd == d and max(0, min(e, r1) - max(s, r0)) / (max(e, r1) - min(s, r0)) > 0.5 for r0, r1, rd in rf[c])
            b0, b1 = s // RES, min(e // RES, chr_bin_number[c] - 1)
            rows.append((c, s, e, d, e - s, sup, frac[c][b0:b1 + 1].mean(), set().union(*kind[c][b0:b1 + 1])))
    print('\ncalls vs reference, IoU > 0.5:')
    for grp, name in [(True, 'matched'), (False, 'unmatched')]:
        g = [r for r in rows if r[5] == grp]
        mean = f'{np.mean([r[6] for r in g]):.2f}' if g else '-'
        print(f'  {name:9s}: n={len(g):2d} | <2 Mb: {sum(r[4] < 2e6 for r in g):2d} | mean fraction in problem regions '
              f'{mean} | calls >= {args.min_frac:.0%} in problem regions: {sum(r[6] >= args.min_frac for r in g)}')
    print('\nunmatched calls (chrom, start Mb, end Mb, dir, size Mb, fraction in problem regions, region types):')
    for r in sorted([r for r in rows if not r[5]], key=lambda x: (x[0], x[1])):
        print(f'  {r[0]:5s} {r[1] / 1e6:7.2f} {r[2] / 1e6:7.2f} {r[3]:4s} {r[4] / 1e6:6.2f}  {r[6]:.2f}  {",".join(sorted(r[7])) or "-"}')
    lr = load_bin_dir_to_ratio(args.bins)
    ext = [(c, b) for c in lr for b in range(len(lr[c])) if abs(lr[c][b]) >= args.extreme]
    print(f'\nbins with |log2 ratio| >= {args.extreme}: {len(ext)} | >= {args.min_frac:.0%} in problem regions: '
          f'{sum(frac[c][b] >= args.min_frac for c, b in ext)}')
    print('  by chromosome:', dict(collections.Counter(c for c, _ in ext).most_common(10)))


if __name__ == '__main__':
    main()
