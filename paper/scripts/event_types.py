"""Which reference CNV events are matched, by event type (Discussion): large low-amplitude events vs focal events.

Reference events = reference segments with |log2| >= 0.3 (K562: WGS BIC-seq2 lambda = 4; PANC-1 and Caki2: DepMap WGS
copy-number segments, log2(SegmentMean) floored at -5, as in score_reference.py), not merged. Event classes:
  large_lowamp  size >= 20 Mb and 0.3 <= |log2| < 0.6
  focal_lt5     size < 5 Mb
An event is matched by a Hi-C segmentation if a Hi-C call (|log2| >= 0.3) of the same direction has IoU > 0.5 with it
(the matching rule of the paper). Reference coordinates are 1-based (shifted by 1 bp), Hi-C segment coordinates 0-based.
Reference events on Caki2 chr20 are excluded, because chr20 is absent from every mean-filter Caki2 segmentation.

Settings scored (complete Hi-C, BIC-seq2): unsmoothed at lambda 0.5, 1, 1.5, 3; mean filter w = 3 and w = 9 at lambda 1, 3.

DEPENDENCY: this script reads existing segmentations only. Run the K562 workflow (workflow/10_run_k562.sh) and the
independent-line workflow (workflow/20_run_independent.sh) first; the segmentation files are found through the `seg` column
of their results tables (--k562-results, --panc1-results, --caki2-results; relative paths are resolved against
--seg-root unless they exist as given). For PANC-1 and Caki2 the problem-region excluded run is used (rows with excluded == True when that column
exists).

Run from paper/:
  python scripts/event_types.py --out output/analyses/event_types.tsv
"""
import argparse
import csv
import math
import os

import pandas as pd

THR = 0.3
CONFIGS = [('none', 0, 0.5), ('none', 0, 1), ('none', 0, 1.5), ('none', 0, 3),
           ('mean', 3, 1), ('mean', 3, 3), ('mean', 9, 1), ('mean', 9, 3)]


def tag(method, w, lam):
    lam = f'{lam:g}'
    return f'none_l{lam}' if method == 'none' else f'mean{w}_l{lam}'


def chrom_len(ref_dir, genome):
    d = {}
    for line in open(os.path.join(ref_dir, genome, f'{genome}.len')):
        c, n = line.split('\t')[:2]
        if c in ('chrM', 'chrY', 'chrMT'):
            continue
        d[c] = int(n)
    return d


def read_seg(path, clen, offset):
    rows = []
    with open(path) as f:
        f.readline()
        for line in f:
            p = line.split()
            c, s, e, v = p[0], int(p[1]), int(p[2]), float(p[6])
            if c not in clen:
                continue
            if offset:
                s, e = max(s - 1, 0), min(e - 1, clen[c])
            else:
                s, e = max(s, 0), min(e, clen[c])
            rows.append((c, s, e, v))
    return rows


def read_depmap_csv(path, clen):
    """DepMap CN segments (ProfileID,Chromosome,Start,End,SegmentMean,...) -> the tuples read_seg gives for the seg file
    written by score_reference.depmap_to_seg (log2 floored at -5; that file's 1-bp offset removed again)."""
    rows = []
    for r in csv.DictReader(open(path)):
        c = 'chr' + r['Chromosome']
        if c not in clen:
            continue
        m = float(r['SegmentMean'])
        v = max(math.log2(m) if m > 0 else -5.0, -5.0)
        rows.append((c, max(int(r['Start']), 0), min(int(r['End']), clen[c]), v))
    return rows


def read_reference(path, clen):
    return read_depmap_csv(path, clen) if path.endswith('.csv') else read_seg(path, clen, offset=True)


def calls_by_chrom(segs, clen):
    d = {c: [] for c in clen}
    for c, s, e, v in segs:
        if v >= THR:
            d[c].append((s, e, 'gain', v))
        elif v <= -THR:
            d[c].append((s, e, 'loss', v))
    return d


def matched_ref(trg, ref):
    hit = set()
    for s, e, d, _ in trg:
        for j, (rs, re_, rd, _) in enumerate(ref):
            if rd == d and max(0, min(e, re_) - max(s, rs)) / (max(e, re_) - min(s, rs)) > 0.5:
                hit.add(j)
    return hit


def event_class(size, amp):
    if size >= 20e6 and amp < 0.6:
        return 'large_lowamp'
    if size < 5e6:
        return 'focal_lt5'
    return None


def seg_lookup(results, method, w, lam, seg_root):
    d = pd.read_csv(results, sep='\t')
    if 'segmenter' in d:
        d = d[d.segmenter == 'bicseq2']
    if 'excluded' in d:
        d = d[d.excluded.astype(str) == 'True']
    x = d[(d.setting == 'complete') & (d.method == method) & (d.param.astype(float) == w) & (d['lambda'].astype(float) == lam)]
    if len(x) != 1:
        raise SystemExit(f'{results}: expected one segmentation for {method} w={w} lambda={lam}, found {len(x)}')
    seg = x.seg.iloc[0]
    return seg if os.path.isabs(seg) or os.path.exists(seg) else os.path.join(seg_root, seg)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--k562-results', default='output/k562/results_MboI.tsv')
    ap.add_argument('--panc1-results', default='output/panc1/excl/results.tsv')
    ap.add_argument('--caki2-results', default='output/caki2/excl/results.tsv')
    ap.add_argument('--k562-ref', default='data/k562/wgs/bicseq2/K562_WGS_CNV_lambda4.txt')
    ap.add_argument('--seg-root', default='output', help="relative 'seg' paths not found from the working directory are resolved against this directory")
    ap.add_argument('--panc1-ref', default='data/panc1/depmap/PANC1_PR-Tm9gPE_CNsegments.csv', help='DepMap CN segments (.csv) or a .seg file')
    ap.add_argument('--caki2-ref', default='data/caki2/depmap/Caki2_PR-g3gtC7_CNsegments.csv', help='DepMap CN segments (.csv) or a .seg file')
    ap.add_argument('--ref-dir', default='data/ref_genome')
    ap.add_argument('--exclude', nargs='*', default=['Caki2:chr20'], help='LINE:CHROM pairs whose reference events are dropped')
    ap.add_argument('--out', default='output/analyses/event_types.tsv')
    args = ap.parse_args()

    lines = [('K562', 'hg19', args.k562_ref, args.k562_results),
             ('PANC1', 'hg38', args.panc1_ref, args.panc1_results),
             ('Caki2', 'hg38', args.caki2_ref, args.caki2_results)]
    drop = {tuple(x.split(':')) for x in args.exclude}
    tags = [tag(*c) for c in CONFIGS]
    events = []
    for line, genome, ref_file, results in lines:
        clen = chrom_len(args.ref_dir, genome)
        ref_calls = calls_by_chrom(read_reference(ref_file, clen), clen)
        flags = {}
        for (m, w, lam), t in zip(CONFIGS, tags):
            seg = seg_lookup(results, m, w, lam, args.seg_root)
            calls = calls_by_chrom(read_seg(seg, clen, offset=False), clen)
            flags[t] = {c: matched_ref(calls[c], ref_calls[c]) for c in clen}
        for c in clen:
            if (line, c) in drop:
                continue
            for j, (s, e, d, v) in enumerate(ref_calls[c]):
                k = event_class(e - s, abs(v))
                if k:
                    events.append(dict(line=line, chrom=c, event_class=k, **{t: int(j in flags[t][c]) for t in tags}))
    E = pd.DataFrame(events)

    rows = []
    for k in ('large_lowamp', 'focal_lt5'):
        for line in ('all', 'K562', 'PANC1', 'Caki2'):
            e = E[E.event_class == k] if line == 'all' else E[(E.event_class == k) & (E.line == line)]
            rows.append(dict(event_class=k, line=line, n_events=len(e), **{t: int(e[t].sum()) for t in tags}))
    T = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(args.out) or '.', exist_ok=True)
    T.to_csv(args.out, sep='\t', index=False)
    print(T.to_string(index=False))
    a = T.set_index(['event_class', 'line'])
    la, fo = a.loc[('large_lowamp', 'all')], a.loc[('focal_lt5', 'all')]
    print(f"\nlarge low-amplitude (>= 20 Mb, |log2| 0.3-0.6) matched at lambda 3: unsmoothed {la['none_l3']}, w = 3 {la['mean3_l3']}, "
          f"w = 9 {la['mean9_l3']} of {la['n_events']}; unsmoothed lambda 1: {la['none_l1']}")
    print(f"focal (< 5 Mb) matched at lambda 3: unsmoothed {fo['none_l3']}, w = 3 {fo['mean3_l3']}, w = 9 {fo['mean9_l3']} of {fo['n_events']}")
    print('wrote', args.out)


if __name__ == '__main__':
    main()
