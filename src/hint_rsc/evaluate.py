"""Agreement of a Hi-C segmentation with a reference copy-number profile (the paper's metrics).

Segment level: Pearson and Spearman correlation between the per-bin log2 ratios of the two segmentations.
Bin level: the same correlations between the Hi-C bin log2 ratios (no segmentation) and the reference.
Calls: segments with |log2| >= 0.3; a Hi-C call matches a reference call of the same direction when their
intersection over union exceeds 0.5. Precision = matched Hi-C calls / Hi-C calls; recall = matched / reference calls.
"""
import math
import os

import numpy as np
from scipy import stats

SIZE_EDGES = (1e6, 2e6, 5e6, 1e7, 2e7, 5e7)


def read_reference(path, genome, fmt='auto'):
    """Reference segments as (chrom, start, end, log2) in 1-based seg coordinates.

    fmt 'seg': a header line, then BIC-seq2-style columns (chrom, start, end, binNum, observed, expected, log2[, pvalue])
    with 1-based coordinates (as in BIC-seq2 WGS output). fmt 'bed': chrom, start, end, log2 (0-based, half-open;
    '#' lines skipped). 'auto' picks 'bed' for .bed/.bed.gz files and 'seg' otherwise.
    """
    if fmt == 'auto':
        fmt = 'bed' if path.endswith(('.bed', '.bed.gz')) else 'seg'
    out = []
    with open(path) as f:
        if fmt == 'seg':
            f.readline()
        for line in f:
            if not line.strip() or line.startswith(('#', 'track')):
                continue
            p = line.split()
            if p[0] not in genome.nbins:
                continue
            if fmt == 'seg':
                out.append((p[0], int(p[1]), int(p[2]), float(p[6])))
            else:
                out.append((p[0], int(p[1]) + 1, int(p[2]) + 1, float(p[3])))
    return out


def read_hic_seg(path, genome):
    """Hi-C BIC-seq2 segments (0-based bin starts) as (chrom, start, end, log2)."""
    out = []
    with open(path) as f:
        f.readline()
        for line in f:
            p = line.split()
            if p[0] in genome.nbins:
                out.append((p[0], int(p[1]), int(p[2]), float(p[6])))
    return out


def seg_to_bins(segments, genome):
    res = genome.resolution
    d = {c: np.zeros(n) for c, n in genome.nbins.items()}
    for c, s, e, v in segments:
        for b in range(math.floor(s / res), math.floor(e / res) + 1):
            if b >= genome.nbins[c]:
                break
            d[c][b] = v
    return d


def bins_to_ratio(bin_dir, genome):
    """Per-bin log2(obs / exp) from a BIC-seq2 bin directory (0 for bins not written or with obs or exp not positive)."""
    res = genome.resolution
    d = {c: np.zeros(n) for c, n in genome.nbins.items()}
    for name in sorted(os.listdir(bin_dir)):
        if name.startswith('.') or not name.endswith('.bin'):
            continue
        c = name.split('.')[0].split('_')[-1]
        with open(os.path.join(bin_dir, name)) as f:
            f.readline()
            for line in f:
                s, _, obs, exp, _ = line.split()
                obs, exp = float(obs), float(exp)
                d[c][int(s) // res] = math.log2(obs / exp) if obs > 0 and exp > 0 else 0
    return d


def correlation(a, b, genome):
    x = np.concatenate([a[c] for c in genome.nbins])
    y = np.concatenate([b[c] for c in genome.nbins])
    return round(stats.pearsonr(x, y)[0], 3), round(stats.spearmanr(x, y)[0], 3)


def calls(segments, genome, threshold=0.3, offset=False):
    """{chrom: [(start, end, 'gain'|'loss')]}; offset=True converts 1-based reference coordinates."""
    d = {c: [] for c in genome.nbins}
    for c, s, e, v in segments:
        L = genome.length[c]
        if offset:
            s, e = max(s - 1, 0), min(e - 1, L)
        else:
            s, e = max(s, 0), min(e, L)
        if v >= threshold:
            d[c].append((s, e, 'gain'))
        elif v <= -threshold:
            d[c].append((s, e, 'loss'))
    return d


def match_calls(hic_calls, ref_calls, iou=0.5):
    """Return (number of Hi-C calls, number of reference calls, sizes of the matched reference calls)."""
    hit_sizes = []
    for c, trg in hic_calls.items():
        for s, e, d in trg:
            for rs, re_, rd in ref_calls[c]:
                inter = max(0, min(e, re_) - max(s, rs))
                if inter / (max(e, re_) - min(s, rs)) > iou and d == rd:
                    hit_sizes.append(re_ - rs)
    return sum(map(len, hic_calls.values())), sum(map(len, ref_calls.values())), hit_sizes


def size_table(hit_sizes, ref_calls):
    edges = (0,) + SIZE_EDGES + (math.inf,)
    ref_sizes = [e - s for v in ref_calls.values() for s, e, _ in v]
    return ' '.join(f'{sum(lo <= x < hi for x in hit_sizes)}/{sum(lo <= x < hi for x in ref_sizes)}' for lo, hi in zip(edges, edges[1:]))


def evaluate(seg_file, bin_dir, reference, genome, threshold=0.3, iou=0.5):
    """Metrics of one Hi-C segmentation (and its bin directory) against reference segments from read_reference()."""
    hic = read_hic_seg(seg_file, genome)
    ref_bins = seg_to_bins(reference, genome)
    row = {}
    row['seg_PCC'], row['seg_SCC'] = correlation(seg_to_bins(hic, genome), ref_bins, genome)
    if bin_dir:
        row['bin_PCC'], row['bin_SCC'] = correlation(bins_to_ratio(bin_dir, genome), ref_bins, genome)
    ref_calls = calls(reference, genome, threshold, offset=True)
    n_trg, n_ref, hits = match_calls(calls(hic, genome, threshold), ref_calls, iou)
    row.update({'calls': n_trg, 'ref_calls': n_ref, 'matched': len(hits),
                'precision': round(len(hits) / n_trg, 3) if n_trg else None,
                'recall': round(len(hits) / n_ref, 3) if n_ref else None,
                'matched_by_ref_size': size_table(hits, ref_calls)})
    return row
