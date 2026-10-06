"""One-dimensional Hi-C coverage (row sums of the raw contact matrix) from a .hic file, as in HiNT-CNV.

rowsum_complete counts intra- and inter-chromosomal contacts; rowsum_intra counts intra-chromosomal contacts only
(about ten times faster to compute, since only the diagonal blocks are read).
"""
import time

import numpy as np


def _hic_chromosomes(hic):
    import hicstraw
    return {c.name: (c.index, c.length) for c in hicstraw.HiCFile(hic).getChromosomes()}


def _records(hic, a, b, resolution):
    import hicstraw
    rec = hicstraw.straw('observed', 'NONE', hic, a, b, 'BP', resolution)
    n = len(rec)
    x = np.fromiter((r.binX for r in rec), dtype=np.int64, count=n) // resolution
    y = np.fromiter((r.binY for r in rec), dtype=np.int64, count=n) // resolution
    v = np.nan_to_num(np.fromiter((r.counts for r in rec), dtype=float, count=n))
    return x, y, v


def coverage_profiles(hic, genome, chroms=None, intra_only=False, log=print):
    """Return {chrom: (complete, intra)} arrays over the genome's bins (complete is None with intra_only).

    straw returns each block once, with binX on the chromosome that comes first in the .hic file, so every
    inter-chromosomal block is read once and added to both chromosomes.
    """
    res = genome.resolution
    info = _hic_chromosomes(hic)

    def name(c):
        if c in info:
            return c
        if c.startswith('chr') and c[3:] in info:
            return c[3:]
        raise SystemExit(f'{hic}: chromosome {c} not found (the .hic file has {", ".join(sorted(info))})')

    for c in genome.chroms:
        if info[name(c)][1] // res + 1 != genome.nbins[c]:
            raise SystemExit(f'{hic}: {c} has length {info[name(c)][1]} in the .hic file but {genome.length[c]} in genome '
                             f'{genome.name}; pass the matching --genome')
    targets = list(chroms) if chroms else genome.chroms
    intra = {c: np.zeros(genome.nbins[c]) for c in targets}
    inter = {c: np.zeros(genome.nbins[c]) for c in targets}
    order = sorted(genome.chroms, key=lambda c: info[name(c)][0])
    for i, a in enumerate(order):
        t = time.time()
        if a in targets:
            x, y, v = _records(hic, name(a), name(a), res)
            off = x != y                                   # straw gives the upper triangle; mirror off-diagonal contacts
            np.add.at(intra[a], x, v)
            np.add.at(intra[a], y[off], v[off])
        if not intra_only:
            for b in order[i + 1:]:
                if a not in targets and b not in targets:
                    continue
                x, y, v = _records(hic, name(a), name(b), res)
                if a in targets:
                    np.add.at(inter[a], x, v)
                if b in targets:
                    np.add.at(inter[b], y, v)
        if a in targets:
            log(f'{a}: {time.time() - t:.0f}s')
    return {c: (None if intra_only else intra[c] + inter[c], intra[c]) for c in targets}


def write_coverage(profiles, out):
    intra_only = next(iter(profiles.values()))[0] is None
    with open(out, 'w') as f:
        f.write('chrom\tbin_num\trowsum_intra\n' if intra_only else 'chrom\tbin_num\trowsum_complete\trowsum_intra\n')
        for chrom, (complete, intra) in profiles.items():
            for i in range(len(intra)):
                f.write(f'{chrom}\t{i}\t{intra[i]}\n' if intra_only else f'{chrom}\t{i}\t{complete[i]}\t{intra[i]}\n')
    return out
