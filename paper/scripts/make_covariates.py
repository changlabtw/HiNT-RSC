"""Build HiNT-style per-bin covariates (restriction sites, GC, mappability) from a genome FASTA.

Definitions follow the HiNT helper functions copied in load.py:
  - GC: mean of 1-kb window GC fractions ((G+C) / window length, N counted in the length) within each bin
  - cutSitesNumber: number of enzyme motif occurrences whose 1-based start falls in the bin
  - mapscore: base-pair-weighted mean of a mappability bedGraph over the bin (0 where not covered)
Output columns: chrom, bin_num, cutSitesNumber, gcper, mapscore

The shipped data/covariates/hg38_HindIII_cov.txt (PANC-1, Caki2) was built from UCSC hg38.fa.gz (goldenPath/hg38/bigZips)
and the Umap k50 track (gbdb/hg38/hoffmanMappability/k50.Umap.MultiTrackMappability.bw, converted with
`bigWigToBedGraph k50.Umap.MultiTrackMappability.bw /dev/stdout | gzip > k50.Umap.bedGraph.gz`), run from paper/:
  python scripts/make_covariates.py --fasta downloads/hg38.fa.gz --enzyme HindIII
         --mappability downloads/k50.Umap.bedGraph.gz -o data/covariates/hg38_HindIII_cov.txt
"""
import argparse
import gzip
import math
import re

import numpy as np

ENZYME_MOTIF = {'HindIII': 'AAGCTT', 'MboI': 'GATC', 'DpnII': 'GATC'}


def read_len(len_file):
    chroms = []
    with open(len_file) as f:
        for line in f:
            name, length = line.split('\t')[:2]
            chroms.append((name, int(length)))
    return chroms


def iter_fasta(path, wanted):
    name, seq = None, []
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                if name in wanted:
                    yield name, ''.join(seq).upper()
                name, seq = line[1:].split()[0], []
            elif name in wanted:
                seq.append(line.strip())
    if name in wanted:
        yield name, ''.join(seq).upper()


def gc_per_bin(seq, length, resolution, window=1000):
    arr = np.frombuffer(seq.encode(), dtype=np.uint8)
    is_gc = ((arr == ord('G')) | (arr == ord('C'))).astype(np.int32)
    n_win = math.ceil(length / window)
    win_gc = np.add.reduceat(is_gc, np.arange(0, length, window)) / np.minimum(window, length - np.arange(0, length, window))
    per = resolution // window
    n_bin = math.floor(length / resolution) + 1
    gc = np.zeros(n_bin)
    for i in range(n_bin):
        chunk = win_gc[i * per:(i + 1) * per]
        gc[i] = chunk.mean() if len(chunk) else 0.0
    assert len(win_gc) == n_win
    return gc


def sites_per_bin(seq, length, resolution, motif):
    n_bin = math.floor(length / resolution) + 1
    counts = np.zeros(n_bin)
    for m in re.finditer(f'(?={motif})', seq):
        counts[min((m.start() + 1) // resolution, n_bin - 1)] += 1  # HiNT site files use 1-based positions
    return counts


def map_per_bin(bedgraph_gz, chroms, resolution):
    lengths = dict(chroms)
    out = {c: np.zeros(math.floor(l / resolution) + 1) for c, l in chroms}
    with gzip.open(bedgraph_gz, 'rt') as f:
        for line in f:
            c, s, e, v = line.split()[:4]
            if c not in out:
                continue
            s, e, v = int(s), min(int(e), lengths[c]), float(v)
            while s < e:
                b = s // resolution
                stop = min(e, (b + 1) * resolution)
                out[c][b] += v * (stop - s)
                s = stop
    for c, l in chroms:
        sizes = np.full(len(out[c]), resolution, dtype=float)
        sizes[-1] = l - (len(out[c]) - 1) * resolution or resolution
        out[c] /= sizes
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--len', default='data/ref_genome/hg38/hg38.len', help='chromosome length file (chrom<TAB>length)')
    ap.add_argument('--enzyme', required=True, choices=sorted(ENZYME_MOTIF))
    ap.add_argument('--mappability', help='bedGraph.gz (chrom start end value); omit to write 0')
    ap.add_argument('--resolution', type=int, default=50000)
    ap.add_argument('-o', '--output', required=True)
    args = ap.parse_args()

    chroms = read_len(args.len)
    wanted = {c for c, _ in chroms}
    mapp = map_per_bin(args.mappability, chroms, args.resolution) if args.mappability else None
    rows = {}
    for name, seq in iter_fasta(args.fasta, wanted):
        length = dict(chroms)[name]
        assert len(seq) == length, f'{name}: fasta length {len(seq)} != {length}'
        gc = gc_per_bin(seq, length, args.resolution)
        sites = sites_per_bin(seq, length, args.resolution, ENZYME_MOTIF[args.enzyme])
        mp = mapp[name] if mapp else np.zeros(len(gc))
        rows[name] = (sites, gc, mp)
        print(f'{name}: {len(gc)} bins', flush=True)
    with open(args.output, 'w') as f:
        f.write('chrom\tbin_num\tcutSitesNumber\tgcper\tmapscore\n')
        for name, _ in chroms:
            sites, gc, mp = rows[name]
            for i in range(len(gc)):
                f.write(f'{name}\t{i}\t{sites[i]}\t{gc[i]}\t{mp[i]}\n')


if __name__ == '__main__':
    main()
