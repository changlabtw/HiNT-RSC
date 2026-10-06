"""Per-bin GAM covariates (restriction sites, GC content, mappability) from a genome FASTA, as in HiNT-CNV.

  cutSitesNumber  occurrences of the enzyme motif(s) whose 1-based start falls in the bin
  gcper           mean of the 1-kb window GC fractions in the bin ((G + C) / window length, N counted in the length)
  mapscore        base-pair-weighted mean of a mappability bedGraph over the bin (0 where not covered)
"""
import gzip
import math
import re

import numpy as np

ENZYME_MOTIFS = {'HindIII': 'AAGCTT', 'MboI': 'GATC', 'DpnII': 'GATC', 'NcoI': 'CCATGG', 'Arima': 'GATC,GANTC'}


def _open(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def iter_fasta(path, wanted):
    name, seq = None, []
    with _open(path) as f:
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
    win_gc = np.add.reduceat(is_gc, np.arange(0, length, window)) / np.minimum(window, length - np.arange(0, length, window))
    per = resolution // window
    gc = np.zeros(math.floor(length / resolution) + 1)
    for i in range(len(gc)):
        chunk = win_gc[i * per:(i + 1) * per]
        gc[i] = chunk.mean() if len(chunk) else 0.0
    return gc


def sites_per_bin(seq, length, resolution, motifs):
    """motifs: comma-separated, N matches any base (e.g. 'GATC,GANTC' for the Arima kit)."""
    counts = np.zeros(math.floor(length / resolution) + 1)
    for motif in motifs.split(','):
        for m in re.finditer('(?=' + motif.upper().replace('N', '[ACGT]') + ')', seq):
            counts[min((m.start() + 1) // resolution, len(counts) - 1)] += 1   # 1-based positions, as in HiNT site files
    return counts


def map_per_bin(bedgraph, genome):
    res = genome.resolution
    out = {c: np.zeros(n) for c, n in genome.nbins.items()}
    with _open(bedgraph) as f:
        for line in f:
            c, s, e, v = line.split()[:4]
            if c not in out:
                continue
            s, e, v = int(s), min(int(e), genome.length[c]), float(v)
            while s < e:
                b = s // res
                stop = min(e, (b + 1) * res)
                out[c][b] += v * (stop - s)
                s = stop
    for c, l in genome.length.items():
        sizes = np.full(len(out[c]), res, dtype=float)
        sizes[-1] = l - (len(out[c]) - 1) * res or res
        out[c] /= sizes
    return out


def build_covariates(fasta, genome, motifs, mappability=None, out=None, log=print):
    mapp = map_per_bin(mappability, genome) if mappability else None
    rows = {}
    for name, seq in iter_fasta(fasta, set(genome.chroms)):
        length = genome.length[name]
        if len(seq) != length:
            raise SystemExit(f'{fasta}: {name} has {len(seq)} bp but {length} in genome {genome.name}')
        gc = gc_per_bin(seq, length, genome.resolution)
        rows[name] = (sites_per_bin(seq, length, genome.resolution, motifs), gc, mapp[name] if mapp else np.zeros(len(gc)))
        log(f'{name}: {len(gc)} bins')
    missing = [c for c in genome.chroms if c not in rows]
    if missing:
        raise SystemExit(f'{fasta}: chromosomes missing: {", ".join(missing)}')
    with open(out, 'w') as f:
        f.write('chrom\tbin_num\tcutSitesNumber\tgcper\tmapscore\n')
        for name in genome.chroms:
            sites, gc, mp = rows[name]
            for i in range(len(gc)):
                f.write(f'{name}\t{i}\t{sites[i]}\t{gc[i]}\t{mp[i]}\n')
    return out
