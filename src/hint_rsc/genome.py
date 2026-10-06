"""Genome definition: chromosome lengths and fixed-size bins, and the files bundled with the package."""
import math
import os
from importlib import resources

EXCLUDED_CHROMS = ('chrM', 'chrY', 'chrMT')
BUILTIN_GENOMES = ('hg19', 'hg38')
BUILTIN_COVARIATES = {('hg19', 'MboI'): 'hg19_MboI.tsv.gz', ('hg38', 'HindIII'): 'hg38_HindIII.tsv.gz'}
PROBLEM_REGION_FILES = ('cytoBand.txt.gz', 'gap.txt.gz', 'blacklist_v2.bed.gz')


def resource(*parts):
    """Path of a file bundled in hint_rsc/resources (or hint_rsc/R when parts[0] == 'R')."""
    root = resources.files('hint_rsc')
    if parts[0] != 'R':
        root = root / 'resources'
    for p in parts:
        root = root / p
    return str(root)


class Genome:
    """Chromosomes (chrY and chrM excluded, as in HiNT-CNV) and their bins at a fixed resolution.

    spec is a bundled genome name (hg19, hg38) or a chromosome-sizes file (chrom<TAB>length; order is kept).
    A chromosome of length L has floor(L / resolution) + 1 bins.
    """

    def __init__(self, spec, resolution=50000):
        if spec in BUILTIN_GENOMES:
            self.name, path = spec, resource('genomes', spec + '.len')
        elif os.path.exists(spec):
            self.name, path = os.path.basename(spec).split('.')[0], spec
        else:
            raise SystemExit(f'genome {spec!r}: not a bundled genome ({", ".join(BUILTIN_GENOMES)}) or a chromosome-sizes file')
        self.resolution = resolution
        self.length = {}
        with open(path) as f:
            for line in f:
                if not line.strip():
                    continue
                chrom, length = line.split('\t')[:2]
                if chrom not in EXCLUDED_CHROMS:
                    self.length[chrom] = int(length)
        self.nbins = {c: math.floor(l / resolution) + 1 for c, l in self.length.items()}

    @property
    def chroms(self):
        return list(self.length)

    def builtin_covariates(self, enzyme):
        f = BUILTIN_COVARIATES.get((self.name, enzyme))
        return resource('covariates', f) if f else None

    def builtin_problem_regions(self):
        if self.name != 'hg38':
            return None
        return [resource('problem_regions', 'hg38', f) for f in PROBLEM_REGION_FILES]
