# HiNT-RSC

HiNT-RSC is a residual-processing module for [HiNT-CNV](https://github.com/parklab/HiNT), the copy-number module of
HiNT (Wang *et al.*, *Genome Biology* 2020). HiNT-CNV turns a Hi-C contact map into a one-dimensional coverage profile,
normalizes it with a generalized additive model (GC content, mappability, restriction sites) and segments the residuals
with BIC-seq2. HiNT-RSC works after that normalization step. It smooths the residuals with a mean filter (window *w*).
It can also subtract a fraction of a control Hi-C residual, or build the profile from intra-chromosomal contacts only.
The residuals then go through HiNT's per-chromosome min-shift transform and BIC-seq2 segmentation, as in HiNT-CNV.
It adds no new statistical model or segmenter, and it is not a standalone CNV caller.

The accompanying paper finds that residual smoothing and the BIC-seq2 penalty λ jointly set segmentation granularity:
- On K562, at the default λ = 3, smoothing (w = 9) added matched calls, but lowering λ without smoothing reached a similar or better balance.
- The K562-selected setting did not transfer to PANC-1.
- The most consistent benefit of smoothing was a higher bin-level agreement with WGS copy-number profiles.

The recommendations are to tune λ, to exclude hg38 problem regions, and to smooth when a bin-level copy-number profile
is the goal.

## Status

| Version | Content |
|---|---|
| **v1.1.0** (this release) | Code, small derived inputs and frozen results that reproduce every table, figure and number of the APBC 2026 paper (`paper/`). |
| v1.2.0 (planned) | A packaged command-line tool for running HiNT-RSC on your own `.hic` data. |

Until v1.2.0, the scripts in `paper/scripts/` show the full pipeline. For a `.hic` file on hg38, the PANC-1/Caki2 path
is the closest template:
- `make_coverage.py` builds the coverage profile.
- `make_covariates.py` builds the covariates.
- `panc1_rerun.py` runs the GAM, smoothing, BIC-seq2 and scoring.

## Repository layout

```
paper/
  scripts/    all code (Python and R); every script is run from paper/
  workflow/   numbered shell workflows (00 download, 10-50 runs, run_all.sh)
  data/       small derived inputs (K562, PANC-1, Caki2, covariates, hg38 annotations, genome lengths)
  expected/   frozen result tables that the workflows compare against
  output/     created by the workflows (not tracked)
env/          conda environments (environment.yml, bicseq2.yml)
```

## Install

You need [conda](https://conda-forge.org/download/). The tested platform is macOS with conda packages for osx-64.

```bash
git clone https://github.com/changlabtw/HiNT-RSC.git && cd HiNT-RSC
conda env create -f env/environment.yml        # Python 3.11, R 4.5 (mgcv, DNAcopy), htslib, UCSC tools, hic-straw
conda env create -f env/bicseq2.yml            # BIC-seq2 v0.7.2 (kept separate: it pins its own R and Perl)
# Apple Silicon: prefix both commands with CONDA_SUBDIR=osx-64
```

## Quick start

```bash
conda activate hint-rsc
export BICSEQ="$(conda run -n bicseq2 which NBICseq-seg.pl)"
cd paper

# Tables from the frozen results (a few minutes, no segmentation)
bash workflow/40_k562_tables_figures.sh                    # Tables 1-6, Supplementary Tables S2-S15, S19, S20, Figure 3
EXPECTED_ONLY=1 bash workflow/50_independent_tables_figures.sh    # Table 7, Supplementary Tables S16-S18

# Full reproduction from the shipped inputs (about 2.5 h with 2 workers)
bash workflow/run_all.sh
```

`run_all.sh` compares every result with `paper/expected/` and stops at the first difference. See
[`paper/README.md`](paper/README.md) for the workflows, the provenance of every input, and a map from each manuscript
table, figure and number to the command and output file. No download is needed. `workflow/00_download_data.sh` fetches
the raw public data only if you want to rebuild the shipped inputs.

## Citation

If you use HiNT-RSC, please cite:

> Wei-Han Chen and Jia-Ming Chang. HiNT-RSC: Residual Smoothing and Segmentation Granularity in Hi-C-based Copy Number
> Variation Detection. APBC 2026 (manuscript under revision).

Please also cite HiNT:

> Su Wang *et al.* HiNT: a computational method for detecting copy number variations and translocations from Hi-C data.
> *Genome Biology* 21, 73 (2020). https://doi.org/10.1186/s13059-020-01986-5

Machine-readable metadata are in [`CITATION.cff`](CITATION.cff).

## License and attribution

HiNT-RSC is released under the [MIT License](LICENSE).

It builds on HiNT ([parklab/HiNT](https://github.com/parklab/HiNT)), also MIT-licensed
(Copyright (c) 2019 Su Wang and contributors). The permission notice is the same as in [`LICENSE`](LICENSE). From HiNT:
- The functions `getmappability`, `getGCpercent` and `getFragmentsNumber` in `paper/scripts/load.py` are taken from the HiNT source code.
- The GAM model (`paper/scripts/GAM.R`) and the residual min-shift transform (`paper/scripts/binlib.py`) follow HiNT-CNV.
- `paper/data/k562/hint_MboI/` and `paper/data/covariates/hg19_MboI_cov.txt` are output and input files of a HiNT-CNV run.

Third-party data in `paper/data/` keep their own terms. The DepMap extracts are from DepMap Public 24Q4 (Broad
Institute, CC BY 4.0). The other sources (ENCODE, UCSC, ENCODE blacklist, GEO GSE63525) are listed in
[`paper/README.md`](paper/README.md).

## Earlier versions

This work is adapted from the 2022 master's thesis of Wei-Han Chen (National Chengchi University). The thesis notebooks
and experiment folders, which are not part of the paper, are no longer on the main branch. They remain available in
the git history and under the tag `v1.0.0` (the repository as it stood when the paper was submitted).
