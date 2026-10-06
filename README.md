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
| **v1.2.0** (this release) | The `hint-rsc` command-line tool for your own `.hic` data (`src/hint_rsc/`), with tests that check it against the paper. |
| v1.1.0 | Code, small derived inputs and frozen results that reproduce every table, figure and number of the APBC 2026 paper (`paper/`). |

## Install

You need [conda](https://conda-forge.org/download/). The tested platform is macOS with conda packages for osx-64.

```bash
git clone https://github.com/changlabtw/HiNT-RSC.git && cd HiNT-RSC
conda env create -f env/environment.yml        # Python 3.11, R 4.5 (mgcv, DNAcopy), htslib, UCSC tools, hic-straw
conda env create -f env/bicseq2.yml            # BIC-seq2 v0.7.2 (kept separate: it pins its own R and Perl)
# Apple Silicon: prefix both commands with CONDA_SUBDIR=osx-64
conda activate hint-rsc
pip install --no-deps -e .                     # installs the hint-rsc command
export BICSEQ="$(conda run -n bicseq2 which NBICseq-seg.pl)"
hint-rsc check                                 # R with mgcv, BIC-seq2, Python modules
```

## Using HiNT-RSC on your own data

The input is a `.hic` file. HiNT-RSC builds the 50-kb coverage profile and fits HiNT-CNV's GAM, then segments with BIC-seq2.
Covariates are bundled for hg19 + MboI and hg38 + HindIII. Problem regions are bundled for hg38: centromeres, assembly gaps
and the ENCODE blacklist v2.

```bash
# 1. HiNT-CNV baseline (no smoothing, lambda = 3); coverage is cached as out/coverage.tsv (a few minutes)
hint-rsc run --hic sample.hic --genome hg38 --enzyme HindIII -o out/

# 2. Tune the segmentation granularity: every window x lambda, with call counts (and agreement with a reference, if any)
hint-rsc sweep --coverage out/coverage.tsv --genome hg38 --enzyme HindIII --workers 4 -o out/sweep
#    optional: --reference wgs_segments.bed  (chrom, start, end, log2 ratio)

# 3. Final run with the chosen setting
hint-rsc run --coverage out/coverage.tsv --genome hg38 --enzyme HindIII --lambda 1 -o out/lambda1
```

`run` writes the following files to the output folder:

| File | Content |
|---|---|
| `segments.seg` | BIC-seq2 segments (chrom, start, end, bins, observed, expected, log2 ratio, p-value) |
| `cnv_calls.bed` | segments with \|log2 ratio\| ≥ 0.3 (`--threshold`), labelled gain or loss |
| `profile.tsv` | bin-level copy-number profile (residual and log2 ratio per GAM-fitted bin), mean-filtered with `--profile-window` (default 9; 0 = none) |
| `metrics.tsv` | with `--reference`: segment- and bin-level Pearson/Spearman correlation, matched calls, precision and recall |
| `run.json` | all parameters and summary counts |

The defaults follow the paper's recommendations:
- **Tune λ.** The default λ = 3 (HiNT-CNV's default) under-segmented K562. The best unsmoothed λ was 0.5–1.5 across lines, so
  explore that range with `sweep`.
- **Exclude hg38 problem regions.** This is on by default for hg38 (`--problem-regions auto`). For other genomes, give BED files or
  `none`.
- **Smooth when a bin-level copy-number profile is the goal.** For example, the profile can be the input to CNV-aware Hi-C
  normalization. Smoothing raised bin-level agreement with WGS, most for shallow libraries. Before segmentation, smoothing
  (`-w/--window`, `--smoother`) is off by default: once λ was tuned, it did not improve CNV calls.

Other options:
- `--contacts intra` uses intra-chromosomal contacts only. It is faster to compute, at some cost in segment-level agreement.
- `--control-hic control.hic --control-fraction 0.3` subtracts a fraction of a control sample's residual.
- `--smoother` selects one of `mean`, `gaussian`, `median`, `savgol`, `haar` or `sym4`.

For another genome or enzyme, build the covariates once and pass them with `--covariates`:

```bash
bigWigToBedGraph k50.Umap.MultiTrackMappability.bw /dev/stdout | gzip > k50.Umap.bedGraph.gz
hint-rsc covariates --fasta hg38.fa.gz --genome hg38 --enzyme DpnII --mappability k50.Umap.bedGraph.gz -o hg38_DpnII.tsv
```

`--genome` also accepts a chromosome-sizes file (chrom, length). `hint-rsc coverage` and `hint-rsc evaluate` run the
coverage and scoring steps on their own. See `hint-rsc <command> --help`.

### Tests

```bash
pip install pytest
pytest tests/          # unit tests; with R (mgcv) and BIC-seq2 available, also the paper regression (about 2 min)
```

The regression tests run `hint-rsc` on the shipped PANC-1 and Caki2 inputs and compare every metric with
`paper/expected/`. With `HINT_RSC_TEST_HIC=/path/to/ENCFF545SGZ.hic`, they also compare the coverage computed from the
`.hic` file with the shipped coverage table. In addition, the following were checked byte for byte against the paper's
outputs:
- PANC-1 and Caki2: coverage computed from the ENCODE `.hic` files, the GAM residuals, the BIC-seq2 input bins and the
  segmentations.
- K562: the intra-only GAM and segmentation, and the control correction when given the paper's covariate values.

## Reproducing the paper

The paper's own scripts are kept unchanged in `paper/`.

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

## Repository layout

```
src/hint_rsc/   the hint-rsc package (CLI, pipeline, bundled genomes, covariates, hg38 problem regions, R/GAM.R)
tests/          unit tests and the regression against paper/expected/
paper/
  scripts/      the paper's code (Python and R); every script is run from paper/
  workflow/     numbered shell workflows (00 download, 10-50 runs, run_all.sh)
  data/         small derived inputs (K562, PANC-1, Caki2, covariates, hg38 annotations, genome lengths)
  expected/     frozen result tables that the workflows compare against
  output/       created by the workflows (not tracked)
env/            conda environments (environment.yml, bicseq2.yml)
```

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
- The GAM model (`paper/scripts/GAM.R`, `src/hint_rsc/R/GAM.R`) and the residual min-shift transform
  (`paper/scripts/binlib.py`, `src/hint_rsc/segment.py`) follow HiNT-CNV.
- `paper/data/k562/hint_MboI/` and `paper/data/covariates/hg19_MboI_cov.txt` (bundled as
  `src/hint_rsc/resources/covariates/hg19_MboI.tsv.gz`) are output and input files of a HiNT-CNV run.

Third-party data in `paper/data/` and `src/hint_rsc/resources/` keep their own terms. The DepMap extracts are from DepMap Public 24Q4 (Broad
Institute, CC BY 4.0). The other sources (ENCODE, UCSC, ENCODE blacklist, GEO GSE63525) are listed in
[`paper/README.md`](paper/README.md).

## Earlier versions

This work is adapted from the 2022 master's thesis of Wei-Han Chen (National Chengchi University). The thesis notebooks
and experiment folders, which are not part of the paper, are no longer on the main branch. They remain available in
the git history and under the tag `v1.0.0` (the repository as it stood when the paper was submitted).
