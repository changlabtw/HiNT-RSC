# Reproducing the HiNT-RSC paper

This folder holds everything needed to reproduce the tables, figures and numbers of *HiNT-RSC: Residual Smoothing and
Segmentation Granularity in Hi-C-based Copy Number Variation Detection* (W.-H. Chen and J.-M. Chang, APBC 2026). It
uses the final numbering of the manuscript:
- main text: Tables 1–7, Figures 1–6;
- supplement: Supplementary Tables S1–S21, Supplementary Figures S1–S6, Supplementary Note S1.

## Conventions

- **Run everything from `paper/`.** The workflows `cd` there themselves, so `bash paper/workflow/run_all.sh` also works
  from the repository root. All paths below are relative to `paper/`.
- `scripts/` holds all code. `data/` holds the small derived inputs (28 MB). `expected/` holds the frozen result tables.
  `output/` is created by the runs and is not tracked (about 1.2 GB after a full run).
- Environment (see the root README): `conda activate hint-rsc`, and BIC-seq2 v0.7.2 from `env/bicseq2.yml`.
- Every workflow reads these optional variables:

  | Variable | Default | Meaning |
  |---|---|---|
  | `PYTHON` | `python` | Python interpreter |
  | `RSCRIPT` | `Rscript` | R interpreter |
  | `BICSEQ` | `NBICseq-seg.pl` on `PATH` | BIC-seq2 segmenter |
  | `WORKERS` | `2` | parallel GAM fits and segmentations |

- The genome is chosen with `HINT_RSC_REF`, which is a folder under `data/ref_genome/`. The workflows set it to hg19 for
  K562 and to hg38 for PANC-1 and Caki2. When you run a K562 script by hand, use `HINT_RSC_REF=data/ref_genome/hg19`.
- The residual transform is HiNT-CNV's per-chromosome min-shift: obs = r − min(r) and exp = −min(r), applied after any
  smoothing or control correction (`scripts/binlib.py`).

## Run order

| Step | Command | What it does | Time* | Built-in check |
|---|---|---|---|---|
| 00 | `bash workflow/00_download_data.sh [groups]` | Optional. Downloads the raw public data behind `data/` (see [Data provenance](#data-provenance)). | network | MD5 of every file |
| 10 | `bash workflow/10_run_k562.sh` | K562: GAM fits, every smoothing, control, intra-only, filter and λ configuration, BIC-seq2 and CBS. Output: `output/k562/results_MboI.tsv`. | ~30 min | Baseline (unsmoothed, λ = 3) identical to HiNT-CNV's own MboI segmentation. Every cell equal to `expected/k562_MboI.tsv`. |
| 20 | `bash workflow/20_run_independent.sh` | PANC-1 pre-registered run without exclusion, DepMap re-score, problem-region diagnosis, then PANC-1 and Caki2 runs with problem regions excluded. | ~30 min | Every column equal to `expected/panc1_*.tsv` and `expected/caki2_excluded_grid.tsv`. |
| 30 | `bash workflow/30_analyses.sh` | Supplementary Table S1 (all columns), depth experiment (Supplementary Table S21, Supplementary Figure S6), event-type counts. Needs steps 10 and 20. | ~40 min | Output byte-identical to `expected/profile_references.tsv`, `depth_per_seed.tsv` and `event_types.tsv`. |
| 40 | `bash workflow/40_k562_tables_figures.sh` | K562 tables and figures. By default it builds from `expected/k562_MboI.tsv`; use `RESULTS=output/k562/results_MboI.tsv` for a fresh run. Figures 2, 4, 5 and S1, S3 need step 10. | minutes | — |
| 50 | `bash workflow/50_independent_tables_figures.sh` | Table 7, Supplementary Tables S16–S18, Figure 6, Supplementary Figures S4 and S5. Needs step 20; with `EXPECTED_ONLY=1` it builds the tables from `expected/` only. | minutes | Tables from the rerun byte-identical to the tables from `expected/`. |
| all | `bash workflow/run_all.sh` | Steps 10 → 50; `STEPS="40 50"` selects steps. Logs go to `output/run_all_step<N>.txt`. | ~2.5 h | Stops at the first failed check. |

\* With 2 workers on a laptop (macOS, osx-64 conda environment).

Other options:
- `SKIP_DEPTH=1` (step 30) reuses `output/depth/`.
- `REBUILD_INPUTS=1` (step 20) rebuilds every shipped PANC-1/Caki2 input from `downloads/` and compares it with `data/`.

## Data provenance

Nothing in `data/` needs to be downloaded to run the workflows. The **Rebuild** column says how each file can be
regenerated from the raw data that `workflow/00_download_data.sh` fetches into `downloads/`.

| Path in `data/` | Content | Source | Rebuild |
|---|---|---|---|
| `ref_genome/hg19/hg19.len`, `ref_genome/hg38/hg38.len` | chromosome lengths, chr1–22 and chrX | UCSC `hg19.chrom.sizes`, `hg38.chrom.sizes` (chrY and chrM dropped) | — |
| `k562/K562_bias_complete.txt` | K562 per-bin coverage (intra + inter contacts) with HiNT's hg19 covariates (GC, mappability, HindIII sites) | Hi-C: GEO [GSE63525](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525) `GSE63525_K562_combined.hic` (Rao *et al.* 2014, MboI). Built in 2022 with `scripts/load.py`. | not scripted (see the TODO in `00_download_data.sh`, group `k562`) |
| `k562/K562_bias_intra.txt` | the same table from intra-chromosomal contacts only | as above | as above |
| `k562/K562_bias_control_target.txt`, `k562/GM12878_bias_control.txt` | K562 and GM12878 tables for the control-correction model | as above, plus `GSE63525_GM12878_insitu_primary.hic` | as above |
| `covariates/hg19_MboI_cov.txt` | MboI restriction-site counts, GC and mappability per hg19 bin. `scripts/k562_rerun.py` uses its `cutSitesNumber` in place of the HindIII column of the K562 tables. | HiNT-CNV's own MboI regression input for K562 | GC and cut sites: `scripts/make_covariates.py` (see the [check](#numbers-without-a-dedicated-script)) |
| `k562/hint_MboI/` | HiNT-CNV's own K562 MboI run: residual bins (`b50000/`) and segmentation (λ = 3). Used as the complete-Hi-C residual and for the QC of step 10. | HiNT-CNV output | — |
| `k562/wgs/bicseq2/` | K562 WGS BIC-seq2 segmentations at λ = 3, 4, 5 and 50-kb bins (λ = 4 is the primary reference) | CCLE K562 WGS (GDC legacy archive), computed in 2022 | — (no accession in the paper; see TODO in `00_download_data.sh`) |
| `k562/wgs/freec/G15509.K-562.2.bam_ratio.txt` | Control-FREEC ratio profile | same WGS data | — |
| `k562/wgs/cnvpytor/` | CNVpytor 50-kb calls | same WGS data | — |
| `covariates/hg38_HindIII_cov.txt` | hg38 HindIII sites, GC, Umap k50 mappability per bin | UCSC `hg38.fa.gz`, Umap k50 (`k50.Umap.MultiTrackMappability.bw`) | `scripts/make_covariates.py` (step 20, `REBUILD_INPUTS=1`) |
| `panc1/`, `caki2/` | Hi-C coverage per bin, PANC-1 array references, DepMap extracts | ENCODE ENCFF545SGZ, ENCFF356GEH, ENCFF392GUL, ENCFF155BFQ; DepMap 24Q4 | see [`data/panc1/README.md`](data/panc1/README.md) and [`data/caki2/README.md`](data/caki2/README.md) |
| `annotations/hg38/` | cytoBand (centromeres), gap, ENCODE blacklist v2 | UCSC, Boyle lab | unchanged copies, see [`data/annotations/README.md`](data/annotations/README.md) |

Accessions, as given in the manuscript's Data Availability statement:

| Data | Accessions |
|---|---|
| K562 and GM12878 Hi-C | GSE63525 |
| PANC-1 Hi-C | ENCSR440CTR / ENCFF545SGZ |
| Caki2 Hi-C | ENCSR401TBQ / ENCFF356GEH |
| PANC-1 array | ENCSR000AZP (ENCFF392GUL, ENCFF155BFQ) |
| DepMap 24Q4 | [doi:10.25452/figshare.plus.27993248.v1](https://doi.org/10.25452/figshare.plus.27993248.v1): PANC-1 ACH-000164 / PR-Tm9gPE; Caki2 ACH-000234 / PR-g3gtC7 |

The download script checks every file against its MD5. ENCODE and figshare publish these checksums. For the other
files, the checksum is that of the copy used for the paper.

## Frozen results (`expected/`)

| File | Content | Produced by |
|---|---|---|
| `k562_MboI.tsv` | one row per K562 configuration (95 rows: setting, filter, parameter, segmenter, λ; all metrics) | step 10 |
| `panc1_noexcl_array.tsv`, `panc1_noexcl_depmap.tsv` | PANC-1 without exclusion: array references (pre-registered) and DepMap | step 20 (`run_validation.py`, `score_reference.py`) |
| `panc1_excluded_grid.tsv`, `caki2_excluded_grid.tsv` | PANC-1 and Caki2 with problem regions excluded: λ × window grid, alternative filters, intra-only; `absent_chromosomes` lists chromosomes missing from a segmentation | step 20 (`panc1_rerun.py`) |
| `profile_references.tsv` | Supplementary Table S1 | step 30 (`profile_references.py`) |
| `depth_per_seed.tsv` | per-depth, per-seed values behind Supplementary Table S21 and Supplementary Figure S6 | step 30 (`depth_downsampling.py`, `depth_summary.py`) |
| `event_types.tsv` | reference events matched by type (Discussion) | step 30 (`event_types.py`) |

The `seg` column holds the path of the segmentation file inside `output/`, so it only resolves after a run.

## Where each manuscript item comes from

Outputs are under `output/`. In the step column, "40 (needs 10)" means that the item comes from step 40 but needs the
outputs of step 10.

### Main text

| Item | Script | Step | Output |
|---|---|---|---|
| Figure 1 (workflow overview) | — (schematic, drawn by hand) | — | — |
| Table 1 (WGS caller agreement) | `table1_wgs_callers.py` | 40 | `tables/k562/table1.tsv` |
| Figure 2 (chr6, smoothing) | `make_figures.py` | 40 (needs 10) | `figures/figure2_chr6_smoothing.png` |
| Table 2 (smoothing, λ = 3) | `make_tables.py` | 40 | `tables/k562/table2.tsv` |
| Figure 3 (K562, λ sweep vs window) | `plot_grid.py --cell k562` | 40 | `figures/figure3_k562_grid.png` |
| Table 3 (smoothing methods) | `make_tables.py` | 40 | `tables/k562/table3.tsv` |
| Table 4 (control correction) | `make_tables.py` | 40 | `tables/k562/table4.tsv` |
| Figure 4 (chr6, control) | `make_figures.py` | 40 (needs 10) | `figures/figure4_chr6_control.png` |
| Table 5 (complete vs intra-only) | `make_tables.py` | 40 | `tables/k562/table5.tsv` (parsing-time column: see below) |
| Figure 5 (chr6, intra-only) | `make_figures.py` | 40 (needs 10) | `figures/figure5_chr6_intra.png` |
| Table 6 (control + smoothing) | `make_tables.py` | 40 | `tables/k562/table6.tsv` |
| Table 7 (PANC-1, Caki2) | `make_independent_tables.py` | 50 | `tables/Table7_independent_lines.tsv` (from `expected/`: `tables_expected/`) |
| Figure 6 (PANC-1 grid) | `plot_grid.py --cell panc1` | 50 (needs 20) | `figures/Figure6_PANC1_grid.png` |

### Supplementary material

| Item | Script | Step | Output |
|---|---|---|---|
| Table S1 (bin level, three WGS references) | `profile_references.py` (all columns); `make_tables.py` (BIC-seq2 columns only) | 30 (needs 10); 40 | `analyses/profile_references.tsv`; `tables/k562/tableS1.tsv` |
| Tables S2–S15 (size strata, consistency, λ × window grid, filters, control, intra-only, combined) | `make_tables.py` | 40 | `tables/k562/tableS2.tsv` … `tableS15.tsv` |
| Table S16 (PANC-1 without exclusion) | `make_independent_tables.py` | 50 | `tables/TableS16_PANC1_no_exclusion.tsv` |
| Table S17 (PANC-1 grids) | `make_independent_tables.py` | 50 | `tables/TableS17_PANC1_grid.tsv` |
| Table S18 (Caki2 grids) | `make_independent_tables.py` | 50 | `tables/TableS18_Caki2_grid.tsv` |
| Table S19 (CBS) | `k562_rerun.py` + `cbs_variant.py` (runs), `make_tables.py` (table) | 10, 40 | `tables/k562/tableS19.tsv` |
| Table S20 (matching rule, reference λ) | `make_tables.py` | 40 | `tables/k562/tableS20.tsv` |
| Table S21 (sequencing depth) | `depth_downsampling.py`, `depth_summary.py` | 30 (needs 10) | `analyses/tableS21_depth.tsv` |
| Figure S1 (chr5, chr6: WGS vs HiNT-CNV) | `make_figures.py` | 40 (needs 10) | `figures/figureS1_chr5_chr6.png` |
| Figure S2 (WGS callers, A–C) | `table1_wgs_callers.py` | 40 | `figures/figureS2a_wgs_callers.png`, `…S2b…`, `…S2c…` |
| Figure S3 (chr6, window sizes) | `make_figures.py` | 40 (needs 10) | `figures/figureS3_chr6_windows.png` |
| Figure S4 (PANC-1 problem regions) | `plot_panc1_problem_regions.py` | 50 (needs 20) | `figures/FigureS4_PANC1_problem_regions.png` |
| Figure S5 (Caki2 grid) | `plot_grid.py --cell caki2` | 50 (needs 20) | `figures/FigureS5_Caki2_grid.png` |
| Figure S6 (sequencing depth) | `figS6_depth.py` | 30 | `figures/figS6_depth.png`, `.pdf` (redraw from `--per-seed expected/depth_per_seed.tsv` without a rerun) |

### Numbers quoted in the text

| Number | Where | Source |
|---|---|---|
| Reproduced HiNT-CNV segmentation exactly (λ = 3) | Section 2.4, Note S1 | step 10: `QC baseline … PASS` in `output/k562/run.log` |
| Covariate builder on hg19: GC identical; MboI r = 1.000000, HindIII r = 0.999995 (99.96 % of bins identical) | Note S1 | commands below |
| 5,160 of 60,630 bins excluded; 16 of 24 unmatched calls (mean 64 %) vs 1 of 16 matched calls (mean 4 %) in problem regions; 1,429 of 2,362 bins with \|log2\| ≥ 1.5; none after exclusion | Section 3.6, Note S1, Discussion | step 20: `diagnose_problem_regions.py` → `output/panc1/diagnosis.txt`, `diagnosis_after_exclusion.txt` |
| Baseline precision against DepMap 0.400 → 0.625 after exclusion | Section 3.6, Discussion | Tables S16 and 7 |
| Reference calls: DepMap PANC-1 441 (197 gains, 244 losses), Caki2 193 (82, 111); array-state 145 (15, 130), replicate 75; 411 of 418 segments lifted | Note S1 | `*_ref` columns of the results; lift-over counts are logged by `build_array_reference.py` (step 20, `REBUILD_INPUTS=1`); gain/loss split: command below |
| Array-state replicate agreement: bin-level SCC 0.764; precision 0.52, recall 0.27 | Section 3.6, Discussion | command below |
| Best unsmoothed λ (PANC-1 0.796 DepMap, 0.670 array; Caki2 0.903); light windows 0.713 (array); intra-only values | Section 3.6, Discussion | Tables S17, S18 |
| Chromosomes absent from segmentations (†; Caki2 chr20) | Section 3.6, Note S1 | `absent_chromosomes` column of `expected/*_excluded_grid.tsv`, Tables S17 and S18 |
| Hi-C residual autocorrelation on chr6 (lag 1): 0.67, and 0.94–1.00 after smoothing | Section 3.7 | `cbs_noise_check.R`, step 40 (needs 10) → `output/tables/k562/cbs_noise_chr6.txt` |
| Event types: w = 9 at λ = 3 matched 25 of 60 large low-amplitude events vs 17 unsmoothed; unsmoothed λ = 1: 25; no extra focal events | Section 2.9, Discussion | `event_types.py`, step 30 → `output/analyses/event_types.tsv` |
| Depth: bin-level SCC gain of w = 9 from 0.018 (full depth) to 0.117 (0.1 % of contacts); best tuned segment-level SCC nearly unchanged down to 0.1 % | Section 2.9, Discussion | Table S21; `expected/depth_per_seed.tsv` |
| Highest recall in the λ × window grid 0.581 (precision ≤ 0.500); 25 of 74 K562 reference calls < 1 Mb | Discussion | Tables S4, S2 |
| Contact parsing time ~1 hour (complete) vs ~6 minutes (intra-only) | Abstract, Section 3.5.2, Table 5 | not recomputed (see below) |
| HiNT's published PANC-1 values (bin-level SCC 0.865; 84 % of calls > 2 Mb supported) | Section 3.6 | Wang *et al.*, *Genome Biology* 2020 (not computed here) |

## Numbers without a dedicated script

**Parsing time (Table 5, Abstract, Section 3.5.2).** The time to extract coverage from the K562 `.hic` file, about
1 hour for all contacts and about 6 minutes for intra-chromosomal contacts only, was measured once, on a single CPU core
with 8 GB RAM, in the 2022 thesis (Chen 2022). `make_tables.py` writes it as fixed text.

**HiNT's PANC-1 results (Section 3.6)** are quoted from the HiNT paper for a qualitative comparison only.

**Covariate builder on hg19 (Supplementary Note S1).** The MboI check is reproducible from the shipped data plus
`hg19.fa.gz`:

```bash
bash workflow/00_download_data.sh hg19
python scripts/make_covariates.py --fasta downloads/hg19.fa.gz --len data/ref_genome/hg19/hg19.len \
    --enzyme MboI -o output/hg19_MboI_rebuilt.txt
python - <<'EOF'
import numpy as np, pandas as pd
a = pd.read_csv('output/hg19_MboI_rebuilt.txt', sep='\t')
b = pd.read_csv('data/covariates/hg19_MboI_cov.txt', sep='\t')        # HiNT's MboI input
m = a.merge(b, on=['chrom', 'bin_num'], suffixes=('_new', '_hint'))
print(len(m), 'bins; cut sites r = %.6f, identical in %.2f %% of bins; GC max |difference| = %.1e' % (
    np.corrcoef(m.cutSitesNumber_new, m.cutSitesNumber_hint)[0, 1],
    100 * (m.cutSitesNumber_new == m.cutSitesNumber_hint).mean(), (m.gcper_new - m.gcper_hint).abs().max()))
EOF
# 60739 bins; cut sites r = 1.000000, identical in 99.96 % of bins; GC max |difference| = 4.1e-09 (floating-point rounding)
```

The HindIII check (r = 0.999995, 99.96 %) runs the same comparison with `--enzyme HindIII` against HiNT's hg19 HindIII
inputs. Those inputs are not shipped in `data/`. They are kept in the `v1.0.0` tree as
`CNV/HiNTcnv_OUTPUT/K562_combined_dataForRegression/regressionInfo_<chrom>_K562_combined.txt`, with columns bin, cut
sites, GC, mappability and coverage.

**PANC-1 array-state replicate agreement (Section 3.6, Discussion):**

```bash
PYTHONPATH=scripts HINT_RSC_REF=data/ref_genome/hg38 python - <<'EOF'
from load import load_seg_to_bin
from metric import calculate_dict_correlation, calculate_precision_racall
a, b = 'data/panc1/array/ENCFF392GUL', 'data/panc1/array/ENCFF155BFQ'
calculate_dict_correlation(load_seg_to_bin(b + '.logR.seg'), load_seg_to_bin(a + '.logR.seg'))  # Spearman 0.764
calculate_precision_racall(b + '.state_calls.seg', groundtruth_file=a + '.state_calls.seg')     # precision 0.52, recall 0.269
EOF
```

**Gain/loss split of the reference calls (Supplementary Note S1),** after step 20:

```bash
for f in output/panc1/excl/depmap_PR-Tm9gPE.seg output/caki2/excl/depmap_PR-g3gtC7.seg \
         data/panc1/array/ENCFF392GUL.state_calls.seg data/panc1/array/ENCFF155BFQ.state_calls.seg; do
  awk -F'\t' -v f="$f" 'NR > 1 { if ($7 >= 0.3) g++; else if ($7 <= -0.3) l++ }
                        END { print f ": " g + l " calls, " g " gains, " l " losses" }' "$f"
done
```
