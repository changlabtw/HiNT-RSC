#!/usr/bin/env bash
# Reproduce every result of the paper from the inputs shipped in data/ (no download needed), in this order:
#   10  K562 runs (Tables 2-6, Supplementary Tables S1-S15, S19, S20)          ~30 min
#   20  PANC-1 and Caki2 runs (Table 7, Supplementary Tables S16-S18)          ~30 min
#   30  Supplementary Table S1 (all columns), depth (S21, Figure S6), event types  ~40 min
#   40  K562 tables and figures, built from the fresh K562 run                  minutes
#   50  PANC-1 / Caki2 tables and figures                                       minutes
# Steps 10, 20, 30 and 50 compare their results with expected/ and exit non-zero on any difference, which stops the
# run. Runtimes are for 2 parallel workers on a laptop. Each step's log is kept in output/run_all_step<N>.txt.
#
# Usage (from anywhere; paths are relative to paper/):
#   bash workflow/run_all.sh
#   STEPS="40 50" bash workflow/run_all.sh        # only some steps (each step needs the outputs of earlier ones)
# Environment (passed to every step):
#   PYTHON   (python)   RSCRIPT (Rscript)   BICSEQ (NBICseq-seg.pl, BIC-seq2 v0.7.2)   WORKERS (2)
#   SKIP_DEPTH=1        reuse output/depth/ in step 30 instead of rerunning the depth experiment
#   REBUILD_INPUTS=1    in step 20, also rebuild the shipped PANC-1/Caki2 inputs from downloads/ and compare them with
#                       data/ (needs workflow/00_download_data.sh first)
set -euo pipefail
cd "$(dirname "$0")/.."

export PYTHON=${PYTHON:-python}
export RSCRIPT=${RSCRIPT:-Rscript}
export BICSEQ=${BICSEQ:-NBICseq-seg.pl}
export WORKERS=${WORKERS:-2}
# the steps use OUTDIR / RESULTS / RUN with different meanings; never let a value from the caller leak into them
unset OUTDIR RESULTS RUN
STEPS=${STEPS:-10 20 30 40 50}

mkdir -p output
for step in $STEPS; do
    case "$step" in
        10) cmd=(bash workflow/10_run_k562.sh) ;;
        20) cmd=(bash workflow/20_run_independent.sh) ;;
        30) cmd=(bash workflow/30_analyses.sh) ;;
        40) cmd=(env RESULTS=output/k562/results_MboI.tsv bash workflow/40_k562_tables_figures.sh) ;;
        50) cmd=(bash workflow/50_independent_tables_figures.sh) ;;
        *) echo "unknown step: $step (use 10 20 30 40 50)" >&2; exit 2 ;;
    esac
    echo "=== step $step: ${cmd[*]}  ($(date '+%Y-%m-%d %H:%M:%S'))"
    start=$SECONDS
    "${cmd[@]}" 2>&1 | tee "output/run_all_step$step.txt"
    echo "=== step $step done in $(( (SECONDS - start) / 60 )) min"
done
echo "all steps passed: tables in output/tables/, figures in output/figures/, analyses in output/analyses/"
