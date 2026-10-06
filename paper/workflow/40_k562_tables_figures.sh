#!/usr/bin/env bash
# K562 tables and figures (final manuscript numbering).
#
# Run from paper/:   bash workflow/40_k562_tables_figures.sh
#   Tables 1-6, Supplementary Tables S1 (columns 2-3), S2-S15, S19, S20 and Figure 3 are built from RESULTS
#   (default expected/k562_MboI.tsv, so no rerun is needed). Figures 2, 4, 5, Supplementary Figures S1, S3 and the
#   CBS noise diagnostics (Results 3.7) need the bins of a K562 run (workflow/10_run_k562.sh) and are skipped otherwise.
#   Supplementary Figure S2 (WGS callers) is built with Table 1.
# Environment (optional):
#   PYTHON   (default: python)        RSCRIPT (default: Rscript)
#   RESULTS  results TSV              (default: expected/k562_MboI.tsv; use output/k562/results_MboI.tsv for a fresh run)
#   RUN      K562 run directory       (default: output/k562)
#   OUTDIR   output root              (default: output; tables -> $OUTDIR/tables/k562, figures -> $OUTDIR/figures)
set -euo pipefail
cd "$(dirname "$0")/.."

PYTHON=${PYTHON:-python}
RSCRIPT=${RSCRIPT:-Rscript}
RESULTS=${RESULTS:-expected/k562_MboI.tsv}
RUN=${RUN:-output/k562}
OUTDIR=${OUTDIR:-output}
export HINT_RSC_REF=data/ref_genome/hg19

# Table 1 + Supplementary Figure S2
"$PYTHON" scripts/table1_wgs_callers.py --outdir "$OUTDIR"

# Tables 2-6, Supplementary Tables S1 (cols 2-3), S2-S15, S19, S20
"$PYTHON" scripts/make_tables.py --results "$RESULTS" --outdir "$OUTDIR/tables/k562"

# Figure 3 (lambda sweep vs smoothing window)
"$PYTHON" scripts/plot_grid.py --cell k562 --results "$RESULTS" -o "$OUTDIR/figures/figure3_k562_grid.png"

# Figures 2, 4, 5, Supplementary Figures S1, S3 and CBS noise diagnostics: need the bins of a K562 run
if [ -f "$RUN/complete/mean_9/_chr6.bin" ] && [ -f "$RUN/intra/none_0/_chr6.bin" ] && [ -f "$RUN/gam/control/residual.txt" ]; then
  "$PYTHON" scripts/make_figures.py --run "$RUN" --outdir "$OUTDIR/figures"
  "$RSCRIPT" scripts/cbs_noise_check.R "$RUN" chr6 | tee "$OUTDIR/tables/k562/cbs_noise_chr6.txt"
else
  echo "skipped Figures 2, 4, 5, S1, S3 and the CBS noise check: no K562 run in $RUN (run workflow/10_run_k562.sh first)" >&2
fi
echo "done: $OUTDIR/tables/k562, $OUTDIR/figures"
