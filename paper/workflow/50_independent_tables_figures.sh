#!/usr/bin/env bash
# Independent cell lines: main-text Table 7, Supplementary Tables S16-S18, Figure 6 (PANC-1 grid),
# Supplementary Figure S4 (PANC-1 problem regions) and Supplementary Figure S5 (Caki2 grid).
#
# Needs the outputs of workflow/20_run_independent.sh. The tables are also built from the frozen results in expected/
# (output/tables_expected/), and the two sets are compared; set EXPECTED_ONLY=1 to build only the tables from expected/
# without any rerun (the figures need the segmentations and are then skipped).
#
# Usage (from anywhere; paths are relative to paper/):
#   bash workflow/50_independent_tables_figures.sh
# Environment: PYTHON (python), EXPECTED_ONLY (0).
set -euo pipefail
cd "$(dirname "$0")/.."

PYTHON=${PYTHON:-python}
export HINT_RSC_REF=data/ref_genome/hg38

"$PYTHON" scripts/make_independent_tables.py --expected --outdir output/tables_expected
if [ "${EXPECTED_ONLY:-0}" = 1 ]; then
    exit 0
fi

# Table 7, Supplementary Tables S16-S18 from the rerun, compared with the tables from expected/
"$PYTHON" scripts/make_independent_tables.py --outdir output/tables
status=0
for f in Table7_independent_lines.tsv TableS16_PANC1_no_exclusion.tsv TableS17_PANC1_grid.tsv TableS18_Caki2_grid.tsv; do
    if cmp -s "output/tables/$f" "output/tables_expected/$f"; then echo "identical to expected: $f"; else echo "DIFFERS from expected: $f"; status=1; fi
done

# Figure 6 (PANC-1) and Supplementary Figure S5 (Caki2): accuracy vs number of calls over the lambda and window grids
mkdir -p output/figures
"$PYTHON" scripts/plot_grid.py --cell panc1 --results output/panc1/excl/results.tsv -o output/figures/Figure6_PANC1_grid.png
"$PYTHON" scripts/plot_grid.py --cell caki2 --results output/caki2/excl/results.tsv -o output/figures/FigureS5_Caki2_grid.png

# Supplementary Figure S4: PANC-1 profile before and after problem-region exclusion
"$PYTHON" scripts/plot_panc1_problem_regions.py -o output/figures/FigureS4_PANC1_problem_regions.png
exit $status
