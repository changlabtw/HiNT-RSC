#!/usr/bin/env bash
# Additional analyses: Supplementary Table S1 (bin-level agreement with WGS BIC-seq2 bins, Control-FREEC and CNVpytor),
# the sequencing-depth experiment (Supplementary Table S21, Supplementary Figure S6) and the event-type counts quoted
# in the Discussion (large low-amplitude vs focal reference events).
#
# Run from paper/:   bash workflow/30_analyses.sh
# Requires the K562 run (workflow/10_run_k562.sh -> output/k562) and, for the event-type counts, the PANC-1 and Caki2
# runs with problem regions excluded (workflow/20_run_independent.sh -> output/{panc1,caki2}/excl); the segmentations
# are found through the 'seg' column of their results tables.
# Environment (optional):
#   PYTHON         python with pandas, numpy, scipy, pywt, matplotlib            (default: python)
#   RSCRIPT        Rscript with mgcv, optparse                                    (default: Rscript)
#   BICSEQ         BIC-seq2 NBICseq-seg.pl, v0.7.2                                 (default: NBICseq-seg.pl on PATH)
#   WORKERS        parallel GAM fits / segmentation jobs                          (default: 2)
#   K562_RESULTS   K562 results table                                             (default: output/k562/results_MboI.tsv)
#   PANC1_RESULTS  PANC-1 results table (problem regions excluded)                (default: output/panc1/excl/results.tsv)
#   CAKI2_RESULTS  Caki2 results table (problem regions excluded)                 (default: output/caki2/excl/results.tsv)
#   SKIP_DEPTH=1   reuse output/depth/results_depth.tsv instead of rerunning the depth experiment
# Output: output/analyses/{profile_references.tsv, depth_per_seed.tsv, tableS21_depth.tsv, event_types.tsv},
#         output/depth/ (thinned tables, GAM fits, segmentations, results_depth.tsv, noise_depth.tsv),
#         output/figures/figS6_depth.{png,pdf}; the three tables with a released copy are compared with expected/.
# Runtime: the depth experiment is about 25 GAM fits and 800 BIC-seq2 runs (roughly 1 h with 2 workers); the rest
# takes a few minutes.
set -euo pipefail
cd "$(dirname "$0")/.."

PYTHON=${PYTHON:-python}
RSCRIPT=${RSCRIPT:-Rscript}
BICSEQ=${BICSEQ:-NBICseq-seg.pl}
WORKERS=${WORKERS:-2}
K562_RESULTS=${K562_RESULTS:-output/k562/results_MboI.tsv}
PANC1_RESULTS=${PANC1_RESULTS:-output/panc1/excl/results.tsv}
CAKI2_RESULTS=${CAKI2_RESULTS:-output/caki2/excl/results.tsv}
export HINT_RSC_REF=data/ref_genome/hg19
mkdir -p output/analyses output/figures

status=0
check() {  # check <output> <expected>
    if cmp -s "$1" "$2"; then echo "PASS  $1 == $2"; else echo "DIFF  $1 != $2"; status=1; fi
}

# 1. Supplementary Table S1: bin-level SCC / PCC of the residual profiles vs three WGS references
"$PYTHON" scripts/profile_references.py --results "$K562_RESULTS" \
    --out output/analyses/profile_references.tsv --long output/analyses/profile_references_all.tsv
check output/analyses/profile_references.tsv expected/profile_references.tsv

# 2. Sequencing depth (Supplementary Table S21, Supplementary Figure S6): binomial thinning of the K562 MboI coverage,
#    p = 1 ... 0.0001, seeds 1 and 2, unsmoothed and mean filter w = 3, 5, 9, lambda 0.5-5
if [ "${SKIP_DEPTH:-0}" != 1 ]; then
    "$PYTHON" scripts/depth_downsampling.py --bias output/k562/gam/complete/bias.txt --outdir output/depth \
        --fractions 1,0.5,0.2,0.1,0.05,0.02,0.01,0.005,0.002,0.001,0.0005,0.0002,0.0001 --seeds 1,2 \
        --windows 3,5,9 --lambdas 0.5,1,1.5,2,2.5,3,4,5 \
        --bicseq "$BICSEQ" --rscript "$RSCRIPT" --workers "$WORKERS"
fi
"$PYTHON" scripts/depth_summary.py --results output/depth/results_depth.tsv --noise output/depth/noise_depth.tsv \
    --per-seed output/analyses/depth_per_seed.tsv --table output/analyses/tableS21_depth.tsv
check output/analyses/depth_per_seed.tsv expected/depth_per_seed.tsv
"$PYTHON" scripts/figS6_depth.py --per-seed output/analyses/depth_per_seed.tsv --out output/figures/figS6_depth

# 3. Discussion: reference events matched by type (K562, PANC-1, Caki2; Caki2 chr20 excluded)
"$PYTHON" scripts/event_types.py --k562-results "$K562_RESULTS" --panc1-results "$PANC1_RESULTS" \
    --caki2-results "$CAKI2_RESULTS" --seg-root output --out output/analyses/event_types.tsv
check output/analyses/event_types.tsv expected/event_types.tsv

[ "$status" = 0 ] && echo "done: all compared tables match expected/" || echo "done: some tables differ from expected/ (see DIFF lines)"
exit "$status"
