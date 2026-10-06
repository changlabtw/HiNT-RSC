#!/usr/bin/env bash
# Full K562 HiNT-RSC run (hg19, MboI covariate, 50-kb bins): every configuration behind Tables 2-6,
# Supplementary Tables S1-S15, S19, S20 and Figure 3.
#
# Run from paper/:   bash workflow/10_run_k562.sh
# Environment (optional):
#   PYTHON   python with pandas, numpy, scipy, pywt, matplotlib      (default: python)
#   RSCRIPT  Rscript with mgcv, optparse, DNAcopy                    (default: Rscript)
#   BICSEQ   BIC-seq2 NBICseq-seg.pl, v0.7.2                          (default: NBICseq-seg.pl on PATH)
#   WORKERS  parallel segmentation jobs                               (default: 2)
#   OUTDIR   output directory                                         (default: output/k562)
# Output: $OUTDIR/results_MboI.tsv (one row per configuration; compare with expected/k562_MboI.tsv),
#         plus cached bins and segmentations under $OUTDIR/{gam,complete,intra,control,combined}.
# Runtime: about 15-20 min with 2 workers.
set -euo pipefail
cd "$(dirname "$0")/.."

PYTHON=${PYTHON:-python}
RSCRIPT=${RSCRIPT:-Rscript}
BICSEQ=${BICSEQ:-NBICseq-seg.pl}
WORKERS=${WORKERS:-2}
OUTDIR=${OUTDIR:-output/k562}
export HINT_RSC_REF=data/ref_genome/hg19
mkdir -p "$OUTDIR"

# 1. GAM fits, smoothing / control / intra-only grids, BIC-seq2 (lambda sweep) and default CBS; QC against HiNT's own output
"$PYTHON" scripts/k562_rerun.py --outdir "$OUTDIR" --bicseq "$BICSEQ" --rscript "$RSCRIPT" --workers "$WORKERS" | tee "$OUTDIR/run.log"
grep -q 'QC baseline .*: PASS' "$OUTDIR/run.log" || { echo 'QC FAILED: baseline differs from HiNT-CNV MboI segmentation' >&2; exit 1; }

# 2. Post-hoc CBS variant (smooth.CNA + sdundo; Supplementary Table S19)
"$PYTHON" scripts/cbs_variant.py --outdir "$OUTDIR" --rscript "$RSCRIPT" --workers "$WORKERS" | tee -a "$OUTDIR/run.log"

# 3. Mean-filter windows at lambda = 0.5 (Supplementary Table S4, Figure 3)
"$PYTHON" scripts/grid_extend.py --cell k562 --outdir "$OUTDIR" --bicseq "$BICSEQ" --workers "$WORKERS" --lambdas 0.5 | tee -a "$OUTDIR/run.log"

# 4. Regression check against the published results (every column except the 'seg' path, row-matched)
"$PYTHON" - "$OUTDIR/results_MboI.tsv" expected/k562_MboI.tsv <<'PY'
import sys
import pandas as pd
key = ['setting', 'method', 'param', 'segmenter', 'lambda']
new, ref = (pd.read_csv(f, sep='\t', dtype=str, keep_default_na=False) for f in sys.argv[1:3])
cols = [c for c in ref.columns if c != 'seg']
m = ref.merge(new, on=key, how='outer', suffixes=('_ref', '_new'), indicator=True)
missing = m[m._merge != 'both']
bad = [(tuple(r[k] for k in key), c, r[c + '_ref'], r[c + '_new']) for _, r in m[m._merge == 'both'].iterrows()
       for c in cols if c not in key and r[c + '_ref'] != r[c + '_new']]
print(f'regression: {len(ref)} expected rows, {len(new)} rerun rows, {len(missing)} unmatched, {len(bad)} differing cells')
for b in bad[:20]:
    print('  differs:', *b)
sys.exit(1 if len(missing) or bad else 0)
PY
echo "done: $OUTDIR/results_MboI.tsv (identical to expected/k562_MboI.tsv)"
