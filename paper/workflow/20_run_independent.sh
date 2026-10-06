#!/usr/bin/env bash
# Independent cell lines PANC-1 and Caki2 (main-text Section 3.6; Table 7; Supplementary Tables S16-S18).
#
#   1. PANC-1, pre-registered run without problem-region exclusion (array-state reference; DepMap WGS added post hoc)
#   2. PANC-1 problem-region diagnosis (Supplementary Note S1: 16 of 24 unmatched calls; 1,429 of 2,362 extreme bins)
#   3. PANC-1 and Caki2 with hg38 problem regions excluded: primary comparison and post hoc grid
#   4. comparison of every numeric column with expected/
#
# Inputs are the small derived files shipped in data/ (coverage, hg38 HindIII covariates, array and DepMap references,
# hg38 annotations), so no download is needed. Set REBUILD_INPUTS=1 to first rebuild those inputs from the raw downloads
# in $DOWNLOADS (default downloads/) into output/rebuilt_inputs/ and compare them with data/.
#
# Usage (from anywhere; paths are relative to paper/):
#   bash workflow/20_run_independent.sh
# Environment: PYTHON (python), RSCRIPT (Rscript), BICSEQ (NBICseq-seg.pl), WORKERS (2 parallel BIC-seq2 jobs),
#              REBUILD_INPUTS (0), DOWNLOADS (downloads), LIFTOVER (liftOver).
set -euo pipefail
cd "$(dirname "$0")/.."

PYTHON=${PYTHON:-python}
RSCRIPT=${RSCRIPT:-Rscript}
BICSEQ=${BICSEQ:-NBICseq-seg.pl}
WORKERS=${WORKERS:-2}
DOWNLOADS=${DOWNLOADS:-downloads}
LIFTOVER=${LIFTOVER:-liftOver}
export HINT_RSC_REF=data/ref_genome/hg38

if [ "${REBUILD_INPUTS:-0}" = 1 ]; then
    # Raw downloads: ENCFF545SGZ.hic (PANC-1), ENCFF356GEH.hic (Caki2), ENCFF392GUL.bed.gz and ENCFF155BFQ.bed.gz
    # (PANC-1 arrays), hg38.fa.gz, k50.Umap.bedGraph.gz, hg19ToHg38.over.chain.gz, OmicsCNSegmentsProfile.csv (DepMap 24Q4).
    RB=output/rebuilt_inputs
    mkdir -p "$RB"
    "$PYTHON" scripts/make_covariates.py --fasta "$DOWNLOADS/hg38.fa.gz" --enzyme HindIII \
        --mappability "$DOWNLOADS/k50.Umap.bedGraph.gz" -o "$RB/hg38_HindIII_cov.txt"
    "$PYTHON" scripts/make_coverage.py --hic "$DOWNLOADS/ENCFF545SGZ.hic" -o "$RB/PANC1_coverage_50kb.txt"
    "$PYTHON" scripts/make_coverage.py --hic "$DOWNLOADS/ENCFF356GEH.hic" -o "$RB/Caki2_coverage_50kb.txt"
    for a in ENCFF392GUL ENCFF155BFQ; do
        "$PYTHON" scripts/build_array_reference.py --bed "$DOWNLOADS/$a.bed.gz" --chain "$DOWNLOADS/hg19ToHg38.over.chain.gz" \
            --liftover "$LIFTOVER" -o "$RB/$a"
    done
    { head -1 "$DOWNLOADS/OmicsCNSegmentsProfile.csv"; grep '^PR-Tm9gPE,' "$DOWNLOADS/OmicsCNSegmentsProfile.csv"; } > "$RB/PANC1_PR-Tm9gPE_CNsegments.csv"
    { head -1 "$DOWNLOADS/OmicsCNSegmentsProfile.csv"; grep '^PR-g3gtC7,' "$DOWNLOADS/OmicsCNSegmentsProfile.csv"; } > "$RB/Caki2_PR-g3gtC7_CNsegments.csv"
    for pair in hg38_HindIII_cov.txt:data/covariates/hg38_HindIII_cov.txt \
                PANC1_coverage_50kb.txt:data/panc1/PANC1_coverage_50kb.txt \
                Caki2_coverage_50kb.txt:data/caki2/Caki2_coverage_50kb.txt \
                ENCFF392GUL.state_calls.seg:data/panc1/array/ENCFF392GUL.state_calls.seg \
                ENCFF392GUL.logR.seg:data/panc1/array/ENCFF392GUL.logR.seg \
                ENCFF155BFQ.state_calls.seg:data/panc1/array/ENCFF155BFQ.state_calls.seg \
                ENCFF155BFQ.logR.seg:data/panc1/array/ENCFF155BFQ.logR.seg \
                PANC1_PR-Tm9gPE_CNsegments.csv:data/panc1/depmap/PANC1_PR-Tm9gPE_CNsegments.csv \
                Caki2_PR-g3gtC7_CNsegments.csv:data/caki2/depmap/Caki2_PR-g3gtC7_CNsegments.csv; do
        if cmp -s "$RB/${pair%%:*}" "${pair#*:}"; then echo "identical: ${pair#*:}"; else echo "DIFFERS: ${pair#*:}"; fi
    done
fi

# 1. PANC-1, pre-registered: no problem-region exclusion (array-state reference), then the DepMap WGS reference
"$PYTHON" scripts/run_validation.py --outdir output/panc1/noexcl --bicseq "$BICSEQ" --rscript "$RSCRIPT" --workers "$WORKERS"
"$PYTHON" scripts/score_reference.py --outdir output/panc1/noexcl \
    --depmap-csv data/panc1/depmap/PANC1_PR-Tm9gPE_CNsegments.csv --profile PR-Tm9gPE

# 2. PANC-1 diagnosis of the unexcluded baseline (complete Hi-C, lambda = 3) against DepMap
"$PYTHON" scripts/diagnose_problem_regions.py --bins output/panc1/noexcl/complete/w0 \
    --ref output/panc1/noexcl/depmap_PR-Tm9gPE.seg | tee output/panc1/diagnosis.txt

# 3. PANC-1 and Caki2 with problem regions excluded (primary comparison and post hoc grid)
for line in panc1 caki2; do
    "$PYTHON" scripts/panc1_rerun.py --line "$line" --outdir "output/$line/excl" --bicseq "$BICSEQ" --rscript "$RSCRIPT" --workers "$WORKERS"
done
# after exclusion, no |log2| >= 1.5 bin of the PANC-1 baseline remains in a problem region
"$PYTHON" scripts/diagnose_problem_regions.py --bins output/panc1/excl/complete/none_0 \
    --seg output/panc1/excl/complete/none_0/bicseq_lambda3.seg --ref output/panc1/excl/depmap_PR-Tm9gPE.seg \
    | tee output/panc1/diagnosis_after_exclusion.txt

# 4. compare with the frozen results
"$PYTHON" - <<'EOF'
import numpy as np
import pandas as pd

pairs = [('output/panc1/noexcl/results.tsv', 'expected/panc1_noexcl_array.tsv', ['setting', 'w', 'lambda']),
         ('output/panc1/noexcl/results_depmap_PR-Tm9gPE.tsv', 'expected/panc1_noexcl_depmap.tsv', ['setting', 'w', 'lambda']),
         ('output/panc1/excl/results.tsv', 'expected/panc1_excluded_grid.tsv', ['setting', 'method', 'param', 'lambda']),
         ('output/caki2/excl/results.tsv', 'expected/caki2_excluded_grid.tsv', ['setting', 'method', 'param', 'lambda'])]
bad = 0
for new, ref, key in pairs:
    a, b = pd.read_csv(new, sep='\t'), pd.read_csv(ref, sep='\t')
    for d in (a, b):
        d['lambda'] = d['lambda'].astype(float)
    m = a.merge(b, on=key, how='outer', suffixes=('', '_exp'), indicator=True)
    cols = [c for c in b.columns if c not in key and c != 'seg']
    diff = [c for c in cols if not (m[c].astype(str) == m[c + '_exp'].astype(str)).all()]
    rows_ok = (m['_merge'] == 'both').all() and len(a) == len(b)
    num = [c for c in cols if pd.api.types.is_numeric_dtype(b[c])]
    print(f'{new}: {len(a)} rows, {len(num)} numeric + {len(cols) - len(num)} text columns; '
          + ('all identical' if rows_ok and not diff else f'DIFFERENT rows_ok={rows_ok} columns={diff}'))
    bad += (not rows_ok) or bool(diff)
raise SystemExit(1 if bad else 0)
EOF
