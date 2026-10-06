# Caki2 inputs (hg38, 50-kb bins)

Small files derived from public data, so that the Caki2 analysis (main-text Section 3.6, Table 7; Supplementary Table S18,
Supplementary Figure S5) runs without the 5.0 GB `.hic` download. The hg38 HindIII covariates shared with PANC-1 are in
`../covariates/hg38_HindIII_cov.txt`. Commands are run from `paper/`; `downloads/` holds the raw files.

| File | Content | Source | Rebuild |
|---|---|---|---|
| `Caki2_coverage_50kb.txt` | Hi-C coverage per bin: `rowsum_complete` (intra + inter), `rowsum_intra` | ENCODE in situ Hi-C, HindIII, experiment ENCSR401TBQ, file [ENCFF356GEH](https://www.encodeproject.org/files/ENCFF356GEH/) (`.hic`, GRCh38, 4,998,439,888 bytes, md5 `5b6b60ff438ad044b39299591175ab6a`) | `python scripts/make_coverage.py --hic downloads/ENCFF356GEH.hic -o data/caki2/Caki2_coverage_50kb.txt` |
| `depmap/Caki2_PR-g3gtC7_CNsegments.csv` | DepMap WGS copy-number segments of Caki2 (ACH-000234, profile PR-g3gtC7): the rows of this profile, unchanged | DepMap Public 24Q4, `OmicsCNSegmentsProfile.csv` ([figshare file 51065333](https://ndownloader.figshare.com/files/51065333), md5 `79fc5cf5c46f215ef085dddaa0e56860`); DepMap, Broad Institute, CC BY 4.0 | `{ head -1 downloads/OmicsCNSegmentsProfile.csv; grep '^PR-g3gtC7,' downloads/OmicsCNSegmentsProfile.csv; } > data/caki2/depmap/Caki2_PR-g3gtC7_CNsegments.csv` |

No array reference is used for Caki2. The DepMap reference seg file (log2(SegmentMean), floored at −5; calls |log2| ≥ 0.3;
193 calls) is written by `panc1_rerun.py --line caki2` into the run folder. Reference rules are given in Supplementary Note S1.
