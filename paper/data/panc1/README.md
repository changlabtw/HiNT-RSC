# PANC-1 inputs (hg38, 50-kb bins)

Small files derived from public data, so that the PANC-1 analysis (main-text Section 3.6, Table 7, Figure 6;
Supplementary Tables S16–S17, Supplementary Figure S4) runs without the 4.9 GB `.hic` download.
The hg38 HindIII covariates shared with Caki2 are in `../covariates/hg38_HindIII_cov.txt`.
Commands are run from `paper/`; `downloads/` holds the raw files.

| File | Content | Source | Rebuild |
|---|---|---|---|
| `PANC1_coverage_50kb.txt` | Hi-C coverage per bin: `rowsum_complete` (intra + inter), `rowsum_intra` | ENCODE in situ Hi-C, HindIII, experiment ENCSR440CTR, file [ENCFF545SGZ](https://www.encodeproject.org/files/ENCFF545SGZ/) (`.hic`, GRCh38, 4,875,835,311 bytes, md5 `bfd3c1813178319e1d839a63189d9918`) | `python scripts/make_coverage.py --hic downloads/ENCFF545SGZ.hic -o data/panc1/PANC1_coverage_50kb.txt` |
| `array/ENCFF392GUL.state_calls.seg`, `array/ENCFF392GUL.logR.seg` | Array-state reference calls (145) and array logR segments, lifted to hg38 (411 of 418 segments) | ENCODE HAIB genotyping array, experiment ENCSR000AZP, file [ENCFF392GUL](https://www.encodeproject.org/files/ENCFF392GUL/) (hg19 bed, md5 `938074031cc7340c7a080c96f1d512c4`); UCSC `hg19ToHg38.over.chain.gz` | `python scripts/build_array_reference.py --bed downloads/ENCFF392GUL.bed.gz --chain downloads/hg19ToHg38.over.chain.gz -o data/panc1/array/ENCFF392GUL` |
| `array/ENCFF155BFQ.state_calls.seg`, `array/ENCFF155BFQ.logR.seg` | Second array segmentation of the same experiment (75 calls; 218 of 224 segments lifted) | ENCODE file [ENCFF155BFQ](https://www.encodeproject.org/files/ENCFF155BFQ/) (hg19 bed, md5 `a9ba7dd87dedcc1ab0130181d45d4d2f`) | as above with `ENCFF155BFQ` |
| `depmap/PANC1_PR-Tm9gPE_CNsegments.csv` | DepMap WGS copy-number segments of PANC-1 (ACH-000164, profile PR-Tm9gPE): the rows of this profile, unchanged | DepMap Public 24Q4, `OmicsCNSegmentsProfile.csv` ([figshare file 51065333](https://ndownloader.figshare.com/files/51065333), md5 `79fc5cf5c46f215ef085dddaa0e56860`); DepMap, Broad Institute, CC BY 4.0 | `{ head -1 downloads/OmicsCNSegmentsProfile.csv; grep '^PR-Tm9gPE,' downloads/OmicsCNSegmentsProfile.csv; } > data/panc1/depmap/PANC1_PR-Tm9gPE_CNsegments.csv` |

The DepMap reference seg file (log2(SegmentMean), floored at −5; calls |log2| ≥ 0.3; 441 calls) is written by
`score_reference.py` and `panc1_rerun.py` (function `score_reference.depmap_to_seg`) into the run folder.
Reference rules are given in Supplementary Note S1.

`workflow/20_run_independent.sh` with `REBUILD_INPUTS=1` rebuilds every file above from `downloads/` and compares it
with the shipped copy.
