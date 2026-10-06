# hg38 problem-region annotations

Used for the problem-region exclusion of PANC-1 and Caki2 (Supplementary Note S1): a 50-kb bin is excluded when at least
25 % of it lies in the union of these regions (5,160 of 60,630 bins). Files are unchanged copies of the sources.

| File | Content | Source |
|---|---|---|
| `hg38/cytoBand.txt.gz` | cytogenetic bands; only `acen` (centromere) bands are used | UCSC, https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz |
| `hg38/gap.txt.gz` | assembly gaps | UCSC, https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gap.txt.gz |
| `hg38/blacklist_v2.bed.gz` | ENCODE blacklist v2 (Amemiya, Kundaje and Boyle, Sci. Rep. 2019) | https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz |
