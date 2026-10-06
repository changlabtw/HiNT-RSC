#!/usr/bin/env bash
# Download the raw public data behind the small derived inputs shipped in data/.
#
# NOT needed to reproduce the paper: every workflow (10-50, run_all.sh) runs on data/ alone. Use this script only to
# rebuild or check the shipped inputs, e.g. REBUILD_INPUTS=1 bash workflow/20_run_independent.sh.
#
# Usage (from anywhere; files go to paper/downloads/, or to $DOWNLOADS relative to paper/):
#   bash workflow/00_download_data.sh                  # default groups: hic arrays depmap hg38 annotations (~12 GB)
#   bash workflow/00_download_data.sh arrays depmap    # only some groups
#   bash workflow/00_download_data.sh all              # every group, including hg19 and k562 (~60 GB)
# Groups:
#   hic          ENCODE in situ Hi-C (HindIII, GRCh38 .hic): PANC-1 ENCFF545SGZ, Caki2 ENCFF356GEH         (~9.9 GB)
#                -> data/panc1/PANC1_coverage_50kb.txt, data/caki2/Caki2_coverage_50kb.txt
#   arrays       ENCODE PANC-1 genotyping-array CNV segments (hg19) ENCFF392GUL, ENCFF155BFQ + UCSC hg19ToHg38 chain
#                -> data/panc1/array/
#   depmap       DepMap Public 24Q4 OmicsCNSegmentsProfile.csv (doi:10.25452/figshare.plus.27993248.v1)     (41 MB)
#                -> data/panc1/depmap/, data/caki2/depmap/
#   hg38         UCSC hg38.fa.gz + Umap k50 multi-read mappability (bigWig, converted to k50.Umap.bedGraph.gz) (~2.3 GB)
#                -> data/covariates/hg38_HindIII_cov.txt
#   annotations  UCSC hg38 cytoBand and gap tables, ENCODE blacklist v2; compared with the copies in data/annotations/hg38/
#   hg19         UCSC hg19.fa.gz, for the covariate-builder check on hg19 (paper/README.md)                    (0.9 GB)
#   k562         GSE63525 (Rao et al. 2014) K562 and GM12878 MboI .hic files                                   (~47 GB)
#                -> source of the coverage in data/k562/*_bias_*.txt (see TODO below)
# Every file is checked against its MD5: from the ENCODE portal and figshare for those sources, and otherwise the MD5 of
# the copy used for the paper (downloaded 2026-10). Existing files with the right MD5 are skipped; partial downloads resume.
# Tools: curl, md5sum or md5, and bigWigToBedGraph for the hg38 group (env/environment.yml).
set -uo pipefail
cd "$(dirname "$0")/.."

DOWNLOADS=${DOWNLOADS:-downloads}
GROUPS_ARG=${*:-hic arrays depmap hg38 annotations}
[ "$GROUPS_ARG" = all ] && GROUPS_ARG="hic arrays depmap hg38 annotations hg19 k562"
mkdir -p "$DOWNLOADS"
status=0

ENCODE=https://www.encodeproject.org/files
UCSC=https://hgdownload.soe.ucsc.edu
GEO=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl

md5_of() { if command -v md5sum >/dev/null; then md5sum "$1" | cut -d' ' -f1; else md5 -q "$1"; fi; }

fetch() {  # fetch <file name under $DOWNLOADS> <url> <md5, or '-' if no checksum is recorded>
    local out="$DOWNLOADS/$1" url=$2 want=$3 got
    mkdir -p "$(dirname "$out")"
    if [ -f "$out" ] && { [ "$want" = - ] || [ "$(md5_of "$out")" = "$want" ]; }; then
        echo "present   $out"; return
    fi
    echo "download  $out  <-  $url"
    if ! curl -fL --retry 3 -C - -o "$out.part" "$url"; then
        echo "FAILED    $out" >&2; status=1; return
    fi
    mv "$out.part" "$out"
    got=$(md5_of "$out")
    if [ "$want" = - ]; then
        echo "md5       $out  $got  (no reference checksum)"
    elif [ "$got" = "$want" ]; then
        echo "md5 OK    $out"
    else
        echo "MD5 MISMATCH $out: got $got, expected $want (the file differs from the one used for the paper)" >&2
        status=1
    fi
}

same_as_shipped() {  # same_as_shipped <downloaded file> <shipped copy in data/>
    if cmp -s "$1" "$2"; then echo "identical to shipped copy: $2"; else echo "DIFFERS from shipped copy: $2" >&2; status=1; fi
}

for group in $GROUPS_ARG; do
    echo "== $group"
    case "$group" in
    hic)
        fetch ENCFF545SGZ.hic "$ENCODE/ENCFF545SGZ/@@download/ENCFF545SGZ.hic" bfd3c1813178319e1d839a63189d9918  # PANC-1, ENCSR440CTR
        fetch ENCFF356GEH.hic "$ENCODE/ENCFF356GEH/@@download/ENCFF356GEH.hic" 5b6b60ff438ad044b39299591175ab6a  # Caki2, ENCSR401TBQ
        ;;
    arrays)
        fetch ENCFF392GUL.bed.gz "$ENCODE/ENCFF392GUL/@@download/ENCFF392GUL.bed.gz" 938074031cc7340c7a080c96f1d512c4  # ENCSR000AZP
        fetch ENCFF155BFQ.bed.gz "$ENCODE/ENCFF155BFQ/@@download/ENCFF155BFQ.bed.gz" a9ba7dd87dedcc1ab0130181d45d4d2f  # replicate
        fetch hg19ToHg38.over.chain.gz "$UCSC/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz" 35887f73fe5e2231656504d1f6430900
        ;;
    depmap)
        fetch OmicsCNSegmentsProfile.csv https://ndownloader.figshare.com/files/51065333 79fc5cf5c46f215ef085dddaa0e56860
        ;;
    hg38)
        fetch hg38.fa.gz "$UCSC/goldenPath/hg38/bigZips/hg38.fa.gz" 1c9dcaddfa41027f17cd8f7a82c7293b
        fetch k50.Umap.MultiTrackMappability.bw "$UCSC/gbdb/hg38/hoffmanMappability/k50.Umap.MultiTrackMappability.bw" \
            722edf38316a518c1dc6e97b021c2028
        if [ -f "$DOWNLOADS/k50.Umap.bedGraph.gz" ]; then
            echo "present   $DOWNLOADS/k50.Umap.bedGraph.gz"
        elif command -v bigWigToBedGraph >/dev/null; then
            echo "convert   $DOWNLOADS/k50.Umap.bedGraph.gz"
            bigWigToBedGraph "$DOWNLOADS/k50.Umap.MultiTrackMappability.bw" /dev/stdout | gzip > "$DOWNLOADS/k50.Umap.bedGraph.gz.part" \
                && mv "$DOWNLOADS/k50.Umap.bedGraph.gz.part" "$DOWNLOADS/k50.Umap.bedGraph.gz" || status=1
        else
            echo "bigWigToBedGraph not found (conda package ucsc-bigwigtobedgraph); k50.Umap.bedGraph.gz not made" >&2; status=1
        fi
        ;;
    annotations)
        fetch annotations/hg38/cytoBand.txt.gz "$UCSC/goldenPath/hg38/database/cytoBand.txt.gz" 6677dd50caf19b0c3626b99ce67ca073
        fetch annotations/hg38/gap.txt.gz "$UCSC/goldenPath/hg38/database/gap.txt.gz" da739b566b6a0db5eb0dd6c8c206adb6
        fetch annotations/hg38/blacklist_v2.bed.gz \
            https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz 83fe6bf8187a64dee8079b80f75ba289
        for f in cytoBand.txt.gz gap.txt.gz blacklist_v2.bed.gz; do
            [ -f "$DOWNLOADS/annotations/hg38/$f" ] && same_as_shipped "$DOWNLOADS/annotations/hg38/$f" "data/annotations/hg38/$f"
        done
        ;;
    hg19)
        fetch hg19.fa.gz "$UCSC/goldenPath/hg19/bigZips/hg19.fa.gz" 806c02398f5ac5da8ffd6da2d1d5d1a9
        ;;
    k562)
        # TODO: no MD5 is recorded for these GEO files, and rebuilding data/k562/*_bias_*.txt from them (done in 2022 with
        # scripts/load.py) is not part of the release workflows; the K562 analyses start from the shipped coverage tables.
        fetch GSE63525_K562_combined.hic "$GEO/GSE63525_K562_combined.hic" -                          # 12,975,001,416 bytes
        fetch GSE63525_GM12878_insitu_primary.hic "$GEO/GSE63525_GM12878_insitu_primary.hic" -        # 34,210,807,977 bytes
        # TODO: the K562 WGS data (CCLE, GDC legacy archive; Control-FREEC sample G15509.K-562.2) are not downloaded: the
        # manuscript gives no accession for them. The WGS-derived profiles are shipped in data/k562/wgs/.
        ;;
    *)
        echo "unknown group: $group (use hic arrays depmap hg38 annotations hg19 k562 all)" >&2; exit 2 ;;
    esac
done

[ "$status" = 0 ] && echo "done: all files in $DOWNLOADS/ verified" || echo "done with problems (see messages above)" >&2
exit "$status"
