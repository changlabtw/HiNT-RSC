# Circular binary segmentation (DNAcopy) of a HiNT-RSC bin directory, as an alternative to BIC-seq2 (Supplementary Table S19).
# Usage (from paper/): Rscript scripts/cbs.R <bin_dir> <out.seg> [alpha] [sdundo]
#   called by k562_rerun.py (default CBS, alpha = 0.01) and cbs_variant.py (post-hoc "sdundo" variant).
# With "sdundo": smooth.CNA() outlier smoothing, then undo.splits = "sdundo", undo.SD = 3 (standard practice for noisy data).
# Input bins: start end obs exp(ected) var (one file per chromosome, _chr*.bin); log2 ratio = log2(obs/exp).
# Output: BIC-seq2-like seg (chrom start end binNum observed expected log2.copyRatio), 0-based bin coordinates.
suppressMessages(library(DNAcopy))
args <- commandArgs(trailingOnly = TRUE)
bin_dir <- args[1]; out <- args[2]
alpha <- if (length(args) >= 3) as.numeric(args[3]) else 0.01
sdundo <- length(args) >= 4 && args[4] == "sdundo"
files <- list.files(bin_dir, pattern = "^_chr.*\\.bin$", full.names = TRUE)
chrom <- c(); pos <- c(); val <- c()
for (f in files) {
  cname <- sub("^_", "", sub("\\.bin$", "", basename(f)))
  d <- read.table(f, header = TRUE)
  expcol <- if ("exp" %in% names(d)) d$exp else d$expected
  keep <- d$obs > 0 & expcol > 0
  chrom <- c(chrom, rep(cname, sum(keep)))
  pos <- c(pos, d$start[keep])
  val <- c(val, log2(d$obs[keep] / expcol[keep]))
}
cna <- CNA(val, chrom, pos, data.type = "logratio", sampleid = "hic")
set.seed(1)
if (sdundo) {
  cna <- smooth.CNA(cna)
  seg <- segment(cna, alpha = alpha, undo.splits = "sdundo", undo.SD = 3, verbose = 0)$output
} else {
  seg <- segment(cna, alpha = alpha, verbose = 0)$output
}
res <- data.frame(chrom = seg$chrom, start = seg$loc.start, end = seg$loc.end + 49999, binNum = seg$num.mark,
                  observed = NA, expected = NA, log2.copyRatio = seg$seg.mean)
write.table(res, out, sep = "\t", quote = FALSE, row.names = FALSE)
