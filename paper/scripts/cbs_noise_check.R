# Noise diagnostics behind the CBS discussion (Results 3.7, Supplementary Table S19): on one chromosome of the
# K562 run, the spread of the Hi-C log2 ratio (SD), the CBS-style noise estimate from first differences
# (mad(diff(x)) / sqrt(2), the quantity DNAcopy's undo.SD uses) and the lag-1 autocorrelation, for the unsmoothed
# residual and the mean-filtered residuals. CBS assumes independent noise between bins, i.e. lag-1 autocorrelation ~ 0.
#
# Usage (from paper/, after workflow/10_run_k562.sh):  Rscript scripts/cbs_noise_check.R [run_dir] [chrom]
#   run_dir  k562_rerun.py output directory (default output/k562)
#   chrom    chromosome (default chr6)
args <- commandArgs(trailingOnly = TRUE)
run <- if (length(args) >= 1) args[1] else "output/k562"
chrom <- if (length(args) >= 2) args[2] else "chr6"
cat(sprintf("%-8s %-6s %8s %24s %12s\n", "input", "chrom", "SD(x)", "mad(diff(x))/sqrt(2)", "lag-1 acf"))
for (cfg in c("none_0", paste0("mean_", c(3, 5, 7, 9, 11, 21, 41)))) {
  d <- read.table(sprintf("%s/complete/%s/_%s.bin", run, cfg, chrom), header = TRUE)
  ok <- d$obs > 0 & d$exp > 0
  x <- log2(d$obs[ok] / d$exp[ok])
  cat(sprintf("%-8s %-6s %8.3f %24.4f %12.2f\n", cfg, chrom, sd(x), mad(diff(x)) / sqrt(2), cor(x[-1], x[-length(x)])))
}
