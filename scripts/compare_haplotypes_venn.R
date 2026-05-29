#!/usr/bin/env Rscript
#
# Render a Venn diagram comparing the haplotype catalogs produced by the
# original (two-pass-union) and new (per-sample adjacency) compute_haplotypes
# algorithms on chr22 (HGDP+1KG, AF >= 0.005 in at least one population, 25 bp
# window). Counts are sourced from `scripts/compare_haplotypes.py` output.
#
# Output: data/analysis/compute_haplotypes/algo_comparison.venn.png

suppressPackageStartupMessages({
  library(eulerr)
})

out_dir <- "data/analysis/compute_haplotypes"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_path <- file.path(out_dir, "algo_comparison.venn.png")

# Counts from compare_haplotypes.py on chr22 fixtures
venn_counts <- c(
  "Original" = 1333,
  "New" = 1680,
  "Original&New" = 29548
)

fit <- euler(venn_counts)

png(out_path, height = 900, width = 1200, res = 150)
plot(
  fit,
  quantities = list(cex = 1.1),
  fills = list(fill = c("#4e79a7", "#f28e2b"), alpha = 0.55),
  labels = list(fontfamily = "sans", cex = 1.2),
  main = "chr22 haplotypes: original vs new algorithm"
)
invisible(dev.off())

cat(sprintf("wrote %s\n", out_path))
