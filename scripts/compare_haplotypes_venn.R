#!/usr/bin/env Rscript
#
# Render a Venn diagram comparing the haplotype catalogs produced by the
# original (two-pass-union) and new (per-sample adjacency) compute_haplotypes
# algorithms on the autosomes (chr1-22) (HGDP+1KG, AF >= 0.005 in at least one population,
# 25 bp window). Counts are read from the summary TSV emitted by
# `scripts/compare_haplotypes.py` so the figure stays in sync with the data.
#
# Output: data/analysis/compute_haplotypes/algo_comparison.venn.png

suppressPackageStartupMessages({
  library(eulerr)
  library(optparse)
})

option_list <- list(
  make_option("--summary",
    type = "character",
    default = "data/analysis/compute_haplotypes/algo_comparison.summary.tsv",
    help = "Path to algo_comparison.summary.tsv [default: %default]"
  ),
  make_option("--output",
    type = "character",
    default = "data/analysis/compute_haplotypes/algo_comparison.venn.png",
    help = "Output PNG path [default: %default]"
  )
)
opts <- parse_args(OptionParser(option_list = option_list))

summary_df <- read.table(opts$summary, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
get_metric <- function(name) {
  row <- summary_df[summary_df$metric == name, "value"]
  if (length(row) != 1) {
    stop(sprintf("expected exactly one row for metric '%s' in %s, got %d",
                 name, opts$summary, length(row)))
  }
  as.integer(row)
}

n_old_only <- get_metric("old_only")
n_new_only <- get_metric("new_only")
n_shared <- get_metric("shared")

venn_counts <- c(
  "Original" = n_old_only,
  "New" = n_new_only,
  "Original&New" = n_shared
)

fit <- euler(venn_counts)

dir.create(dirname(opts$output), showWarnings = FALSE, recursive = TRUE)
png(opts$output, height = 900, width = 1200, res = 150)
plot(
  fit,
  quantities = list(cex = 1.1),
  fills = list(fill = c("#4e79a7", "#f28e2b"), alpha = 0.55),
  labels = list(fontfamily = "sans", cex = 1.2),
  main = "Autosome haplotypes: original vs new algorithm"
)
invisible(dev.off())

cat(sprintf("wrote %s\n", opts$output))
