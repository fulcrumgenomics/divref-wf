#!/usr/bin/env Rscript
#
# Plot the per-category distribution of estimated_gnomad_AF for the original
# vs new compute_haplotypes algorithm comparison on chr22. Reads the TSV
# produced by `scripts/compare_haplotypes.py`. Each haplotype contributes one
# point per category.
#
# For Shared and New-only categories the est_af value comes from the new
# algorithm. For Tiers 1-5 (original-only) the value comes from the original
# algorithm.
#
# Output: data/analysis/compute_haplotypes/est_af_dist_by_tier.png

suppressPackageStartupMessages({
  library(tidyverse)
})

in_path <- "data/analysis/compute_haplotypes/algo_comparison.est_af.tsv"
out_dir <- "data/analysis/compute_haplotypes"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_path <- file.path(out_dir, "est_af_dist_by_tier.png")

af <- read_tsv(in_path, show_col_types = FALSE) %>%
  filter(!is.na(est_af) & est_af > 0)

category_levels <- c(
  "Shared (both algorithms)",
  "New-only (found only by new)",
  "Tier 1 (fully variant-disjoint)",
  "Tier 2 (partially missing variants)",
  "Tier 3 (split across new)",
  "Tier 4 (subset, not contiguous)",
  "Tier 5 (sub-fragment of new)"
)
af <- af %>%
  mutate(category = factor(category, levels = category_levels))

# Per-category n, for legend annotations.
ns <- af %>%
  count(category) %>%
  mutate(label = sprintf("%s (n=%d)", category, n))
label_map <- setNames(ns$label, ns$category)
af <- af %>%
  mutate(category_with_n = factor(label_map[as.character(category)],
                                  levels = label_map[category_levels]))

palette <- c(
  "Shared (both algorithms)"             = "#999999",
  "New-only (found only by new)"         = "#009E73",
  "Tier 1 (fully variant-disjoint)"      = "#D55E00",
  "Tier 2 (partially missing variants)"  = "#E69F00",
  "Tier 3 (split across new)"            = "#CC79A7",
  "Tier 4 (subset, not contiguous)"      = "#56B4E9",
  "Tier 5 (sub-fragment of new)"         = "#0072B2"
)
# Re-map palette names to the n-annotated labels.
palette_with_n <- setNames(palette[category_levels], label_map[category_levels])

p <- ggplot(af, aes(x = category_with_n, y = est_af, fill = category_with_n)) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.55,
    width = 0.55,
    color = "grey25"
  ) +
  # Overlay individual points for the small-n tiers (3 and 4) so they are
  # legible; jitter for visual separation. Larger categories use boxplots only.
  geom_jitter(
    data = af %>% filter(grepl("Tier 3|Tier 4", as.character(category))),
    aes(color = category_with_n),
    width = 0.15,
    height = 0,
    size = 2,
    alpha = 0.8
  ) +
  geom_hline(
    yintercept = 0.005,
    linetype = "dashed",
    color = "grey30"
  ) +
  annotate(
    "text",
    x = 0.6,
    y = 0.005,
    label = "0.5% threshold",
    vjust = -0.5,
    hjust = 0,
    size = 3.5,
    color = "grey30"
  ) +
  scale_y_log10(
    breaks = c(0.005, 0.01, 0.05, 0.1, 0.5, 1.0),
    labels = scales::percent_format(accuracy = 0.1)
  ) +
  scale_fill_manual(values = palette_with_n) +
  scale_color_manual(values = palette_with_n) +
  labs(
    title = "chr22 estimated gnomAD AF distribution, by category",
    subtitle = paste(
      "est_af source: new algorithm for Shared and New-only;",
      "original algorithm for Tiers 1-5. Boxes show median + IQR; whiskers extend",
      "to 1.5*IQR.",
      sep = " "
    ),
    x = NULL,
    y = "estimated gnomAD AF (log scale)",
    fill = NULL,
    color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 25, hjust = 1, size = 10),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey40", size = 10),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 10, r = 20, b = 20, l = 60)
  )

ggsave(out_path, plot = p, width = 11, height = 6.5, dpi = 150)
cat(sprintf("wrote %s\n", out_path))
