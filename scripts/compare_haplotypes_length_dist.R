#!/usr/bin/env Rscript
#
# Plot the per-tier length distribution (as proportions within each tier) of
# original-only haplotypes from the compute_haplotypes algorithm comparison
# on chr22. Shared is included as a reference baseline. Tiers 3 and 4 are
# omitted from the plot — both have n < 20 and don't read meaningfully as
# proportions.
#
# Counts are taken from `scripts/compare_haplotypes.py` output. Re-run that
# script and update the data block below if the underlying outputs change.
#
# Output: data/analysis/compute_haplotypes/length_dist_by_tier.png

suppressPackageStartupMessages({
  library(tidyverse)
})

out_dir <- "data/analysis/compute_haplotypes"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_path <- file.path(out_dir, "length_dist_by_tier.png")

# Long-format counts: one row per (category, length, count).
# Lengths >= 11 are collapsed into a single "11+" bin to keep the x-axis tidy.
counts <- tribble(
  ~category, ~length, ~count,
  "Tier 1 (fully variant-disjoint)",          2, 854,
  "Tier 1 (fully variant-disjoint)",          3,  64,
  "Tier 1 (fully variant-disjoint)",          4,   5,
  "Tier 1 (fully variant-disjoint)",          5,   2,
  "Tier 1 (fully variant-disjoint)",          6,   4,
  "Tier 1 (fully variant-disjoint)",          7,   1,
  "Tier 2 (partially missing variants)",      2, 266,
  "Tier 2 (partially missing variants)",      3, 105,
  "Tier 2 (partially missing variants)",      4,  24,
  "Tier 2 (partially missing variants)",      5,  11,
  "Tier 2 (partially missing variants)",      6,   2,
  "Tier 2 (partially missing variants)",      7,   4,
  "Tier 2 (partially missing variants)",      8,   2,
  "Tier 2 (partially missing variants)",     10,   1,
  "Tier 5 (sub-fragment of new)",     2, 911,
  "Tier 5 (sub-fragment of new)",     3, 200,
  "Tier 5 (sub-fragment of new)",     4,  65,
  "Tier 5 (sub-fragment of new)",     5,  31,
  "Tier 5 (sub-fragment of new)",     6,  22,
  "Tier 5 (sub-fragment of new)",     7,   9,
  "Tier 5 (sub-fragment of new)",     8,   6,
  "Tier 5 (sub-fragment of new)",     9,   4,
  "Tier 5 (sub-fragment of new)",    10,   1,
  "Tier 5 (sub-fragment of new)",    11,   1,
  "Shared (both algorithms)",                 2, 24103,
  "Shared (both algorithms)",                 3,  3137,
  "Shared (both algorithms)",                 4,   684,
  "Shared (both algorithms)",                 5,   208,
  "Shared (both algorithms)",                 6,    68,
  "Shared (both algorithms)",                 7,    34,
  "Shared (both algorithms)",                 8,    14,
  "Shared (both algorithms)",                 9,     7,
  "Shared (both algorithms)",                10,     5,
  "Shared (both algorithms)",                11,     4,
  "Shared (both algorithms)",                13,     1,
  "New-only (found only by new)",             2, 1203,
  "New-only (found only by new)",             3,  371,
  "New-only (found only by new)",             4,  142,
  "New-only (found only by new)",             5,   60,
  "New-only (found only by new)",             6,   33,
  "New-only (found only by new)",             7,   23,
  "New-only (found only by new)",             8,   11,
  "New-only (found only by new)",             9,    9,
  "New-only (found only by new)",            10,    8,
  "New-only (found only by new)",            11,    1,
  "New-only (found only by new)",            12,    7,
  "New-only (found only by new)",            13,    1,
  "New-only (found only by new)",            14,    2,
)

length_bin_levels <- c("2", "3", "4", "5", "6", "7", "8", "9", "10", "11+")

proportions <- counts %>%
  mutate(length_bin = if_else(length >= 11, "11+", as.character(length))) %>%
  group_by(category, length_bin) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  group_by(category) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup() %>%
  complete(
    category, length_bin = length_bin_levels,
    fill = list(count = 0, proportion = 0)
  ) %>%
  mutate(length_bin = factor(length_bin, levels = length_bin_levels))

# Order the legend: Shared and New-only first (the two baseline / outside-the-
# tier-partition categories), then the original-only tiers in numerical order.
category_levels <- c(
  "Shared (both algorithms)",
  "New-only (found only by new)",
  "Tier 1 (fully variant-disjoint)",
  "Tier 2 (partially missing variants)",
  "Tier 5 (sub-fragment of new)"
)
proportions <- proportions %>%
  mutate(category = factor(category, levels = category_levels))

# Okabe-Ito palette — color-blind friendly. Shared is muted grey; tiers and
# the new-only set get distinguishing colors.
palette <- c(
  "Shared (both algorithms)"            = "#999999",
  "New-only (found only by new)"        = "#009E73",
  "Tier 1 (fully variant-disjoint)"     = "#D55E00",
  "Tier 2 (partially missing variants)" = "#E69F00",
  "Tier 5 (sub-fragment of new)"        = "#0072B2"
)

p <- ggplot(
  proportions,
  aes(x = length_bin, y = proportion, color = category, group = category)
) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.4) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_color_manual(values = palette, guide = guide_legend(nrow = 2)) +
  labs(
    title = "chr22 haplotype length distribution, by tier",
    subtitle = "Proportion of haplotypes within each category, by number of variants",
    x = "Haplotype length (number of variants)",
    y = "Proportion within category",
    color = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 11),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey40"),
    panel.grid.minor = element_blank()
  )

ggsave(out_path, plot = p, width = 9, height = 5.5, dpi = 150)
cat(sprintf("wrote %s\n", out_path))
