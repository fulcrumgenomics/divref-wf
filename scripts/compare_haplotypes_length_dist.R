#!/usr/bin/env Rscript
#
# Plot the per-tier length distribution (as proportions within each tier) of
# original-only haplotypes from the compute_haplotypes algorithm comparison
# on chr22. Shared and New-only are included as reference baselines. Tiers
# 1-4 are omitted from the plot — after the min_variant_frequency formula
# correction they're all n < 5 and don't read meaningfully as proportions.
# The dominant original-only story is now Tier 5 (1,327 / 1,333 = 99.5%).
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
  "Tier 5 (sub-fragment of new)",     2, 967,
  "Tier 5 (sub-fragment of new)",     3, 208,
  "Tier 5 (sub-fragment of new)",     4,  73,
  "Tier 5 (sub-fragment of new)",     5,  32,
  "Tier 5 (sub-fragment of new)",     6,  22,
  "Tier 5 (sub-fragment of new)",     7,  11,
  "Tier 5 (sub-fragment of new)",     8,   8,
  "Tier 5 (sub-fragment of new)",     9,   4,
  "Tier 5 (sub-fragment of new)",    10,   1,
  "Tier 5 (sub-fragment of new)",    11,   1,
  "Shared (both algorithms)",                 2, 25168,
  "Shared (both algorithms)",                 3,  3306,
  "Shared (both algorithms)",                 4,   708,
  "Shared (both algorithms)",                 5,   221,
  "Shared (both algorithms)",                 6,    74,
  "Shared (both algorithms)",                 7,    37,
  "Shared (both algorithms)",                 8,    15,
  "Shared (both algorithms)",                 9,     8,
  "Shared (both algorithms)",                10,     6,
  "Shared (both algorithms)",                11,     4,
  "Shared (both algorithms)",                13,     1,
  "New-only (found only by new)",             2, 1067,
  "New-only (found only by new)",             3,  333,
  "New-only (found only by new)",             4,  128,
  "New-only (found only by new)",             5,   59,
  "New-only (found only by new)",             6,   30,
  "New-only (found only by new)",             7,   22,
  "New-only (found only by new)",             8,   12,
  "New-only (found only by new)",             9,   10,
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
  "Tier 5 (sub-fragment of new)"
)
proportions <- proportions %>%
  mutate(category = factor(category, levels = category_levels))

# Okabe-Ito palette — color-blind friendly. Shared is muted grey; tiers and
# the new-only set get distinguishing colors.
palette <- c(
  "Shared (both algorithms)"            = "#999999",
  "New-only (found only by new)"        = "#009E73",
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
