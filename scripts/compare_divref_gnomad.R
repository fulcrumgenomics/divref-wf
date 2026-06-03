#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(duckdb)
  library(duckplyr)
  library(eulerr)
  library(logger)
  library(optparse)
  library(tidyverse)
})

option_list <- list(
  make_option("--contig",
    type = "character", default = "chr22",
    help = "Contig to filter from DivRef index [default: %default]"
  ),
  make_option("--divref_duckdb",
    type = "character",
    default = "data/resources/DivRef-v1.1.haplotypes_gnomad_merge.index.duckdb",
    help = "Path to DivRef DuckDB index [default: %default]"
  ),
  make_option("--gnomad_tsv",
    type = "character",
    help = "Path to gnomAD allele frequency TSV"
  ),
  make_option("--gnomad_label",
    type = "character",
    help = "Label for gnomAD data in plot axes and titles"
  ),
  make_option("--output_base",
    type = "character",
    help = "Base path for output files (suffixes added per file)"
  )
)

opts <- parse_args(OptionParser(option_list = option_list))

populations <- c("afr", "amr", "eas", "sas", "nfe")

# Load DivRef data ----

con <- dbConnect(duckdb(), read_only = TRUE, dbdir = opts$divref_duckdb)

divref <- dbGetQuery(
  con,
  "SELECT * FROM sequences WHERE contig = ? AND source = 'gnomAD_variant'",
  params = list(opts$contig)
) %>%
  select(-c(sequence, sequence_length, sequence_id, n_variants, source, contig)) %>%
  rename(divref_maxpop = max_pop) %>%
  mutate(
    popmax_empirical_AN = if_else(
      !is.na(popmax_empirical_AF) & popmax_empirical_AF > 0,
      ceiling(popmax_empirical_AC / popmax_empirical_AF),
      NA
    ),
    gnomAD_AF_afr = as.numeric(gnomAD_AF_afr),
    gnomAD_AF_amr = as.numeric(gnomAD_AF_amr),
    gnomAD_AF_eas = as.numeric(gnomAD_AF_eas),
    gnomAD_AF_sas = as.numeric(gnomAD_AF_sas),
    gnomAD_AF_nfe = as.numeric(gnomAD_AF_nfe),
  ) %>%
  complete(fill = list(
    "gnomAD_AF_afr" = 0.0, "gnomAD_AF_amr" = 0.0, "gnomAD_AF_eas" = 0.0,
    "gnomAD_AF_sas" = 0.0, "gnomAD_AF_nfe" = 0.0
  ))

dbDisconnect(con)

log_info("Loaded ", nrow(divref), " DivRef variants for ", opts$contig)

# Load gnomAD data ----

gnomad <- read_tsv(opts$gnomad_tsv, show_col_types = FALSE) %>%
  rename(gnomad_maxpop = maxpop) %>%
  mutate(
    afr = as.numeric(afr),
    amr = as.numeric(amr),
    eas = as.numeric(eas),
    sas = as.numeric(sas),
    nfe = as.numeric(nfe)
  ) %>%
  complete(fill = list("afr" = 0.0, "amr" = 0.0, "eas" = 0.0, "sas" = 0.0, "nfe" = 0.0))

log_info("Loaded ", nrow(gnomad), " gnomAD variants for ", opts$contig)

# Join and split ----

divref_merged_with_gnomad <- divref %>%
  left_join(gnomad, by = join_by(variants == variant))

divref_in_gnomad <- divref_merged_with_gnomad %>%
  filter(!if_all(all_of(populations), is.na))

divref_not_in_gnomad <- divref_merged_with_gnomad %>%
  filter(if_all(all_of(populations), is.na)) %>%
  select(-all_of(populations))

divref_in_gnomad %>% write_tsv(paste0(opts$output_base, ".divref_in_gnomad.tsv"))

# Plot: Venn diagram of DivRef vs gnomAD variant overlap ----

# Count by distinct variant key, not by joined-row count: a duplicate `variant` in the gnomAD
# table would otherwise inflate the overlap and could push n_gnomad_only negative (which euler()
# rejects). On clean inputs these equal the row counts. The joined frames above are still used for
# the per-variant AF-difference plots.
divref_keys <- unique(divref$variants)
gnomad_keys <- unique(gnomad$variant)
n_both <- length(intersect(divref_keys, gnomad_keys))
n_divref_only <- length(setdiff(divref_keys, gnomad_keys))
n_gnomad_only <- length(setdiff(gnomad_keys, divref_keys))

log_info(n_both, " DivRef variants found in gnomAD")
log_info(n_divref_only, " DivRef variants not found in gnomAD")
log_info(n_gnomad_only, " gnomAD variants not found in DivRef")

venn_counts <- c(n_divref_only, n_gnomad_only, n_both)
names(venn_counts) <- c("DivRef 1.1", opts$gnomad_label, paste0("DivRef 1.1&", opts$gnomad_label))
fit <- euler(venn_counts)

png(paste0(opts$output_base, ".venn.png"), height = 800, width = 800)
plot(fit, quantities = TRUE)
invisible(dev.off())

# Plot: AF differences for variants found in gnomAD ----

divref_in_gnomad_with_af_diffs <- divref_in_gnomad %>%
  mutate(
    diff_afr = afr - gnomAD_AF_afr,
    diff_amr = amr - gnomAD_AF_amr,
    diff_eas = eas - gnomAD_AF_eas,
    diff_sas = sas - gnomAD_AF_sas,
    diff_nfe = nfe - gnomAD_AF_nfe,
  )

p <- divref_in_gnomad_with_af_diffs %>%
  select(variants, diff_afr, diff_amr, diff_eas, diff_sas, diff_nfe) %>%
  pivot_longer(
    cols = c(diff_afr, diff_amr, diff_eas, diff_sas, diff_nfe),
    names_to = "population", values_to = "diff_freq", names_prefix = "diff_"
  ) %>%
  ggplot(aes(x = diff_freq)) +
  geom_histogram() +
  facet_wrap(~population, nrow = 5) +
  scale_y_log10() +
  theme_bw() +
  xlab(paste0(opts$gnomad_label, " AF - DivRef 1.1 AF")) +
  ylab("Variants")

ggsave(paste0(opts$output_base, ".af_diffs.png"), p, height = 10, width = 6)

p <- divref_in_gnomad_with_af_diffs %>%
  select(variants, diff_afr, diff_amr, diff_eas, diff_sas, diff_nfe) %>%
  pivot_longer(
    cols = c(diff_afr, diff_amr, diff_eas, diff_sas, diff_nfe),
    names_to = "population", values_to = "diff_freq", names_prefix = "diff_"
  ) %>%
  ggplot(aes(x = diff_freq)) +
  geom_histogram() +
  scale_y_log10() +
  theme_bw() +
  xlab(paste0(opts$gnomad_label, " AF - DivRef 1.1 AF")) +
  ylab("Variants")

ggsave(paste0(opts$output_base, ".af_diffs_all.png"), p, height = 6, width = 6)

# Count: differences in maxpop

divref_in_gnomad %>%
  filter(!(gnomad_maxpop %in% populations)) %>%
  write_tsv(paste0(opts$output_base, ".max_pop_not_in_populations.tsv"))

divref_in_gnomad %>%
  filter(gnomad_maxpop %in% populations) %>%
  filter(divref_maxpop != gnomad_maxpop) %>%
  write_tsv(paste0(opts$output_base, ".max_pop_diffs.tsv"))

n_gnomad_maxpop_not_in_pops <- divref_in_gnomad %>%
  filter(!(gnomad_maxpop %in% populations)) %>%
  nrow

n_divref_in_gnomad_diff_maxpop <- divref_in_gnomad %>%
  filter(gnomad_maxpop %in% populations) %>%
  filter(divref_maxpop != gnomad_maxpop) %>%
  nrow

log_info(n_gnomad_maxpop_not_in_pops, " gnomAD variants in DivRef with no gnomAD maxpop in populations")
log_info(n_divref_in_gnomad_diff_maxpop, " DivRef variants with different maxpop as gnomAD")

# Plot: AN differences for variants found in gnomAD ----

divref_in_gnomad_with_an_diffs <- divref_in_gnomad %>%
  filter(gnomad_maxpop %in% populations) %>%
  filter(divref_maxpop == gnomad_maxpop) %>%
  mutate(diff_popmax_AN = popmax_AN - popmax_empirical_AN)

p <- divref_in_gnomad_with_an_diffs %>%
  ggplot(aes(x = diff_popmax_AN)) +
  geom_histogram() +
  scale_y_log10() +
  theme_bw() +
  xlab(paste0(opts$gnomad_label, " popmax AN - DivRef 1.1 popmax AN")) +
  ylab("Variants")

ggsave(paste0(opts$output_base, ".popmax_an_diffs.png"), p, height = 6, width = 6)

# Count: variants with large AF differences ----

large_af_diff <- divref_in_gnomad_with_af_diffs %>%
  select(variants, diff_afr, diff_amr, diff_eas, diff_sas, diff_nfe) %>%
  pivot_longer(
    cols = c(diff_afr, diff_amr, diff_eas, diff_sas, diff_nfe),
    names_to = "population", values_to = "diff_freq", names_prefix = "diff_"
  ) %>%
  mutate(abs_diff_freq = abs(diff_freq)) %>%
  filter(abs_diff_freq >= 0.001)

large_af_diff %>% write_tsv(paste0(opts$output_base, ".large_pop_af_diff.tsv"))

n_large_af_diff <- large_af_diff %>% distinct(variants) %>% nrow

log_info(n_large_af_diff, " DivRef variants found in gnomAD with |AF diff| >= 0.001 in any population")

if (nrow(divref_not_in_gnomad) == 0) {
  log_info("No DivRef variants missing from gnomAD, exiting")
  quit(save = "no")
}

# Plot: DivRef AF distribution for variants not found in gnomAD ----

p <- divref_not_in_gnomad %>%
  select(variants, gnomAD_AF_afr, gnomAD_AF_amr, gnomAD_AF_eas, gnomAD_AF_sas, gnomAD_AF_nfe) %>%
  pivot_longer(
    cols = c(gnomAD_AF_afr, gnomAD_AF_amr, gnomAD_AF_eas, gnomAD_AF_sas, gnomAD_AF_nfe),
    names_to = "population", values_to = "freq", names_prefix = "gnomAD_AF_"
  ) %>%
  ggplot(aes(x = freq)) +
  geom_histogram() +
  facet_wrap(~population, nrow = 5) +
  scale_y_log10() +
  theme_bw() +
  xlab("DivRef 1.1 AF") +
  ylab("Variants") +
  labs(title = paste0(
    "DivRef 1.1 'gnomAD_variant' variants not found in ", opts$gnomad_label
  ))

ggsave(paste0(opts$output_base, ".not_in_gnomad_afs.png"), p, height = 10, width = 6)

# Write not-in-gnomAD variant list ----

divref_not_in_gnomad %>%
  select(variants) %>%
  write_tsv(paste0(opts$output_base, ".divref_not_in_gnomad.tsv"))
