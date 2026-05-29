# Reproducing `docs/blog_draft.md`

This document lists every numeric claim, table cell, and figure in [`blog_draft.md`](blog_draft.md) and points to the command, output file, and log line that produces it.
The aim is that someone running the workflows on their own machine can verify each value against current data without re-reading the analysis scripts.

## Prerequisite data fixtures

The blog draws on two independent sets of pre-computed Hail / DuckDB artefacts.
All paths in this document are relative to the repository root.

### From the main `generate_divref` workflow (chr22 standard run)

Produced by:

```bash
pixi run snakemake \
    -j 1
    -s workflows/generate_divref.smk \
    --configfile workflows/config/config_gcs.yml \
    --config 'chromosomes=["chr22"]'
```

Outputs:

- `data/work/inputs/hgdp_1kg.phased_genotypes.chr22.vcf.gz` (+ `.tbi` index)
- `data/work/inputs/hgdp_1kg.sites.chr22.ht` (HGDP+1KG per-pop AFs)
- `data/work/inputs/hgdp_1kg.sample_metadata.ht`
- `data/work/inputs/gnomad.sites.chr22.ht` (gnomAD single-variant track for the DuckDB build)
- `data/work/inputs/Homo_sapiens_assembly38.fasta` (+ `.fai` index)
- `data/work/inputs/table_pairs.tsv`
- `data/work/haplotypes/hgdp_1kg.haplotypes.chr22.ht` (new algorithm haplotype output)
- `data/work/output/hgdp_1kg.haplotypes_gnomad_merge.index.duckdb` (final DivRef-style DuckDB)
- `data/work/output/hgdp_1kg.haplotypes_gnomad_merge.chr22.fasta` (per-chromosome FASTA)

### From the low-haplotype-AF rebuild (for the six-case deep dive only)

The standard workflow filters out the six original-only haplotypes the blog inspects (their new-algorithm `estimated_gnomad_AF` falls below the 0.005 threshold) and deletes the per-chromosome intermediate Hail tables that the inspector script needs.
Reproduce the intermediates with a direct `divref compute-haplotypes` invocation that uses a lower threshold and writes outside `data/work/`:

```bash
mkdir -p data/analysis/compute_haplotypes/test_data_new_low_af
pixi run divref compute-haplotypes \
    --vcfs-path data/work/inputs/hgdp_1kg.phased_genotypes.chr22.vcf.gz \
    --gnomad-va-file data/work/inputs/hgdp_1kg.sites.chr22.ht \
    --gnomad-sa-file data/work/inputs/hgdp_1kg.sample_metadata.ht \
    --window-size 25 \
    --variant-freq-threshold 0.005 \
    --haplotype-freq-threshold 0.002 \
    --output-base data/analysis/compute_haplotypes/test_data_new_low_af/hgdp_1kg.haplotypes.chr22 \
    --spark-driver-memory-gb 16 --spark-executor-memory-gb 16
```

Produces:

- `data/analysis/compute_haplotypes/test_data_new_low_af/hgdp_1kg.haplotypes.chr22.variants.ht`
- `data/analysis/compute_haplotypes/test_data_new_low_af/hgdp_1kg.haplotypes.chr22.hap_ac.ht`
- `data/analysis/compute_haplotypes/test_data_new_low_af/hgdp_1kg.haplotypes.chr22.parents.ht`
- `data/analysis/compute_haplotypes/test_data_new_low_af/hgdp_1kg.haplotypes.chr22.ht`

## Section: "DivRef 1.1 doesn't contain gnomAD 4.1 variants"

The 3 by 4 table at blog lines 35-39.

**Command**:

```bash
pixi run snakemake -j1 -s workflows/compare_divref_gnomad.smk
```

For each `gnomad_version` in `{hgdp_1kg_312, genomes_312, joint_41}` the workflow produces a TSV of all gnomAD variants on chr22 above the 0.5% threshold, runs the R comparator against the DivRef 1.1 single-variant track, and writes a log file.
The workflow also downloads the published DivRef 1.1 DuckDB index to `data/analysis/input/DivRef-v1.1.haplotypes_gnomad_merge.index.duckdb`; the haplotype-comparison section below reads it as `--old-duckdb`.
The blog's columns map to log/TSV values as follows.

| Blog column | Source |
|---|---|
| gnomAD variants | `wc -l < data/analysis/compare_divref_gnomad/chr22.<v>.tsv` minus the header row (equivalently the `Loaded N gnomAD variants for chr22` line in `chr22.<v>.log`) |
| In DivRef 1.1 | `N DivRef variants found in gnomAD` in `chr22.<v>.log` |
| DivRef 1.1-only | `N DivRef variants not found in gnomAD` in `chr22.<v>.log` |
| gnomAD-only | "gnomAD variants" column minus "In DivRef 1.1"; also written as `n_gnomad_only` on the `N gnomAD variants not found in DivRef` log line |

Row mapping:

| Blog row | Log file |
|---|---|
| 3.1.2 HGDP+1KG subset (4K) | `chr22.hgdp_1kg_312.log` |
| 3.1.2 genomes (76K) | `chr22.genomes_312.log` |
| 4.1 joint (76K genomes, 730K exomes) | `chr22.joint_41.log` |

## Section: "Comparing against the original algorithm"

### Venn figure (blog line 119)

`data/analysis/compute_haplotypes/algo_comparison.venn.png`

**Commands**:

```bash
pixi run python scripts/compare_haplotypes.py
pixi run Rscript scripts/compare_haplotypes_venn.R
```

The Python script writes `data/analysis/compute_haplotypes/algo_comparison.summary.tsv`; the R script reads the `shared` / `old_only` / `new_only` rows from there and renders the figure.
No counts are hardcoded in the R script.

### 1,327 of 1,333 + the 1,197 / 130 split + the new-only 1,590 / 54 / 36 / 21 / 8 numbers (blog lines 121-126)

**Command**: `pixi run python scripts/compare_haplotypes.py`

| Blog claim | Summary TSV row |
|---|---|
| 1,333 original-only | `old_only` |
| 1,680 new-only | `new_only` |
| 29,548 shared | `shared` |
| 1,327 original-only are proper sub-fragments of some new | `old_only_subfragment_of_new` |
| 1,197 of those: both algorithms find the longer fragment | `old_only_subfragment_both_have_longer` |
| 130 of those: only new finds the longer fragment | `old_only_subfragment_only_new_has_longer` |
| 1,590 new-only with all variants present in old | `new_only_all_variants_in_old` |
| 54 new-only with a mix of shared and novel variants | `new_only_mixed_variants` |
| 36 new-only with all-novel variants | `new_only_all_novel` |
| 21 of those 36 are 2-variant | `new_only_all_novel_length_2` |
| max length 8 | `new_only_all_novel_max_length` |

The remaining 6 (`old_only` minus `old_only_subfragment_of_new`) are the deep-dive cases analysed by `inspect_threshold_edge_haplotypes.py`.

## Section: "A closer look at the six original-only haplotypes"

The 6-row case table (blog lines 134-150) and the 4-row Case 6 chromosome breakdown (blog lines 173-178).

**Command**:

```bash
pixi run python scripts/inspect_threshold_edge_haplotypes.py \
    --base data/analysis/compute_haplotypes/test_data_new_low_af/hgdp_1kg.haplotypes.chr22
```

For each case (`Adjacent #1`, `VNTR`, `Short-gap`, `Intermediate`, `Triple+skip`, `Non-contig`) the script prints two lines that map directly to the blog columns:

```
  OLD DUCKDB: max_pop=... popmax_AC=... popmax_AF=... fp=... est_af=...
  NEW HT:     max_pop=... popmax_AC=... popmax_AF=... min_var_freq=... fp=... est_af=...
```

| Blog case # | Script case label |
|---|---|
| 1 | Adjacent #1 |
| 2 | VNTR |
| 3 | Short-gap |
| 4 | Intermediate |
| 5 | Triple+skip |
| 6 | Non-contig |

The "Mechanism 1 / Mechanism 2" classification in the prose (blog lines 157-168) is derived by inspection from the per-case OLD/NEW print lines:

- Mechanism 1 ("`fp = 1.0000` exactly in the original") applies when the OLD DUCKDB line reports `fp=1.0000` and the NEW HT line reports a smaller `fp` at the same `max_pop` (cases 1, 4, 5).
- Mechanism 2 ("`max_pop` shifts under containment counting") applies when the two `max_pop=` values differ (cases 2, 3, 6).

### Case 6 chromosome breakdown (blog lines 173-178)

In the same script output, the `Non-contig` case prints a `PARENT BLOCK PATTERNS touching at least one case variant` section.
The four parent-block patterns there map to the blog rows:

| Blog row | Parent-block pattern (positions) | Count column |
|---|---|---|
| 84 — `(v1, v2, v3)`, contiguous | `[24626868:C>T, 24626884:T>C, 24626892:A>G]` `[PURE]` | sum of per-pop counts |
| 45 — `(v1, v3)`, pure (no alt between) | `[24626868:C>T, 24626892:A>G]` `[PURE]` | sum of per-pop counts |
| 2 — `(v1, i1, v3)`, also alt at intermediate `i1=24626883` | `[24626868:C>T, 24626883:..., 24626892:A>G]` `[non-contig]` | sum of per-pop counts |
| 37 — singleton `v1`, dropped (length < 2) | `[24626868:C>T]` | sum of per-pop counts |

The "i1=24626883" position falls inside the ±100 bp scan window of the script's `variants.ht` filter; the exact ref/alt for `i1` appears in the script's per-case `variants` print line.

### AC counts table (blog lines 188-191)

**Command**: `pixi run python scripts/compare_haplotypes.py` (same run as the Venn / decomposition counts above).

| Blog row | Summary TSV row |
|---|---|
| Same popmax AC (25,718) | `shared_same_ac` |
| New > old (3,785) | `shared_new_higher_ac` |
| Old > new (45) | `shared_old_higher_ac` |
