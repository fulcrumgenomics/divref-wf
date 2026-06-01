# Reproducing `docs/blog_draft.md`

This document lists every numeric claim, table cell, and figure in [`blog_draft.md`](blog_draft.md) and points to the command, output file, and log line that produces it.
The aim is that someone running the workflows on their own machine can verify each value against current data without re-reading the analysis scripts.

## Prerequisite: standard chr22 workflow run

All paths in this document are relative to the repository root.
Every section below depends on the standard `generate_divref` workflow run for chr22:

```bash
pixi run snakemake \
    -j 1 \
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

## Section: "DivRef 1.1 doesn't contain gnomAD 4.1 variants"

The 3-by-4 comparison table.

**Command**:

```bash
pixi run snakemake -j1 -s workflows/compare_divref_gnomad.smk
```

For each `gnomad_version` in `{hgdp_1kg_312, genomes_312, joint_41}` the workflow produces a TSV of all gnomAD variants on chr22 above the 0.5% threshold, runs the R comparator against the DivRef 1.1 single-variant track, and writes a log file.
The workflow also downloads the published DivRef 1.1 DuckDB index to `data/analysis/input/DivRef-v1.1.haplotypes_gnomad_merge.index.duckdb`; the haplotype-comparison section below reads it as `--old-duckdb`.

The `joint_41` extract applies an AC0-tolerant intersection filter: genome filter must be PASS, and exome filter must be PASS or contain only `AC0`.
The filter logic lives in `_apply_filters` in `divref/divref/tools/extract_gnomad_single_afs.py`; the test for it is `test_apply_filters_joint41_keeps_when_genome_pass_and_exome_pass_or_only_AC0` in `divref/tests/tools/test_extract_gnomad_single_afs.py`.

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

### Why DivRef-only variants are missing from gnomAD v4.1 joint

`scripts/explain_divref_only_v41.py` looks each DivRef-only variant up in the v4.1 joint sites Hail table and bins by reason.

**Command**:

```bash
pixi run python scripts/explain_divref_only_v41.py \
    --variants-tsv data/analysis/compare_divref_gnomad/chr22.joint_41.divref_not_in_gnomad.tsv \
    &> logs/explain_divref_only_v41.chr22.log
```

The script prints per-bucket counts (`absent_from_v41`, `below_af_threshold`, `exome_filter_only`, `genome_filter_only`, `both_filters_nonempty`, and an `unexpected_pass` sanity-check bucket) and writes a per-variant TSV alongside the input at `data/analysis/compare_divref_gnomad/chr22.joint_41.divref_not_in_gnomad.explained.tsv`.

### Whole-genome extrapolation

The blog's whole-genome extrapolations from chr22 use an empirical multiplier derived from DivRef 1.1's own DuckDB index:

```bash
pixi run duckdb -readonly data/analysis/input/DivRef-v1.1.haplotypes_gnomad_merge.index.duckdb -c "
SELECT
  COUNT(*) AS total_gnomad_variant_rows,
  SUM(CASE WHEN contig = 'chr22' THEN 1 ELSE 0 END) AS chr22_rows,
  CAST(COUNT(*) AS DOUBLE) / SUM(CASE WHEN contig = 'chr22' THEN 1 ELSE 0 END) AS genome_to_chr22_ratio
FROM sequences
WHERE source = 'gnomAD_variant';
"
```

## Section: "Comparing against the original algorithm"

**Prerequisite**: a whole-genome `generate_divref` run (autosomes chr1-22) that produced the new DuckDB index `data/work/output/hgdp_1kg.haplotypes_gnomad_merge.index.duckdb`, plus the original DivRef 1.1 DuckDB at `data/analysis/input/DivRef-v1.1.haplotypes_gnomad_merge.index.duckdb`.
This section compares the two DuckDBs directly, so the per-chromosome `.ht` files are not needed here (only the chr22 deep-dive below still uses the chr22 workflow outputs).

### Venn figure

`data/analysis/compute_haplotypes/algo_comparison.venn.png`

**Commands**:

```bash
mkdir -p logs
pixi run python scripts/compare_haplotypes.py &> logs/compare_haplotypes.log
pixi run Rscript scripts/compare_haplotypes_venn.R &> logs/compare_haplotypes_venn.log
```

The Python script defaults to comparing the original DivRef 1.1 DuckDB (`--old-duckdb`) against the whole-genome DuckDB (`--new-duckdb`) over the autosomes (`--contigs`), and writes `data/analysis/compute_haplotypes/algo_comparison.summary.tsv`; the R script reads the `shared` / `old_only` / `new_only` rows from there and renders the figure.
No counts are hardcoded in the R script.

### Sub-fragment counts and new-only decomposition

**Command**: `pixi run python scripts/compare_haplotypes.py &> logs/compare_haplotypes.log` (same invocation as the Venn section above).

| Blog claim | Summary TSV row |
|---|---|
| original-only count | `old_only` |
| new-only count | `new_only` |
| shared count | `shared` |
| original-only proper sub-fragments of some new | `old_only_subfragment_of_new` |
| of those: both algorithms find the longer fragment | `old_only_subfragment_both_have_longer` |
| of those: only new finds the longer fragment | `old_only_subfragment_only_new_has_longer` |
| new-only with all variants present in old | `new_only_all_variants_in_old` |
| new-only with a mix of shared and novel variants | `new_only_mixed_variants` |
| new-only with all-novel variants | `new_only_all_novel` |
| all-novel new-only of length 2 | `new_only_all_novel_length_2` |
| max length among all-novel new-only | `new_only_all_novel_max_length` |
| original-only that are NOT a sub-fragment of any new (residual) | `old_only_not_subfragment` |
| chrX haplotypes the new algorithm adds (coverage note, not part of the autosome comparison) | `new_chrX_haplotypes` |

The `old_only_not_subfragment` residual haplotypes are the deep-dive cases analysed by `inspect_threshold_edge_haplotypes.py` (restricted to chr22 below).

## Section: "A closer look at the six original-only haplotypes"

The 6-row case table and the surrounding Mechanism prose.
This is a chr22 zoom-in: `pixi run python scripts/compare_haplotypes.py --contigs chr22` restricts the comparison to chr22 and reproduces its residual count of six (the genome-wide `old_only_not_subfragment` is larger).

The standard workflow filters out the six original-only haplotypes the blog inspects (their new-algorithm `estimated_gnomad_AF` falls below the 0.005 threshold) and deletes the per-chromosome intermediate Hail tables that the inspector script needs.
Two steps reconstruct just what's needed.

### Step 1: rebuild per-case intermediates with no AF cutoff

Subset the existing chr22 phased VCF to ±500 bp around each case window and run `compute-haplotypes` on that subset with threshold 0.
The case regions sum to ~6 kb (0.01% of chr22), so the run takes a couple of minutes instead of ~20.

```bash
mkdir -p data/analysis/compute_haplotypes/test_data_new_low_af logs

# Subset the chr22 phased VCF to ±500 bp around the six case windows.
# Regions cover the case variants in scripts/inspect_threshold_edge_haplotypes.py CASES.
pixi run bcftools view \
    -r chr22:19950635-19951656,chr22:40456973-40457977,chr22:32402062-32403068,chr22:22626850-22627866,chr22:24600487-24601522,chr22:24626368-24627406 \
    -O z -o data/analysis/compute_haplotypes/test_data_new_low_af/hgdp_1kg.phased_genotypes.chr22.cases.vcf.gz \
    data/work/inputs/hgdp_1kg.phased_genotypes.chr22.vcf.gz
pixi run bcftools index --tbi \
    data/analysis/compute_haplotypes/test_data_new_low_af/hgdp_1kg.phased_genotypes.chr22.cases.vcf.gz

# Run compute_haplotypes on the subset with no haplotype-AF filter so every case is emitted.
pixi run divref compute-haplotypes \
    --vcfs-path data/analysis/compute_haplotypes/test_data_new_low_af/hgdp_1kg.phased_genotypes.chr22.cases.vcf.gz \
    --gnomad-va-file data/work/inputs/hgdp_1kg.sites.chr22.ht \
    --gnomad-sa-file data/work/inputs/hgdp_1kg.sample_metadata.ht \
    --window-size 25 \
    --variant-freq-threshold 0.005 \
    --haplotype-freq-threshold 0 \
    --output-base data/analysis/compute_haplotypes/test_data_new_low_af/hgdp_1kg.haplotypes.chr22 \
    --spark-driver-memory-gb 16 --spark-executor-memory-gb 16 \
    &> logs/compute_haplotypes.chr22.cases.log
```

Produces:

- `data/analysis/compute_haplotypes/test_data_new_low_af/hgdp_1kg.haplotypes.chr22.variants.ht`
- `data/analysis/compute_haplotypes/test_data_new_low_af/hgdp_1kg.haplotypes.chr22.hap_ac.ht`
- `data/analysis/compute_haplotypes/test_data_new_low_af/hgdp_1kg.haplotypes.chr22.parents.ht`
- `data/analysis/compute_haplotypes/test_data_new_low_af/hgdp_1kg.haplotypes.chr22.ht`

The output Hail tables are restricted to variants in the subset regions, but the per-variant call_stats and per-sample parent-block formation use the same logic as the full workflow, so the case-specific `popmax_AC`, `fp`, and `est_af` values match what the full-chr22 run with the same threshold would produce.

### Step 2: run the inspector against the case intermediates

```bash
mkdir -p logs
pixi run python scripts/inspect_threshold_edge_haplotypes.py \
    --base data/analysis/compute_haplotypes/test_data_new_low_af/hgdp_1kg.haplotypes.chr22 \
    &> logs/inspect_threshold_edge_haplotypes.chr22.log
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

The "Mechanism 1 / Mechanism 2" classification in the prose is derived by inspection from the per-case OLD/NEW print lines:

- Mechanism 1 ("`fp = 1.0000` exactly in the original") applies when the OLD DUCKDB line reports `fp=1.0000` and the NEW HT line reports a smaller `fp` at the same `max_pop` (cases 1, 4, 5).
- Mechanism 2 ("`max_pop` shifts under containment counting") applies when the two `max_pop=` values differ (cases 2, 3, 6).

### AC counts table

**Command**: `pixi run python scripts/compare_haplotypes.py &> logs/compare_haplotypes.chr22.log` (same invocation as the Venn / decomposition counts above).

| Blog row | Summary TSV row |
|---|---|
| Same popmax AC | `shared_same_ac` |
| New > old | `shared_new_higher_ac` |
| Old > new | `shared_old_higher_ac` |
