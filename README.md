[![CI](https://github.com/fg-labs/divref-wf/actions/workflows/python_package.yml/badge.svg?branch=main)](https://github.com/fg-labs/divref-wf/actions/workflows/python_package.yml?query=branch%3Amain)
[![Python Versions](https://img.shields.io/badge/python-3.12_|_3.13-blue)](https://github.com/fg-labs/divref-wf)
[![MyPy Checked](http://www.mypy-lang.org/static/mypy_badge.svg)](http://mypy-lang.org/)
[![uv](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/uv/main/assets/badge/v0.json)](https://github.com/astral-sh/uv)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://docs.astral.sh/ruff/)
[![Pixi](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/prefix-dev/pixi/main/assets/badge/v0.json)](https://pixi.sh)

# Snakemake workflow implementation to create DivRef-style resource

This workflow is inspired by the [DivRef](https://github.com/e9genomics/human-diversity-reference) repository which is used to generate a bundle of FASTA sequences and a corresponding DuckDB index of common human variation.

The original implementation is via a set of standalone Python scripts and a Makefile.

This implementation:

1. Wraps the Python scripts in a toolkit with added typing, improved parameterization, and added unit testing.
2. Adds a Snakemake workflow and associated configuration to drive the resource generation process.

## Set up Environment

The environment for this analysis is managed using `pixi`.
Follow the developer [instructions](https://pixi.sh/latest/installation/) to install `pixi`.

The environment and dependencies are automatically created and installed by calling `pixi install` or when calling `pixi run` for the first time.

By default the workflow reads gnomAD inputs from GCS (`gs://gcp-public-data--gnomad/`).
To install the GCS connector for Hail/Spark, run

```bash
pixi run setup-gcs
```

Log in before running any Hail-dependent tools:

```bash
gcloud auth application-default login
```

To use the AWS Open Data S3 mirror (`s3a://gnomad-public-us-east-1/`) instead — relevant
when running compute in `us-east-1` — run

```bash
pixi run setup-s3
```

…and configure AWS credentials via the standard AWS chain (`aws configure`, environment
variables, or an IAM role).

To install both connectors at once, run `pixi run setup-cloud`.

### Running the workflow

The workflow does not bundle a default `configfile:` — pass one explicitly with
`--configfile`. Two ready-made configs are provided under `workflows/config/`:

- `config_gcs.yml` — reads all cloud inputs from GCS (`gs://gcp-public-data--gnomad/`,
  `gs://hail-common/`).
- `config_aws.yml` — reads all cloud inputs from the AWS Open Data S3 mirror
  (`s3://gnomad-public-us-east-1/`, `s3://broad-references/`).

Run on GCS (after `pixi run setup-gcs` and `gcloud auth application-default login`):

```bash
pixi run snakemake -j1 -s workflows/generate_divref.smk \
    --configfile workflows/config/config_gcs.yml
```

Run on AWS (after `pixi run setup-s3` and configuring AWS credentials):

```bash
pixi run snakemake -j1 -s workflows/generate_divref.smk \
    --configfile workflows/config/config_aws.yml
```

To run multi-threaded, set `-j` to be greater than `1`.

To override individual settings (e.g. `chromosomes`, `version`, output paths) without
editing the shipped configs, append `--config key=value …` after the `--configfile`
argument.

## Resource Description

The below statements are from using the [default parameters](workflows/config/config_schema.yml).

### Populations

AFR, AMR, EAS, SAS, NFE

### HGDP_haplotype

Haplotypes are derived from the [gnomAD 3.1.2 HGDP+1KG individual-level phased genotypes](https://gnomad.broadinstitute.org/news/2021-10-gnomad-v3-1-2-minor-release/).

- Individuals are annotated with continental ancestry using [gnomAD labels](https://gnomad.broadinstitute.org/data).
- Only variants in the HGDP+1KG subset of gnomAD 3.1.2 are considered for inclusion, as variants only present in the full genomes dataset do not have associated phased genotypes.
- Variants with less than 0.5%AF in the full gnomAD 3.1.2 genomes (n=76,156) dataset in all of the populations are removed.
- Sets of phased alleles are computed over 100-base-pair windows of the genome, with two passes offset by 50bp so that haplotypes are not split at fixed window boundaries.
- Duplicate haplotypes produced by the overlapping passes are removed when the DuckDB index is built (after sub-haplotype splitting, see below).
- For each haplotype, the population with the highest empirical AF is recorded as `max_pop`. The phase ratio is the empirical AF of the haplotype in `max_pop` divided by the empirical AF of the rarest component variant in `max_pop`.
- The estimated gnomAD haplotype frequency is the phase ratio multiplied by the smallest gnomAD AF among the component variants in `max_pop`.
- Haplotypes with estimated gnomAD frequency < 0.5% are removed.
- Haplotypes spanning variants further than or equal to 25bp apart are broken into sub-haplotypes at those gaps. Sub-haplotypes with fewer than two variants are discarded.

**Note** Given the sample size of the HGDP+1KG phased genotypes dataset, there is limited power to detect haplotypes under 1%.
We have included haplotypes discovered between 0.5% and 1%, but would expect to find more haplotypes in that frequency range with a larger dataset.

### gnomAD_variant

gnomAD variants are derived from the [gnomAD 4.1 joint exomes+genomes sites](https://gnomad.broadinstitute.org/news/2024-04-gnomad-v4-1/).

- Genotypes at variants with less than 0.5% AF in all of the populations are removed.

### Index database and FASTA files

The `HGDP_haplotype` and `gnomAD_variant` data is merged and position-sorted.
The 25bp sequence context around each haplotype/variant is obtained from the GRCh38 reference genome and exported to FASTA.
An accompanying DuckDB index database contains all of the input haplotypes/variants, their population AFs, and for haplotypes, the empirical AF and AC.

The FASTA files are intended to be used with tools for guide design and off-target nomination that already accept reference sequences in FASTA format.

#### `sequences` table columns

The primary table in the DuckDB index is the `sequences` table, which has one row per haplotype or single variant.

| Column | Description |
|---|---|
| `sequence_id` | Unique identifier of the form `DR-{version}-{index}`. `index` is the global row number assigned in genomic-position order across all contigs. |
| `sequence` | Sequence string for the haplotype/variant, with `window_size` bp of flanking reference context on each side. |
| `sequence_length` | Length of `sequence` in bases. |
| `n_variants` | Number of component variants in this row (`1` for `gnomAD_variant`, `>= 2` for `HGDP_haplotype`). |
| `contig` | Reference contig (e.g. `chr1`). |
| `start` | 0-based inclusive start of `sequence` on `contig`: `(first_variant_position − 1) − window_size`. |
| `end` | 0-based exclusive end of `sequence` on `contig`: `(last_variant_position − 1 + length(ref_allele)) + window_size`. |
| `source` | `HGDP_haplotype` or `gnomAD_variant`. |
| `max_pop` | Population code with the highest empirical AF for this haplotype/variant. |
| `popmax_empirical_AF` | For `HGDP_haplotype`: empirical AF of the haplotype in `max_pop` from observed phased genotypes. For `gnomAD_variant`: the gnomAD AF in `max_pop`. |
| `popmax_empirical_AC` | Allele count corresponding to `popmax_empirical_AF`. |
| `popmax_estimated_gnomad_AF` | For `HGDP_haplotype`: `popmax_fraction_phased × min(gnomAD_AF[component, max_pop])` over the component variants — i.e. the phase ratio applied to the rarest component's gnomAD AF in `max_pop`. For `gnomAD_variant`: the gnomAD AF in `max_pop`. Equivalent to `estimated_gnomad_AF_{max_pop}`. |
| `popmax_fraction_phased` | For `HGDP_haplotype`: phase ratio, `popmax_empirical_AF / min(empirical_AF[component, max_pop])` — empirical haplotype AF over empirical AF of the rarest component variant in the same population. `1.0` for `gnomAD_variant`. Equivalent to `fraction_phased_{max_pop}`. |
| `variants` | Comma-separated list of component variants in `chr:pos:ref:alt` form, in the order they appear in the haplotype. |
| `gnomAD_AF_{POP}` | One column per configured population (e.g. `gnomAD_AF_AFR`). Comma-separated per-component gnomAD AFs in that population, in the same order as `variants`, formatted to 5 decimal places. |
| `empirical_AC_{POP}` | One column per pop in `joint_pops_legend`. For `HGDP_haplotype`: count of chromosomes carrying the haplotype in that pop (haplotype-level). For `gnomAD_variant`: the gnomAD AC for that variant in that pop (variant-level). Missing if the pop has no source data for this row. |
| `empirical_AF_{POP}` | For `HGDP_haplotype`: empirical AF of the haplotype in that pop, derived from `empirical_AC_{POP}` and the min AN over component variants. For `gnomAD_variant`: the gnomAD AF for that variant in that pop. Missing if AN is 0 or the pop has no source data. |
| `fraction_phased_{POP}` | Per-pop phase ratio: `empirical_AF_{POP} / min(local_call_stats_AF[component, POP])`. `1.0` for `gnomAD_variant` rows. Uses *that pop's own* denominators (not `max_pop`'s). |
| `estimated_gnomad_AF_{POP}` | Per-pop projection: `fraction_phased_{POP} × min(gnomAD_AF[component, POP])`. Equals `empirical_AF_{POP}` for `gnomAD_variant` rows. |

The DuckDB file also contains three single-row metadata tables: `window_size` (the flanking context size used), `pops_legend` (JSON-encoded ordered population list), and `VERSION`.

## Analysis

### Additional Environment Requirements

To install R packages not available as conda-forge builds for all platforms (duckdb, duckplyr), run

```bash
pixi run -e analysis setup-r-packages
```

### Compare DivRef 1.1 against different gnomAD releases

[DivRef 1.1](https://zenodo.org/records/14802613) states that:

> DivRef is constructed by computing empirical phased haplotypes within 25 BPs over 0.5% allele frequency from the Human Genome Diversity Panel (HGDP) using the phased Hail dataset provided by the gnomAD team at the Broad Institute, merged with single variants over 0.5% AF from the gnomAD v4.1.0 summary release.

Some gnomAD v4.1.0 variants that we expected to see represented in DivRef 1.1 were not.
We checked all 'gnomAD_variant' variants on chr22 from DivRef 1.1 against:

- gnomAD 3.1.2 HGDP+1KG subset
- gnomAD 3.1.2 genomes (~76K genomes)
- gnomAD 4.1 joint (~730K exomes and ~76K genomes)

gnomAD 3.1.2 HGDP+1KG subset is the source used for the haplotypes present in DivRef 1.1.

```bash
pixi run -e analysis snakemake -j1 -s workflows/compare_divref_gnomad.smk
```

**DivRef 1.1 variants present in gnomAD datasets**

We found that all DivRef 1.1 'gnomAD_variant' variants were present in the gnomAD 3.1.2 HGDP+1KG subset and in the gnomAD 3.1.2 genomes set, while 16 were missing from the gnomAD 4.1 joint set.

We further compared the allele frequencies for the 5 populations recorded in the DivRef 1.1 DuckDB index against the frequencies for those populations in the two gnomAD sets, using the Hail tables as input.

- gnomAD 3.1.2 HGDP+1KG subset: within a rounding error of 5e-6, for all variants.
- gnomAD 3.1.2 genomes: 1,400 variants with an AF difference >= 0.001, all of which had lower AF in the genomes set compared to HGDP+1KG, indicating more stringent filtering in the genomes set
- gnomAD 4.1 joint: 60,424 variants with an AF difference >= 0.001, evenly distributed both lower and higher AF

There are 506,983 DivRef 1.1 'gnomAD_variant' variants. When we run the original DivRef script [extract_gnomad_afs.py](https://github.com/e9genomics/human-diversity-reference/blob/main/scripts/extract_gnomad_afs.py) (uses gnomAD 3.1.2 HGDP+1KG sites as input, hard-coded) followed by lines 156-177 of [create_fasta_and_index.py](https://github.com/e9genomics/human-diversity-reference/blob/main/scripts/create_fasta_and_index.py) (lines 156-177) using the parameters specified in the [Makefile](https://github.com/e9genomics/human-diversity-reference/blob/main/Makefile), we also get 506,983 variants. 

We concluded that the DivRef 1.1 documentation was incorrect, and that the actual source of the gnomAD variants in the dataset was the gnomAD 3.1.2 HGDP+1KG subset, the same as for the haplotypes.

### Updated haplotype computation algorithm

#### Old algorithm

**Stage 1: MT prep** (same as new — kept verbatim in the rewrite)
- Read VCF → MatrixTable, filter to phased biallelic variants, attach `row_idx` / `col_idx`.
- Filter entries to alt-carriers; annotate `frequencies_by_pop` from gnomAD.
- Checkpoint to `{output_base}.variants.ht` for downstream row-index → variant lookup.

**Stage 2: two fixed-bin window passes**

The genome is partitioned into non-overlapping bins of `window_size` bp. To avoid systematic edge artefacts (a true haplotype split across a bin boundary would be missed), the tool runs **two passes offset by `window_size / 2`**:
- Pass 1: bins `[0, W)`, `[W, 2W)`, …
- Pass 2: bins `[W/2, 3W/2)`, `[3W/2, 5W/2)`, …

Each pass independently:
1. Assigns every alt-carrier entry to its bin.
2. Per (sample, strand, bin), `agg.collect`s the carried variants into a haplotype block.
3. Drops blocks of length < 2.
4. Writes intermediate `{output_base}.1.ht` / `{output_base}.2.ht`.

**Stage 3: cross-sample collapse** (`collapse_haplos_across_samples`)
- For each pass, group by the haplotype (row-index tuple) and aggregate `empirical_AC` per population: count of distinct (sample, strand) tuples that emitted exactly this block.
- This counts only *exact full-block matches across samples* — two samples whose blocks differed by even one variant produced two separate rows.

**Stage 4: union** (`get_haplotype_summary`)
- `union` the two passes' output tables.
- If the same haplotype appears in both passes (e.g., a block fully contained in one bin of pass 1 and also one bin of pass 2), it is *kept twice* — one row per pass — unless a later `distinct()` collapses them. This is one source of double-counting.

**Stage 5: per-fragment metrics**
- Same as new: derive `max_pop`, `max_empirical_AF`, `max_empirical_AC`, `min_variant_frequency`, `fraction_phased`, `estimated_gnomad_AF`. **Note:** in the old algorithm `min_variant_frequency` was derived from the gnomAD-sites AF rather than the local HGDP+1KG `call_stats.AF[1]`; with the gnomAD AF on both sides of the formula the gnomAD multiplier cancels and `estimated_gnomad_AF` collapses to `max_empirical_AF`. This was corrected as part of the rewrite — the column definitions in the table above describe the new (correct) semantic.

**Stage 6: filter**
- `min_variant_frequency > 0` and `estimated_gnomad_AF >= haplotype_freq_threshold`.

**Stage 7 (in `create_duckdb_index`, not in old `compute_haplotypes`): `split_haplotypes` + `distinct()`**
- Re-split each haplotype at gaps ≥ `window_size`, emitting all sub-blocks of length ≥ 2.
- `distinct()` on `haplotype` — drops duplicate row-index tuples *but discards their AC counts*, so AC is not summed across the duplicates this step creates. This is the AC-underestimation bug that motivated the rewrite.


#### New algorithm

Enumerate every adjacency-contiguous sub-fragment of length ≥ 2 across all observed per-sample parent blocks at a single window W, with each sub-fragment's empirical_AC per population equal to the count of parent blocks containing it as adjacency-contiguous (full or proper sub).

Drop any sub-fragment that is contained in a larger sub-fragment with the same per-population AC vector — the larger fragment carries strictly more information at the same observation count, so the smaller is redundant.

Net: the index contains every distinct full block (always survives, since no larger fragment contains it) plus every sub-fragment whose containment AC genuinely exceeds what its enclosing fragments provide.

Concretely, for parents [V1,V2,V3] (sample A) and [V0,V1,V2] (sample B), the algorithm emits [V1,V2,V3] AC=1, [V0,V1,V2] AC=1, and [V1,V2] AC=2. It does not emit [V2,V3] or [V0,V1] because each is contained in a larger fragment (its parent's full block) with the same AC=1.

#### Limitations of old algorithm

1. **Bin-edge blindness for both passes.** A variant pair with gap < W can still fall on opposite sides of *both* bin grids if the variants happen to straddle position `kW` and `(k+0.5)W` simultaneously. Per-sample adjacency in the new algorithm uses gap-based cuts, so any pair within W bp on the same strand of the same sample is always grouped.
2. **No cross-context AC summation.** Sample A emits `[V1,V2,V3]` and sample B emits `[V0,V1,V2]`. Both contain `[V1,V2]`, but neither emits `[V1,V2]` as its full block, so the old algorithm never produces `[V1,V2]` at all. The new algorithm enumerates all sub-fragments per parent and aggregates AC across containing parents, yielding `[V1,V2]` AC=2.
3. **Double-counting from the union of two passes**, partially offset by `distinct()` in `create_duckdb_index` *dropping* AC information — so the bias goes both directions and is hard to reason about.

#### Comparison of old algorithm to new algorithm on chr22

The full pipeline was run on chr22 with both algorithms.
The old algorithm relied on a two stage process in both `compute_haplotypes.py` and `create_duckdb_index.py`, so we use the DuckDB as the final set of haplotypes.
The new algorithm is a single stage only in `compute_haplotypes.py`, so we can use the Hail table of haplotypes directly.

Old:

```bash
git checkout bd81b8ee55de4f52ec6efcb32a5f4b92b647b17d
pixi run snakemake -j1 -s workflows/generate_divref.smk --config chromosomes=["chr22"]
mkdir -p data/analysis/compute_haplotypes/test_data_old
mv data/work/output/hgdp_1kg.haplotypes_gnomad_merge.index.duckdb data/analysis/compute_haplotypes/test_data_old/
```

New:

```bash
git checkout main
pixi run snakemake -j1 -s workflows/generate_divref.smk --config chromosomes=["chr22"]
mkdir -p data/analysis/compute_haplotypes/test_data_new
mv data/work/haplotypes/hgdp_1kg.haplotypes.chr22.ht data/analysis/compute_haplotypes/test_data_new/
```

Comparison

```bash
pixi run python scripts/compare_haplotypes.py > data/analysis/compute_haplotypes/chr22_old_duckdb_vs_new_compute_haplotypes.txt
```

**Counts**

- Old: 30,881
- New: 30,136
- Shared: 28,265 (≈92% of old, ≈94% of new)
- Old-only: 2,616
- New-only: 1,871

**Length distribution**

| Length | Old | New | Shared | Old-only | New-only |
|---|---|---|---|---|---|
| 2  | 26,140 | 25,306 | 24,103 | 2,037 | 1,203 |
| 3  | 3,515  | 3,508  | 3,137  | 378   | 371   |
| 4  | 781    | 826    | 684    | 97    | 142   |
| 5  | 253    | 268    | 208    | 45    | 60    |
| 6  | 96     | 101    | 68     | 28    | 33    |
| 7  | 48     | 57     | 34     | 14    | 23    |
| 8  | 23     | 25     | 14     | 9     | 11    |
| 9  | 12     | 16     | 7      | 5     | 9     |
| 10 | 7      | 13     | 5      | 2     | 8     |
| 11 | 5      | 5      | 4      | 1     | 1     |
| 12 | 0      | 7      | 0      | 0     | 7     |
| 13 | 1      | 2      | 1      | 0     | 1     |
| 14 | 0      | 2      | 0      | 0     | 2     |

The new algorithm finds long haplotypes (lengths 12, 13, 14) the old one missed entirely. Old required some single sample's bin haplotype to match a length-k sequence exactly to surface it; new enumerates contiguous sub-fragments across all parent blocks and aggregates their containment AC, surfacing long fragments that no single sample produced as its full bin haplotype.

**Old-only — mutually exclusive breakdown** (most missing → least missing)

| Tier | Count | Meaning |
|---|---|---|
| 1. fully variant-disjoint | 930 | No variant in this haplotype appears in any new haplotype. |
| 2. partially missing variants | 415 | At least one variant absent from every new haplotype; others do appear. |
| 3. variants split across new haplotypes | 14 | Every variant is in some new haplotype, but no single new haplotype contains all of them. |
| 4. subset of some new haplotype but not contiguous | 7 | All variants fit inside one new haplotype's variant set, but this old key is not a contiguous sub-array of any new haplotype. |
| 5. proper contiguous sub-fragment of some new haplotype | 1,250 | Would have been emitted but was dropped by the new algorithm's containment dedup (same per-pop AC as an enclosing fragment). |

Tier 5 is the largest bucket (≈48% of old-only) and is the containment-dedup case the new algorithm explicitly drops. Tiers 1–2 (1,345 combined) are dominated by haplotypes whose `estimated_gnomad_AF` cleared the threshold in old only because two-pass-union double-counting inflated metrics; correctly single-pass-counted in new, they fall below the threshold. Tiers 3–4 (21) are residual ordering / non-contiguous selections that old's bin-emission permitted but new's adjacency-contiguous enumeration cannot.

**AC comparison for shared haplotypes (28,265)**

| | Count |
|---|---|
| Same AC          | 24,447 |
| New > old        | 3,773  |
| Old > new        | 45     |

`new > old` is the AC-undercount fix: old's stage-7 `distinct()` collapsed sub-fragments produced from different bin parents and discarded the AC sum; new sums containment AC across all parent blocks. The 45 `old > new` cases are the residual two-pass-union double-count that survived `distinct()` arbitrarily picking the inflated row.

**Summary**

The two algorithms agree on 91% of old haplotypes and 94% of new. The disagreement is structurally accounted for by the three rewrite goals:

1. Containment dedup drops redundant sub-fragments (tier 5: 1,250 old-only).
2. Containment AC summation surfaces sub-fragments old could not emit and corrects under-counted shared rows (1,871 new-only, plus 3,773 shared with new > old).
3. Single-pass evaluation removes two-pass-union double-counting that was passing some haplotypes through the frequency threshold on inflated AC (tiers 1–2: ≈1,345 old-only, plus the 45 shared old > new).
