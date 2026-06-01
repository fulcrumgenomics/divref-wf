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

No AWS credentials are required: the workflow's S3 reads are configured to use
`AnonymousAWSCredentialsProvider` because all inputs (`gnomad-public-us-east-1`,
`broad-references`) are on public Open Data buckets.

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

Run on AWS (after `pixi run setup-s3`; no AWS credentials needed):

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

> [!IMPORTANT]
> **Sample-input assumption: sex-aneuploid samples must not have an assigned population.**
> The haplotype computation assumes all samples that survive the population filter have either an `XX` or `XY` karyotype, so that the chrX non-PAR ploidy correction (males treated as haploid) is well-defined. This invariant is enforced indirectly: samples are kept only when their `gnomad_population_inference.pop` is in the configured set of populations, and the gnomAD 3.1.2 HGDP+1KG sample-meta source declines to assign a population to any sex-aneuploid sample (`X`, `XXY`, `XYY`, `ambiguous`) — their PCA-based pop inference is unreliable on non-diploid sex karyotypes. Aneuploid samples therefore have `pop = None` upstream and get dropped by the population filter before the chrX correction runs.
>
> If you point the workflow at a different sample-metadata source that *does* assign populations to aneuploid samples, those samples will pass the population filter and their genotypes will be fed into the chrX non-PAR correction with incorrect ploidy assumptions. The pipeline does not currently check for this; users supplying custom sample metadata are responsible for filtering aneuploid samples upstream.

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
| `popmax_estimated_gnomad_AF` | For `HGDP_haplotype`: `popmax_fraction_phased × min(gnomAD_AF[component, max_pop])` over the component variants — i.e. the phase ratio applied to the rarest component's gnomAD AF in `max_pop`. For `gnomAD_variant`: the gnomAD AF in `max_pop`. Equivalent to `estimated_gnomAD_haplotype_AF_{max_pop}`. |
| `popmax_fraction_phased` | For `HGDP_haplotype`: phase ratio, `popmax_empirical_AF / min(empirical_AF[component, max_pop])` — empirical haplotype AF over empirical AF of the rarest component variant in the same population. `1.0` for `gnomAD_variant`. Equivalent to `fraction_phased_{max_pop}`. |
| `variants` | Comma-separated list of component variants in `chr:pos:ref:alt` form, in the order they appear in the haplotype. |
| `gnomAD_AF_{POP}` | One column per configured population (e.g. `gnomAD_AF_AFR`). Comma-separated per-component gnomAD AFs in that population, in the same order as `variants`, formatted to 5 decimal places. |
| `empirical_AC_{POP}` | One column per pop in `joint_pops_legend`. For `HGDP_haplotype`: count of chromosomes carrying the haplotype in that pop (haplotype-level). For `gnomAD_variant`: the gnomAD AC for that variant in that pop (variant-level). Missing if the pop has no source data for this row. |
| `empirical_AF_{POP}` | For `HGDP_haplotype`: empirical AF of the haplotype in that pop, derived from `empirical_AC_{POP}` and the min AN over component variants. For `gnomAD_variant`: the gnomAD AF for that variant in that pop. Missing if AN is 0 or the pop has no source data. |
| `fraction_phased_{POP}` | Per-pop phase ratio: `empirical_AF_{POP} / min(local_call_stats_AF[component, POP])`. `1.0` for `gnomAD_variant` rows. Uses *that pop's own* denominators (not `max_pop`'s). |
| `estimated_gnomAD_haplotype_AF_{POP}` | Per-pop projection: `fraction_phased_{POP} × min(gnomAD_AF[component, POP])`. Equals `empirical_AF_{POP}` for `gnomAD_variant` rows. |
| `haplotype_filter` | VCF-style compatibility flag. `PASS` for `gnomAD_variant` rows and for `HGDP_haplotype` rows whose component variants do not overlap; otherwise the `;`-joined reason(s) the haplotype is incompatible (e.g. `snp_in_deletion`, `overlapping_deletions`, `same_position`, `indel_in_deletion`). Incompatible haplotypes are component variants that cannot co-occur on one chromosome — upstream phasing artifacts at tandem repeats — and are flagged rather than dropped. |

The DuckDB file also contains three single-row metadata tables: `window_size` (the flanking context size used), `pops_legend` (JSON-encoded ordered population list), and `VERSION`.

## Analysis

The analyses behind this workflow are written up in the project blog post, [`docs/blog.md`](docs/blog.md):

- why DivRef 1.1's `gnomAD_variant` sites come from gnomAD 3.1.2 HGDP+1KG rather than the documented gnomAD 4.1 release, and how the two releases differ;
- the rewrite of the haplotype computation algorithm from two-pass binning to per-sample adjacency with containment counting;
- a whole-genome comparison of the original and new algorithms over the original DivRef 1.1 DuckDB index versus the new one.

Every numeric claim, table, and figure in the blog maps to a reproducible command in [`docs/blog_reproduction.md`](docs/blog_reproduction.md), which also documents the additional `analysis` pixi environment those steps need.

