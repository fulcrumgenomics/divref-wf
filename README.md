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

The resource that is generated is described in [`docs/resource_description.md`](docs/resource_description.md).

## Analysis

The analyses behind this workflow are written up in the project blog post, [`docs/blog.md`](docs/blog.md):

- why DivRef 1.1's `gnomAD_variant` sites come from gnomAD 3.1.2 HGDP+1KG rather than the documented gnomAD 4.1 release, and how the two releases differ;
- the rewrite of the haplotype computation algorithm from two-pass binning to per-sample adjacency with containment counting;
- a whole-genome comparison of the original and new algorithms over the original DivRef 1.1 DuckDB index versus the new one.

Every numeric claim, table, and figure in the blog maps to a reproducible command in [`docs/blog_reproduction.md`](docs/blog_reproduction.md), which also documents the additional `analysis` pixi environment those steps need.

