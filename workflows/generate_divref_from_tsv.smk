####################################################################################################
# Generates a DivRef-format resource from one or more wide single-variant TSVs (pure Python; no
# Hail/Spark).
#
# Final output is, per source: a DuckDB index and a set of per-chromosome FASTA files.
####################################################################################################

import os
from pathlib import Path
from snakemake.exceptions import WorkflowError
from snakemake.utils import validate

####################################################################################################
# Inputs
####################################################################################################
#
# No default `configfile:` is provided; pass one explicitly, e.g.
# `--configfile workflows/config/config_tsv_test.yml`.
#
# `config_tsv_gcs.yml` / `config_tsv_aws.yml` set only the reference-genome URI; layer one AFTER a
# sources-bearing config, e.g. `--configfile my_sources.yml --configfile config_tsv_aws.yml`.


validate(config, os.path.join(workflow.basedir, "config", "config_tsv_schema.yml"))

WORK_DIR: Path = Path(config["work_dir"])
TMP_DIR: Path = Path(config["tmp_dir"])

CHROMS: list[str] = config["chromosomes"]

REFERENCE_GENOME: str = config["reference_genome_base_name"]
REFERENCE_GENOME_URI: str = config["reference_genome_uri"]

SEQUENCE_WINDOW_SIZE: int = config["sequence_window_size"]
POLARS_CHUNK_SIZE: int = config["polars_chunk_size"]

SOURCE_BY_NAME: dict[str, dict[str, str]] = {source["name"]: source for source in config["sources"]}
# The schema's `uniqueItems` only rejects fully-identical entries; two sources sharing a `name`
# would otherwise collapse silently and one would be dropped from the build.
if len(SOURCE_BY_NAME) != len(config["sources"]):
    raise WorkflowError("Each source in `sources` must have a unique `name`.")
SOURCE_NAMES: list[str] = list(SOURCE_BY_NAME.keys())

# Run every shell rule under strict bash so a failing command aborts the rule instead of being
# masked. `run:` directives are unaffected.
shell.prefix("set -euo pipefail; ")

####################################################################################################
# Rules
####################################################################################################


rule all:
    input:
        expand(f"{WORK_DIR}/output/{{source}}.index.duckdb", source=SOURCE_NAMES),
        expand(
            f"{WORK_DIR}/output/{{source}}.{{chrom}}.fasta",
            source=SOURCE_NAMES,
            chrom=CHROMS,
        ),


####################################################################################################
# Downloads and unzips the reference genome.
####################################################################################################
rule download_reference_genome:
    output:
        fasta=f"{WORK_DIR}/inputs/{REFERENCE_GENOME}.fasta",
    log:
        "logs/generate_divref_from_tsv/download_reference_genome.log",
    params:
        fasta_uri=REFERENCE_GENOME_URI,
    shell:
        """
        (
            uri="{params.fasta_uri}"
            # Download to a generic path; gzip is detected via magic bytes after fetch.
            dl="{output.fasta}.download"
            case "$uri" in
                s3://*|s3a://*)
                    # `aws s3` uses the s3:// scheme; rewrite a leading s3a:// if present.
                    s3_uri="${{uri/#s3a:/s3:}}"
                    # Try authenticated first; fall back to --no-sign-request for public
                    # Open Data buckets when no AWS credentials are configured.
                    aws s3 cp "$s3_uri" "$dl" \
                        || aws s3 cp --no-sign-request "$s3_uri" "$dl"
                    ;;
                gs://*)
                    gsutil -m cp "$uri" "$dl"
                    ;;
                *)
                    echo "Unsupported reference_genome_uri scheme: $uri" >&2
                    exit 1
                    ;;
            esac
            # Detect gzip from the file's magic bytes (\\x1f\\x8b) rather than the URI suffix,
            # so misnamed objects are handled correctly.
            magic=$(head -c 2 "$dl" | od -An -tx1 | tr -d ' \\n')
            if [[ "$magic" == "1f8b" ]]; then
                mv "$dl" "{output.fasta}.gz"
                gunzip "{output.fasta}.gz"
            else
                mv "$dl" "{output.fasta}"
            fi
        ) &> {log}
        """


####################################################################################################
# Indexes the reference genome.
####################################################################################################
rule index_reference_genome:
    input:
        fasta=f"{WORK_DIR}/inputs/{REFERENCE_GENOME}.fasta",
    output:
        fai=f"{WORK_DIR}/inputs/{REFERENCE_GENOME}.fai",
    log:
        "logs/generate_divref_from_tsv/index_reference_genome.log",
    shell:
        # Force the `.fai` name (not samtools' default `<fasta>.fasta.fai`) to match the
        # `reference_fasta.with_suffix(".fai")` lookup in create_duckdb_from_tsv.
        """
        (
            samtools faidx \
                {input.fasta} \
                --output {output.fai}
        ) &> {log}
        """


####################################################################################################
# Build one source's DivRef DuckDB index directly from its wide single-variant TSV (no Hail).
####################################################################################################
rule create_tsv_index:
    input:
        variants=lambda wc: SOURCE_BY_NAME[wc.source]["variants_tsv"],
        source_meta=lambda wc: SOURCE_BY_NAME[wc.source]["source_meta"],
        fasta=f"{WORK_DIR}/inputs/{REFERENCE_GENOME}.fasta",
        fai=f"{WORK_DIR}/inputs/{REFERENCE_GENOME}.fai",
    output:
        duckdb=f"{WORK_DIR}/output/{{source}}.index.duckdb",
    log:
        "logs/generate_divref_from_tsv/create_tsv_index.{source}.log",
    params:
        output_base=f"{WORK_DIR}/output/{{source}}",
        window_size=SEQUENCE_WINDOW_SIZE,
        polars_chunk_size=POLARS_CHUNK_SIZE,
        contigs=" ".join(CHROMS),
        tmp_dir=TMP_DIR,
    shell:
        # The tool requires --tmp-dir to already exist; Snakemake only auto-creates output dirs.
        """
        mkdir -p {params.tmp_dir}
        divref create-duckdb-from-tsv \
            --variants-tsv {input.variants} \
            --source-meta {input.source_meta} \
            --output-base {params.output_base} \
            --reference-fasta {input.fasta} \
            --window-size {params.window_size} \
            --contigs {params.contigs} \
            --polars-chunk-size {params.polars_chunk_size} \
            --tmp-dir {params.tmp_dir} \
            --force \
            &> {log}
        """


####################################################################################################
# Write one source's per-chromosome FASTA files from its DivRef DuckDB index.
####################################################################################################
rule create_tsv_fasta:
    input:
        duckdb=f"{WORK_DIR}/output/{{source}}.index.duckdb",
    output:
        fasta=f"{WORK_DIR}/output/{{source}}.{{chrom}}.fasta",
    log:
        "logs/generate_divref_from_tsv/create_tsv_fasta.{source}.{chrom}.log",
    wildcard_constraints:
        chrom=r"chr([1-9]|1[0-9]|2[0-2]|X|Y)",
    params:
        output_base=f"{WORK_DIR}/output/{{source}}",
    shell:
        """
        divref create-divref-fasta \
            --duckdb-path {input.duckdb} \
            --output-base {params.output_base} \
            --contigs {wildcards.chrom} \
            &> {log}
        """
