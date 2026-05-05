####################################################################################################
# Generates a DivRef-format resource of human haplotypes.
#
# Final output is a set of per-chromosome FASTA files and a DuckDB index.
####################################################################################################

import os
from pathlib import Path
from snakemake.utils import validate

####################################################################################################
# Inputs
####################################################################################################


configfile: os.path.join(workflow.basedir, "config", "config.yml")


validate(config, os.path.join(workflow.basedir, "config", "config_schema.yml"))

VERSION: str = config["version"]

WORK_DIR: Path = Path(config["work_dir"])
TMP_DIR: Path = Path(config["tmp_dir"])

CHROMS: list[str] = config["chromosomes"]
POPS: list[str] = config["populations"]

REFERENCE_GENOME: str = config["reference_genome_base_name"]
REFERENCE_GENOME_URI: str = config["reference_genome_uri"]

# The HGDP+1KG phased BCF files are at
# "{HGDP_1KG_PHASED_BCF_PREFIX}.{chrom}.{HGDP_1KG_PHASED_BCF_SUFFIX}"
HGDP_1KG_PHASED_BCF_PREFIX: str = config["hgdp_1kg_phased_bcf_prefix"]
HGDP_1KG_PHASED_BCF_SUFFIX: str = config["hgdp_1kg_phased_bcf_suffix"]
HGDP_1KG_VARIANT_ANNOTATION_HAIL_TABLE: str = config["hgdp_1kg_variant_annotation_hail_table"]
HGDP_1KG_SAMPLE_METADATA_HAIL_TABLE: str = config["hgdp_1kg_sample_metadata_hail_table"]
HGDP_1KG_MIN_POP_VARIANT_AF: float = config["hgdp_1kg_min_pop_variant_allele_freq"]
HGDP_1KG_MIN_POP_HAPLOTYPE_AF: float = config["hgdp_1kg_min_estimated_gnomad_haplotype_allele_freq"]
HGDP_1KG_HAPLOTYPE_WINDOW_SIZE: int = config["hgdp_1kg_haplotype_window_size"]

# gnomAD variants can be from a different source than the haplotypes
GNOMAD_VARIANT_ANNOTATION_SOURCE: str = config["gnomad_variant_annotation_source"]
GNOMAD_VARIANT_ANNOTATION_CLOUD: str = config["gnomad_variant_annotation_cloud"]
GNOMAD_VARIANT_MIN_POP_VARIANT_AF: float = config["gnomad_variant_min_pop_variant_allele_freq"]

SEQUENCE_WINDOW_SIZE: int = config["sequence_window_size"]
POLARS_CHUNK_SIZE: int = config["polars_chunk_size"]

SPARK_DRIVER_MEMORY_GB: int = config["spark_driver_memory_gb"]
SPARK_EXECUTOR_MEMORY_GB: int = config["spark_executor_memory_gb"]

VCF_EXTS: list[str] = [".vcf.gz", ".vcf.gz.tbi"]

####################################################################################################
# Rules
####################################################################################################


rule all:
    input:
        f"{WORK_DIR}/output/hgdp_1kg.haplotypes_gnomad_merge.index.duckdb",
        expand(
            f"{WORK_DIR}/output/hgdp_1kg.haplotypes_gnomad_merge.{{chrom}}.fasta",
            chrom=CHROMS,
        ),


####################################################################################################
# Extracts the phased genotypes for all HGDP+1KG samples in the specified locus.
#
# Removes the INFO field, since this is not required for haplotype computation (allele frequencies
# are re-annotated from the sites table), and it inflates the size on disk and subsequently the time
# for Hail to load and parse the VCF with the `divref compute-haplotypes` tool.
####################################################################################################
rule subset_phased_genotypes:
    output:
        vcf=f"{WORK_DIR}/inputs/hgdp_1kg.phased_genotypes.{{chrom}}.vcf.gz",
        tbi=f"{WORK_DIR}/inputs/hgdp_1kg.phased_genotypes.{{chrom}}.vcf.gz.tbi",
    log:
        "logs/generate_divref/subset_phased_genotypes.{chrom}.log",
    params:
        bcf=f"{HGDP_1KG_PHASED_BCF_PREFIX}{{chrom}}{HGDP_1KG_PHASED_BCF_SUFFIX}",
    shell:
        """
        (
            bcftools annotate \
                --remove INFO \
                --output-type z \
                --output {output.vcf} \
                --write-index=tbi \
                {params.bcf}
        ) &> {log}
        """


####################################################################################################
# Extracts allele frequencies from the HGDP+1KG gnomAD subset for the given populations and subsets
# to sites over the specified minimum allele frequency in at least one population.
####################################################################################################
rule extract_gnomad_afs:
    output:
        variant_ht=directory(f"{WORK_DIR}/inputs/hgdp_1kg.sites.{{chrom}}.ht"),
    log:
        "logs/generate_divref/extract_gnomad_afs.{chrom}.log",
    params:
        variant_ht=HGDP_1KG_VARIANT_ANNOTATION_HAIL_TABLE,
        freq_threshold=HGDP_1KG_MIN_POP_VARIANT_AF,
        populations=" ".join(POPS),
        spark_driver_memory_gb=SPARK_DRIVER_MEMORY_GB,
        spark_executor_memory_gb=SPARK_EXECUTOR_MEMORY_GB,
    shell:
        """
        (
            divref extract-gnomad-afs \
                --in-gnomad-sites-table {params.variant_ht} \
                --out-variant-annotation-table {output.variant_ht} \
                --contig {wildcards.chrom} \
                --freq-threshold {params.freq_threshold} \
                --populations {params.populations} \
                --spark-driver-memory-gb {params.spark_driver_memory_gb} \
                --spark-executor-memory-gb {params.spark_executor_memory_gb}
        ) &> {log}
        """


####################################################################################################
# Extracts selected fields from HGDP+1KG sample metadata.
####################################################################################################
rule extract_sample_metadata:
    output:
        sample_ht=directory(f"{WORK_DIR}/inputs/hgdp_1kg.sample_metadata.ht"),
    log:
        "logs/generate_divref/extract_sample_metadata.log",
    params:
        sample_ht=HGDP_1KG_SAMPLE_METADATA_HAIL_TABLE,
        spark_driver_memory_gb=SPARK_DRIVER_MEMORY_GB,
        spark_executor_memory_gb=SPARK_EXECUTOR_MEMORY_GB,
    shell:
        """
        (
            divref extract-sample-metadata \
                --in-gnomad-hgdp-sample-data {params.sample_ht} \
                --out-sample-metadata {output.sample_ht} \
                --spark-driver-memory-gb {params.spark_driver_memory_gb} \
                --spark-executor-memory-gb {params.spark_executor_memory_gb}
        ) &> {log}
        """


####################################################################################################
# Compute haplotypes from the HGDP+1KG filtered sites, sample metadata, and phased genotypes.
####################################################################################################
rule compute_haplotypes:
    input:
        vcf=f"{WORK_DIR}/inputs/hgdp_1kg.phased_genotypes.{{chrom}}.vcf.gz",
        tbi=f"{WORK_DIR}/inputs/hgdp_1kg.phased_genotypes.{{chrom}}.vcf.gz.tbi",
        variant_ht=f"{WORK_DIR}/inputs/hgdp_1kg.sites.{{chrom}}.ht",
        sample_ht=f"{WORK_DIR}/inputs/hgdp_1kg.sample_metadata.ht",
    output:
        haplotypes_ht=directory(f"{WORK_DIR}/haplotypes/hgdp_1kg.haplotypes.{{chrom}}.ht"),
    log:
        "logs/generate_divref/compute_haplotypes.{chrom}.log",
    params:
        window_size=HGDP_1KG_HAPLOTYPE_WINDOW_SIZE,
        variant_freq_threshold=HGDP_1KG_MIN_POP_VARIANT_AF,
        haplotype_freq_threshold=HGDP_1KG_MIN_POP_HAPLOTYPE_AF,
        output_base=f"{WORK_DIR}/haplotypes/hgdp_1kg.haplotypes.{{chrom}}",
        spark_driver_memory_gb=SPARK_DRIVER_MEMORY_GB,
        spark_executor_memory_gb=SPARK_EXECUTOR_MEMORY_GB,
    shell:
        """
        (
            divref compute-haplotypes \
                --vcfs-path {input.vcf} \
                --gnomad-va-file {input.variant_ht} \
                --gnomad-sa-file {input.sample_ht} \
                --window-size {params.window_size} \
                --variant-freq-threshold {params.variant_freq_threshold} \
                --haplotype-freq-threshold {params.haplotype_freq_threshold} \
                --output-base {params.output_base} \
                --spark-driver-memory-gb {params.spark_driver_memory_gb} \
                --spark-executor-memory-gb {params.spark_executor_memory_gb}

            # remove intermediate files
            rm -r {params.output_base}.[12].ht {params.output_base}.variants.ht
        ) &> {log}
        """


####################################################################################################
# Extracts allele frequencies from the gnomAD sites table for the given populations and subsets to
# sites over the specified minimum allele frequency in at least one population.
####################################################################################################
rule extract_gnomad_variant_afs:
    output:
        variant_ht=directory(f"{WORK_DIR}/inputs/gnomad.sites.{{chrom}}.ht"),
    log:
        "logs/generate_divref/extract_gnomad_variant_afs.{chrom}.log",
    params:
        gnomad_source=GNOMAD_VARIANT_ANNOTATION_SOURCE,
        gnomad_cloud=GNOMAD_VARIANT_ANNOTATION_CLOUD,
        freq_threshold=GNOMAD_VARIANT_MIN_POP_VARIANT_AF,
        populations=" ".join(POPS),
        spark_driver_memory_gb=SPARK_DRIVER_MEMORY_GB,
        spark_executor_memory_gb=SPARK_EXECUTOR_MEMORY_GB,
    shell:
        """
        (
            divref extract-gnomad-single-afs \
                --gnomad-version {params.gnomad_source} \
                --gnomad-cloud {params.gnomad_cloud} \
                --contig {wildcards.chrom} \
                --freq-threshold {params.freq_threshold} \
                --out-sites-hail-table {output.variant_ht} \
                --populations {params.populations} \
                --spark-driver-memory-gb {params.spark_driver_memory_gb} \
                --spark-executor-memory-gb {params.spark_executor_memory_gb}
        ) &> {log}
        """


####################################################################################################
# Downloads and unzips the reference genome.
####################################################################################################
rule download_reference_genome:
    output:
        fasta=f"{WORK_DIR}/inputs/{REFERENCE_GENOME}.fasta",
    log:
        "logs/generate_divref/download_reference_genome.log",
    params:
        fasta_uri=REFERENCE_GENOME_URI,
    shell:
        """
        (
            set -euo pipefail
            uri="{params.fasta_uri}"
            # Strip a trailing .gz to determine whether the source is gzipped.
            if [[ "$uri" == *.gz ]]; then
                dl="{output.fasta}.gz"
            else
                dl="{output.fasta}"
            fi
            case "$uri" in
                s3://*|s3a://*)
                    aws s3 cp "${{uri/s3a:\\/\\//s3:\\/\\/}}" "$dl"
                    ;;
                gs://*)
                    gsutil -m cp "$uri" "$dl"
                    ;;
                *)
                    echo "Unsupported reference_genome_uri scheme: $uri" >&2
                    exit 1
                    ;;
            esac
            if [[ "$dl" == *.gz ]]; then
                gunzip "$dl"
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
        "logs/generate_divref/index_reference_genome.log",
    shell:
        """
        (
            samtools faidx \
                {input.fasta} \
                --output {output.fai}
        ) &> {log}
        """


####################################################################################################
# Writes a TSV listing per-chromosome haplotype and gnomAD sites Hail tables for the index builder.
####################################################################################################
rule create_table_pairs_tsv:
    input:
        haplotypes_hts=expand(
            f"{WORK_DIR}/haplotypes/hgdp_1kg.haplotypes.{{chrom}}.ht",
            chrom=CHROMS,
        ),
        sites_hts=expand(
            f"{WORK_DIR}/inputs/gnomad.sites.{{chrom}}.ht",
            chrom=CHROMS,
        ),
    output:
        tsv=f"{WORK_DIR}/inputs/table_pairs.tsv",
    run:
        with open(output.tsv, "w") as f:
            f.write("contig\thaplotype_table_path\tsites_table_path\n")
            for chrom in CHROMS:
                haplotype_ht = f"{WORK_DIR}/haplotypes/hgdp_1kg.haplotypes.{chrom}.ht"
                sites_ht = f"{WORK_DIR}/inputs/gnomad.sites.{chrom}.ht"
                f.write(f"{chrom}\t{haplotype_ht}\t{sites_ht}\n")


####################################################################################################
# Build the DivRef DuckDB index from all per-chromosome haplotype and gnomAD sites Hail tables.
####################################################################################################
rule create_divref_index:
    input:
        table_pairs_tsv=f"{WORK_DIR}/inputs/table_pairs.tsv",
        fasta=f"{WORK_DIR}/inputs/{REFERENCE_GENOME}.fasta",
        fai=f"{WORK_DIR}/inputs/{REFERENCE_GENOME}.fai",
    output:
        duckdb=f"{WORK_DIR}/output/hgdp_1kg.haplotypes_gnomad_merge.index.duckdb",
    log:
        "logs/generate_divref/create_divref_index.log",
    params:
        window_size=SEQUENCE_WINDOW_SIZE,
        output_base=f"{WORK_DIR}/output/hgdp_1kg",
        version=VERSION,
        polars_chunk_size=POLARS_CHUNK_SIZE,
        tmp_dir=TMP_DIR,
        spark_driver_memory_gb=SPARK_DRIVER_MEMORY_GB,
        spark_executor_memory_gb=SPARK_EXECUTOR_MEMORY_GB,
    shell:
        """
        (
            divref create-duckdb-index \
                --in-table-pairs-tsv {input.table_pairs_tsv} \
                --reference-fasta {input.fasta} \
                --window-size {params.window_size} \
                --output-base {params.output_base} \
                --version {params.version} \
                --polars-chunk-size {params.polars_chunk_size} \
                --tmp-dir {params.tmp_dir} \
                --spark-driver-memory-gb {params.spark_driver_memory_gb} \
                --spark-executor-memory-gb {params.spark_executor_memory_gb}
        ) &> {log}
        """


####################################################################################################
# Write per-chromosome FASTA files from the DivRef DuckDB index.
####################################################################################################
rule create_divref_fasta:
    input:
        duckdb=f"{WORK_DIR}/output/hgdp_1kg.haplotypes_gnomad_merge.index.duckdb",
    output:
        fastas=expand(
            f"{WORK_DIR}/output/hgdp_1kg.haplotypes_gnomad_merge.{{chrom}}.fasta",
            chrom=CHROMS,
        ),
    log:
        "logs/generate_divref/create_divref_fasta.log",
    params:
        output_base=f"{WORK_DIR}/output/hgdp_1kg.haplotypes_gnomad_merge",
        contigs=" ".join(CHROMS),
    shell:
        """
        (
            divref create-divref-fasta \
                --duckdb-path {input.duckdb} \
                --output-base {params.output_base} \
                --contigs {params.contigs}
        ) &> {log}
        """
