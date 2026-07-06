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
#
# No default `configfile:` is provided; pass one explicitly with
# `--configfile workflows/config/config_gcs.yml` (GCS sources) or
# `--configfile workflows/config/config_aws.yml` (S3 / AWS Open Data sources).


validate(config, os.path.join(workflow.basedir, "config", "config_schema.yml"))

CLOUD: str = config["cloud"]
_CLOUD_SCHEMES: dict[str, tuple[str, ...]] = {
    "GCS": ("gs://",),
    "AWS": ("s3://", "s3a://"),
}
# `cloud` uses GCS/AWS; the gnomAD single-AFs tool's `gnomad_cloud` enum uses GCS/S3.
_GNOMAD_CLOUD_FOR: dict[str, str] = {"GCS": "GCS", "AWS": "S3"}


def _validate_cloud_uri(field: str, uri: str) -> None:
    expected = _CLOUD_SCHEMES[CLOUD]
    if not any(uri.startswith(scheme) for scheme in expected):
        raise ValueError(
            f"Config field '{field}' has URI {uri!r} which does not match cloud "
            f"{CLOUD!r} (expected scheme one of {expected})."
        )


for _field in (
    "reference_genome_uri",
    "hgdp_1kg_phased_bcf_prefix",
    "hgdp_1kg_variant_annotation_hail_table",
    "hgdp_1kg_sample_metadata_hail_table",
):
    _validate_cloud_uri(_field, config[_field])

for _chrx_part in ("par1", "non_par", "par2"):
    _validate_cloud_uri(
        f"hgdp_1kg_phased_bcf_chrX.{_chrx_part}",
        config["hgdp_1kg_phased_bcf_chrX"][_chrx_part],
    )


VERSION: str = config["version"]

WORK_DIR: Path = Path(config["work_dir"])
TMP_DIR: Path = Path(config["tmp_dir"])

CHROMS: list[str] = config["chromosomes"]
# Haplotypes are computed for autosomes + chrX (chrX non-PAR uses the haploid-male ploidy
# correction in `divref compute-haplotypes`); chrY contributes single gnomAD variants only.
_HAPLOTYPE_CONTIGS: frozenset[str] = frozenset(f"chr{n}" for n in range(1, 23)) | {"chrX"}
HAPLOTYPE_CHROMS: list[str] = [c for c in CHROMS if c in _HAPLOTYPE_CONTIGS]

# Constrains `{chrom}` wildcards to valid GRCh38 main-contig tokens so they cannot greedily match
# unintended path segments (e.g. an intermediate filename suffix that embeds a contig).
_CONTIG_WILDCARD_REGEX: str = r"chr(\d+|X|Y)"

REFERENCE_GENOME: str = config["reference_genome_base_name"]
REFERENCE_GENOME_URI: str = config["reference_genome_uri"]

# The HGDP+1KG phased BCF files are at
# "{HGDP_1KG_PHASED_BCF_PREFIX}.{chrom}.{HGDP_1KG_PHASED_BCF_SUFFIX}"
HGDP_1KG_PHASED_BCF_PREFIX: str = config["hgdp_1kg_phased_bcf_prefix"]
HGDP_1KG_PHASED_BCF_SUFFIX: str = config["hgdp_1kg_phased_bcf_suffix"]
# chrX uses three BCFs with non-uniform naming (PAR1, non-PAR, PAR2); see
# `subset_phased_genotypes_chrX` below.
HGDP_1KG_PHASED_BCF_CHRX: dict[str, str] = config["hgdp_1kg_phased_bcf_chrX"]
HGDP_1KG_VARIANT_ANNOTATION_HAIL_TABLE: str = config["hgdp_1kg_variant_annotation_hail_table"]
HGDP_1KG_SAMPLE_METADATA_HAIL_TABLE: str = config["hgdp_1kg_sample_metadata_hail_table"]
HGDP_1KG_POPS: list[str] = config["hgdp_1kg_populations"]
HGDP_1KG_MIN_POP_VARIANT_AF: float = config["hgdp_1kg_min_pop_variant_allele_freq"]
HGDP_1KG_MIN_POP_HAPLOTYPE_AF: float = config["hgdp_1kg_min_estimated_gnomad_haplotype_allele_freq"]

# gnomAD variants can be from a different source than the haplotypes; the cloud is
# derived from the workflow-level `cloud` so all inputs come from the same provider.
GNOMAD_VARIANT_ANNOTATION_SOURCE: str = config["gnomad_variant_annotation_source"]
GNOMAD_VARIANT_ANNOTATION_CLOUD: str = _GNOMAD_CLOUD_FOR[CLOUD]
GNOMAD_VARIANT_POPS: list[str] = config["gnomad_variant_populations"]
GNOMAD_VARIANT_MIN_POP_VARIANT_AF: float = config["gnomad_variant_min_pop_variant_allele_freq"]

SEQUENCE_WINDOW_SIZE: int = config["sequence_window_size"]
POLARS_CHUNK_SIZE: int = config["polars_chunk_size"]

SPARK_DRIVER_MEMORY_GB: int = config["spark_driver_memory_gb"]
SPARK_EXECUTOR_MEMORY_GB: int = config["spark_executor_memory_gb"]

# Memory footprint and a core-count hint for the Hail/Spark rules, so `-j` / `--resources mem_mb=...`
# can bound how many of these heavyweight JVMs run concurrently (each launches a driver and an
# executor JVM at the configured per-JVM heap, plus ~2GB for the Python process and off-heap
# buffers). Spark runs in local mode and grabs cores irrespective of Snakemake, so `threads` is a
# scheduling hint that keeps the scheduler from packing several of these jobs onto one node.
_HAIL_RULE_MEM_MB: int = (SPARK_DRIVER_MEMORY_GB + SPARK_EXECUTOR_MEMORY_GB) * 1024 + 2048
_HAIL_RULE_THREADS: int = 4

# Built once so each Hail rule can append it to its CLI invocation. Empty when on AWS,
# since the GCS connector is not loaded and the credentials path is ignored.
GCS_CREDENTIALS_FLAG: str = (
    f"--gcs-credentials-path '{config['gcs_credentials_path']}'" if CLOUD == "GCS" else ""
)

VCF_EXTS: list[str] = [".vcf.gz", ".vcf.gz.tbi"]

# Run every shell rule under strict bash so a failing command aborts the rule instead of being
# masked. Without this the create_divref_index per-contig append loop would continue past a failed
# contig (and still run finalize), and compute_haplotypes would run its intermediate-file cleanup
# even if the tool failed. Subshells `( ... )` inherit these flags. `run:` directives are unaffected.
shell.prefix("set -euo pipefail; ")

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
#
# chrX uses three separate phased BCFs (PAR1, non-PAR, PAR2) — see `subset_phased_genotypes_chrX`.
####################################################################################################
rule subset_phased_genotypes:
    output:
        vcf=f"{WORK_DIR}/inputs/hgdp_1kg.phased_genotypes.{{chrom}}.vcf.gz",
        tbi=f"{WORK_DIR}/inputs/hgdp_1kg.phased_genotypes.{{chrom}}.vcf.gz.tbi",
    log:
        "logs/generate_divref/subset_phased_genotypes.{chrom}.log",
    wildcard_constraints:
        # Autosomes + chrY only — chrX has its own subset rule below.
        chrom=r"chr(\d+|Y)",
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
# Extracts and concatenates the chrX phased genotypes from the three HGDP+1KG BCFs (PAR1, non-PAR,
# PAR2). The three regions are disjoint on GRCh38 — PAR1 ends well before non-PAR starts, and PAR2
# starts after non-PAR ends — so plain `bcftools concat` stitches them in genomic order without
# needing `--allow-overlaps`. The INFO field is dropped on output, matching the autosome subset rule.
####################################################################################################
rule subset_phased_genotypes_chrX:
    output:
        vcf=f"{WORK_DIR}/inputs/hgdp_1kg.phased_genotypes.chrX.vcf.gz",
        tbi=f"{WORK_DIR}/inputs/hgdp_1kg.phased_genotypes.chrX.vcf.gz.tbi",
    log:
        "logs/generate_divref/subset_phased_genotypes.chrX.log",
    params:
        par1=HGDP_1KG_PHASED_BCF_CHRX["par1"],
        non_par=HGDP_1KG_PHASED_BCF_CHRX["non_par"],
        par2=HGDP_1KG_PHASED_BCF_CHRX["par2"],
    shell:
        """
        (
            bcftools concat \
                --output-type u \
                {params.par1} {params.non_par} {params.par2} \
            | bcftools annotate \
                --remove INFO \
                --output-type z \
                --output {output.vcf} \
                --write-index=tbi \
                -
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
    wildcard_constraints:
        chrom=_CONTIG_WILDCARD_REGEX,
    threads: _HAIL_RULE_THREADS
    resources:
        mem_mb=_HAIL_RULE_MEM_MB,
    params:
        variant_ht=HGDP_1KG_VARIANT_ANNOTATION_HAIL_TABLE,
        freq_threshold=HGDP_1KG_MIN_POP_VARIANT_AF,
        populations=" ".join(HGDP_1KG_POPS),
        spark_driver_memory_gb=SPARK_DRIVER_MEMORY_GB,
        spark_executor_memory_gb=SPARK_EXECUTOR_MEMORY_GB,
        use_s3_flag="--use-s3" if CLOUD == "AWS" else "--no-use-s3",
        gcs_credentials_flag=GCS_CREDENTIALS_FLAG,
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
                --spark-executor-memory-gb {params.spark_executor_memory_gb} \
                {params.use_s3_flag} {params.gcs_credentials_flag}
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
    threads: _HAIL_RULE_THREADS
    resources:
        mem_mb=_HAIL_RULE_MEM_MB,
    params:
        sample_ht=HGDP_1KG_SAMPLE_METADATA_HAIL_TABLE,
        spark_driver_memory_gb=SPARK_DRIVER_MEMORY_GB,
        spark_executor_memory_gb=SPARK_EXECUTOR_MEMORY_GB,
        use_s3_flag="--use-s3" if CLOUD == "AWS" else "--no-use-s3",
        gcs_credentials_flag=GCS_CREDENTIALS_FLAG,
    shell:
        """
        (
            divref extract-sample-metadata \
                --in-gnomad-hgdp-sample-data {params.sample_ht} \
                --out-sample-metadata {output.sample_ht} \
                --spark-driver-memory-gb {params.spark_driver_memory_gb} \
                --spark-executor-memory-gb {params.spark_executor_memory_gb} \
                {params.use_s3_flag} {params.gcs_credentials_flag}
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
    wildcard_constraints:
        chrom=_CONTIG_WILDCARD_REGEX,
    threads: _HAIL_RULE_THREADS
    resources:
        mem_mb=_HAIL_RULE_MEM_MB,
    params:
        window_size=SEQUENCE_WINDOW_SIZE,
        variant_freq_threshold=HGDP_1KG_MIN_POP_VARIANT_AF,
        haplotype_freq_threshold=HGDP_1KG_MIN_POP_HAPLOTYPE_AF,
        output_base=f"{WORK_DIR}/haplotypes/hgdp_1kg.haplotypes.{{chrom}}",
        spark_driver_memory_gb=SPARK_DRIVER_MEMORY_GB,
        spark_executor_memory_gb=SPARK_EXECUTOR_MEMORY_GB,
        min_partitions=config["haplotype_min_partitions"],
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
                --spark-executor-memory-gb {params.spark_executor_memory_gb} \
                --min-partitions {params.min_partitions}

            # remove intermediate files
            rm -r {params.output_base}.variants.ht \
                  {params.output_base}.blocks.ht \
                  {params.output_base}.parents.ht \
                  {params.output_base}.hap_ac.ht
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
    wildcard_constraints:
        chrom=_CONTIG_WILDCARD_REGEX,
    threads: _HAIL_RULE_THREADS
    resources:
        mem_mb=_HAIL_RULE_MEM_MB,
    params:
        gnomad_source=GNOMAD_VARIANT_ANNOTATION_SOURCE,
        gnomad_cloud=GNOMAD_VARIANT_ANNOTATION_CLOUD,
        freq_threshold=GNOMAD_VARIANT_MIN_POP_VARIANT_AF,
        populations=" ".join(GNOMAD_VARIANT_POPS),
        spark_driver_memory_gb=SPARK_DRIVER_MEMORY_GB,
        spark_executor_memory_gb=SPARK_EXECUTOR_MEMORY_GB,
        gcs_credentials_flag=GCS_CREDENTIALS_FLAG,
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
                --spark-executor-memory-gb {params.spark_executor_memory_gb} \
                {params.gcs_credentials_flag}
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
        "logs/generate_divref/index_reference_genome.log",
    shell:
        # Force the `.fai` name (not samtools' default `<fasta>.fasta.fai`) to match the
        # `reference_fasta.with_suffix(".fai")` lookup in append_contig_to_duckdb_index.
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
            chrom=HAPLOTYPE_CHROMS,
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
                sites_ht = f"{WORK_DIR}/inputs/gnomad.sites.{chrom}.ht"
                if chrom in HAPLOTYPE_CHROMS:
                    haplotype_ht = f"{WORK_DIR}/haplotypes/hgdp_1kg.haplotypes.{chrom}.ht"
                else:
                    haplotype_ht = ""
                f.write(f"{chrom}\t{haplotype_ht}\t{sites_ht}\n")


####################################################################################################
# Build the DivRef DuckDB index from all per-chromosome haplotype and gnomAD sites Hail tables.
#
# Built one contig at a time in a serial loop: `init-duckdb-index` creates the DuckDB and writes the
# population-legend + version metadata, each `append-contig-to-duckdb-index` adds one contig's rows
# in a fresh JVM (so file descriptors do not accumulate across contigs, which exhausted the
# per-process limit when the whole genome was built in one long-lived process), and
# `finalize-duckdb-index` creates the sequence_id index. Contigs run in `CHROMS` order so sequence
# IDs are assigned deterministically; DuckDB is single-writer, so the loop must stay serial.
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
    threads: _HAIL_RULE_THREADS
    resources:
        mem_mb=_HAIL_RULE_MEM_MB,
    params:
        window_size=SEQUENCE_WINDOW_SIZE,
        output_base=f"{WORK_DIR}/output/hgdp_1kg",
        version=VERSION,
        polars_chunk_size=POLARS_CHUNK_SIZE,
        tmp_dir=TMP_DIR,
        spark_driver_memory_gb=SPARK_DRIVER_MEMORY_GB,
        spark_executor_memory_gb=SPARK_EXECUTOR_MEMORY_GB,
        contigs=" ".join(CHROMS),
    shell:
        """
        (
            # Hail's FASTAReader stages a ~3GB copy of the reference genome per JVM, and Spark a
            # blockmgr scratch dir; both default to $TMPDIR and leak there when a per-contig JVM
            # exits non-cleanly (a crash or kill mid-run). Confine them (via $TMPDIR), along with
            # Hail's own tmp_dir and the per-contig TSV (via --tmp-dir below), to a per-run temp dir
            # and delete it on exit (success or failure); also clear any dir a previously killed run
            # left behind (safe: this index build is serial and single-writer).
            rm -rf "{params.tmp_dir}"/divref_index_tmp.* 2>/dev/null || true
            run_tmp=$(mktemp -d "{params.tmp_dir}"/divref_index_tmp.XXXXXX)
            export TMPDIR="$run_tmp"
            trap 'rm -rf "$run_tmp"' EXIT

            divref init-duckdb-index \
                --in-table-pairs-tsv {input.table_pairs_tsv} \
                --output-base {params.output_base} \
                --version {params.version} \
                --window-size {params.window_size} \
                --force

            for contig in {params.contigs}; do
                divref append-contig-to-duckdb-index \
                    --in-table-pairs-tsv {input.table_pairs_tsv} \
                    --contig "$contig" \
                    --output-base {params.output_base} \
                    --reference-fasta {input.fasta} \
                    --window-size {params.window_size} \
                    --version {params.version} \
                    --polars-chunk-size {params.polars_chunk_size} \
                    --tmp-dir "$run_tmp" \
                    --spark-driver-memory-gb {params.spark_driver_memory_gb} \
                    --spark-executor-memory-gb {params.spark_executor_memory_gb}
            done

            divref finalize-duckdb-index --output-base {params.output_base}
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
