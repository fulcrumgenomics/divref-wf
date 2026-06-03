####################################################################################################
# Compares DivRef 1.1 single-variant sites against allele frequencies from different gnomAD
# releases for chr22.
####################################################################################################

from pathlib import Path

####################################################################################################
# Inputs / constants
####################################################################################################

OUTPUT_DIR: Path = Path("data/analysis")
COMPARISON_NAME: str = "compare_divref_gnomad"
CONTIG: str = "chr22"
FREQUENCY_THRESHOLD: float = 0.005
DIVREF_DUCKDB_URL: str = (
    "https://zenodo.org/records/14802613/files/" "DivRef-v1.1.haplotypes_gnomad_merge.index.duckdb"
)
# Zenodo-published md5 for the DuckDB above (record 14802613). Verified after download so a
# truncated or corrupted fetch fails the rule instead of silently feeding a bad index downstream.
DIVREF_DUCKDB_MD5: str = "6066d92f2d0269e4620602f4ded60b2b"
GNOMAD_VERSIONS: list[str] = ["joint_41", "genomes_312", "hgdp_1kg_312"]

# Maps filename wildcard → plot label used by compare_divref_gnomad.R
GNOMAD_LABEL: dict[str, str] = {
    "joint_41": "gnomAD 4.1 joint",
    "genomes_312": "gnomAD 3.1.2 genomes",
    "hgdp_1kg_312": "gnomAD 3.1.2 HGDP+1KG",
}

OUT_FILE_EXTS: list[str] = [
    ".af_diffs.png",
    ".af_diffs_all.png",
    ".venn.png",
    ".not_in_gnomad_afs.png",
    ".divref_not_in_gnomad.tsv",
    ".log",
]

####################################################################################################
# Rules
####################################################################################################


# Constrain gnomad_version to the known tokens so extract_gnomad_single_afs's
# `{contig}.{gnomad_version}.tsv` cannot greedily match compare_divref_gnomad's
# `{contig}.{gnomad_version}.divref_not_in_gnomad.tsv` (which previously forced a ruleorder).
wildcard_constraints:
    gnomad_version="|".join(GNOMAD_VERSIONS),


rule all:
    input:
        expand(
            f"{OUTPUT_DIR}/{COMPARISON_NAME}/{CONTIG}.{{gnomad_version}}{{ext}}",
            gnomad_version=GNOMAD_VERSIONS,
            ext=OUT_FILE_EXTS,
        ),


####################################################################################################
# Downloads the DivRef 1.1 DuckDB index from Zenodo.
####################################################################################################
rule download_divref_index:
    output:
        duckdb=f"{OUTPUT_DIR}/input/DivRef-v1.1.haplotypes_gnomad_merge.index.duckdb",
    log:
        f"logs/{COMPARISON_NAME}/download_divref_index.log",
    params:
        url=DIVREF_DUCKDB_URL,
        expected_md5=DIVREF_DUCKDB_MD5,
    shell:
        """
        (
            wget --no-verbose -O {output.duckdb} {params.url}
            # Verify against the Zenodo-published md5 (md5sum on Linux, md5 on macOS).
            if command -v md5sum >/dev/null 2>&1; then
                actual=$(md5sum "{output.duckdb}" | awk '{{print $1}}')
            else
                actual=$(md5 -q "{output.duckdb}")
            fi
            if [ "$actual" != "{params.expected_md5}" ]; then
                echo "Checksum mismatch for {output.duckdb}: expected {params.expected_md5}, got $actual" >&2
                exit 1
            fi
        ) &> {log}
        """


####################################################################################################
# Extracts chr22 allele frequencies from a gnomAD sites table and writes a Hail table and a flat
# TSV with one AF column per population.
####################################################################################################
rule extract_gnomad_single_afs:
    output:
        tsv=f"{OUTPUT_DIR}/{COMPARISON_NAME}/{CONTIG}.{{gnomad_version}}.tsv",
    log:
        f"logs/{COMPARISON_NAME}/extract_gnomad_single_afs.{{gnomad_version}}.log",
    params:
        contig=CONTIG,
        freq_threshold=FREQUENCY_THRESHOLD,
        gnomad_version=lambda wildcards: wildcards.gnomad_version.upper(),
        gcs_credentials_path="~/.config/gcloud/application_default_credentials.json",
        spark_driver_memory_gb=16,
        spark_executor_memory_gb=16,
    shell:
        """
        (
            divref extract-gnomad-single-afs \
                --contig {params.contig} \
                --freq-threshold {params.freq_threshold} \
                --gnomad-version {params.gnomad_version} \
                --out-sites-tsv {output.tsv} \
                --gcs-credentials-path '{params.gcs_credentials_path}' \
                --spark-driver-memory-gb {params.spark_driver_memory_gb} \
                --spark-executor-memory-gb {params.spark_executor_memory_gb}
        ) &> {log}
        """


####################################################################################################
# Compares DivRef 1.1 gnomAD_variant sites to the extracted gnomAD allele frequencies.
# If all DivRef variants are found in gnomAD the R script exits early; touch ensures Snakemake's
# output checks are satisfied even in that case.
####################################################################################################
rule compare_divref_gnomad:
    input:
        duckdb=f"{OUTPUT_DIR}/input/DivRef-v1.1.haplotypes_gnomad_merge.index.duckdb",
        tsv=f"{OUTPUT_DIR}/{COMPARISON_NAME}/{CONTIG}.{{gnomad_version}}.tsv",
    output:
        af_diffs_png=f"{OUTPUT_DIR}/{COMPARISON_NAME}/{CONTIG}.{{gnomad_version}}.af_diffs.png",
        af_diffs_all_png=f"{OUTPUT_DIR}/{COMPARISON_NAME}/{CONTIG}.{{gnomad_version}}.af_diffs_all.png",
        venn_png=f"{OUTPUT_DIR}/{COMPARISON_NAME}/{CONTIG}.{{gnomad_version}}.venn.png",
        not_in_gnomad_png=f"{OUTPUT_DIR}/{COMPARISON_NAME}/{CONTIG}.{{gnomad_version}}.not_in_gnomad_afs.png",
        not_in_gnomad_tsv=f"{OUTPUT_DIR}/{COMPARISON_NAME}/{CONTIG}.{{gnomad_version}}.divref_not_in_gnomad.tsv",
    log:
        f"{OUTPUT_DIR}/{COMPARISON_NAME}/{CONTIG}.{{gnomad_version}}.log",
    params:
        contig=CONTIG,
        gnomad_label=lambda wildcards: GNOMAD_LABEL[wildcards.gnomad_version],
        output_base=f"{OUTPUT_DIR}/{COMPARISON_NAME}/{CONTIG}.{{gnomad_version}}",
    shell:
        """
        (
            Rscript scripts/compare_divref_gnomad.R \
                --contig {params.contig} \
                --divref_duckdb {input.duckdb} \
                --gnomad_tsv {input.tsv} \
                --gnomad_label '{params.gnomad_label}' \
                --output_base {params.output_base}

            # The R script exits early when no DivRef variants are absent from gnomAD,
            # so touch ensures these outputs exist for Snakemake's file checks.
            touch {output.not_in_gnomad_png} {output.not_in_gnomad_tsv}
        ) &> {log}
        """
