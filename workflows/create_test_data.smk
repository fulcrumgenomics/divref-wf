####################################################################################################
# Creates the test data in divref/tests/data.
####################################################################################################

from pathlib import Path

####################################################################################################
# Inputs
####################################################################################################

OUTPUT_DIR: Path = Path("divref/tests/data")
LOCUS_CHROM: str = "chr1"
LOCUS: str = "chr1:100001-200000"
LOCUS_FILENAME: str = "chr1_100001_200000"
# chrX non-PAR locus used to exercise the haploid-male ploidy correction in
# `compute_haplotypes`. The window is well inside non-PAR
# (PAR1 ends at 2,781,479; PAR2 starts at 155,701,383 on GRCh38).
CHRX_LOCUS_CHROM: str = "chrX"
CHRX_LOCUS: str = "chrX:50000000-50025000"
CHRX_LOCUS_FILENAME: str = "chrX_50000000_50025000"
CHRX_BCF_NONPAR: str = (
    "gs://gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2/"
    "hgdp1kgp_chrX_non_par.full.shapeit5_rare.bcf"
)
MIN_POP_AF_EXTRACT_GNOMAD_AFS: float = 0.001
MIN_POP_AF_COMPUTE_HAPLOTYPES: float = 0.005
MIN_POPMAX_AF_CREATE_GNOMAD_SITES_VCF: float = 0.01
WINDOW_SIZE_COMPUTE_HAPLOTYPES: int = 25
# Hail-using divref tools require a GCS credentials path when reading from local-only Hail
# tables, because hail_init currently sets `use_s3=False` by default and asserts the path is
# present. Threaded through every rule below for consistency.
GCS_CREDENTIALS_PATH: str = "~/.config/gcloud/application_default_credentials.json"

####################################################################################################
# Rules
####################################################################################################


rule all:
    input:
        f"{OUTPUT_DIR}/{LOCUS_FILENAME}.ht",
        f"{OUTPUT_DIR}/hgdp_1kg_sample_metadata.ht",
        f"{OUTPUT_DIR}/samples.txt",
        f"{OUTPUT_DIR}/{LOCUS_FILENAME}.vcf.gz",
        f"{OUTPUT_DIR}/{LOCUS_FILENAME}.vcf.gz.tbi",
        f"{OUTPUT_DIR}/{LOCUS_FILENAME}.gnomad_afs.ht",
        f"{OUTPUT_DIR}/hgdp_1kg_sample_metadata.extract.ht",
        f"{OUTPUT_DIR}/{LOCUS_FILENAME}_haplotypes.ht",
        f"{OUTPUT_DIR}/{LOCUS_FILENAME}.gnomad_sites.vcf.bgz",
        f"{OUTPUT_DIR}/{CHRX_LOCUS_FILENAME}.ht",
        f"{OUTPUT_DIR}/{CHRX_LOCUS_FILENAME}.vcf.gz",
        f"{OUTPUT_DIR}/{CHRX_LOCUS_FILENAME}.vcf.gz.tbi",
        f"{OUTPUT_DIR}/{CHRX_LOCUS_FILENAME}.gnomad_afs.ht",


####################################################################################################
# Extracts all gnomAD HGDP+1KG variants in the specified locus from
# gs://gcp-public-data--gnomad/release/3.1.2/ht/genomes/gnomad.genomes.v3.1.2.hgdp_1kg_subset_variant_annotations.ht.
#
# Extracts selected fields from gnomAD sample metadata required by downstream tools from
# gs://gcp-public-data--gnomad/release/3.1.2/ht/genomes/gnomad.genomes.v3.1.2.hgdp_1kg_subset_sample_meta.ht.
####################################################################################################
rule subset_gnomad_hail_tables:
    output:
        chr1_variant_ht=directory(f"{OUTPUT_DIR}/{LOCUS_FILENAME}.ht"),
        chrx_variant_ht=directory(f"{OUTPUT_DIR}/{CHRX_LOCUS_FILENAME}.ht"),
        sample_ht=directory(f"{OUTPUT_DIR}/hgdp_1kg_sample_metadata.ht"),
        samples_txt=f"{OUTPUT_DIR}/samples.txt",
    log:
        f"logs/create_test_data/subset_gnomad_hail_tables.log",
    params:
        chr1_locus=LOCUS,
        chrx_locus=CHRX_LOCUS,
    shell:
        """
        (
            divref gnomad-hail-table-test-data \
                --loci {params.chr1_locus} {params.chrx_locus} \
                --out-variant-annotation-dir {OUTPUT_DIR} \
                --out-sample-metadata {output.sample_ht} \
                --out-samples-txt {output.samples_txt}
        ) &> {log}
        """


####################################################################################################
# Extracts the phased genotypes for all HGDP+1KG samples in the specified locus.
####################################################################################################
rule subset_phased_genotypes:
    input:
        samples_txt=f"{OUTPUT_DIR}/samples.txt",
    output:
        vcf=f"{OUTPUT_DIR}/{LOCUS_FILENAME}.vcf.gz",
        tbi=f"{OUTPUT_DIR}/{LOCUS_FILENAME}.vcf.gz.tbi",
    log:
        f"logs/create_test_data/subset_phased_genotypes.{LOCUS_FILENAME}.log",
    params:
        locus=LOCUS,
        bcf=f"gs://gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2/hgdp1kgp_{LOCUS_CHROM}.filtered.SNV_INDEL.phased.shapeit5.bcf",
    shell:
        """
        (
            bcftools view \
                --regions {params.locus} \
                --samples-file {input.samples_txt} \
                --force-samples \
                --output-type z \
                --output {output.vcf} \
                --write-index=tbi \
                {params.bcf}
        ) &> {log}
        """


####################################################################################################
# Extracts allele frequencies for the default populations and subsets to sites over the specified
# minimum frequency.
####################################################################################################
rule extract_gnomad_afs:
    input:
        variant_ht=f"{OUTPUT_DIR}/{LOCUS_FILENAME}.ht",
    output:
        variant_ht=directory(f"{OUTPUT_DIR}/{LOCUS_FILENAME}.gnomad_afs.ht"),
    log:
        f"logs/create_test_data/extract_gnomad_afs.{LOCUS_FILENAME}.log",
    params:
        contig=LOCUS_CHROM,
        freq_threshold=MIN_POP_AF_EXTRACT_GNOMAD_AFS,
        gcs_credentials_path=GCS_CREDENTIALS_PATH,
    shell:
        """
        (
            divref extract-gnomad-afs \
                --in-gnomad-sites-table {input.variant_ht} \
                --out-variant-annotation-table {output.variant_ht} \
                --contig {params.contig} \
                --freq-threshold {params.freq_threshold} \
                --gcs-credentials-path {params.gcs_credentials_path}
        ) &> {log}
        """


####################################################################################################
# Extracts selected fields from sample metadata.
####################################################################################################
rule extract_sample_metadata:
    input:
        sample_ht=f"{OUTPUT_DIR}/hgdp_1kg_sample_metadata.ht",
    output:
        sample_ht=directory(f"{OUTPUT_DIR}/hgdp_1kg_sample_metadata.extract.ht"),
    log:
        f"logs/create_test_data/extract_sample_metadata.{LOCUS_FILENAME}.log",
    params:
        gcs_credentials_path=GCS_CREDENTIALS_PATH,
    shell:
        """
        (
            divref extract-sample-metadata \
                --in-gnomad-hgdp-sample-data {input.sample_ht} \
                --out-sample-metadata {output.sample_ht} \
                --gcs-credentials-path {params.gcs_credentials_path}
        ) &> {log}
        """


####################################################################################################
# Compute haplotypes from the sites and phased genotypes.
####################################################################################################
rule compute_haplotypes:
    input:
        vcf=f"{OUTPUT_DIR}/{LOCUS_FILENAME}.vcf.gz",
        tbi=f"{OUTPUT_DIR}/{LOCUS_FILENAME}.vcf.gz.tbi",
        variant_ht=f"{OUTPUT_DIR}/{LOCUS_FILENAME}.gnomad_afs.ht",
        sample_ht=f"{OUTPUT_DIR}/hgdp_1kg_sample_metadata.extract.ht",
    output:
        haplotypes_ht=directory(f"{OUTPUT_DIR}/{LOCUS_FILENAME}_haplotypes.ht"),
    log:
        f"logs/create_test_data/compute_haplotypes.{LOCUS_FILENAME}.log",
    params:
        window_size=WINDOW_SIZE_COMPUTE_HAPLOTYPES,
        freq_threshold=MIN_POP_AF_COMPUTE_HAPLOTYPES,
        output_base=f"{OUTPUT_DIR}/{LOCUS_FILENAME}_haplotypes",
    shell:
        """
        (
            divref compute-haplotypes \
                --vcfs-path {input.vcf} \
                --gnomad-va-file {input.variant_ht} \
                --gnomad-sa-file {input.sample_ht} \
                --window-size {params.window_size} \
                --variant-freq-threshold {params.freq_threshold} \
                --haplotype-freq-threshold {params.freq_threshold} \
                --output-base {params.output_base}
        ) &> {log}
        """


####################################################################################################
# Create a sites VCF from the gnomAD Hail table.
####################################################################################################
rule create_gnomad_sites_vcf:
    input:
        variant_ht=f"{OUTPUT_DIR}/{LOCUS_FILENAME}.gnomad_afs.ht",
    output:
        vcf=f"{OUTPUT_DIR}/{LOCUS_FILENAME}.gnomad_sites.vcf.bgz",
    log:
        f"logs/create_test_data/create_gnomad_sites_vcf.{LOCUS_FILENAME}.log",
    params:
        min_popmax=MIN_POPMAX_AF_CREATE_GNOMAD_SITES_VCF,
    shell:
        """
        (
            divref create-gnomad-sites-vcf \
                --sites-table-path {input.variant_ht} \
                --output-vcf-path {output.vcf} \
                --min-popmax {params.min_popmax}
        ) &> {log}
        """


####################################################################################################
# Extracts phased genotypes for the chrX non-PAR test locus directly from the non-PAR BCF.
####################################################################################################
rule subset_phased_genotypes_chrX:
    input:
        samples_txt=f"{OUTPUT_DIR}/samples.txt",
    output:
        vcf=f"{OUTPUT_DIR}/{CHRX_LOCUS_FILENAME}.vcf.gz",
        tbi=f"{OUTPUT_DIR}/{CHRX_LOCUS_FILENAME}.vcf.gz.tbi",
    log:
        f"logs/create_test_data/subset_phased_genotypes.{CHRX_LOCUS_FILENAME}.log",
    params:
        locus=CHRX_LOCUS,
        bcf=CHRX_BCF_NONPAR,
    shell:
        """
        (
            bcftools view \
                --regions {params.locus} \
                --samples-file {input.samples_txt} \
                --force-samples \
                --output-type z \
                --output {output.vcf} \
                --write-index=tbi \
                {params.bcf}
        ) &> {log}
        """


####################################################################################################
# Extract gnomAD allele frequencies for the chrX non-PAR test locus.
####################################################################################################
rule extract_gnomad_afs_chrX:
    input:
        variant_ht=f"{OUTPUT_DIR}/{CHRX_LOCUS_FILENAME}.ht",
    output:
        variant_ht=directory(f"{OUTPUT_DIR}/{CHRX_LOCUS_FILENAME}.gnomad_afs.ht"),
    log:
        f"logs/create_test_data/extract_gnomad_afs.{CHRX_LOCUS_FILENAME}.log",
    params:
        contig=CHRX_LOCUS_CHROM,
        freq_threshold=MIN_POP_AF_EXTRACT_GNOMAD_AFS,
        gcs_credentials_path=GCS_CREDENTIALS_PATH,
    shell:
        """
        (
            divref extract-gnomad-afs \
                --in-gnomad-sites-table {input.variant_ht} \
                --out-variant-annotation-table {output.variant_ht} \
                --contig {params.contig} \
                --freq-threshold {params.freq_threshold} \
                --gcs-credentials-path {params.gcs_credentials_path}
        ) &> {log}
        """
