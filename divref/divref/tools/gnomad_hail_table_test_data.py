"""Tool to fetch a small subset of a gnomAD sites hail table and HGDP sample information."""

import logging
from pathlib import Path

import hail as hl
from fgpyo.io import assert_path_is_writable

from divref import defaults
from divref.alias import HailPath
from divref.hail import hail_init

logger = logging.getLogger(__name__)


def gnomad_hail_table_test_data(
    *,
    in_gnomad_hgdp_variant_annotation_table: HailPath = defaults.GNOMAD_HGDP_1KG_VARIANT_ANNOTATION_HAIL_TABLE,  # noqa: E501
    in_gnomad_hgdp_sample_metadata: HailPath = defaults.GNOMAD_HGDP_1KG_SAMPLE_METADATA_HAIL_TABLE,  # noqa: E501
    out_variant_annotation_table: Path,
    out_sample_metadata: Path,
    locus: str = "chr1:100001-200000",
    gcs_credentials_path: Path = Path("~/.config/gcloud/application_default_credentials.json"),
    spark_driver_memory_gb: int = 1,
    spark_executor_memory_gb: int = 1,
) -> None:
    """
    Extract subsets of gnomAD HGDP/1KG variant annotations and sample metadata for testing.

    Args:
        in_gnomad_hgdp_variant_annotation_table: Path to the gnomAD HGDP/1KG variant annotation
            Hail table.
        in_gnomad_hgdp_sample_metadata: Path to the gnomAD HGDP/1KG sample metadata Hail table.
        out_variant_annotation_table: Output path for the subset variant annotation Hail table.
        out_sample_metadata: Output path for the sample metadata Hail table, stripped to key,
            `gnomad_population_inference`, and `gnomad_sex_imputation`.
        locus: Locus interval for variant filtering.
        gcs_credentials_path: Path to GCS default credentials JSON file.
        spark_driver_memory_gb: Memory in GB to allocate to the Spark driver.
        spark_executor_memory_gb: Memory in GB to allocate to the Spark executor.
    """
    assert_path_is_writable(out_variant_annotation_table)
    assert_path_is_writable(out_sample_metadata)

    hail_init(
        gcs_credentials_path.expanduser(),
        spark_driver_memory_gb=spark_driver_memory_gb,
        spark_executor_memory_gb=spark_executor_memory_gb,
    )

    va = hl.read_table(in_gnomad_hgdp_variant_annotation_table)
    roi = [hl.parse_locus_interval(locus, reference_genome=defaults.REFERENCE_GENOME)]

    logger.info(f"Filtering to {locus}.")
    va_subset = hl.filter_intervals(va, roi)

    logger.info(f"Writing {va_subset.count()} variants to {out_variant_annotation_table}.")
    va_subset.write(str(out_variant_annotation_table), overwrite=True)

    sa = hl.read_table(in_gnomad_hgdp_sample_metadata)
    sa = sa.select("gnomad_population_inference", "gnomad_sex_imputation").select_globals()

    logger.info(f"Writing {sa.count()} samples to {out_sample_metadata}.")
    sa.naive_coalesce(1).write(str(out_sample_metadata), overwrite=True)
