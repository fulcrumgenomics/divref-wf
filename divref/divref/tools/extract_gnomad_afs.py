"""Tool to extract gnomAD variant and sample frequency data for the DivRef pipeline."""

from pathlib import Path

import hail as hl
from fgpyo.io import assert_path_is_writable

from divref import defaults
from divref.alias import HailPath
from divref.hail import hail_init
from divref.haplotype import to_hashable_items


def extract_gnomad_afs(
    *,
    in_gnomad_sites_table: HailPath,
    out_variant_annotation_table: Path,
    contig: str,
    freq_threshold: float = 0.001,
    populations: list[str] = defaults.POPULATIONS,
    reference_genome: str = defaults.REFERENCE_GENOME,
    gcs_credentials_path: Path | None = None,
    spark_driver_memory_gb: int = 1,
    spark_executor_memory_gb: int = 1,
    use_s3: bool = False,
) -> None:
    """
    Extract gnomAD variant and sample frequency data for downstream pipeline tools.

    Reads the gnomAD v3.1.2 HGDP/1KG subset, extracts per-population allele frequencies
    for the specified populations, filters to variants above the frequency threshold in at
    least one population, and writes compact variant annotation and sample metadata tables.

    Args:
        in_gnomad_sites_table: Path to the gnomAD HGDP/1KG sites table.
        out_variant_annotation_table: Output path for the variant annotation Hail table.
        contig: Contig to extract sites.
        freq_threshold: Minimum allele frequency in any population to retain a variant.
        populations: List of population codes to extract frequencies for.
        reference_genome: Reference genome to use. Defaults to "GRCh38".
        gcs_credentials_path: Path to GCS default credentials JSON file. Required
            when `use_s3` is `False`; ignored otherwise.
        spark_driver_memory_gb: Memory in GB to allocate to the Spark driver.
        spark_executor_memory_gb: Memory in GB to allocate to the Spark executor.
        use_s3: If `True`, initialize Hail with the S3A connector instead of the GCS
            connector. Set this when `in_gnomad_sites_table` is an `s3a://` URI.
    """
    assert_path_is_writable(out_variant_annotation_table)

    hail_init(
        gcs_credentials_path.expanduser() if gcs_credentials_path is not None else None,
        spark_driver_memory_gb=spark_driver_memory_gb,
        spark_executor_memory_gb=spark_executor_memory_gb,
        use_s3=use_s3,
    )

    va_all = hl.read_table(in_gnomad_sites_table)
    interval = hl.parse_locus_interval(contig, reference_genome=reference_genome)
    va = hl.filter_intervals(va_all, [interval])

    freq_meta = va.globals.gnomad_freq_meta.collect()[0]
    map_to_index = {to_hashable_items(x): i for i, x in enumerate(freq_meta)}

    pop_indices = []
    for pop in populations:
        idx = map_to_index.get(to_hashable_items({"group": "adj", "pop": pop}))
        if idx is None:
            raise ValueError(f"Population {pop!r} not found in gnomAD frequency metadata")
        pop_indices.append(idx)

    # Some filter sets are {} and some are NA; treat NA as passing.
    va = va.filter(hl.coalesce(hl.len(va.filters) == 0, True))

    va = va.select_globals(pops=populations)
    va = va.select(pop_freqs=hl.literal(pop_indices).map(lambda i: va.gnomad_freq[i]))
    va = va.filter(hl.any(lambda x: x.AF >= freq_threshold, va.pop_freqs))
    va.naive_coalesce(64).write(str(out_variant_annotation_table), overwrite=True)
