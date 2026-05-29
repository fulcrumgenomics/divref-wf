"""Tool to extract gnomAD variant and sample frequency data for the DivRef pipeline."""

import operator
from dataclasses import dataclass
from enum import StrEnum
from enum import unique
from pathlib import Path

import hail as hl
from fgpyo.io import assert_path_is_writable

from divref import defaults
from divref.hail import hail_init
from divref.haplotype import to_hashable_items


@unique
class GnomadVersion(StrEnum):
    """gnomAD release labels for the supported sites tables."""

    JOINT_41 = "JOINT_41"
    GENOMES_312 = "GENOMES_312"
    HGDP_1KG_312 = "HGDP_1KG_312"


@unique
class GnomadCloud(StrEnum):
    """Cloud provider hosting the gnomAD sites table."""

    S3 = "S3"
    GCS = "GCS"


_GNOMAD_TABLE_URI: dict[tuple[GnomadVersion, GnomadCloud], str] = {
    (GnomadVersion.JOINT_41, GnomadCloud.S3): (
        "s3a://gnomad-public-us-east-1/release/4.1/ht/joint/gnomad.joint.v4.1.sites.ht"
    ),
    (GnomadVersion.JOINT_41, GnomadCloud.GCS): (
        "gs://gcp-public-data--gnomad/release/4.1/ht/joint/gnomad.joint.v4.1.sites.ht"
    ),
    (GnomadVersion.GENOMES_312, GnomadCloud.S3): (
        "s3a://gnomad-public-us-east-1/release/3.1.2/ht/genomes/gnomad.genomes.v3.1.2.sites.ht"
    ),
    (GnomadVersion.GENOMES_312, GnomadCloud.GCS): (
        "gs://gcp-public-data--gnomad/release/3.1.2/ht/genomes/gnomad.genomes.v3.1.2.sites.ht"
    ),
    (GnomadVersion.HGDP_1KG_312, GnomadCloud.S3): (
        "s3a://gnomad-public-us-east-1/release/3.1.2/ht/genomes/"
        "gnomad.genomes.v3.1.2.hgdp_1kg_subset_variant_annotations.ht"
    ),
    (GnomadVersion.HGDP_1KG_312, GnomadCloud.GCS): (
        "gs://gcp-public-data--gnomad/release/3.1.2/ht/genomes/"
        "gnomad.genomes.v3.1.2.hgdp_1kg_subset_variant_annotations.ht"
    ),
}


@dataclass(frozen=True)
class _GnomadSchema:
    """Per-version schema config for reading population AFs from a gnomAD sites table."""

    freq_meta_path: str  # dot-path into table globals, e.g. "joint_globals.freq_meta"
    row_freq_field: str  # dot-path on table row, e.g. "joint.freq"
    pop_key: str  # key in freq_meta entries for the population code, e.g. "gen_anc" or "pop"
    popmax_AC_key: str  # dot-path on table row, e.g. "joint.grpmax.AC"  # noqa: N815
    popmax_AF_key: str  # dot-path on table row, e.g. "joint.grpmax.AF"  # noqa: N815
    popmax_AN_key: str  # dot-path on table row, e.g. "joint.grpmax.AN"  # noqa: N815
    popmax_maxpop_key: str  # dot-path on table row, e.g. "joint.grpmax.gen_anc"


_GNOMAD_SCHEMA: dict[GnomadVersion, _GnomadSchema] = {
    GnomadVersion.JOINT_41: _GnomadSchema(
        freq_meta_path="joint_globals.freq_meta",
        row_freq_field="joint.freq",
        pop_key="gen_anc",
        popmax_AC_key="joint.grpmax.AC",
        popmax_AF_key="joint.grpmax.AF",
        popmax_AN_key="joint.grpmax.AN",
        popmax_maxpop_key="joint.grpmax.gen_anc",
    ),
    GnomadVersion.GENOMES_312: _GnomadSchema(
        freq_meta_path="freq_meta",
        row_freq_field="freq",
        pop_key="pop",
        popmax_AC_key="popmax.AC",
        popmax_AF_key="popmax.AF",
        popmax_AN_key="popmax.AN",
        popmax_maxpop_key="popmax.pop",
    ),
    GnomadVersion.HGDP_1KG_312: _GnomadSchema(
        freq_meta_path="gnomad_freq_meta",
        row_freq_field="gnomad_freq",
        pop_key="pop",
        popmax_AC_key="gnomad_popmax.AC",
        popmax_AF_key="gnomad_popmax.AF",
        popmax_AN_key="gnomad_popmax.AN",
        popmax_maxpop_key="gnomad_popmax.pop",
    ),
}


def _apply_filters(va: hl.Table, gnomad_version: GnomadVersion) -> hl.Table:
    """
    Apply gnomAD variant filters, keeping only variants that pass all filters.

    For gnomAD 4.1 joint, exome and genome filter sets are maintained separately; a variant is
    kept only when both pass (the filter set is empty in each). For gnomAD 3.1.2 tables, a
    single filter set is used. Entries with missing filter sets are treated as passing.

    Args:
        va: gnomAD sites Hail table filtered to a single contig.
        gnomad_version: gnomAD version determining which filter fields to apply.

    Returns:
        Filtered Hail table.
    """
    if gnomad_version is GnomadVersion.JOINT_41:
        return va.filter(
            hl.coalesce(hl.len(va.exomes.filters) == 0, True)
            & hl.coalesce(hl.len(va.genomes.filters) == 0, True)
        )
    else:
        return va.filter(hl.coalesce(hl.len(va.filters) == 0, True))


def extract_gnomad_single_afs(
    *,
    gnomad_version: GnomadVersion,
    contig: str,
    freq_threshold: float = 0,
    no_apply_filters: bool = False,
    populations: list[str] = defaults.POPULATIONS,
    reference_genome: str = defaults.REFERENCE_GENOME,
    out_sites_hail_table: Path | None = None,
    out_sites_tsv: Path | None = None,
    gnomad_cloud: GnomadCloud = GnomadCloud.GCS,
    gcs_credentials_path: Path | None = None,
    spark_driver_memory_gb: int = 1,
    spark_executor_memory_gb: int = 1,
) -> None:
    """
    Extract gnomAD variant and sample frequency data for downstream pipeline tools.

    Reads a gnomAD sites table and filters to variants above the frequency threshold in at least one
    population. Writes up to two outputs: a Hail table at `out_sites_hail_table` for downstream
    pipeline tools, and a flat TSV at `out_sites_tsv` with columns `variant` (contig:pos:ref:alt),
    one allele-frequency column per population, `popmax_A[CFN]`, and `maxpop`.

    At least one of `out_sites_hail_table` or `out_sites_tsv` must be defined.

    Args:
        gnomad_version: gnomAD sites table to use - the schema varies per version.
        contig: Contig to extract sites.
        freq_threshold: Minimum allele frequency in any population to retain a variant. Defaults 0,
            all variants returned.
        no_apply_filters: If set, don't apply any of the variant filters (e.g. VQSR, AC0)
        populations: List of population codes to extract frequencies for.
        reference_genome: Reference genome to use. Defaults to "GRCh38".
        out_sites_hail_table: Output path for the Hail table. Optional.
        out_sites_tsv: Output path for the TSV file. Optional.
        gnomad_cloud: Cloud provider hosting the gnomAD sites table. Defaults to GCS
            (`gs://gcp-public-data--gnomad`); set to S3 for the
            `s3a://gnomad-public-us-east-1` mirror.
        gcs_credentials_path: Path to GCS default credentials JSON file. Required
            when `gnomad_cloud` is `GCS`; ignored otherwise.
        spark_driver_memory_gb: Memory in GB to allocate to the Spark driver.
        spark_executor_memory_gb: Memory in GB to allocate to the Spark executor.
    """
    if out_sites_hail_table is None and out_sites_tsv is None:
        raise ValueError("At least one of out_sites_hail_table or out_sites_tsv must be provided")

    # validate output paths before starting Hail
    if out_sites_tsv is not None:
        assert_path_is_writable(out_sites_tsv)

    hail_init(
        gcs_credentials_path.expanduser() if gcs_credentials_path is not None else None,
        spark_driver_memory_gb=spark_driver_memory_gb,
        spark_executor_memory_gb=spark_executor_memory_gb,
        use_s3=(gnomad_cloud is GnomadCloud.S3),
    )

    schema = _GNOMAD_SCHEMA[gnomad_version]
    table_uri = _GNOMAD_TABLE_URI[(gnomad_version, gnomad_cloud)]

    va_all = hl.read_table(table_uri)
    interval = hl.parse_locus_interval(contig, reference_genome=reference_genome)
    va = hl.filter_intervals(va_all, [interval])

    freq_meta = operator.attrgetter(schema.freq_meta_path)(va.globals).collect()[0]
    map_to_index = {to_hashable_items(x): i for i, x in enumerate(freq_meta)}

    pop_indices = []
    for pop in populations:
        idx = map_to_index.get(to_hashable_items({"group": "adj", schema.pop_key: pop}))
        if idx is None:
            raise ValueError(f"Population {pop!r} not found in gnomAD frequency metadata")
        pop_indices.append(idx)

    if not no_apply_filters:
        va = _apply_filters(va, gnomad_version)

    va = va.select_globals(pops=populations)
    row_freq = operator.attrgetter(schema.row_freq_field)(va)
    va = va.select(
        pop_freqs=hl.literal(pop_indices).map(lambda i: row_freq[i]),
        popmax_AC=operator.attrgetter(schema.popmax_AC_key)(va),
        popmax_AF=operator.attrgetter(schema.popmax_AF_key)(va),
        popmax_AN=operator.attrgetter(schema.popmax_AN_key)(va),
        maxpop=operator.attrgetter(schema.popmax_maxpop_key)(va),
    )
    if freq_threshold > 0:
        va = va.filter(hl.any(lambda x: x.AF >= freq_threshold, va.pop_freqs))
    va = va.key_by()
    va = va.select(
        "locus",
        "alleles",
        "pop_freqs",
        "popmax_AC",
        "popmax_AF",
        "popmax_AN",
        "maxpop",
    )

    if out_sites_hail_table is not None:
        va.naive_coalesce(64).write(str(out_sites_hail_table), overwrite=True)

    if out_sites_tsv is not None:
        va.select(
            variant=hl.format(
                "%s:%s:%s:%s",
                va.locus.contig,
                hl.str(va.locus.position),
                va.alleles[0],
                va.alleles[1],
            ),
            **{pop: va.pop_freqs[i].AF for i, pop in enumerate(populations)},
            popmax_AC=va.popmax_AC,
            popmax_AF=va.popmax_AF,
            popmax_AN=va.popmax_AN,
            maxpop=va.maxpop,
        ).export(str(out_sites_tsv))
