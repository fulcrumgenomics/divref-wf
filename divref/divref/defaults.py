from typing import Final

from divref.alias import HailPath

GNOMAD_HGDP_1KG_VARIANT_ANNOTATION_HAIL_TABLE: Final[HailPath] = (
    "gs://gcp-public-data--gnomad/release/3.1.2/ht/genomes/gnomad.genomes.v3.1.2.hgdp_1kg_subset_variant_annotations.ht"  # noqa: E501
)
"""HGDP+1KG individual level genotypes (GCS)."""

GNOMAD_HGDP_1KG_VARIANT_ANNOTATION_HAIL_TABLE_S3: Final[HailPath] = (
    "s3a://gnomad-public-us-east-1/release/3.1.2/ht/genomes/gnomad.genomes.v3.1.2.hgdp_1kg_subset_variant_annotations.ht"  # noqa: E501
)
"""HGDP+1KG individual level genotypes (AWS S3 alternative)."""

GNOMAD_HGDP_1KG_SAMPLE_METADATA_HAIL_TABLE: Final[HailPath] = (
    "gs://gcp-public-data--gnomad/release/3.1.2/ht/genomes/gnomad.genomes.v3.1.2.hgdp_1kg_subset_sample_meta.ht"  # noqa: E501
)
"""HGDP+1KG sample metadata (GCS)."""

GNOMAD_HGDP_1KG_SAMPLE_METADATA_HAIL_TABLE_S3: Final[HailPath] = (
    "s3a://gnomad-public-us-east-1/release/3.1.2/ht/genomes/gnomad.genomes.v3.1.2.hgdp_1kg_subset_sample_meta.ht"  # noqa: E501
)
"""HGDP+1KG sample metadata (AWS S3 alternative)."""

POPULATIONS: list[str] = ["afr", "amr", "eas", "sas", "nfe"]
"""Default HGDP+1KG populations."""

REFERENCE_GENOME: str = "GRCh38"
"""Default reference genome assembly."""
