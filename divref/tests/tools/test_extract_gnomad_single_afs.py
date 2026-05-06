"""Tests for the extract_gnomad_single_afs tool."""

from inspect import signature

import pytest

from divref.tools.extract_gnomad_single_afs import _GNOMAD_TABLE_URI
from divref.tools.extract_gnomad_single_afs import GnomadCloud
from divref.tools.extract_gnomad_single_afs import GnomadVersion
from divref.tools.extract_gnomad_single_afs import extract_gnomad_single_afs


@pytest.mark.parametrize(
    "version,cloud,expected_prefix,expected_suffix",
    [
        (
            GnomadVersion.JOINT_41,
            GnomadCloud.S3,
            "s3a://gnomad-public-us-east-1/",
            "release/4.1/ht/joint/gnomad.joint.v4.1.sites.ht",
        ),
        (
            GnomadVersion.JOINT_41,
            GnomadCloud.GCS,
            "gs://gcp-public-data--gnomad/",
            "release/4.1/ht/joint/gnomad.joint.v4.1.sites.ht",
        ),
        (
            GnomadVersion.GENOMES_312,
            GnomadCloud.S3,
            "s3a://gnomad-public-us-east-1/",
            "release/3.1.2/ht/genomes/gnomad.genomes.v3.1.2.sites.ht",
        ),
        (
            GnomadVersion.GENOMES_312,
            GnomadCloud.GCS,
            "gs://gcp-public-data--gnomad/",
            "release/3.1.2/ht/genomes/gnomad.genomes.v3.1.2.sites.ht",
        ),
        (
            GnomadVersion.HGDP_1KG_312,
            GnomadCloud.S3,
            "s3a://gnomad-public-us-east-1/",
            "gnomad.genomes.v3.1.2.hgdp_1kg_subset_variant_annotations.ht",
        ),
        (
            GnomadVersion.HGDP_1KG_312,
            GnomadCloud.GCS,
            "gs://gcp-public-data--gnomad/",
            "gnomad.genomes.v3.1.2.hgdp_1kg_subset_variant_annotations.ht",
        ),
    ],
)
def test_gnomad_table_uri_lookup(
    version: GnomadVersion,
    cloud: GnomadCloud,
    expected_prefix: str,
    expected_suffix: str,
) -> None:
    """Each (version, cloud) maps to a URI on the matching cloud, ending in the right table."""
    uri = _GNOMAD_TABLE_URI[(version, cloud)]
    assert uri.startswith(expected_prefix)
    assert uri.endswith(expected_suffix)


def test_gnomad_table_uri_table_path_matches_across_clouds() -> None:
    """For each version, the S3 and GCS URIs reference the same table path under the bucket."""
    for version in GnomadVersion:
        s3_uri = _GNOMAD_TABLE_URI[(version, GnomadCloud.S3)]
        gcs_uri = _GNOMAD_TABLE_URI[(version, GnomadCloud.GCS)]
        s3_tail = s3_uri.removeprefix("s3a://gnomad-public-us-east-1/")
        gcs_tail = gcs_uri.removeprefix("gs://gcp-public-data--gnomad/")
        assert s3_tail == gcs_tail


def test_gnomad_cloud_default_is_gcs() -> None:
    """The CLI / function default for `gnomad_cloud` is `GnomadCloud.GCS`."""
    sig = signature(extract_gnomad_single_afs)
    assert sig.parameters["gnomad_cloud"].default is GnomadCloud.GCS
