"""Tests for the extract_gnomad_single_afs tool."""

from inspect import signature
from pathlib import Path
from typing import Any

import hail as hl
import pytest

from divref.tools.extract_gnomad_single_afs import _GNOMAD_TABLE_URI
from divref.tools.extract_gnomad_single_afs import GnomadCloud
from divref.tools.extract_gnomad_single_afs import GnomadVersion
from divref.tools.extract_gnomad_single_afs import _apply_filters
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


class _StopReadTableError(Exception):
    """Sentinel raised by the stubbed ``hl.read_table`` to short-circuit the tool."""


@pytest.mark.parametrize(
    "cloud,expected_use_s3,expected_prefix",
    [
        (GnomadCloud.S3, True, "s3a://"),
        (GnomadCloud.GCS, False, "gs://"),
    ],
)
def test_extract_gnomad_single_afs_dispatches_to_correct_cloud(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    cloud: GnomadCloud,
    expected_use_s3: bool,
    expected_prefix: str,
) -> None:
    """`hail_init` receives the matching `use_s3` and `hl.read_table` gets the resolved URI."""
    captured: dict[str, Any] = {}

    def fake_hail_init(*_args: Any, **kwargs: Any) -> None:
        captured["use_s3"] = kwargs.get("use_s3")

    def fake_read_table(uri: str) -> Any:
        captured["uri"] = uri
        raise _StopReadTableError

    monkeypatch.setattr("divref.tools.extract_gnomad_single_afs.hail_init", fake_hail_init)
    monkeypatch.setattr("divref.tools.extract_gnomad_single_afs.hl.read_table", fake_read_table)

    with pytest.raises(_StopReadTableError):
        extract_gnomad_single_afs(
            gnomad_version=GnomadVersion.JOINT_41,
            contig="chr22",
            gnomad_cloud=cloud,
            out_sites_hail_table=tmp_path / "out.ht",
        )

    assert captured["use_s3"] is expected_use_s3
    assert captured["uri"].startswith(expected_prefix)
    assert captured["uri"].endswith("release/4.1/ht/joint/gnomad.joint.v4.1.sites.ht")


def _make_joint41_filter_ht(
    rows: list[tuple[int, list[str] | None, list[str] | None]],
) -> hl.Table:
    """
    Build a synthetic JOINT_41-shaped Hail table for `_apply_filters` testing.

    Each input row is `(position, exomes_filters, genomes_filters)`. A `None` filter list
    encodes a missing filter set (the `hl.coalesce(..., True)` branch in `_apply_filters`,
    which should be treated as passing). A list (possibly empty) encodes a present filter set.
    """
    schema = hl.tstruct(
        locus=hl.tstruct(contig=hl.tstr, position=hl.tint32),
        alleles=hl.tarray(hl.tstr),
        exomes=hl.tstruct(filters=hl.tarray(hl.tstr)),
        genomes=hl.tstruct(filters=hl.tarray(hl.tstr)),
    )
    table_rows = []
    for pos, exome_filters, genome_filters in rows:
        table_rows.append({
            "locus": {"contig": "chr22", "position": pos},
            "alleles": ["A", "T"],
            "exomes": {"filters": exome_filters},
            "genomes": {"filters": genome_filters},
        })
    return hl.Table.parallelize(table_rows, schema=schema)


def test_apply_filters_joint41_keeps_only_variants_passing_both(
    hail_context: None,  # noqa: ARG001
) -> None:
    """JOINT_41 filter requires both exome and genome filter sets to be empty (or missing)."""
    ht = _make_joint41_filter_ht([
        (100, [], []),  # both empty → kept
        (200, [], ["AC0"]),  # exome empty, genome non-empty → dropped
        (300, ["AC0"], []),  # exome non-empty, genome empty → dropped
        (400, ["AC0"], ["AC0"]),  # both non-empty → dropped
        (500, None, []),  # exome missing, genome empty → kept (missing coalesces to pass)
        (600, [], None),  # exome empty, genome missing → kept
        (700, None, None),  # both missing → kept
        (800, ["AC0"], None),  # exome non-empty, genome missing → dropped (exome fails)
    ])

    survivors = sorted(
        r.locus.position for r in _apply_filters(ht, GnomadVersion.JOINT_41).collect()
    )
    assert survivors == [100, 500, 600, 700]


def test_extract_gnomad_single_afs_propagates_hail_init_failure(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """If `hail_init` raises, the exception propagates out of `extract_gnomad_single_afs`."""

    def boom(*_args: Any, **_kwargs: Any) -> None:
        raise RuntimeError("hail_init failed")

    monkeypatch.setattr("divref.tools.extract_gnomad_single_afs.hail_init", boom)

    with pytest.raises(RuntimeError, match="hail_init failed"):
        extract_gnomad_single_afs(
            gnomad_version=GnomadVersion.JOINT_41,
            contig="chr22",
            gnomad_cloud=GnomadCloud.GCS,
            out_sites_hail_table=tmp_path / "out.ht",
        )
