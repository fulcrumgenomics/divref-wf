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

    Each input row is `(position, exomes_filters, genomes_filters)`. A `None` filter set encodes
    a missing field (the `hl.coalesce(..., True)` branch in `_apply_filters`, which should be
    treated as passing). A list (possibly empty) is loaded into a `tset[tstr]` field to mirror
    the real gnomAD joint HT schema.
    """
    schema = hl.tstruct(
        locus=hl.tstruct(contig=hl.tstr, position=hl.tint32),
        alleles=hl.tarray(hl.tstr),
        exomes=hl.tstruct(filters=hl.tset(hl.tstr)),
        genomes=hl.tstruct(filters=hl.tset(hl.tstr)),
    )
    table_rows = []
    for pos, exome_filters, genome_filters in rows:
        table_rows.append({
            "locus": {"contig": "chr22", "position": pos},
            "alleles": ["A", "T"],
            "exomes": {"filters": set(exome_filters) if exome_filters is not None else None},
            "genomes": {"filters": set(genome_filters) if genome_filters is not None else None},
        })
    return hl.Table.parallelize(table_rows, schema=schema)


def test_apply_filters_joint41_keeps_when_genome_pass_and_exome_pass_or_only_ac0(
    hail_context: None,  # noqa: ARG001
) -> None:
    """
    JOINT_41 filter: genome must be empty/missing, exome must be empty/missing or only AC0.

    `AC0` on the exome side is dominated by exome capture-region absence (gnomAD assigns AC0
    when no sample meets the high-quality genotype threshold; for a whole-genome-sourced
    catalog like DivRef this is overwhelmingly explained by the position being outside the
    exome capture footprint). Treating `AC0`-only as PASS on the exome side recovers good
    genome-supported variants that the strict both-sides-empty intersection would discard.
    """
    ht = _make_joint41_filter_ht([
        (100, [], []),  # both PASS → kept
        (200, ["AC0"], []),  # exome only AC0, genome PASS → kept
        (300, ["AC0", "AS_VQSR"], []),  # exome has AC0+AS_VQSR, genome PASS → dropped
        (400, ["AS_VQSR"], []),  # exome AS_VQSR alone, genome PASS → dropped
        (500, [], ["AS_VQSR"]),  # exome PASS, genome non-empty → dropped
        (600, ["AC0"], ["AS_VQSR"]),  # exome AC0 OK but genome non-empty → dropped
        (700, None, []),  # exome missing, genome PASS → kept
        (800, [], None),  # exome PASS, genome missing → kept
        (900, None, None),  # both missing → kept
        (1000, ["AC0"], None),  # exome only AC0, genome missing → kept
        (1100, ["AS_VQSR"], None),  # exome AS_VQSR, genome missing → dropped
    ])

    survivors = sorted(
        r.locus.position for r in _apply_filters(ht, GnomadVersion.JOINT_41).collect()
    )
    assert survivors == [100, 200, 700, 800, 900, 1000]


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


def _make_joint41_sites_ht() -> hl.Table:
    """
    Build a two-variant JOINT_41-shaped sites table with the fields the tool reads.

    Both variants pass the AF threshold; the second fails the genome filter, so it survives only
    when the filters are not applied.
    """
    freq = hl.tarray(
        hl.tstruct(AC=hl.tint32, AF=hl.tfloat64, AN=hl.tint32, homozygote_count=hl.tint32)
    )
    grpmax = hl.tstruct(AC=hl.tint32, AF=hl.tfloat64, AN=hl.tint32, gen_anc=hl.tstr)
    schema = hl.tstruct(
        locus=hl.tlocus("GRCh38"),
        alleles=hl.tarray(hl.tstr),
        joint=hl.tstruct(freq=freq, grpmax=grpmax),
        exomes=hl.tstruct(filters=hl.tset(hl.tstr)),
        genomes=hl.tstruct(filters=hl.tset(hl.tstr)),
    )
    pop_freqs = [
        {"AC": 10, "AF": 0.1, "AN": 100, "homozygote_count": 0},
        {"AC": 20, "AF": 0.2, "AN": 100, "homozygote_count": 0},
    ]
    rows = [
        {
            "locus": hl.Locus("chr22", position, reference_genome="GRCh38"),
            "alleles": ["A", "T"],
            "joint": {
                "freq": pop_freqs,
                "grpmax": {"AC": 20, "AF": 0.2, "AN": 100, "gen_anc": "nfe"},
            },
            "exomes": {"filters": set()},
            "genomes": {"filters": genome_filters},
        }
        for position, genome_filters in [(20_000_000, set()), (20_000_100, {"AS_VQSR"})]
    ]
    ht = hl.Table.parallelize(rows, schema=schema, key=["locus", "alleles"])
    return ht.annotate_globals(
        joint_globals=hl.struct(
            freq_meta=[{"group": "adj", "gen_anc": "afr"}, {"group": "adj", "gen_anc": "nfe"}]
        )
    )


@pytest.mark.parametrize(
    ("no_apply_filters", "expected_rows"),
    [
        pytest.param(False, 1, id="filters_applied_and_recorded"),
        pytest.param(True, 2, id="filters_skipped_and_recorded"),
    ],
)
def test_extract_gnomad_single_afs_records_build_parameters(
    hail_context: None,  # noqa: ARG001
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    no_apply_filters: bool,
    expected_rows: int,
) -> None:
    """The written sites table records the run's settings, and `apply_filters` matches the rows."""
    monkeypatch.setattr("divref.tools.extract_gnomad_single_afs.hail_init", lambda **_: None)
    monkeypatch.setattr(
        "divref.tools.extract_gnomad_single_afs.hl.read_table",
        lambda _uri: _make_joint41_sites_ht(),
    )
    out = tmp_path / "sites.ht"
    extract_gnomad_single_afs(
        gnomad_version=GnomadVersion.JOINT_41,
        contig="chr22",
        freq_threshold=0.05,
        no_apply_filters=no_apply_filters,
        populations=["afr", "nfe"],
        out_sites_hail_table=out,
    )
    # The stub replaced `hl.read_table` on the shared hail module; restore it to read the output.
    monkeypatch.undo()
    written = hl.read_table(str(out))
    assert hl.eval(written.index_globals().build_parameters) == hl.Struct(
        gnomad_version="JOINT_41", freq_threshold=0.05, apply_filters=not no_apply_filters
    )
    assert written.count() == expected_rows
