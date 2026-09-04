"""Tests for the source-agnostic sequences-TSV ingest helpers in duckdb_index."""

from pathlib import Path

import duckdb
import polars
import pytest

from divref.duckdb_index import _stream_tsv_into_sequences
from divref.duckdb_index import sequences_tsv_columns
from divref.duckdb_index import with_compatibility_flag


@pytest.mark.parametrize(
    ("af_prefix", "popmax_estimated_col", "expected"),
    [
        pytest.param(
            "gnomAD",
            "popmax_estimated_gnomad_AF",
            [
                "sequence",
                "sequence_length",
                "sequence_id",
                "n_variants",
                "contig",
                "start",
                "end",
                "popmax_empirical_AF",
                "popmax_empirical_AC",
                "source",
                "popmax_estimated_gnomad_AF",
                "popmax_fraction_phased",
                "max_pop",
                "variants",
                "gnomAD_AF_afr",
                "gnomAD_AF_eas",
                "empirical_AC_afr",
                "empirical_AF_afr",
                "fraction_phased_afr",
                "estimated_gnomAD_haplotype_AF_afr",
                "empirical_AC_eas",
                "empirical_AF_eas",
                "fraction_phased_eas",
                "estimated_gnomAD_haplotype_AF_eas",
            ],
            id="gnomad-defaults",
        ),
        pytest.param(
            "src",
            "popmax_estimated_src_AF",
            [
                "sequence",
                "sequence_length",
                "sequence_id",
                "n_variants",
                "contig",
                "start",
                "end",
                "popmax_empirical_AF",
                "popmax_empirical_AC",
                "source",
                "popmax_estimated_src_AF",
                "popmax_fraction_phased",
                "max_pop",
                "variants",
                "src_AF_afr",
                "src_AF_eas",
                "empirical_AC_afr",
                "empirical_AF_afr",
                "fraction_phased_afr",
                "estimated_src_haplotype_AF_afr",
                "empirical_AC_eas",
                "empirical_AF_eas",
                "fraction_phased_eas",
                "estimated_src_haplotype_AF_eas",
            ],
            id="source-prefixed",
        ),
    ],
)
def test_sequences_tsv_columns(
    af_prefix: str, popmax_estimated_col: str, expected: list[str]
) -> None:
    """Column order: 14 scalars, then per-pop AF, then per-pop families grouped by pop."""
    assert (
        sequences_tsv_columns(
            ["afr", "eas"], af_prefix=af_prefix, popmax_estimated_col=popmax_estimated_col
        )
        == expected
    )


def test_stream_empty_tsv_creates_sequences_table(tmp_path: Path) -> None:
    """
    A contig that yields zero rows still leaves a valid (empty) sequences table behind.

    Guards the finalize step: it must always find a `sequences` table, even if the first appended
    contig produced no rows. No Hail needed — this exercises the polars/DuckDB streaming path
    directly with a header-only TSV.
    """
    joint_pops_legend = ["afr"]
    # Mirror the real export header (export_sequences_table_to_tsv), including the `variants`,
    # `source`, and `n_variants` columns the compatibility-flag step reads.
    columns = [
        "sequence",
        "sequence_length",
        "sequence_id",
        "n_variants",
        "contig",
        "start",
        "end",
        "popmax_empirical_AF",
        "popmax_empirical_AC",
        "source",
        "popmax_estimated_gnomad_AF",
        "popmax_fraction_phased",
        "max_pop",
        "variants",
        "gnomAD_AF_afr",
        "empirical_AC_afr",
        "empirical_AF_afr",
        "fraction_phased_afr",
        "estimated_gnomAD_haplotype_AF_afr",
    ]
    header_only_tsv = tmp_path / "empty.sequences.tsv"
    header_only_tsv.write_text("\t".join(columns) + "\n")

    db = tmp_path / "idx.duckdb"
    with duckdb.connect(str(db)) as conn:
        appended = _stream_tsv_into_sequences(
            conn,
            tsv=header_only_tsv,
            joint_pops_legend=joint_pops_legend,
            chunk_size=100,
        )
        assert appended == 0
        exists = conn.execute(
            "SELECT 1 FROM information_schema.tables WHERE table_name = 'sequences'"
        ).fetchone()
        assert exists is not None
        count = conn.execute("SELECT COUNT(*) FROM sequences").fetchone()
        assert count is not None
        assert count[0] == 0


def test_stream_tsv_flags_overlapping_haplotype(tmp_path: Path) -> None:
    """
    Streaming a TSV with an overlapping haplotype persists its incompatibility reason.

    The committed e2e test data has no overlapping haplotypes, so this injects the case at the
    TSV -> DuckDB boundary where the flag is applied: an overlapping HGDP haplotype (a SNP at a
    deleted base) is flagged `snp_in_deletion`, while a clean haplotype and a `gnomAD_variant`
    row stay `PASS`.
    """
    columns = [
        "sequence",
        "sequence_length",
        "sequence_id",
        "n_variants",
        "contig",
        "start",
        "end",
        "popmax_empirical_AF",
        "popmax_empirical_AC",
        "source",
        "popmax_estimated_gnomad_AF",
        "popmax_fraction_phased",
        "max_pop",
        "variants",
        "gnomAD_AF_afr",
        "empirical_AC_afr",
        "empirical_AF_afr",
        "fraction_phased_afr",
        "estimated_gnomAD_haplotype_AF_afr",
    ]
    rows = [
        # overlapping: a SNP at a base the deletion removes -> snp_in_deletion
        [
            "ACGT",
            "4",
            "DR-9.9-0",
            "2",
            "chr1",
            "274",
            "326",
            "0.05",
            "7",
            "HGDP_haplotype",
            "0.004",
            "0.5",
            "afr",
            "chr1:300:AT:A,chr1:301:T:A",
            "0.1,0.2",
            "7",
            "0.05",
            "0.5",
            "0.004",
        ],
        # clean two-variant haplotype -> PASS
        [
            "ACGT",
            "4",
            "DR-9.9-1",
            "2",
            "chr1",
            "174",
            "235",
            "0.05",
            "5",
            "HGDP_haplotype",
            "0.004",
            "0.5",
            "afr",
            "chr1:200:A:T,chr1:210:C:G",
            "0.1,0.2",
            "5",
            "0.05",
            "0.5",
            "0.004",
        ],
        # single gnomAD variant -> PASS
        [
            "A",
            "1",
            "DR-9.9-2",
            "1",
            "chr1",
            "24",
            "76",
            "0.1",
            "100",
            "gnomAD_variant",
            "0.1",
            "1.0",
            "afr",
            "chr1:50:A:T",
            "0.1",
            "100",
            "0.1",
            "1.0",
            "0.1",
        ],
    ]
    tsv = tmp_path / "with_rows.sequences.tsv"
    tsv.write_text("\t".join(columns) + "\n" + "\n".join("\t".join(r) for r in rows) + "\n")

    db = tmp_path / "idx.duckdb"
    with duckdb.connect(str(db)) as conn:
        appended = _stream_tsv_into_sequences(
            conn, tsv=tsv, joint_pops_legend=["afr"], chunk_size=100
        )
        assert appended == 3
        flags = dict(conn.execute("SELECT sequence_id, haplotype_filter FROM sequences").fetchall())
    assert flags == {
        "DR-9.9-0": "snp_in_deletion",
        "DR-9.9-1": "PASS",
        "DR-9.9-2": "PASS",
    }


def test_with_compatibility_flag_values() -> None:
    """PASS for gnomAD/single/clean rows; the incompatibility reason for an overlapping pair."""
    df = polars.DataFrame({
        "variants": [
            "chr1:300:AT:A,chr1:301:T:A",  # SNP at a deleted base
            "chr1:200:A:T,chr1:210:C:G",  # clean haplotype
            "chr1:50:A:T",  # gnomAD single variant
        ],
        "source": ["HGDP_haplotype", "HGDP_haplotype", "gnomAD_variant"],
        "n_variants": [2, 2, 1],
    })
    out = with_compatibility_flag(df)
    assert out["haplotype_filter"].to_list() == ["snp_in_deletion", "PASS", "PASS"]


def test_with_compatibility_flag_empty_frame() -> None:
    """An empty frame (used to create the table) still gains a String haplotype_filter column."""
    df = polars.DataFrame(
        {"variants": [], "source": [], "n_variants": []},
        schema={
            "variants": polars.String,
            "source": polars.String,
            "n_variants": polars.Int64,
        },
    )
    out = with_compatibility_flag(df)
    assert out.height == 0
    assert out.schema["haplotype_filter"] == polars.String
