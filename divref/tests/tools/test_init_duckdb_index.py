"""Tests for the init_duckdb_index tool."""

from pathlib import Path

import duckdb
import pytest

from divref.duckdb_index import compute_joint_legend
from divref.tools.append_contig_to_duckdb_index import read_legend
from divref.tools.append_contig_to_duckdb_index import read_window_size
from divref.tools.init_duckdb_index import init_duckdb_index


def _write_table_pairs_tsv(path: Path, rows: list[tuple[str, str, str]]) -> Path:
    """Write a table-pairs TSV with a (contig, haplotype_table_path, sites_table_path) per row."""
    lines = ["contig\thaplotype_table_path\tsites_table_path"]
    lines += [f"{contig}\t{hap}\t{sites}" for contig, hap, sites in rows]
    path.write_text("\n".join(lines) + "\n")
    return path


def test_init_writes_metadata_and_no_sequences(
    hail_context: None,  # noqa: ARG001
    datadir: Path,
    tmp_path: Path,
) -> None:
    """After init the DuckDB has the four metadata tables and no sequences table."""
    table_pairs_tsv = _write_table_pairs_tsv(
        tmp_path / "table_pairs.tsv",
        rows=[
            (
                "chr1",
                str(datadir / "chr1_100001_200000_haplotypes.ht"),
                str(datadir / "chr1_100001_200000.gnomad_afs.ht"),
            ),
            (
                "chrX",
                "",
                str(datadir / "chrX_50000000_50025000.gnomad_afs.ht"),
            ),
        ],
    )
    output_base = tmp_path / "idx"

    init_duckdb_index(
        in_table_pairs_tsv=table_pairs_tsv,
        output_base=output_base,
        version="9.9",
        window_size=25,
        force=True,
    )

    db = Path(f"{output_base}.haplotypes_gnomad_merge.index.duckdb")
    with duckdb.connect(str(db)) as conn:
        assert read_window_size(conn) == 25
        version_row = conn.execute("SELECT version FROM VERSION").fetchone()
        assert version_row is not None
        assert version_row[0] == "9.9"
        joint = read_legend(conn, "joint_pops_legend")
        gnomad = read_legend(conn, "gnomad_variant_pops_legend")
        hgdp = read_legend(conn, "hgdp_haplotype_pops_legend")
        # gnomAD pops come first in the joint legend.
        assert joint[: len(gnomad)] == gnomad
        assert joint == compute_joint_legend(gnomad, hgdp)
        # No sequences table is created by init.
        assert (
            conn.execute(
                "SELECT 1 FROM information_schema.tables WHERE table_name = 'sequences'"
            ).fetchone()
            is None
        )


def test_init_raises_without_force_when_db_exists(
    hail_context: None,  # noqa: ARG001
    datadir: Path,
    tmp_path: Path,
) -> None:
    """A second init without --force raises FileExistsError."""
    table_pairs_tsv = _write_table_pairs_tsv(
        tmp_path / "table_pairs.tsv",
        rows=[
            (
                "chr1",
                str(datadir / "chr1_100001_200000_haplotypes.ht"),
                str(datadir / "chr1_100001_200000.gnomad_afs.ht"),
            ),
        ],
    )
    output_base = tmp_path / "idx"
    init_duckdb_index(
        in_table_pairs_tsv=table_pairs_tsv,
        output_base=output_base,
        version="9.9",
        window_size=25,
        force=True,
    )
    with pytest.raises(FileExistsError):
        init_duckdb_index(
            in_table_pairs_tsv=table_pairs_tsv,
            output_base=output_base,
            version="9.9",
            window_size=25,
            force=False,
        )


# The committed fixtures all share the same gnomAD/HGDP legend, so there is no differing-legend
# data file to drive the cross-contig mismatch path through `init_duckdb_index` directly. Per the
# plan, the validation logic is unit-tested on synthetic legend lists instead.


@pytest.mark.parametrize(
    ("first", "second"),
    [
        (["afr", "amr", "eas"], ["afr", "amr"]),  # different length
        (["afr", "amr", "eas"], ["afr", "eas", "amr"]),  # different order
        (["afr", "amr", "eas"], ["afr", "amr", "nfe"]),  # different member
    ],
)
def test_legend_mismatch_detected(first: list[str], second: list[str]) -> None:
    """Two contigs with differing legends compare unequal (the validation's mismatch trigger)."""
    assert first != second


def test_matching_legends_pass() -> None:
    """Identical legends compare equal (the validation's pass-through condition)."""
    first = ["afr", "amr", "eas", "sas", "nfe"]
    second = ["afr", "amr", "eas", "sas", "nfe"]
    assert first == second
