"""Tests for the finalize_duckdb_index tool."""

import logging
from pathlib import Path

import duckdb
import pytest

from divref.tools.append_contig_to_duckdb_index import append_contig_to_duckdb_index
from divref.tools.finalize_duckdb_index import finalize_duckdb_index
from divref.tools.init_duckdb_index import init_duckdb_index


def test_finalize_warns_on_empty_sequences(
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Finalizing an index whose `sequences` table has no rows logs a warning (pure DuckDB)."""
    output_base = tmp_path / "idx"
    with duckdb.connect(str(Path(f"{output_base}.haplotypes_gnomad_merge.index.duckdb"))) as conn:
        conn.execute("CREATE TABLE sequences (sequence_id VARCHAR)")

    with caplog.at_level(logging.WARNING):
        finalize_duckdb_index(output_base=output_base)

    assert "empty 'sequences' table" in caplog.text


def _write_table_pairs_tsv(path: Path, rows: list[tuple[str, str, str]]) -> Path:
    """Write a table-pairs TSV with a (contig, haplotype_table_path, sites_table_path) per row."""
    lines = ["contig\thaplotype_table_path\tsites_table_path"]
    lines += [f"{contig}\t{hap}\t{sites}" for contig, hap, sites in rows]
    path.write_text("\n".join(lines) + "\n")
    return path


def _table_pairs(datadir: Path, tmp_path: Path) -> Path:
    """Build the chr1 table-pairs TSV (full haplotype + gnomAD sites pair)."""
    return _write_table_pairs_tsv(
        tmp_path / "table_pairs.tsv",
        rows=[
            (
                "chr1",
                str(datadir / "chr1_100001_200000_haplotypes.ht"),
                str(datadir / "chr1_100001_200000.gnomad_afs.ht"),
            ),
        ],
    )


def _reference_fasta(datadir: Path) -> Path:
    """Path to the committed masked mini-reference FASTA."""
    return datadir / "test_reference.chr1_chrX.fa.gz"


def _db_path(output_base: Path) -> Path:
    """The DuckDB index path for a given output base."""
    return Path(f"{output_base}.haplotypes_gnomad_merge.index.duckdb")


def test_finalize_creates_index(
    hail_context: None,  # noqa: ARG001
    datadir: Path,
    tmp_path: Path,
) -> None:
    """After init + one append, finalize creates the idx_sequence_id index."""
    table_pairs_tsv = _table_pairs(datadir, tmp_path)
    output_base = tmp_path / "idx"
    init_duckdb_index(
        in_table_pairs_tsv=table_pairs_tsv,
        output_base=output_base,
        version="9.9",
        window_size=25,
        force=True,
    )
    append_contig_to_duckdb_index(
        in_table_pairs_tsv=table_pairs_tsv,
        contig="chr1",
        output_base=output_base,
        reference_fasta=_reference_fasta(datadir),
        window_size=25,
        version="9.9",
    )

    finalize_duckdb_index(output_base=output_base)

    with duckdb.connect(str(_db_path(output_base))) as conn:
        idx = conn.execute(
            "SELECT index_name FROM duckdb_indexes() WHERE index_name = 'idx_sequence_id'"
        ).fetchone()
        assert idx is not None


def test_finalize_is_idempotent(
    hail_context: None,  # noqa: ARG001
    datadir: Path,
    tmp_path: Path,
) -> None:
    """Running finalize twice is a no-op the second time rather than raising on the index."""
    table_pairs_tsv = _table_pairs(datadir, tmp_path)
    output_base = tmp_path / "idx"
    init_duckdb_index(
        in_table_pairs_tsv=table_pairs_tsv,
        output_base=output_base,
        version="9.9",
        window_size=25,
        force=True,
    )
    append_contig_to_duckdb_index(
        in_table_pairs_tsv=table_pairs_tsv,
        contig="chr1",
        output_base=output_base,
        reference_fasta=_reference_fasta(datadir),
        window_size=25,
        version="9.9",
    )

    finalize_duckdb_index(output_base=output_base)
    finalize_duckdb_index(output_base=output_base)

    with duckdb.connect(str(_db_path(output_base))) as conn:
        idx = conn.execute(
            "SELECT index_name FROM duckdb_indexes() WHERE index_name = 'idx_sequence_id'"
        ).fetchone()
        assert idx is not None


def test_finalize_without_sequences_raises(
    hail_context: None,  # noqa: ARG001
    datadir: Path,
    tmp_path: Path,
) -> None:
    """Finalize on an init-only DB (no sequences table) raises ValueError."""
    table_pairs_tsv = _table_pairs(datadir, tmp_path)
    output_base = tmp_path / "idx"
    init_duckdb_index(
        in_table_pairs_tsv=table_pairs_tsv,
        output_base=output_base,
        version="9.9",
        window_size=25,
        force=True,
    )

    with pytest.raises(ValueError):
        finalize_duckdb_index(output_base=output_base)
