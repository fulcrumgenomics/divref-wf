"""Tests for the create_divref_fasta tool."""

import logging
from pathlib import Path

import duckdb
import pytest

from divref.tools.create_divref_fasta import create_divref_fasta


def _build_index(db_path: Path, rows: list[tuple[str, str, str]]) -> None:
    """
    Build a minimal DuckDB index with just the columns create_divref_fasta reads.

    Args:
        db_path: Path to the DuckDB file to create.
        rows: ``(sequence_id, sequence, contig)`` rows to insert into ``sequences``.
    """
    with duckdb.connect(str(db_path)) as conn:
        conn.execute(
            "CREATE TABLE sequences (sequence_id VARCHAR, sequence VARCHAR, contig VARCHAR)"
        )
        conn.executemany("INSERT INTO sequences VALUES (?, ?, ?)", rows)


def _parse_fasta(path: Path) -> list[tuple[str, str]]:
    """
    Parse a FASTA file into ``(id, sequence)`` records, asserting the two-line-per-record format.

    Args:
        path: Path to the FASTA file.

    Returns:
        One ``(sequence_id, sequence)`` tuple per record.
    """
    lines = path.read_text().splitlines()
    assert len(lines) % 2 == 0, "FASTA should have an even number of lines (header + sequence)"
    records: list[tuple[str, str]] = []
    for header, seq in zip(lines[0::2], lines[1::2], strict=True):
        assert header.startswith(">"), f"expected a header line, got {header!r}"
        records.append((header[1:], seq))
    return records


def test_writes_only_the_requested_contig(tmp_path: Path) -> None:
    """A FASTA is written for the requested contig with that contig's rows only."""
    db_path = tmp_path / "index.duckdb"
    _build_index(
        db_path,
        rows=[
            ("DR-1-0", "ACGT", "chr1"),
            ("DR-1-1", "TTTT", "chr1"),
            ("DR-1-2", "GGGG", "chr2"),
        ],
    )
    output_base = tmp_path / "out"

    create_divref_fasta(duckdb_path=db_path, output_base=output_base, contigs=["chr1"])

    chr1_fasta = Path(f"{output_base}.chr1.fasta")
    assert set(_parse_fasta(chr1_fasta)) == {("DR-1-0", "ACGT"), ("DR-1-1", "TTTT")}
    # chr2 was not requested, so its FASTA is not written.
    assert not Path(f"{output_base}.chr2.fasta").exists()


def test_writes_one_file_per_requested_contig(tmp_path: Path) -> None:
    """Each requested contig gets its own FASTA with the matching rows."""
    db_path = tmp_path / "index.duckdb"
    _build_index(
        db_path,
        rows=[
            ("DR-1-0", "ACGT", "chr1"),
            ("DR-1-1", "GGGG", "chr2"),
        ],
    )
    output_base = tmp_path / "out"

    create_divref_fasta(duckdb_path=db_path, output_base=output_base, contigs=["chr1", "chr2"])

    assert _parse_fasta(Path(f"{output_base}.chr1.fasta")) == [("DR-1-0", "ACGT")]
    assert _parse_fasta(Path(f"{output_base}.chr2.fasta")) == [("DR-1-1", "GGGG")]


def test_contig_with_no_rows_writes_empty_file_and_warns(
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """A requested contig absent from the index yields an empty FASTA plus a warning."""
    db_path = tmp_path / "index.duckdb"
    _build_index(db_path, rows=[("DR-1-0", "ACGT", "chr1")])
    output_base = tmp_path / "out"

    with caplog.at_level(logging.WARNING):
        create_divref_fasta(duckdb_path=db_path, output_base=output_base, contigs=["chrX"])

    out_path = Path(f"{output_base}.chrX.fasta")
    assert out_path.exists()
    assert out_path.read_text() == ""
    assert "No sequences found for contig chrX" in caplog.text


def test_empty_contigs_raises_value_error(tmp_path: Path) -> None:
    """An empty contig list is rejected."""
    db_path = tmp_path / "index.duckdb"
    _build_index(db_path, rows=[("DR-1-0", "ACGT", "chr1")])

    with pytest.raises(ValueError, match="[Cc]ontig"):
        create_divref_fasta(duckdb_path=db_path, output_base=tmp_path / "out", contigs=[])


def test_streams_all_rows_across_chunk_boundaries(tmp_path: Path) -> None:
    """All rows are written even when the result spans multiple read chunks."""
    db_path = tmp_path / "index.duckdb"
    rows = [(f"DR-1-{i}", "A" * (i + 1), "chr1") for i in range(5)]
    _build_index(db_path, rows=rows)
    output_base = tmp_path / "out"

    create_divref_fasta(
        duckdb_path=db_path,
        output_base=output_base,
        contigs=["chr1"],
        polars_chunk_size=2,
    )

    written = _parse_fasta(Path(f"{output_base}.chr1.fasta"))
    assert set(written) == {(sid, seq) for sid, seq, _ in rows}
    assert len(written) == 5
