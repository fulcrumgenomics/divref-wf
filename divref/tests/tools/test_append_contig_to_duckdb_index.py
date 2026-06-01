"""Tests for the append_contig_to_duckdb_index tool."""

import shutil
from pathlib import Path

import duckdb
import pytest

from divref.tools.append_contig_to_duckdb_index import append_contig_to_duckdb_index
from divref.tools.init_duckdb_index import init_duckdb_index


def _write_table_pairs_tsv(path: Path, rows: list[tuple[str, str, str]]) -> Path:
    """Write a table-pairs TSV with a (contig, haplotype_table_path, sites_table_path) per row."""
    lines = ["contig\thaplotype_table_path\tsites_table_path"]
    lines += [f"{contig}\t{hap}\t{sites}" for contig, hap, sites in rows]
    path.write_text("\n".join(lines) + "\n")
    return path


def _table_pairs(datadir: Path, tmp_path: Path) -> Path:
    """Build the two-contig table-pairs TSV (chr1 = full pair, chrX = sites-only)."""
    return _write_table_pairs_tsv(
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


def _reference_fasta(datadir: Path) -> Path:
    """Path to the committed masked mini-reference FASTA."""
    return datadir / "test_reference.chr1_chrX.fa.gz"


def _db_path(output_base: Path) -> Path:
    """The DuckDB index path for a given output base."""
    return Path(f"{output_base}.haplotypes_gnomad_merge.index.duckdb")


def test_append_single_contig_creates_sequences(
    hail_context: None,  # noqa: ARG001
    datadir: Path,
    tmp_path: Path,
) -> None:
    """Appending one contig creates `sequences` with rows whose ids start at DR-9.9-0."""
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

    with duckdb.connect(str(_db_path(output_base))) as conn:
        count_row = conn.execute("SELECT COUNT(*) FROM sequences").fetchone()
        assert count_row is not None
        assert count_row[0] > 0
        first_row = conn.execute(
            "SELECT sequence_id FROM sequences ORDER BY sequence_id LIMIT 1"
        ).fetchone()
        assert first_row is not None
        assert first_row[0] == "DR-9.9-0"


def test_append_two_contigs_contiguous_ids(
    hail_context: None,  # noqa: ARG001
    datadir: Path,
    tmp_path: Path,
) -> None:
    """Appending chr1 then chrX yields contiguous DR-9.9-0..N-1 ids continuing the offset."""
    table_pairs_tsv = _table_pairs(datadir, tmp_path)
    output_base = tmp_path / "idx"
    reference_fasta = _reference_fasta(datadir)
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
        reference_fasta=reference_fasta,
        window_size=25,
        version="9.9",
    )
    with duckdb.connect(str(_db_path(output_base))) as conn:
        chr1_count_row = conn.execute("SELECT COUNT(*) FROM sequences").fetchone()
        assert chr1_count_row is not None
        chr1_count: int = chr1_count_row[0]

    append_contig_to_duckdb_index(
        in_table_pairs_tsv=table_pairs_tsv,
        contig="chrX",
        output_base=output_base,
        reference_fasta=reference_fasta,
        window_size=25,
        version="9.9",
    )

    with duckdb.connect(str(_db_path(output_base))) as conn:
        total_row = conn.execute("SELECT COUNT(*) FROM sequences").fetchone()
        assert total_row is not None
        total: int = total_row[0]
        assert total > chr1_count

        ids = [
            row[0]
            for row in conn.execute(
                "SELECT sequence_id FROM sequences ORDER BY CAST(SPLIT_PART("
                "sequence_id, '-', 3) AS INTEGER)"
            ).fetchall()
        ]
        # The ids form the contiguous set DR-9.9-0 .. DR-9.9-(total-1).
        assert ids == [f"DR-9.9-{i}" for i in range(total)]

        # chrX's rows are exactly those beyond chr1's count; its first id continues the offset.
        chrx_first_row = conn.execute(
            "SELECT sequence_id FROM sequences WHERE contig = 'chrX' "
            "ORDER BY CAST(SPLIT_PART(sequence_id, '-', 3) AS INTEGER) LIMIT 1"
        ).fetchone()
        assert chrx_first_row is not None
        assert chrx_first_row[0] == f"DR-9.9-{chr1_count}"


def test_append_legend_mismatch_raises(
    hail_context: None,  # noqa: ARG001
    datadir: Path,
    tmp_path: Path,
) -> None:
    """A stored legend that disagrees with the contig's source legend raises ValueError."""
    table_pairs_tsv = _table_pairs(datadir, tmp_path)
    output_base = tmp_path / "idx"
    init_duckdb_index(
        in_table_pairs_tsv=table_pairs_tsv,
        output_base=output_base,
        version="9.9",
        window_size=25,
        force=True,
    )

    # Corrupt a stored legend in a copy of the DB so the chr1 source legends no longer match.
    db = _db_path(output_base)
    corrupt_output_base = tmp_path / "corrupt"
    corrupt_db = _db_path(corrupt_output_base)
    shutil.copy(db, corrupt_db)
    with duckdb.connect(str(corrupt_db)) as conn:
        conn.execute("DROP TABLE gnomad_variant_pops_legend")
        conn.execute(
            "CREATE TABLE gnomad_variant_pops_legend AS SELECT ? AS pops_legend",
            ['["totally", "wrong"]'],
        )

    with pytest.raises(ValueError):
        append_contig_to_duckdb_index(
            in_table_pairs_tsv=table_pairs_tsv,
            contig="chr1",
            output_base=corrupt_output_base,
            reference_fasta=_reference_fasta(datadir),
            window_size=25,
            version="9.9",
        )
