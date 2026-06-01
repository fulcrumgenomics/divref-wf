"""
End-to-end equivalence test for the per-chromosome DuckDB index flow.

Runs the full init -> append(chr1) -> append(chrX) -> finalize flow over the committed
test fixtures and asserts the exported `sequences` table is identical, row for row, to a
golden TSV produced by the previous monolithic `create_duckdb_index` tool. This proves the
split per-chromosome tools reproduce the old tool's output exactly.
"""

from pathlib import Path

import duckdb

from divref.tools.append_contig_to_duckdb_index import append_contig_to_duckdb_index
from divref.tools.finalize_duckdb_index import finalize_duckdb_index
from divref.tools.init_duckdb_index import init_duckdb_index


def _write_table_pairs_tsv(path: Path, rows: list[tuple[str, str, str]]) -> Path:
    """Write a table-pairs TSV with a (contig, haplotype_table_path, sites_table_path) per row."""
    lines = ["contig\thaplotype_table_path\tsites_table_path"]
    lines += [f"{contig}\t{hap}\t{sites}" for contig, hap, sites in rows]
    path.write_text("\n".join(lines) + "\n")
    return path


def _db_path(output_base: Path) -> Path:
    """The DuckDB index path for a given output base."""
    return Path(f"{output_base}.haplotypes_gnomad_merge.index.duckdb")


def _assert_files_identical(actual: Path, expected: Path) -> None:
    """Assert two text files are line-for-line identical, reporting the first divergence."""
    actual_lines = actual.read_text().splitlines()
    expected_lines = expected.read_text().splitlines()

    for line_number, (got, want) in enumerate(
        zip(actual_lines, expected_lines, strict=False), start=1
    ):
        assert got == want, (
            f"Mismatch at line {line_number}:\n  golden:   {want!r}\n  produced: {got!r}"
        )

    common = min(len(actual_lines), len(expected_lines))
    if len(actual_lines) != len(expected_lines):
        longer = actual_lines if len(actual_lines) > len(expected_lines) else expected_lines
        which = "produced" if len(actual_lines) > len(expected_lines) else "golden"
        raise AssertionError(
            f"Line count differs: produced {len(actual_lines)} lines, "
            f"golden has {len(expected_lines)} lines. "
            f"First extra line (line {common + 1}, from {which}):\n  {longer[common]!r}"
        )


def test_new_flow_matches_golden(
    hail_context: None,  # noqa: ARG001
    datadir: Path,
    tmp_path: Path,
) -> None:
    """The split init->append->append->finalize flow reproduces the golden sequences TSV."""
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
    reference_fasta = datadir / "test_reference.chr1_chrX.fa.gz"
    output_base = tmp_path / "idx"

    init_duckdb_index(
        in_table_pairs_tsv=table_pairs_tsv,
        output_base=output_base,
        version="9.9",
        window_size=25,
        force=True,
    )
    # chr1 BEFORE chrX: the append order sets the sequence_id numbering and must match the golden.
    append_contig_to_duckdb_index(
        in_table_pairs_tsv=table_pairs_tsv,
        contig="chr1",
        output_base=output_base,
        reference_fasta=reference_fasta,
        window_size=25,
        version="9.9",
    )
    append_contig_to_duckdb_index(
        in_table_pairs_tsv=table_pairs_tsv,
        contig="chrX",
        output_base=output_base,
        reference_fasta=reference_fasta,
        window_size=25,
        version="9.9",
    )
    finalize_duckdb_index(output_base=output_base)

    produced_tsv = tmp_path / "produced.sequences.tsv"
    with duckdb.connect(str(_db_path(output_base))) as conn:
        conn.execute(
            "COPY (SELECT * FROM sequences ORDER BY sequence_id) "
            f"TO '{produced_tsv}' (HEADER, DELIMITER E'\\t')"
        )

    golden_tsv = datadir / "duckdb_index_golden" / "sequences.chr1_chrX.tsv"
    _assert_files_identical(produced_tsv, golden_tsv)
