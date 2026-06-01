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


def _read_rows_sorted_by_sequence_number(tsv: Path) -> tuple[str, list[str]]:
    """
    Read a sequences TSV into its header and data rows sorted by numeric sequence number.

    Rows are keyed on the integer suffix of `sequence_id` (`DR-{version}-{n}`) rather than on the
    raw string, so the comparison does not depend on lexical-vs-numeric ordering of `sequence_id`
    in either file (e.g. lexical order interleaves `-10` between `-1` and `-2`).

    Args:
        tsv: Path to a sequences TSV with a header row.

    Returns:
        A tuple of (header line, data rows sorted by sequence number).
    """
    lines = tsv.read_text().splitlines()
    header = lines[0]
    sequence_id_col = header.split("\t").index("sequence_id")

    def sequence_number(row: str) -> int:
        return int(row.split("\t")[sequence_id_col].rsplit("-", 1)[1])

    return header, sorted(lines[1:], key=sequence_number)


def _assert_sequences_equivalent(actual: Path, expected: Path) -> None:
    """Assert two sequences TSVs hold identical rows, comparing in numeric sequence-id order."""
    actual_header, actual_rows = _read_rows_sorted_by_sequence_number(actual)
    expected_header, expected_rows = _read_rows_sorted_by_sequence_number(expected)

    assert actual_header == expected_header, (
        f"Header mismatch:\n  golden:   {expected_header!r}\n  produced: {actual_header!r}"
    )

    for got, want in zip(actual_rows, expected_rows, strict=False):
        assert got == want, f"Row mismatch:\n  golden:   {want!r}\n  produced: {got!r}"

    if len(actual_rows) != len(expected_rows):
        common = min(len(actual_rows), len(expected_rows))
        longer = actual_rows if len(actual_rows) > len(expected_rows) else expected_rows
        which = "produced" if len(actual_rows) > len(expected_rows) else "golden"
        raise AssertionError(
            f"Row count differs: produced {len(actual_rows)} rows, "
            f"golden has {len(expected_rows)} rows. "
            f"First extra row (from {which}):\n  {longer[common]!r}"
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
    _assert_sequences_equivalent(produced_tsv, golden_tsv)
