from pathlib import Path

import duckdb
import polars as pl

from divref.tools.create_duckdb_from_tsv import create_duckdb_from_tsv


def test_tsv_index_e2e(datadir: Path, tmp_path: Path) -> None:
    """Golden end-to-end check: the tool's `sequences` output matches a committed dump."""
    out_base = tmp_path / "test_cohort"
    create_duckdb_from_tsv(
        variants_tsv=datadir / "tsv_source" / "variants.tsv",
        source_meta=datadir / "tsv_source" / "source_meta.yml",
        output_base=out_base,
        reference_fasta=datadir / "test_reference.chr1_chrX.fa.gz",
        window_size=25,
        contigs=["chr1"],
        force=True,
    )

    conn = duckdb.connect(str(f"{out_base}.index.duckdb"), read_only=True)
    got = conn.execute("SELECT * FROM sequences ORDER BY sequence_id").pl()
    conn.close()

    exp = pl.read_csv(
        datadir / "tsv_source" / "tsv_index_golden.sequences.tsv",
        separator="\t",
        infer_schema_length=0,
    )
    assert (
        got.select(sorted(got.columns)).cast(pl.String).to_dicts()
        == exp.select(sorted(exp.columns)).cast(pl.String).to_dicts()
    )

    # Every row is a single-variant, PASS-filtered row from the "test_cohort" source, and the
    # per-pop annotation columns are prefixed with the source name.
    assert got["n_variants"].to_list() == [1] * got.height
    assert got["source"].to_list() == ["test_cohort"] * got.height
    assert got["haplotype_filter"].to_list() == ["PASS"] * got.height
    assert "test_cohort_AF_afr" in got.columns
    assert "test_cohort_AF_eas" in got.columns


def test_tsv_index_all_empty_contigs_yields_valid_empty_index(
    datadir: Path, tmp_path: Path
) -> None:
    """Requesting only contigs absent from the TSV finalizes a valid, empty index (no crash)."""
    out_base = tmp_path / "empty"
    create_duckdb_from_tsv(
        variants_tsv=datadir / "tsv_source" / "variants.tsv",  # chr1-only fixture
        source_meta=datadir / "tsv_source" / "source_meta.yml",
        output_base=out_base,
        reference_fasta=datadir / "test_reference.chr1_chrX.fa.gz",
        window_size=25,
        contigs=["chr2"],  # absent from the TSV, so every requested contig is empty
        force=True,
    )

    conn = duckdb.connect(str(f"{out_base}.index.duckdb"), read_only=True)
    try:
        count = conn.execute("SELECT COUNT(*) FROM sequences").fetchone()
    finally:
        conn.close()
    assert count is not None and count[0] == 0
