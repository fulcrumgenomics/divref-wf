"""Shared DuckDB helpers for the DivRef index build."""

import json
import logging
from collections.abc import Iterator
from pathlib import Path

import duckdb
import polars

from divref.haplotype_compat import compatibility_flag

logger = logging.getLogger(__name__)


def sequences_table_exists(conn: duckdb.DuckDBPyConnection) -> bool:
    """
    Return whether the `sequences` table exists in the connected DuckDB index.

    Args:
        conn: Open DuckDB connection to a DivRef index.

    Returns:
        True if a `sequences` table is present.
    """
    return (
        conn.execute(
            "SELECT 1 FROM information_schema.tables WHERE table_name = 'sequences'"
        ).fetchone()
        is not None
    )


def write_metadata_tables(
    conn: duckdb.DuckDBPyConnection,
    *,
    window_size: int,
    haplotype_pops_legend: list[str],
    variant_pops_legend: list[str],
    joint_pops_legend: list[str],
    annotation_af_prefix: str,
    version: str,
) -> None:
    """
    Write the window_size, three *_pops_legend, annotation_af_prefix, and VERSION metadata tables.

    Args:
        conn: Open DuckDB connection to the index being initialized.
        window_size: Flanking reference-context size stored in the `window_size` table.
        haplotype_pops_legend: Haplotype-source population codes stored as JSON in
            `haplotype_pops_legend`.
        variant_pops_legend: Variant-annotation-source population codes stored as JSON in
            `variant_pops_legend`.
        joint_pops_legend: Joint population codes stored as JSON in `joint_pops_legend`.
        annotation_af_prefix: Annotation-AF column prefix stored in `annotation_af_prefix`.
        version: Version identifier stored in the `VERSION` table.
    """
    # Write the six tables in one transaction so an interrupted init leaves no partially
    # populated index (which a later append/finalize would then read as corrupt metadata).
    conn.execute("BEGIN TRANSACTION")
    try:
        conn.execute("CREATE TABLE window_size AS SELECT ? AS window_size", [window_size])
        conn.execute(
            "CREATE TABLE haplotype_pops_legend AS SELECT ? AS pops_legend",
            [json.dumps(haplotype_pops_legend)],
        )
        conn.execute(
            "CREATE TABLE variant_pops_legend AS SELECT ? AS pops_legend",
            [json.dumps(variant_pops_legend)],
        )
        conn.execute(
            "CREATE TABLE joint_pops_legend AS SELECT ? AS pops_legend",
            [json.dumps(joint_pops_legend)],
        )
        conn.execute(
            "CREATE TABLE annotation_af_prefix AS SELECT ? AS af_prefix",
            [annotation_af_prefix],
        )
        conn.execute("CREATE TABLE VERSION AS SELECT ? AS version", [version])
        conn.execute("COMMIT")
    except BaseException:
        conn.execute("ROLLBACK")
        raise


_LEGEND_TABLES = frozenset({
    "joint_pops_legend",
    "variant_pops_legend",
    "haplotype_pops_legend",
})


def read_legend(conn: duckdb.DuckDBPyConnection, table: str) -> list[str]:
    """Read a stored *_pops_legend table back into a list of pop codes."""
    if table not in _LEGEND_TABLES:
        raise ValueError(
            f"Unknown legend table {table!r}; expected one of {sorted(_LEGEND_TABLES)}."
        )
    # `table` is checked against the allowlist above, so the interpolation is safe; bandit (S608)
    # still flags the f-string statically, hence the suppression.
    row = conn.execute(f"SELECT pops_legend FROM {table}").fetchone()  # noqa: S608
    if row is None:
        raise ValueError(f"Metadata table {table} is empty; run init_duckdb_index first.")
    return list(json.loads(row[0]))


def read_window_size(conn: duckdb.DuckDBPyConnection) -> int:
    """Read the stored window_size metadata value."""
    row = conn.execute("SELECT window_size FROM window_size").fetchone()
    if row is None:
        raise ValueError("Metadata table window_size is empty; run init_duckdb_index first.")
    return int(row[0])


def sequences_row_count(conn: duckdb.DuckDBPyConnection) -> int:
    """Current number of rows in `sequences`, or 0 if the table does not exist yet."""
    if not sequences_table_exists(conn):
        return 0
    row = conn.execute("SELECT COUNT(*) FROM sequences").fetchone()
    if row is None:
        raise ValueError("COUNT(*) on sequences returned no row.")
    return int(row[0])


def contig_already_appended(conn: duckdb.DuckDBPyConnection, contig: str) -> bool:
    """Whether the `sequences` table already holds any rows for `contig` (False if no table yet)."""
    if not sequences_table_exists(conn):
        return False
    row = conn.execute("SELECT 1 FROM sequences WHERE contig = ? LIMIT 1", [contig]).fetchone()
    return row is not None


def create_sequence_id_index(conn: duckdb.DuckDBPyConnection) -> None:
    """
    Create the idempotent `idx_sequence_id` index on `sequences`.

    `IF NOT EXISTS` makes this a no-op on an already-finalized index (a re-triggered Snakemake rule
    or manual retry) rather than a raw DuckDB "index already exists" error.
    """
    conn.execute("CREATE INDEX IF NOT EXISTS idx_sequence_id ON sequences(sequence_id)")


def scan_sequences_tsv(
    tsv: Path,
    joint_pops_legend: list[str],
    *,
    af_prefix: str = "gnomAD",
    popmax_estimated_col: str = "popmax_estimated_gnomad_AF",
) -> polars.LazyFrame:
    """
    Build the typed polars LazyFrame for a sequences TSV.

    The `sequence_id` and `{af_prefix}_AF_*` columns are explicitly typed as strings so that
    schema inference cannot misread comma-delimited per-variant AFs as floats.

    Every numeric column is pinned explicitly as well. Each contig is appended in its own process
    with independent polars inference, and the `sequences` table's schema is fixed by whichever
    contig is appended first (`CREATE TABLE ... AS SELECT`). A per-pop column that is all-`NA` for
    one contig (e.g. an HGDP-only pop on a sites-only contig like chrX) would otherwise infer a
    different dtype than the same column on a contig where it is populated, and the later
    `INSERT INTO` would fail or coerce. Pinning the integer and float columns makes the schema
    contig-order-independent.

    Args:
        tsv: Path to the sequences TSV (plain or bgz-compressed).
        joint_pops_legend: Ordered list of population codes used to name `{af_prefix}_AF_{pop}`
            columns; must match what `export_sequences_table_to_tsv` wrote.
        af_prefix: Source prefix for the per-variant AF and estimated-haplotype-AF columns.
            Defaults to the gnomAD annotation-source naming.
        popmax_estimated_col: Name of the popmax estimated-AF column. Defaults to the gnomAD
            annotation-source naming.

    Returns:
        A LazyFrame over `tsv` with the full sequences schema applied.
    """
    schema_overrides: dict[str, type[polars.DataType]] = {
        "sequence_id": polars.String,
        "sequence_length": polars.Int64,
        "n_variants": polars.Int64,
        "start": polars.Int64,
        "end": polars.Int64,
        "popmax_empirical_AF": polars.Float64,
        "popmax_empirical_AC": polars.Int64,
        popmax_estimated_col: polars.Float64,
        "popmax_fraction_phased": polars.Float64,
        **{f"{af_prefix}_AF_{pop}": polars.String for pop in joint_pops_legend},
        **{f"empirical_AC_{pop}": polars.Int64 for pop in joint_pops_legend},
        **{f"empirical_AF_{pop}": polars.Float64 for pop in joint_pops_legend},
        **{f"fraction_phased_{pop}": polars.Float64 for pop in joint_pops_legend},
        **{
            f"estimated_{af_prefix}_haplotype_AF_{pop}": polars.Float64 for pop in joint_pops_legend
        },
    }
    # Hail's TSV export emits "NA" for missing scalar fields; "null" is included for
    # robustness against other writers.
    lf = polars.scan_csv(
        tsv,
        separator="\t",
        schema_overrides=schema_overrides,
        null_values=["NA", "null"],
    )
    # `null_values` applies globally and can convert a bare "NA" cell to null even though
    # the column is declared as String in `schema_overrides`. Restore "NA" so downstream
    # consumers (e.g. `remap_divref.Haplotype`, which types `gnomad_afs` as `dict[str, str]`)
    # always see a string. This matters mostly for single-variant rows where the per-pop
    # cell can degenerate to a bare "NA".
    lf = lf.with_columns([
        polars.col(f"{af_prefix}_AF_{pop}").fill_null("NA").cast(polars.String)
        for pop in joint_pops_legend
    ])
    return lf


def iter_dataframe_chunks(
    *,
    tsv: Path,
    joint_pops_legend: list[str],
    chunk_size: int,
    af_prefix: str = "gnomAD",
    popmax_estimated_col: str = "popmax_estimated_gnomad_AF",
) -> Iterator[polars.DataFrame]:
    """
    Yield non-empty polars DataFrames of up to `chunk_size` rows from a sequences TSV.

    Args:
        tsv: Path to the sequences TSV (plain or bgz-compressed).
        joint_pops_legend: Ordered list of population codes used to name `{af_prefix}_AF_{pop}`
            columns; must match what `export_sequences_table_to_tsv` wrote.
        chunk_size: Maximum rows per yielded DataFrame.
        af_prefix: Source prefix for the per-variant AF and estimated-haplotype-AF columns.
            Defaults to the gnomAD annotation-source naming.
        popmax_estimated_col: Name of the popmax estimated-AF column. Defaults to the gnomAD
            annotation-source naming.

    Yields:
        Polars DataFrame batches read from `tsv`.
    """
    lf = scan_sequences_tsv(
        tsv,
        joint_pops_legend,
        af_prefix=af_prefix,
        popmax_estimated_col=popmax_estimated_col,
    )
    for df in lf.collect_batches(chunk_size=chunk_size):
        if df.height > 0:
            yield df


def with_compatibility_flag(df: polars.DataFrame) -> polars.DataFrame:
    """
    Append the `haplotype_filter` column to a sequences DataFrame.

    The value is VCF-style: `PASS` for `gnomAD_variant` rows, single-variant rows, and
    haplotypes whose component variants do not overlap; otherwise the `;`-joined incompatibility
    reason(s) (e.g. `snp_in_deletion`). Computed with `divref.haplotype_compat.compatibility_flag`
    over the row's `variants` string; the classifier runs only on multi-variant HGDP rows.

    Args:
        df: A sequences DataFrame with `variants`, `source`, and `n_variants` columns. May be empty.

    Returns:
        The DataFrame with a `haplotype_filter` String column appended.
    """
    # A plain Python loop with a cheap source/n_variants short-circuit, deliberately not a polars
    # when/then + map_elements: the latter evaluates the classifier branch for every row, whereas
    # only multi-variant HGDP rows can be incompatible. The bulk (single gnomAD_variant rows) thus
    # skips parsing/classification entirely. Keep this structure rather than "vectorising" it.
    flags = [
        "PASS" if (source == "gnomAD_variant" or n_variants < 2) else compatibility_flag(variants)
        for variants, source, n_variants in zip(
            df["variants"].to_list(),
            df["source"].to_list(),
            df["n_variants"].to_list(),
            strict=True,
        )
    ]
    return df.with_columns(polars.Series("haplotype_filter", flags, dtype=polars.String))


def sequences_tsv_columns(
    joint_pops_legend: list[str], *, af_prefix: str, popmax_estimated_col: str
) -> list[str]:
    """
    The sequences-TSV column order that `export_sequences_table_to_tsv` selects by.

    The pure-Python TSV builder consumes the same list, so the two column layouts cannot drift.

    Args:
        joint_pops_legend: Ordered list of population codes used to name the per-pop columns.
        af_prefix: Source prefix for the per-variant AF and estimated-haplotype-AF columns.
        popmax_estimated_col: Name of the popmax estimated-AF scalar column.

    Returns:
        Column names in the order `export_sequences_table_to_tsv` writes them.
    """
    scalar_columns = [
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
        popmax_estimated_col,
        "popmax_fraction_phased",
        "max_pop",
        "variants",
    ]
    af_columns = [f"{af_prefix}_AF_{pop}" for pop in joint_pops_legend]
    per_pop_columns = [
        column
        for pop in joint_pops_legend
        for column in (
            f"empirical_AC_{pop}",
            f"empirical_AF_{pop}",
            f"fraction_phased_{pop}",
            f"estimated_{af_prefix}_haplotype_AF_{pop}",
        )
    ]
    return scalar_columns + af_columns + per_pop_columns


def _stream_tsv_into_sequences(
    conn: duckdb.DuckDBPyConnection,
    *,
    tsv: Path,
    joint_pops_legend: list[str],
    chunk_size: int,
    af_prefix: str = "gnomAD",
    popmax_estimated_col: str = "popmax_estimated_gnomad_AF",
) -> int:
    """
    Stream a per-contig sequences TSV into the DuckDB `sequences` table.

    Creates the `sequences` table from the typed schema if it does not exist yet, then appends each
    batch of `chunk_size` rows. Creating the table up front (from a zero-row typed frame rather than
    the first non-empty batch) means a contig that yields no rows still leaves a valid `sequences`
    table behind, so `finalize_duckdb_index` never depends on some contig having had data.

    Args:
        conn: Open connection to the DuckDB index.
        tsv: Path to the per-contig sequences TSV produced by `export_sequences_table_to_tsv`.
        joint_pops_legend: Ordered joint pop legend used to type the `{af_prefix}_AF_*` columns.
        chunk_size: Maximum number of rows per polars read batch.
        af_prefix: Source prefix for the per-variant AF and estimated-haplotype-AF columns.
        popmax_estimated_col: Name of the popmax estimated-AF scalar column.

    Returns:
        The number of rows appended for this contig.
    """
    if not sequences_table_exists(conn):
        # `haplotype_filter` is appended last in both the schema-defining empty frame and every
        # inserted chunk, so column order stays consistent across the CREATE and the INSERTs.
        empty_df = with_compatibility_flag(
            scan_sequences_tsv(
                tsv,
                joint_pops_legend,
                af_prefix=af_prefix,
                popmax_estimated_col=popmax_estimated_col,
            )
            .limit(0)
            .collect()
        )
        conn.register("empty_df", empty_df)
        conn.execute("CREATE TABLE sequences AS SELECT * FROM empty_df")
        conn.unregister("empty_df")

    appended_rows: int = 0
    for chunk in iter_dataframe_chunks(
        tsv=tsv,
        joint_pops_legend=joint_pops_legend,
        chunk_size=chunk_size,
        af_prefix=af_prefix,
        popmax_estimated_col=popmax_estimated_col,
    ):
        chunk_df = with_compatibility_flag(chunk)
        conn.register("chunk_df", chunk_df)
        # BY NAME matches columns by name, not position, so a future change to the export/scan
        # column order surfaces as a bind error rather than silently misaligning columns.
        conn.execute("INSERT INTO sequences BY NAME SELECT * FROM chunk_df")
        conn.unregister("chunk_df")
        appended_rows += chunk_df.height
    return appended_rows


def stream_sequences_tsv_into_duckdb(
    conn: duckdb.DuckDBPyConnection,
    *,
    tsv: Path,
    joint_pops_legend: list[str],
    chunk_size: int,
    af_prefix: str = "gnomAD",
    popmax_estimated_col: str = "popmax_estimated_gnomad_AF",
) -> int:
    """
    Stream a sequences TSV into `sequences` within one transaction; return rows appended.

    Args:
        conn: Open connection to the DuckDB index.
        tsv: Path to the per-contig sequences TSV produced by `export_sequences_table_to_tsv`.
        joint_pops_legend: Ordered joint pop legend used to type the `{af_prefix}_AF_*` columns.
        chunk_size: Maximum number of rows per polars read batch.
        af_prefix: Source prefix for the per-variant AF and estimated-haplotype-AF columns.
            Defaults to the gnomAD annotation-source naming.
        popmax_estimated_col: Name of the popmax estimated-AF scalar column. Defaults to the
            gnomAD annotation-source naming.

    Returns:
        The number of rows appended for this contig.
    """
    conn.execute("BEGIN TRANSACTION")
    try:
        rows = _stream_tsv_into_sequences(
            conn,
            tsv=tsv,
            joint_pops_legend=joint_pops_legend,
            chunk_size=chunk_size,
            af_prefix=af_prefix,
            popmax_estimated_col=popmax_estimated_col,
        )
        conn.execute("COMMIT")
    except BaseException:
        conn.execute("ROLLBACK")
        raise
    return rows
