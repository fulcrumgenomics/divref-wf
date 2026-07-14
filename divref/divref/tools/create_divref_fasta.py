"""Tool to write per-chromosome FASTA files from a DivRef DuckDB index."""

import logging
from collections.abc import Iterator
from pathlib import Path

import duckdb
import polars
from fgpyo.io import assert_path_is_readable
from fgpyo.io import assert_path_is_writable

from divref.duckdb_index import sequences_table_exists

logger = logging.getLogger(__name__)


def create_divref_fasta(
    *,
    duckdb_path: Path,
    output_base: Path,
    contigs: list[str],
    polars_chunk_size: int = 100_000,
) -> None:
    """
    Write per-chromosome FASTA files from a DivRef DuckDB index.

    Streams sequences for each contig out of the DuckDB sequences table in batches of
    ``polars_chunk_size`` rows and writes them directly to ``{output_base}.{contig}.fasta`` without
    materialising the per-contig result as a single in-process DataFrame. If a requested contig
    has no rows in the index, an empty FASTA is written and a warning is logged.

    Args:
        duckdb_path: Path to an existing DivRef DuckDB index.
        output_base: Base path for output FASTA files; chromosome name is appended as a suffix.
        contigs: List of contigs to write FASTA files for.
        polars_chunk_size: Maximum number of rows per polars DataFrame batch read from DuckDB.

    Raises:
        ValueError: If ``contigs`` is empty, or the index has no ``sequences`` table.
    """
    if not contigs:
        raise ValueError("Contig list must be provided.")

    assert_path_is_readable(duckdb_path)

    # Validate all output paths before processing
    out_paths: dict[str, Path] = {
        contig: Path(f"{output_base}.{contig}.fasta") for contig in contigs
    }
    for out_path in out_paths.values():
        assert_path_is_writable(out_path)

    with duckdb.connect(str(duckdb_path), read_only=True) as conn:
        if not sequences_table_exists(conn):
            raise ValueError(
                f"DuckDB index {duckdb_path} has no 'sequences' table; "
                f"run append_contig_to_duckdb_index and finalize_duckdb_index first."
            )
        for contig in contigs:
            logger.info(f"Creating FASTA for chromosome {contig} at {out_paths[contig]}")
            rows_written: int = 0
            with out_paths[contig].open("w") as fh:
                for df in iter_sequence_chunks(
                    conn=conn, contig=contig, chunk_size=polars_chunk_size
                ):
                    for sequence_id, sequence in df.iter_rows():
                        fh.write(f">{sequence_id}\n{sequence}\n")
                    rows_written += df.height
            if rows_written == 0:
                logger.warning(
                    f"No sequences found for contig {contig}; wrote empty FASTA at "
                    f"{out_paths[contig]}"
                )
            else:
                logger.info(f"Wrote {rows_written} sequences for {contig}")


def iter_sequence_chunks(
    *,
    conn: duckdb.DuckDBPyConnection,
    contig: str,
    chunk_size: int,
) -> Iterator[polars.DataFrame]:
    """
    Yield polars DataFrames of ``(sequence_id, sequence)`` rows for one contig, in genomic order.

    Rows are ordered by genomic ``start`` (then ``sequence_id`` to break ties) so the FASTA output
    is deterministic across runs. Streams the DuckDB result set as Arrow record batches via
    ``to_arrow_reader`` and converts each batch to polars, bounding in-process memory by
    ``chunk_size`` rows.

    Args:
        conn: Open DuckDB connection to the DivRef index.
        contig: Contig to filter on.
        chunk_size: Maximum rows per yielded DataFrame.

    Yields:
        Polars DataFrame batches with ``sequence_id`` and ``sequence`` columns.
    """
    # ORDER BY genomic start (contig is fixed by the filter) so FASTA output is deterministic
    # across runs; sequence_id is a unique tiebreak for rows sharing a start position.
    result = conn.execute(
        "SELECT sequence_id, sequence FROM sequences WHERE contig = $contig "
        "ORDER BY start, sequence_id",
        {"contig": contig},
    )
    for batch in result.to_arrow_reader(chunk_size):
        df = polars.from_arrow(batch)
        # from_arrow on a RecordBatch always returns a DataFrame; assert to narrow the type and
        # fail loudly rather than silently dropping rows should that ever change.
        assert isinstance(df, polars.DataFrame)
        if df.height > 0:
            yield df
