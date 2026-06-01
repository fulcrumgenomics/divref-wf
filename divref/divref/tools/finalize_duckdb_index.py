"""Tool to finalize the DivRef DuckDB index by creating its sequence-id index."""

import logging
from pathlib import Path

import duckdb
from fgpyo.io import assert_path_is_readable

logger = logging.getLogger(__name__)


def finalize_duckdb_index(
    *,
    output_base: Path,
) -> None:
    """
    Create the `idx_sequence_id` index on the DivRef DuckDB `sequences` table.

    This is pure DuckDB work: it does not initialize Hail or Spark. The `sequences` table must
    already exist (created by `append_contig_to_duckdb_index`); otherwise a `ValueError` is raised.

    Args:
        output_base: Base path; reads `{output_base}.haplotypes_gnomad_merge.index.duckdb`.

    Raises:
        ValueError: If the DuckDB index has no `sequences` table.
    """
    out_duckdb_file: Path = Path(f"{str(output_base)}.haplotypes_gnomad_merge.index.duckdb")
    assert_path_is_readable(out_duckdb_file)

    with duckdb.connect(str(out_duckdb_file)) as conn:
        sequences_exists = conn.execute(
            "SELECT 1 FROM information_schema.tables WHERE table_name = 'sequences'"
        ).fetchone()
        if sequences_exists is None:
            raise ValueError(
                f"DuckDB index {out_duckdb_file} has no 'sequences' table; "
                f"run append_contig_to_duckdb_index before finalizing."
            )
        conn.execute("CREATE INDEX idx_sequence_id ON sequences(sequence_id)")

    logger.info(f"Created idx_sequence_id index on {out_duckdb_file}.")
