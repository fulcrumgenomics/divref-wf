"""Shared helpers for building the DivRef DuckDB index from haplotype Hail tables."""

import json
import logging

import duckdb

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
    hgdp_pops_legend: list[str],
    gnomad_pops_legend: list[str],
    joint_pops_legend: list[str],
    version: str,
) -> None:
    """
    Write the window_size, three *_pops_legend, and VERSION metadata tables.

    Args:
        conn: Open DuckDB connection to the index being initialized.
        window_size: Flanking reference-context size stored in the `window_size` table.
        hgdp_pops_legend: HGDP-source population codes stored as JSON in
            `hgdp_haplotype_pops_legend`.
        gnomad_pops_legend: gnomAD-source population codes stored as JSON in
            `gnomad_variant_pops_legend`.
        joint_pops_legend: Joint population codes stored as JSON in `joint_pops_legend`.
        version: Version identifier stored in the `VERSION` table.
    """
    # Write the five tables in one transaction so an interrupted init leaves no partially
    # populated index (which a later append/finalize would then read as corrupt metadata).
    conn.execute("BEGIN TRANSACTION")
    try:
        conn.execute("CREATE TABLE window_size AS SELECT ? AS window_size", [window_size])
        conn.execute(
            "CREATE TABLE hgdp_haplotype_pops_legend AS SELECT ? AS pops_legend",
            [json.dumps(hgdp_pops_legend)],
        )
        conn.execute(
            "CREATE TABLE gnomad_variant_pops_legend AS SELECT ? AS pops_legend",
            [json.dumps(gnomad_pops_legend)],
        )
        conn.execute(
            "CREATE TABLE joint_pops_legend AS SELECT ? AS pops_legend",
            [json.dumps(joint_pops_legend)],
        )
        conn.execute("CREATE TABLE VERSION AS SELECT ? AS version", [version])
        conn.execute("COMMIT")
    except BaseException:
        conn.execute("ROLLBACK")
        raise
