"""Shared helpers for building the DivRef DuckDB index from haplotype Hail tables."""

import json
import logging
from pathlib import Path

import duckdb
import hail as hl
from fgmetric import Metric

logger = logging.getLogger(__name__)


class TablePair(Metric):
    """
    Helper class to link a pair of tables for the same contig.

    `haplotype_table_path` may be `None` (empty TSV cell) for contigs that contribute only gnomAD
    single variants (e.g. chrX/chrY in the divref workflow). `sites_table_path` is always required.

    Attributes:
        contig: Contig name.
        haplotype_table_path: HGDP haplotypes Hail table, or `None` if this contig has no haplotype
            track.
        sites_table_path: gnomAD variant Hail table.
    """

    contig: str
    haplotype_table_path: Path | None
    sites_table_path: Path


def read_pops_legend(table_path: Path) -> list[str]:
    """
    Read a Hail table's population legend from its globals.

    `pops` is stored as a global (`globals.pops`), so this reads only the table's globals
    file and does not scan rows.

    Args:
        table_path: Path to a Hail table with a `pops` global.

    Returns:
        The ordered population codes.
    """
    return list(hl.eval(hl.read_table(str(table_path)).index_globals().pops))


def read_and_validate_pops_legends(table_pairs: list[TablePair]) -> tuple[list[str], list[str]]:
    """
    Read and cross-contig-validate the gnomAD and HGDP population legends.

    Reads only `globals.pops` of each input table (no row scan). The gnomAD legend is taken from
    the first pair; the HGDP legend bootstraps from the first pair that has a haplotype table
    (`[]` if none). Every other pair must share the same gnomAD legend, and every pair with a
    haplotype table must share the same HGDP legend.

    Args:
        table_pairs: The per-contig table pairs read from the input TSV.

    Returns:
        A tuple of `(gnomad_pops_legend, hgdp_pops_legend)`.

    Raises:
        ValueError: If any contig's gnomAD or HGDP legend disagrees with the bootstrapped legend.
    """
    first_with_hap: TablePair | None = next(
        (tp for tp in table_pairs if tp.haplotype_table_path is not None), None
    )
    hgdp_pops_legend: list[str] = []
    if first_with_hap is not None:
        assert first_with_hap.haplotype_table_path is not None  # narrowed by the next() predicate
        hgdp_pops_legend = read_pops_legend(first_with_hap.haplotype_table_path)
    gnomad_pops_legend: list[str] = read_pops_legend(table_pairs[0].sites_table_path)

    # All pairs must share the same pops legends so a single remap into the joint legend is valid
    # for every contig; otherwise the exported gnomAD_AF_* columns would be misaligned. Rows
    # without a haplotype table are skipped on the haplotype-side check.
    for tp in table_pairs[1:]:
        tp_gnomad_pops: list[str] = read_pops_legend(tp.sites_table_path)
        if tp_gnomad_pops != gnomad_pops_legend:
            raise ValueError(
                f"gnomAD pops legend mismatch for contig {tp.contig}: "
                f"{tp_gnomad_pops} vs {gnomad_pops_legend}."
            )
    for tp in table_pairs:
        if tp is first_with_hap or tp.haplotype_table_path is None:
            continue
        tp_hgdp_pops: list[str] = read_pops_legend(tp.haplotype_table_path)
        if tp_hgdp_pops != hgdp_pops_legend:
            raise ValueError(
                f"HGDP haplotype pops legend mismatch for contig {tp.contig}: "
                f"{tp_hgdp_pops} vs {hgdp_pops_legend}."
            )

    return gnomad_pops_legend, hgdp_pops_legend


def compute_joint_legend(gnomad_pops: list[str], hgdp_pops: list[str]) -> list[str]:
    """Joint legend: gnomAD pops in original order, then HGDP-only pops appended."""
    return list(gnomad_pops) + [p for p in hgdp_pops if p not in gnomad_pops]


def write_metadata_tables(
    conn: duckdb.DuckDBPyConnection,
    *,
    window_size: int,
    hgdp_pops_legend: list[str],
    gnomad_pops_legend: list[str],
    joint_pops_legend: list[str],
    version: str,
) -> None:
    """Write the window_size, three *_pops_legend, and VERSION metadata tables."""
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
