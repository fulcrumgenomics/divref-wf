"""Tool to initialize the DivRef DuckDB index file and write its metadata tables."""

import logging
import os
from pathlib import Path

import duckdb
import hail as hl
from fgpyo.io import assert_directory_exists
from fgpyo.io import assert_path_is_readable
from fgpyo.io import assert_path_is_writable
from hail.context import Env

from divref.duckdb_index import TablePair
from divref.duckdb_index import compute_joint_legend
from divref.duckdb_index import read_and_validate_pops_legends
from divref.duckdb_index import write_metadata_tables

logger = logging.getLogger(__name__)


def init_duckdb_index(
    *,
    in_table_pairs_tsv: Path,
    output_base: Path,
    version: str,
    window_size: int,
    force: bool = False,
) -> None:
    """
    Create the DuckDB index file and write its population-legend + version metadata.

    Reads only the `globals.pops` of each input Hail table (no row scan), validates that every
    contig shares the same gnomAD and HGDP population legends, computes the joint legend, and
    writes the `window_size`, `hgdp_haplotype_pops_legend`, `gnomad_variant_pops_legend`,
    `joint_pops_legend`, and `VERSION` tables. Does not create `sequences` — the first
    `append_contig_to_duckdb_index` does that.

    Args:
        in_table_pairs_tsv: TSV with 'contig', 'haplotype_table_path' (optional), and
            'sites_table_path' columns.
        output_base: Base path; writes `{output_base}.haplotypes_gnomad_merge.index.duckdb`.
        version: Version identifier embedded later in sequence IDs.
        window_size: Flanking reference-context size stored in the index.
        force: Overwrite an existing DuckDB; otherwise raise FileExistsError.
    """
    assert_path_is_readable(in_table_pairs_tsv)

    out_duckdb_file: Path = Path(f"{str(output_base)}.haplotypes_gnomad_merge.index.duckdb")
    if out_duckdb_file.exists():
        if not force:
            raise FileExistsError(
                f"DuckDB output already exists at {out_duckdb_file}. Pass --force to overwrite."
            )
        out_duckdb_file.unlink()
    assert_path_is_writable(out_duckdb_file)

    table_pairs: list[TablePair] = list(TablePair.read(in_table_pairs_tsv))
    if not table_pairs:
        raise ValueError(f"No table pairs found in {in_table_pairs_tsv}.")

    # fail fast on input Hail tables; haplotype_table_path is optional per row
    for table_pair in table_pairs:
        if table_pair.haplotype_table_path is not None:
            assert_directory_exists(table_pair.haplotype_table_path)
        assert_directory_exists(table_pair.sites_table_path)

    # Light Hail init for the globals-only legend reads. Skip if a context already exists (e.g. a
    # shared test-session context) so this stays idempotent within a process.
    if Env._hc is None:
        os.environ["PYSPARK_SUBMIT_ARGS"] = "--driver-memory 1g --executor-memory 1g pyspark-shell"
        hl.init()

    # Read each table's globals-only pop legend and validate cross-contig consistency.
    gnomad_pops_legend, hgdp_pops_legend = read_and_validate_pops_legends(table_pairs)
    joint_pops_legend: list[str] = compute_joint_legend(gnomad_pops_legend, hgdp_pops_legend)

    with duckdb.connect(str(out_duckdb_file)) as conn:
        write_metadata_tables(
            conn,
            window_size=window_size,
            hgdp_pops_legend=hgdp_pops_legend,
            gnomad_pops_legend=gnomad_pops_legend,
            joint_pops_legend=joint_pops_legend,
            version=version,
        )

    logger.info(
        f"Initialized DuckDB index {out_duckdb_file} "
        f"(joint legend: {joint_pops_legend}, window_size: {window_size}, version: {version})."
    )
