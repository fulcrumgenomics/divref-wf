"""Tool to initialize the DivRef DuckDB index file and write its metadata tables."""

import logging
import os
from collections.abc import Mapping
from dataclasses import fields
from pathlib import Path

import duckdb
import hail as hl
from fgpyo.io import assert_directory_exists
from fgpyo.io import assert_path_is_readable
from fgpyo.io import assert_path_is_writable
from hail.context import Env

from divref.duckdb_index import VariantBuildParameters
from divref.duckdb_index import write_metadata_tables
from divref.gnomad_index_source import TablePair
from divref.gnomad_index_source import compute_joint_legend
from divref.gnomad_index_source import read_and_validate_pops_legends
from divref.gnomad_index_source import read_haplotype_build_parameters
from divref.gnomad_index_source import read_variant_build_parameters

logger = logging.getLogger(__name__)


def _check_consistent_variant_build_parameters(
    variant_build_parameters: Mapping[str, VariantBuildParameters],
) -> None:
    """
    Check that every contig's single-variant track came from the same source and settings.

    The index has one `annotation_af_prefix` and one variant legend, so the single-variant AF
    columns must mean the same thing on every contig. A NULL (unrecorded) value never conflicts.

    Args:
        variant_build_parameters: Per-contig `extract_gnomad_single_afs` parameters.

    Raises:
        ValueError: If two contigs record different non-NULL values for a parameter.
    """
    for field in fields(VariantBuildParameters):
        values = {
            contig: getattr(parameters, field.name)
            for contig, parameters in variant_build_parameters.items()
            if getattr(parameters, field.name) is not None
        }
        if len(set(values.values())) > 1:
            raise ValueError(f"Contigs disagree on {field.name}: {values}.")


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

    Reads only the globals of each input Hail table (no row scan), validates that every contig
    shares the same gnomAD and HGDP population legends, computes the joint legend, and writes the
    `window_size`, `haplotype_pops_legend`, `variant_pops_legend`, `joint_pops_legend`,
    `annotation_af_prefix`, and `VERSION` tables. It also writes `haplotype_build_parameters`, one
    row per haplotype contig, from each haplotype table's `build_parameters` global (NULLs, with a
    warning, for a table built before that global existed). It also writes
    `variant_build_parameters`, one row per contig, from each sites table's `build_parameters`
    global, with the same NULL handling; those values must agree across contigs. Does not create
    `sequences` — the first `append_contig_to_duckdb_index` does that.

    Args:
        in_table_pairs_tsv: TSV with 'contig', 'haplotype_table_path' (optional), and
            'sites_table_path' columns.
        output_base: Base path; writes `{output_base}.haplotypes_gnomad_merge.index.duckdb`.
        version: Version identifier embedded later in sequence IDs.
        window_size: Flanking reference-context size stored in the index.
        force: Overwrite an existing DuckDB; otherwise raise FileExistsError.

    Raises:
        FileExistsError: If the output DuckDB already exists and `force` is False.
        ValueError: If `in_table_pairs_tsv` contains no table pairs or lists a contig twice, if
            the contigs' gnomAD or HGDP population legends disagree, or if their recorded
            single-variant build parameters disagree.
    """
    assert_path_is_readable(in_table_pairs_tsv)

    out_duckdb_file: Path = Path(f"{output_base}.haplotypes_gnomad_merge.index.duckdb")
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
    seen_contigs: set[str] = set()
    for table_pair in table_pairs:
        if table_pair.contig in seen_contigs:
            raise ValueError(f"Duplicate contig {table_pair.contig} in {in_table_pairs_tsv}.")
        seen_contigs.add(table_pair.contig)

    # fail fast on input Hail tables; haplotype_table_path is optional per row
    for table_pair in table_pairs:
        if table_pair.haplotype_table_path is not None:
            assert_directory_exists(table_pair.haplotype_table_path)
        assert_directory_exists(table_pair.sites_table_path)

    # Light Hail init for the globals-only legend reads. Skip if a context already exists (e.g. a
    # shared test-session context) so this stays idempotent within a process. No `tmp_dir` is
    # passed and 1g memory is hardcoded: this reads only `globals.pops` (no row scan, no
    # checkpoints), unlike `append_contig_to_duckdb_index` which needs both.
    if Env._hc is None:
        os.environ["PYSPARK_SUBMIT_ARGS"] = "--driver-memory 1g --executor-memory 1g pyspark-shell"
        hl.init()

    # Read each table's globals-only pop legend and validate cross-contig consistency.
    gnomad_pops_legend, hgdp_pops_legend = read_and_validate_pops_legends(table_pairs)
    joint_pops_legend: list[str] = compute_joint_legend(gnomad_pops_legend, hgdp_pops_legend)
    haplotype_build_parameters = {
        table_pair.contig: read_haplotype_build_parameters(table_pair.haplotype_table_path)
        for table_pair in table_pairs
        if table_pair.haplotype_table_path is not None
    }
    variant_build_parameters = {
        table_pair.contig: read_variant_build_parameters(table_pair.sites_table_path)
        for table_pair in table_pairs
    }
    _check_consistent_variant_build_parameters(variant_build_parameters)

    with duckdb.connect(str(out_duckdb_file)) as conn:
        write_metadata_tables(
            conn,
            window_size=window_size,
            haplotype_pops_legend=hgdp_pops_legend,
            variant_pops_legend=gnomad_pops_legend,
            joint_pops_legend=joint_pops_legend,
            annotation_af_prefix="gnomAD",
            version=version,
            haplotype_build_parameters=haplotype_build_parameters,
            variant_build_parameters=variant_build_parameters,
        )

    logger.info(
        f"Initialized DuckDB index {out_duckdb_file} "
        f"(joint legend: {joint_pops_legend}, window_size: {window_size}, version: {version})."
    )
