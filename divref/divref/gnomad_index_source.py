"""Hail/gnomAD-specific per-contig legend inputs for the DuckDB index build."""

import logging
from dataclasses import fields
from pathlib import Path

import hail as hl
from fgmetric import Metric

from divref.duckdb_index import HaplotypeBuildParameters
from divref.duckdb_index import VariantBuildParameters

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


def _read_build_parameters[ParametersT: (HaplotypeBuildParameters, VariantBuildParameters)](
    *,
    table_path: Path,
    parameters_type: type[ParametersT],
    table_kind: str,
    rerun_tools: str,
) -> ParametersT:
    """
    Read a table's `build_parameters` global into `parameters_type`.

    Reads only the table's globals file. A table built before the global existed yields
    all-None parameters and a warning. A `build_parameters` without a field (built before that
    field existed) yields None for it and a warning.

    Args:
        table_path: Path to a Hail table.
        parameters_type: The build-parameters dataclass to fill.
        table_kind: Table description for the warning (e.g. "Haplotype").
        rerun_tools: Tools the warning says to re-run to record the parameters.

    Returns:
        The recorded parameters: all None when the table has no `build_parameters` global, and
        None for each field absent from it.
    """
    table_globals = hl.read_table(str(table_path)).index_globals()
    names = [field.name for field in fields(parameters_type)]
    if "build_parameters" not in table_globals:
        logger.warning(
            "%s table %s has no build_parameters global, so its build parameters are unknown "
            "(NULL). Re-run %s to record them.",
            table_kind,
            table_path,
            rerun_tools,
        )
        return parameters_type(**dict.fromkeys(names))
    recorded = hl.eval(table_globals.build_parameters)
    missing = [name for name in names if name not in recorded]
    if missing:
        logger.warning(
            "%s table %s build_parameters has no %s, so those parameters are unknown (NULL). "
            "Re-run %s to record them.",
            table_kind,
            table_path,
            ", ".join(missing),
            rerun_tools,
        )
    return parameters_type(**{name: recorded.get(name) for name in names})


def read_haplotype_build_parameters(haplotype_table_path: Path) -> HaplotypeBuildParameters:
    """
    Read the `build_parameters` global that `compute_haplotypes` writes on its output table.

    Args:
        haplotype_table_path: Path to a haplotype Hail table.

    Returns:
        The recorded parameters: all None when the table has no `build_parameters` global, and
        `sites_freq_threshold` None when that field is absent.
    """
    return _read_build_parameters(
        table_path=haplotype_table_path,
        parameters_type=HaplotypeBuildParameters,
        table_kind="Haplotype",
        rerun_tools="extract_gnomad_afs and compute_haplotypes",
    )


def read_variant_build_parameters(sites_table_path: Path) -> VariantBuildParameters:
    """
    Read the `build_parameters` global that `extract_gnomad_single_afs` writes on its output table.

    Args:
        sites_table_path: Path to a gnomAD sites Hail table.

    Returns:
        The recorded parameters, or all None when the table has no `build_parameters` global.
    """
    return _read_build_parameters(
        table_path=sites_table_path,
        parameters_type=VariantBuildParameters,
        table_kind="Sites",
        rerun_tools="extract_gnomad_single_afs",
    )


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
        ValueError: If `table_pairs` is empty, or if any contig's gnomAD or HGDP legend disagrees
            with the bootstrapped legend.
    """
    if not table_pairs:
        raise ValueError("read_and_validate_pops_legends requires at least one table pair.")

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
    """
    Compute the joint population legend across both variation sources.

    Args:
        gnomad_pops: gnomAD-source population codes, in their original order.
        hgdp_pops: HGDP-source population codes.

    Returns:
        The joint legend: every gnomAD population in its original order, followed by the HGDP
        populations not already present (e.g. `["afr", "nfe"]` + `["afr", "oth"]` ->
        `["afr", "nfe", "oth"]`).
    """
    return list(gnomad_pops) + [p for p in hgdp_pops if p not in gnomad_pops]
