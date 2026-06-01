"""Tool to append a single contig's sequences to an existing DivRef DuckDB index."""

import logging
import os
from pathlib import Path

import duckdb
import hail as hl
from fgpyo.io import assert_directory_exists
from fgpyo.io import assert_path_is_readable
from hail.context import Env

from divref import defaults
from divref.duckdb_index import TablePair
from divref.duckdb_index import at_joint
from divref.duckdb_index import build_contig_sequences_table
from divref.duckdb_index import export_sequences_table_to_tsv
from divref.duckdb_index import iter_dataframe_chunks
from divref.duckdb_index import read_and_validate_pops_legends
from divref.duckdb_index import read_legend
from divref.duckdb_index import read_window_size
from divref.duckdb_index import sequences_row_count
from divref.duckdb_index import to_joint

logger = logging.getLogger(__name__)


class _RemapArrays:
    """
    The four pop-legend remap arrays used to build a contig's sequences table.

    Attributes:
        hgdp_to_joint: For each haplotype-source pop index, its index in the joint legend.
        gnomad_to_joint: For each gnomAD-source pop index, its index in the joint legend.
        hgdp_at_joint: For each joint pop index, the haplotype-source index or -1 if absent.
        gnomad_at_joint: For each joint pop index, the gnomAD-source index or -1 if absent.
    """

    def __init__(
        self,
        *,
        hgdp_to_joint: list[int],
        gnomad_to_joint: list[int],
        hgdp_at_joint: list[int],
        gnomad_at_joint: list[int],
    ) -> None:
        self.hgdp_to_joint = hgdp_to_joint
        self.gnomad_to_joint = gnomad_to_joint
        self.hgdp_at_joint = hgdp_at_joint
        self.gnomad_at_joint = gnomad_at_joint


def _resolve_legends_and_remaps(
    conn: duckdb.DuckDBPyConnection,
    table_pair: TablePair,
) -> tuple[list[str], _RemapArrays]:
    """
    Read the stored legends, validate the contig's source legends, and build the remap arrays.

    Reads back the three stored population legends from the DuckDB, re-reads this contig's source
    legends via `read_and_validate_pops_legends` (bootstrapping the HGDP legend from this contig
    only, so a sites-only contig yields an empty HGDP legend), and checks them against the stored
    gnomAD and HGDP legends. The HGDP check is skipped when this contig has no haplotype table.

    Args:
        conn: Open connection to the DuckDB index initialized by `init_duckdb_index`.
        table_pair: The single contig's haplotype + gnomAD sites table pair.

    Returns:
        A tuple of the joint pop legend and the four remap arrays into that joint legend.

    Raises:
        ValueError: If this contig's gnomAD or HGDP source legend disagrees with the stored legend.
    """
    stored_gnomad_pops: list[str] = read_legend(conn, "gnomad_variant_pops_legend")
    stored_hgdp_pops: list[str] = read_legend(conn, "hgdp_haplotype_pops_legend")
    joint_pops_legend: list[str] = read_legend(conn, "joint_pops_legend")

    # Re-read this contig's source legends with the shared validator. Passing a single-element list
    # bootstraps the HGDP legend from this contig alone (`[]` for a sites-only contig).
    contig_gnomad_pops, contig_hgdp_pops = read_and_validate_pops_legends([table_pair])

    if contig_gnomad_pops != stored_gnomad_pops:
        raise ValueError(
            f"gnomAD pops legend mismatch for contig {table_pair.contig}: "
            f"{contig_gnomad_pops} vs stored {stored_gnomad_pops}."
        )
    if table_pair.haplotype_table_path is not None and contig_hgdp_pops != stored_hgdp_pops:
        raise ValueError(
            f"HGDP haplotype pops legend mismatch for contig {table_pair.contig}: "
            f"{contig_hgdp_pops} vs stored {stored_hgdp_pops}."
        )

    remaps = _RemapArrays(
        hgdp_to_joint=to_joint(stored_hgdp_pops, joint_pops_legend),
        gnomad_to_joint=to_joint(stored_gnomad_pops, joint_pops_legend),
        hgdp_at_joint=at_joint(stored_hgdp_pops, joint_pops_legend),
        gnomad_at_joint=at_joint(stored_gnomad_pops, joint_pops_legend),
    )
    return joint_pops_legend, remaps


def _stream_tsv_into_sequences(
    conn: duckdb.DuckDBPyConnection,
    *,
    tsv: Path,
    joint_pops_legend: list[str],
    chunk_size: int,
) -> int:
    """
    Stream a per-contig sequences TSV into the DuckDB `sequences` table.

    Reads the TSV in batches of `chunk_size` rows. The first batch (when `sequences` has no rows
    yet) creates the table; every later batch appends to it.

    Args:
        conn: Open connection to the DuckDB index.
        tsv: Path to the per-contig sequences TSV produced by `export_sequences_table_to_tsv`.
        joint_pops_legend: Ordered joint pop legend used to type the `gnomAD_AF_*` columns.
        chunk_size: Maximum number of rows per polars read batch.

    Returns:
        The number of rows appended for this contig.
    """
    appended_rows: int = 0
    for df in iter_dataframe_chunks(
        tsv=tsv,
        joint_pops_legend=joint_pops_legend,
        chunk_size=chunk_size,
    ):
        if sequences_row_count(conn) == 0:
            conn.execute("CREATE TABLE sequences AS SELECT * FROM df")
        else:
            conn.execute("INSERT INTO sequences SELECT * FROM df")
        appended_rows += df.height
    return appended_rows


def append_contig_to_duckdb_index(
    *,
    in_table_pairs_tsv: Path,
    contig: str,
    output_base: Path,
    reference_fasta: Path,
    window_size: int,
    version: str,
    reference_genome: str = defaults.REFERENCE_GENOME,
    tmp_dir: Path = Path("/tmp"),
    polars_chunk_size: int = 100_000,
    retain_per_contig_tsvs: bool = False,
    spark_driver_memory_gb: int = 1,
    spark_executor_memory_gb: int = 1,
) -> None:
    """
    Append one contig's sequences to an existing DivRef DuckDB index.

    Opens the DuckDB created by `init_duckdb_index`, validates this contig's source population
    legends against the stored gnomAD and HGDP legends, builds the contig's sequences Hail table,
    and streams it into the `sequences` table. The first append creates `sequences`; subsequent
    appends `INSERT INTO` it. Sequence IDs continue the global offset given by the current
    `sequences` row count, so running contigs in canonical order reproduces a contiguous
    `DR-{version}-N` numbering across processes.

    Args:
        in_table_pairs_tsv: TSV with 'contig', 'haplotype_table_path' (optional), and
            'sites_table_path' columns. Only the row matching `contig` is processed.
        contig: The contig to append. Must match a row in the table-pairs TSV.
        output_base: Base path; reads/writes `{output_base}.haplotypes_gnomad_merge.index.duckdb`.
        reference_fasta: Path to the indexed reference FASTA for sequence extraction. Its FASTA
            index is read from the `.fai` sibling path.
        window_size: Flanking reference context size around each haplotype/variant. Must match the
            value stored in the index by `init_duckdb_index`.
        version: Version identifier embedded in sequence IDs (e.g. "1.0").
        reference_genome: Reference genome to use. Defaults to "GRCh38".
        tmp_dir: Temporary directory for Hail checkpoints and (when not retained) the per-contig
            intermediate TSV.
        polars_chunk_size: Maximum number of rows per polars read batch when streaming the
            per-contig TSV into DuckDB.
        retain_per_contig_tsvs: If True, write the per-contig TSV alongside the DuckDB output
            rather than into `tmp_dir`, and do not delete it.
        spark_driver_memory_gb: Memory in GB to allocate to the Spark driver.
        spark_executor_memory_gb: Memory in GB to allocate to the Spark executor.

    Raises:
        ValueError: If Spark memory is below 1GB, the contig is not in the TSV, the stored window
            size disagrees with `window_size`, or this contig's source legends disagree with the
            stored legends.
    """
    assert_path_is_readable(in_table_pairs_tsv)
    assert_path_is_readable(reference_fasta)
    reference_fai: Path = reference_fasta.with_suffix(".fai")
    assert_path_is_readable(reference_fai)
    assert_directory_exists(tmp_dir)

    if spark_driver_memory_gb < 1:
        raise ValueError(
            f"Spark driver memory must be at least 1GB. Saw {spark_driver_memory_gb}GB."
        )
    if spark_executor_memory_gb < 1:
        raise ValueError(
            f"Spark executor memory must be at least 1GB. Saw {spark_executor_memory_gb}GB."
        )

    out_duckdb_file: Path = Path(f"{str(output_base)}.haplotypes_gnomad_merge.index.duckdb")
    assert_path_is_readable(out_duckdb_file)

    table_pairs: list[TablePair] = list(TablePair.read(in_table_pairs_tsv))
    table_pair: TablePair | None = next((tp for tp in table_pairs if tp.contig == contig), None)
    if table_pair is None:
        raise ValueError(f"Contig {contig} not found in {in_table_pairs_tsv}.")

    if table_pair.haplotype_table_path is not None:
        assert_directory_exists(table_pair.haplotype_table_path)
    assert_directory_exists(table_pair.sites_table_path)

    # Same Spark-memory env setup + Hail init as the old tool. Skip if a context already exists
    # (e.g. a shared test-session context) so this stays idempotent within a process.
    if Env._hc is None:
        os.environ["PYSPARK_SUBMIT_ARGS"] = (
            f"--driver-memory {spark_driver_memory_gb}g "
            f"--executor-memory {spark_executor_memory_gb}g "
            "pyspark-shell"
        )
        hl.init(tmp_dir=str(tmp_dir))

    with duckdb.connect(str(out_duckdb_file)) as conn:
        stored_window_size: int = read_window_size(conn)
        if stored_window_size != window_size:
            raise ValueError(
                f"Stored window_size {stored_window_size} does not match requested window_size "
                f"{window_size}."
            )

        joint_pops_legend, remaps = _resolve_legends_and_remaps(conn, table_pair)

        # In production each contig runs in a fresh JVM, so the reference has no sequence yet.
        # Guard against re-adding when a context is reused within a process (e.g. the shared
        # pytest-session Hail context), since add_sequence raises if a sequence is already set.
        reference = hl.get_reference(reference_genome)
        if not reference.has_sequence():
            reference.add_sequence(str(reference_fasta), str(reference_fai))

        sequence_id_offset: int = sequences_row_count(conn)
        contig_seq_ht = build_contig_sequences_table(
            table_pair=table_pair,
            window_size=window_size,
            version=version,
            sequence_id_offset=sequence_id_offset,
            hgdp_to_joint=remaps.hgdp_to_joint,
            gnomad_to_joint=remaps.gnomad_to_joint,
            hgdp_at_joint=remaps.hgdp_at_joint,
            gnomad_at_joint=remaps.gnomad_at_joint,
        )

        per_contig_tsv_dir: Path = output_base.parent if retain_per_contig_tsvs else tmp_dir
        contig_tsv: Path = (
            per_contig_tsv_dir / f"{output_base.name}.haplotypes_gnomad_merge.{contig}.tsv.bgz"
        )
        export_sequences_table_to_tsv(
            ht=contig_seq_ht,
            out_file=contig_tsv,
            joint_pops_legend=joint_pops_legend,
        )

        contig_rows: int = _stream_tsv_into_sequences(
            conn,
            tsv=contig_tsv,
            joint_pops_legend=joint_pops_legend,
            chunk_size=polars_chunk_size,
        )

        if not retain_per_contig_tsvs and contig_tsv.exists():
            contig_tsv.unlink()

    logger.info(
        f"Appended {contig_rows} rows for contig {contig} "
        f"(starting at sequence_id offset {sequence_id_offset})."
    )
