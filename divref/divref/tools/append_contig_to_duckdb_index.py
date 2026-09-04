"""Tool to append a single contig's sequences to an existing DivRef DuckDB index."""

import logging
import os
from dataclasses import dataclass
from pathlib import Path

import duckdb
import hail as hl
from fgpyo.io import assert_directory_exists
from fgpyo.io import assert_path_is_readable
from hail.context import Env

from divref import defaults
from divref.duckdb_index import contig_already_appended
from divref.duckdb_index import read_legend
from divref.duckdb_index import read_window_size
from divref.duckdb_index import sequences_row_count
from divref.duckdb_index import sequences_tsv_columns
from divref.duckdb_index import stream_sequences_tsv_into_duckdb
from divref.gnomad_index_source import TablePair
from divref.gnomad_index_source import read_and_validate_pops_legends
from divref.haplotype import get_haplo_sequence
from divref.haplotype import haplo_coordinates

logger = logging.getLogger(__name__)


def to_joint(source_pops: list[str], joint_pops: list[str]) -> list[int]:
    """For each index i in source_pops, its index in joint_pops."""
    return [joint_pops.index(p) for p in source_pops]


def at_joint(source_pops: list[str], joint_pops: list[str]) -> list[int]:
    """For each index j in joint_pops, the source_pops index or -1 if absent."""
    return [source_pops.index(p) if p in source_pops else -1 for p in joint_pops]


def build_hgdp_haplotype_table_entries(
    haplotypes_table_path: Path,
    hgdp_to_joint: list[int],
    hgdp_at_joint: list[int],
) -> hl.Table:
    """
    Build HGDP_haplotype entries for the "sequences" table.

    Reads the haplotype table and annotates with source and population frequencies. The haplotypes
    are already at the desired granularity from `compute_haplotypes`; this function performs no
    algorithmic transformation. `max_pop` and `all_pop_freqs[*].pop` integer indices are remapped
    from the haplotype table's native pop ordering into the joint pop legend, and each row's
    `gnomad_freqs` inner array is reshuffled and padded to the joint legend's length so it indexes
    positionally by the joint legend (missing pops become a missing struct).

    Args:
        haplotypes_table_path: Path to the computed haplotypes Hail table.
        hgdp_to_joint: For each index `i` in the haplotype table's pop legend, the corresponding
            index in the joint pop legend.
        hgdp_at_joint: For each index `j` in the joint pop legend, the corresponding index in the
            haplotype table's pop legend, or `-1` if that pop is not present on the haplotype side.

    Returns:
        Hail table with source and per-pop frequency annotations.
    """
    ht = hl.read_table(str(haplotypes_table_path)).key_by().drop("haplotype")
    hgdp_remap = hl.literal(hgdp_to_joint)
    hgdp_at_joint_lit = hl.literal(hgdp_at_joint)
    inner_struct_type = ht.gnomad_freqs.dtype.element_type.element_type
    ht = ht.annotate(
        source="HGDP_haplotype",
        max_pop=hgdp_remap[ht.max_pop],
        all_pop_freqs=ht.all_pop_freqs.map(
            lambda x: hl.struct(
                pop=hgdp_remap[x.pop],
                empirical_AC=x.empirical_AC,
                empirical_AF=x.empirical_AF,
                fraction_phased=x.fraction_phased,
                estimated_gnomad_AF=x.estimated_gnomad_AF,
            )
        ),
        gnomad_freqs=ht.gnomad_freqs.map(
            lambda inner: hl.range(hl.len(hgdp_at_joint_lit)).map(
                lambda j: hl.if_else(
                    hgdp_at_joint_lit[j] >= 0,
                    inner[hgdp_at_joint_lit[j]],
                    hl.missing(inner_struct_type),
                )
            )
        ),
    )
    return ht


def build_gnomad_variant_table_entries(
    sites_table_path: Path,
    gnomad_to_joint: list[int],
    gnomad_at_joint: list[int],
) -> hl.Table:
    """
    Build gnomAD_variant entries for the "sequences" table.

    Reads the gnomAD table and annotates entries to match the HGDP_haplotype entries. `max_pop`
    and `all_pop_freqs[*].pop` integer indices are remapped from the gnomAD-source legend into the
    joint pop legend, and the per-variant `gnomad_freqs` inner array is reshuffled and padded to
    the joint legend's length so it indexes positionally by the joint legend (missing pops become
    a missing struct).

    Args:
        sites_table_path: Path to the gnomAD variant annotations Hail table.
        gnomad_to_joint: For each index `i` in the gnomAD sites table's pop legend, the
            corresponding index in the joint pop legend.
        gnomad_at_joint: For each index `j` in the joint pop legend, the corresponding index in
            the gnomAD sites table's pop legend, or `-1` if that pop is not present on the gnomAD
            side.

    Returns:
        Hail table of gnomAD single-variant entries annotated to match the HGDP_haplotype schema.
    """
    va = hl.read_table(str(sites_table_path))
    va = va.rename({"pop_freqs": "gnomad_freqs"})
    va = va.key_by()
    argmax_pop = hl.argmax(va.gnomad_freqs.map(lambda x: x.AF))
    gnomad_remap = hl.literal(gnomad_to_joint)
    gnomad_at_joint_lit = hl.literal(gnomad_at_joint)
    inner_struct_type = va.gnomad_freqs.dtype.element_type
    gnomad_freqs_joint = hl.range(hl.len(gnomad_at_joint_lit)).map(
        lambda j: hl.if_else(
            gnomad_at_joint_lit[j] >= 0,
            va.gnomad_freqs[gnomad_at_joint_lit[j]],
            hl.missing(inner_struct_type),
        )
    )
    va = va.select(
        max_pop=gnomad_remap[argmax_pop],
        max_empirical_AF=va.gnomad_freqs[argmax_pop].AF,
        fraction_phased=1.0,
        estimated_gnomad_AF=va.gnomad_freqs[argmax_pop].AF,
        max_empirical_AC=va.gnomad_freqs[argmax_pop].AC,
        all_pop_freqs=hl.range(hl.len(va.gnomad_freqs)).map(
            lambda i: hl.struct(
                pop=gnomad_remap[i],
                empirical_AC=va.gnomad_freqs[i].AC,
                empirical_AF=va.gnomad_freqs[i].AF,
                # For single-variant rows, the haplotype is a single allele, so phasing is
                # trivially complete and the "estimated" gnomAD AF is just the gnomAD AF in
                # the pop — matching the scalar `fraction_phased=1.0` and
                # `estimated_gnomad_AF=va.gnomad_freqs[argmax_pop].AF` convention above.
                fraction_phased=1.0,
                estimated_gnomad_AF=va.gnomad_freqs[i].AF,
            )
        ),
        source="gnomAD_variant",
        variants=[hl.struct(locus=va.locus, alleles=va.alleles)],
        gnomad_freqs=[gnomad_freqs_joint],
    )
    return va


def build_contig_sequences_table(
    *,
    table_pair: TablePair,
    window_size: int,
    version: str,
    sequence_id_offset: int,
    hgdp_to_joint: list[int],
    gnomad_to_joint: list[int],
    hgdp_at_joint: list[int],
    gnomad_at_joint: list[int],
) -> hl.Table:
    """
    Build the per-contig sequences hail table with sequences, coordinates, and IDs.

    Reads the HGDP haplotype + gnomAD sites tables for one contig, unions them, sorts by genomic
    position, and applies the same per-row annotations as the cross-contig table. Sequence IDs are
    offset by `sequence_id_offset` so they remain unique across contigs. When
    `table_pair.haplotype_table_path` is `None`, the haplotype side is skipped and only the gnomAD
    sites table contributes rows for this contig.

    Args:
        table_pair: Per-contig pair of haplotype + gnomAD sites table paths. The haplotype side
            may be `None`.
        window_size: Flanking reference context size for sequence construction.
        version: Version identifier for sequence IDs.
        sequence_id_offset: Number of rows already written for prior contigs; added to this
            contig's local index to produce a globally unique sequence ID.
        hgdp_to_joint: Remap from the haplotype-source pop legend into the joint pop legend.
        gnomad_to_joint: Remap from the gnomAD-source pop legend into the joint pop legend.
        hgdp_at_joint: Inverse remap: for each joint index, the haplotype-source index or -1.
        gnomad_at_joint: Inverse remap: for each joint index, the gnomAD-source index or -1.

    Returns:
        Hail table with sequences, coordinates, and variant strings annotated.
    """
    gnomad_variants_ht: hl.Table = build_gnomad_variant_table_entries(
        sites_table_path=table_pair.sites_table_path,
        gnomad_to_joint=gnomad_to_joint,
        gnomad_at_joint=gnomad_at_joint,
    )
    seq_ht: hl.Table
    if table_pair.haplotype_table_path is None:
        seq_ht = gnomad_variants_ht
    else:
        hgdp_haplotypes_ht: hl.Table = build_hgdp_haplotype_table_entries(
            haplotypes_table_path=table_pair.haplotype_table_path,
            hgdp_to_joint=hgdp_to_joint,
            hgdp_at_joint=hgdp_at_joint,
        )
        seq_ht = hgdp_haplotypes_ht.union(gnomad_variants_ht, unify=True)

    seq_ht = seq_ht.rename({
        "max_empirical_AF": "popmax_empirical_AF",
        "max_empirical_AC": "popmax_empirical_AC",
    })

    seq_ht = seq_ht.annotate(variant_strs=seq_ht.variants.map(lambda x: hl.variant_str(x)))
    # Total ordering so the row index — and therefore the assigned sequence_id — is determined by
    # the data, not by Spark partitioning. Primary key is the genomic start; ties (e.g. a single
    # variant co-located with a haplotype's leftmost variant) break on source then the joined
    # variant string, which together uniquely identify a row.
    seq_ht = seq_ht.annotate(
        min_pos=hl.sorted(seq_ht.variants, key=lambda v: v.locus.position)[0].locus.position,
        seq_sort_key=hl.delimit(seq_ht.variant_strs, ","),
    )
    seq_ht = seq_ht.order_by(seq_ht.min_pos, seq_ht.source, seq_ht.seq_sort_key).drop(
        "min_pos", "seq_sort_key"
    )
    seq_ht = seq_ht.add_index()
    coords = haplo_coordinates(window_size, seq_ht.variants)
    seq_ht = seq_ht.annotate(
        sequence=get_haplo_sequence(window_size, seq_ht.variants),
        contig=seq_ht.variants[0].locus.contig,
        # `start`/`end` are 0-based half-open (from `haplo_coordinates`; variant loci are 1-based).
        start=coords.start,
        end=coords.end,
    )
    seq_ht = seq_ht.annotate(
        sequence_length=hl.len(seq_ht.sequence),
        sequence_id=hl.str(f"DR-{version}-") + hl.str(seq_ht.idx + sequence_id_offset),
        n_variants=hl.len(seq_ht.variants),
    ).drop("idx")

    return seq_ht


def export_sequences_table_to_tsv(
    ht: hl.Table,
    out_file: Path,
    joint_pops_legend: list[str],
) -> None:
    """
    Export the sequences Hail table to a single bgz-compressed TSV.

    One `gnomAD_AF_{pop}` column is emitted per pop in `joint_pops_legend`, in order. Each row's
    `gnomad_freqs` inner array is already reshuffled to the joint legend at source-table
    construction time (with missing-padding for pops absent from a source), so a uniform
    positional lookup is safe regardless of which source the row came from.

    Four further per-pop columns are emitted per joint legend entry: `empirical_AC_{pop}`,
    `empirical_AF_{pop}`, `fraction_phased_{pop}`, `estimated_gnomAD_haplotype_AF_{pop}`. Values
    come from `all_pop_freqs` by joint-pop-index dict lookup; pops absent from a row's source
    are emitted as missing.

    The scalar columns `popmax_fraction_phased` and `popmax_estimated_gnomad_AF` are renamed
    from `fraction_phased` / `estimated_gnomad_AF` to make the max-pop semantic explicit
    alongside `popmax_empirical_AF` / `popmax_empirical_AC`.

    The emitted column order is driven by `divref.duckdb_index.sequences_tsv_columns`.

    Args:
        ht: Annotated haplotype/variant table with sequences and variant strings.
        out_file: Path for the output TSV file.
        joint_pops_legend: Ordered list of all population codes across both input sources; used to
            resolve `max_pop` integer indices to labels and to name `gnomAD_AF_{pop}` columns.
    """
    # Per-joint-pop dict lookup over `all_pop_freqs`. After build_*_table_entries the entries'
    # `pop` field is already in the joint legend index space; pops absent from this row's source
    # have no entry, and `.get(i, missing_struct)` returns missing for those.
    pop_freq_value_type = ht.all_pop_freqs.dtype.element_type
    ht = ht.annotate(
        _pop_lookup=hl.dict(ht.all_pop_freqs.map(lambda x: (x.pop, x))),
    )
    missing_pop_struct = hl.missing(pop_freq_value_type)

    expr_by_name: dict[str, hl.Expression] = {
        "sequence": ht.sequence,
        "sequence_length": ht.sequence_length,
        "sequence_id": ht.sequence_id,
        "n_variants": ht.n_variants,
        "contig": ht.contig,
        "start": ht.start,
        "end": ht.end,
        "popmax_empirical_AF": ht.popmax_empirical_AF,
        "popmax_empirical_AC": ht.popmax_empirical_AC,
        "source": ht.source,
        "popmax_estimated_gnomad_AF": ht.estimated_gnomad_AF,
        "popmax_fraction_phased": ht.fraction_phased,
        "max_pop": hl.literal(joint_pops_legend)[ht.max_pop],
        "variants": hl.delimit(ht.variant_strs, ","),
    }
    for i, pop in enumerate(joint_pops_legend):
        entry = ht._pop_lookup.get(i, missing_pop_struct)
        # Substitute "NA" per element for variants where this pop's AF is missing (e.g. an
        # HGDP_haplotype row at a joint pop that isn't in the HGDP source legend). Without the
        # substitution Hail collapses the whole all-missing array to a single missing token,
        # which polars then reads as a SQL NULL and trips the downstream Haplotype model.
        # Always emitting a comma-delimited string of length `n_variants` keeps the column
        # shape consistent regardless of which source emitted the row.
        expr_by_name[f"gnomAD_AF_{pop}"] = hl.delimit(
            ht.gnomad_freqs.map(
                lambda x, _i=i: hl.if_else(
                    hl.is_defined(x[_i].AF),
                    hl.format("%.5f", x[_i].AF),
                    hl.literal("NA"),
                )
            ),
            ",",
        )
        expr_by_name[f"empirical_AC_{pop}"] = entry.empirical_AC
        expr_by_name[f"empirical_AF_{pop}"] = entry.empirical_AF
        expr_by_name[f"fraction_phased_{pop}"] = entry.fraction_phased
        expr_by_name[f"estimated_gnomAD_haplotype_AF_{pop}"] = entry.estimated_gnomad_AF

    cols = sequences_tsv_columns(
        joint_pops_legend, af_prefix="gnomAD", popmax_estimated_col="popmax_estimated_gnomad_AF"
    )
    assert sorted(expr_by_name) == sorted(cols), (
        "sequences_tsv_columns and the built expressions have drifted apart."
    )
    ht.select(**{name: expr_by_name[name] for name in cols}).export(str(out_file))


# eq=False so the frozen dataclass keeps a (default, identity-based) __hash__ despite its list
# fields; without it, frozen+eq would synthesize a __hash__ that raises TypeError on unhashable
# lists. The struct is only ever field-accessed, never compared or hashed by value.
@dataclass(frozen=True, kw_only=True, eq=False)
class _RemapArrays:
    """
    The four pop-legend remap arrays used to build a contig's sequences table.

    Attributes:
        hgdp_to_joint: For each haplotype-source pop index, its index in the joint legend.
        gnomad_to_joint: For each gnomAD-source pop index, its index in the joint legend.
        hgdp_at_joint: For each joint pop index, the haplotype-source index or -1 if absent.
        gnomad_at_joint: For each joint pop index, the gnomAD-source index or -1 if absent.
    """

    hgdp_to_joint: list[int]
    gnomad_to_joint: list[int]
    hgdp_at_joint: list[int]
    gnomad_at_joint: list[int]


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


def _export_and_stream_contig(
    conn: duckdb.DuckDBPyConnection,
    *,
    seq_ht: hl.Table,
    contig_tsv: Path,
    joint_pops_legend: list[str],
    chunk_size: int,
    retain_tsv: bool,
) -> int:
    """
    Export a contig's sequences to a TSV and stream it into `sequences` within one transaction.

    The per-contig TSV (which can be ~GB) is deleted on any failure unless `retain_tsv`, and the
    DuckDB writes (the `sequences` CREATE on the first contig plus all chunk INSERTs) run inside a
    single transaction, so a crash mid-stream rolls back to leave the index exactly as it was
    rather than a half-written contig.

    Args:
        conn: Open connection to the DuckDB index.
        seq_ht: The contig's annotated sequences Hail table.
        contig_tsv: Path the per-contig TSV is written to.
        joint_pops_legend: Ordered joint pop legend used to type the streamed columns.
        chunk_size: Maximum number of rows per polars read batch.
        retain_tsv: If True, leave the per-contig TSV in place after streaming.

    Returns:
        The number of rows appended for this contig.
    """
    try:
        export_sequences_table_to_tsv(
            ht=seq_ht,
            out_file=contig_tsv,
            joint_pops_legend=joint_pops_legend,
        )
        contig_rows = stream_sequences_tsv_into_duckdb(
            conn,
            tsv=contig_tsv,
            joint_pops_legend=joint_pops_legend,
            chunk_size=chunk_size,
        )
    finally:
        if not retain_tsv and contig_tsv.exists():
            contig_tsv.unlink()
    return contig_rows


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

    Each contig is written in a single transaction, so a crash mid-stream rolls back and leaves the
    index unchanged. Appends are expected to run once per contig, in order; re-appending a contig
    whose rows are already present is rejected (rebuild via `init_duckdb_index --force` instead),
    so the tool is safe to invoke standalone and not only inside the ordered workflow loop.

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
            size disagrees with `window_size`, this contig's source legends disagree with the
            stored legends, or the contig already has rows in the index.
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

    out_duckdb_file: Path = Path(f"{output_base}.haplotypes_gnomad_merge.index.duckdb")
    # readable (not writable): enforces the index file already exists, i.e. init_duckdb_index ran.
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

        # Appends run once per contig, in canonical order, after init_duckdb_index. Re-appending
        # onto an index that already holds this contig's rows would duplicate them with fresh
        # sequence IDs, so fail loudly instead; to rebuild, re-run init_duckdb_index with --force
        # (how the workflow redoes the whole index). This keeps the tool safe to invoke standalone.
        if contig_already_appended(conn, contig):
            raise ValueError(
                f"Contig {contig} already has rows in {out_duckdb_file}. Appends run once per "
                f"contig; re-run init_duckdb_index with --force to rebuild from scratch."
            )

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
        contig_rows: int = _export_and_stream_contig(
            conn,
            seq_ht=contig_seq_ht,
            contig_tsv=contig_tsv,
            joint_pops_legend=joint_pops_legend,
            chunk_size=polars_chunk_size,
            retain_tsv=retain_per_contig_tsvs,
        )

    logger.info(
        f"Appended {contig_rows} rows for contig {contig} "
        f"(starting at sequence_id offset {sequence_id_offset})."
    )
