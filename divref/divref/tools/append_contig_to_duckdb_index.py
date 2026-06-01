"""Tool to append a single contig's sequences to an existing DivRef DuckDB index."""

import json
import logging
import os
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

import duckdb
import hail as hl
import polars
from fgpyo.io import assert_directory_exists
from fgpyo.io import assert_path_is_readable
from hail.context import Env

from divref import defaults
from divref.duckdb_index import TablePair
from divref.duckdb_index import read_and_validate_pops_legends
from divref.haplotype import get_haplo_sequence
from divref.haplotype import haplo_coordinates

logger = logging.getLogger(__name__)


def to_joint(source_pops: list[str], joint_pops: list[str]) -> list[int]:
    """For each index i in source_pops, its index in joint_pops."""
    return [joint_pops.index(p) for p in source_pops]


def at_joint(source_pops: list[str], joint_pops: list[str]) -> list[int]:
    """For each index j in joint_pops, the source_pops index or -1 if absent."""
    return [source_pops.index(p) if p in source_pops else -1 for p in joint_pops]


_LEGEND_TABLES = frozenset({
    "joint_pops_legend",
    "gnomad_variant_pops_legend",
    "hgdp_haplotype_pops_legend",
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
    exists = conn.execute(
        "SELECT 1 FROM information_schema.tables WHERE table_name = 'sequences'"
    ).fetchone()
    if exists is None:
        return 0
    row = conn.execute("SELECT COUNT(*) FROM sequences").fetchone()
    if row is None:
        raise ValueError("COUNT(*) on sequences returned no row.")
    return int(row[0])


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
    count_orig: int = va.count()
    logger.info(f"Variant table {sites_table_path} contains {count_orig} variants.")

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

    seq_ht = seq_ht.annotate(
        min_pos=hl.sorted(seq_ht.variants, key=lambda v: v.locus.position)[0].locus.position
    )
    seq_ht = seq_ht.order_by(seq_ht.min_pos).drop("min_pos")
    seq_ht = seq_ht.add_index()
    coords = haplo_coordinates(window_size, seq_ht.variants)
    seq_ht = seq_ht.annotate(
        sequence=get_haplo_sequence(window_size, seq_ht.variants),
        contig=seq_ht.variants[0].locus.contig,
        start=coords.start,
        end=coords.end,
    )
    seq_ht = seq_ht.annotate(variant_strs=seq_ht.variants.map(lambda x: hl.variant_str(x)))
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
    per_pop_columns: dict[str, hl.Expression] = {}
    for i, pop in enumerate(joint_pops_legend):
        entry = ht._pop_lookup.get(i, missing_pop_struct)
        per_pop_columns[f"empirical_AC_{pop}"] = entry.empirical_AC
        per_pop_columns[f"empirical_AF_{pop}"] = entry.empirical_AF
        per_pop_columns[f"fraction_phased_{pop}"] = entry.fraction_phased
        per_pop_columns[f"estimated_gnomAD_haplotype_AF_{pop}"] = entry.estimated_gnomad_AF

    ht.select(
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
        popmax_estimated_gnomad_AF=ht.estimated_gnomad_AF,
        popmax_fraction_phased=ht.fraction_phased,
        max_pop=hl.literal(joint_pops_legend)[ht.max_pop],
        variants=hl.delimit(ht.variant_strs, ","),
        # Substitute "NA" per element for variants where this pop's AF is missing (e.g. an
        # HGDP_haplotype row at a joint pop that isn't in the HGDP source legend). Without the
        # substitution Hail collapses the whole all-missing array to a single missing token,
        # which polars then reads as a SQL NULL and trips the downstream Haplotype model.
        # Always emitting a comma-delimited string of length `n_variants` keeps the column
        # shape consistent regardless of which source emitted the row.
        **{
            f"gnomAD_AF_{pop}": hl.delimit(
                ht.gnomad_freqs.map(
                    lambda x, _i=i: hl.if_else(
                        hl.is_defined(x[_i].AF),
                        hl.format("%.5f", x[_i].AF),
                        hl.literal("NA"),
                    )
                ),
                ",",
            )
            for i, pop in enumerate(joint_pops_legend)
        },
        **per_pop_columns,
    ).export(str(out_file))


def _scan_sequences_tsv(tsv: Path, joint_pops_legend: list[str]) -> polars.LazyFrame:
    """
    Build the typed polars LazyFrame for a sequences TSV.

    The `sequence_id` and `gnomAD_AF_*` columns are explicitly typed as strings so that
    schema inference cannot misread comma-delimited per-variant AFs as floats.

    Every numeric column is pinned explicitly as well. Each contig is appended in its own process
    with independent polars inference, and the `sequences` table's schema is fixed by whichever
    contig is appended first (`CREATE TABLE ... AS SELECT`). A per-pop column that is all-`NA` for
    one contig (e.g. an HGDP-only pop on a sites-only contig like chrX) would otherwise infer a
    different dtype than the same column on a contig where it is populated, and the later
    `INSERT INTO` would fail or coerce. Pinning the integer and float columns makes the schema
    contig-order-independent.

    Args:
        tsv: Path to the sequences TSV (bgz-compressed).
        joint_pops_legend: Ordered list of population codes used to name `gnomAD_AF_{pop}`
            columns; must match what `export_sequences_table_to_tsv` wrote.

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
        "popmax_estimated_gnomad_AF": polars.Float64,
        "popmax_fraction_phased": polars.Float64,
        **{f"gnomAD_AF_{pop}": polars.String for pop in joint_pops_legend},
        **{f"empirical_AC_{pop}": polars.Int64 for pop in joint_pops_legend},
        **{f"empirical_AF_{pop}": polars.Float64 for pop in joint_pops_legend},
        **{f"fraction_phased_{pop}": polars.Float64 for pop in joint_pops_legend},
        **{f"estimated_gnomAD_haplotype_AF_{pop}": polars.Float64 for pop in joint_pops_legend},
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
        polars.col(f"gnomAD_AF_{pop}").fill_null("NA").cast(polars.String)
        for pop in joint_pops_legend
    ])
    return lf


def iter_dataframe_chunks(
    *,
    tsv: Path,
    joint_pops_legend: list[str],
    chunk_size: int,
) -> Iterator[polars.DataFrame]:
    """
    Yield non-empty polars DataFrames of up to `chunk_size` rows from a sequences TSV.

    Args:
        tsv: Path to the sequences TSV (bgz-compressed).
        joint_pops_legend: Ordered list of population codes used to name `gnomAD_AF_{pop}`
            columns; must match what `export_sequences_table_to_tsv` wrote.
        chunk_size: Maximum rows per yielded DataFrame.

    Yields:
        Polars DataFrame batches read from `tsv`.
    """
    lf = _scan_sequences_tsv(tsv, joint_pops_legend)
    for df in lf.collect_batches(chunk_size=chunk_size):
        if df.height > 0:
            yield df


@dataclass(frozen=True, kw_only=True)
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


def _stream_tsv_into_sequences(
    conn: duckdb.DuckDBPyConnection,
    *,
    tsv: Path,
    joint_pops_legend: list[str],
    chunk_size: int,
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
        joint_pops_legend: Ordered joint pop legend used to type the `gnomAD_AF_*` columns.
        chunk_size: Maximum number of rows per polars read batch.

    Returns:
        The number of rows appended for this contig.
    """
    table_exists = conn.execute(
        "SELECT 1 FROM information_schema.tables WHERE table_name = 'sequences'"
    ).fetchone()
    if table_exists is None:
        empty_df = _scan_sequences_tsv(tsv, joint_pops_legend).limit(0).collect()  # noqa: F841
        conn.execute("CREATE TABLE sequences AS SELECT * FROM empty_df")

    appended_rows: int = 0
    for df in iter_dataframe_chunks(
        tsv=tsv,
        joint_pops_legend=joint_pops_legend,
        chunk_size=chunk_size,
    ):
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
