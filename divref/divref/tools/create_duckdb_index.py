"""Tool to build DivRef DuckDB index from haplotype Hail tables."""

import json
import logging
import os
from collections.abc import Iterator
from pathlib import Path

import duckdb
import hail as hl
import polars
from fgmetric import Metric
from fgpyo.io import assert_directory_exists
from fgpyo.io import assert_path_is_readable
from fgpyo.io import assert_path_is_writable

from divref import defaults
from divref.haplotype import get_haplo_sequence
from divref.haplotype import haplo_coordinates

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


def create_duckdb_index(  # noqa: C901
    *,
    in_table_pairs_tsv: Path,
    reference_fasta: Path,
    window_size: int,
    output_base: Path,
    version: str,
    reference_genome: str = defaults.REFERENCE_GENOME,
    tmp_dir: Path = Path("/tmp"),
    polars_chunk_size: int = 100_000,
    retain_per_contig_tsvs: bool = False,
    force: bool = False,
    spark_driver_memory_gb: int = 1,
    spark_executor_memory_gb: int = 1,
) -> None:
    """
    Convert per-chr haplotype and gnomAD variant Hail tables into a searchable DuckDB index.

    Streams sequences one contig at a time, then sub-divides each contig at the polars read step:
    each per-contig TSV is read back in batches of `polars_chunk_size` rows, and each batch is
    appended to the DuckDB `sequences` table. The first batch creates the table; subsequent
    batches `INSERT INTO` it. Sequence IDs are assigned with a running offset so they remain
    unique across contigs and batches.

    The two input sources may carry different population legends (e.g. when the gnomAD variant
    track is drawn from a release with more populations than the HGDP haplotype track). Both per-
    source pop-index spaces are remapped into a single `joint_pops_legend` (the union of the two,
    with gnomAD pops first) before union, so the `max_pop` and `all_pop_freqs[*].pop` integers in
    the index are comparable across sources. The DuckDB output stores three legend tables:
    `hgdp_haplotype_pops_legend`, `gnomad_variant_pops_legend`, and `joint_pops_legend`.

    Args:
        in_table_pairs_tsv: Path to a TSV file with fields 'contig', 'haplotype_table_path', and
            'sites_table_path'. The 'haplotype_table_path' cell may be empty for contigs that
            contribute only gnomAD single variants; 'sites_table_path' must be populated on every
            row. A run with zero haplotype rows is supported (the joint pop legend then collapses
            to the gnomAD legend, and the written `hgdp_haplotype_pops_legend` table is `[]`).
        reference_fasta: Path to the indexed reference FASTA for sequence extraction.
        window_size: Flanking reference context size around each haplotype/variant; used when
            constructing sequence strings and stored in the index.
        output_base: Base path for output. Writes
            `{output_base}.haplotypes_gnomad_merge.index.duckdb` and, when
            `retain_per_contig_tsvs` is True, one
            `{output_base}.haplotypes_gnomad_merge.{contig}.tsv.bgz` per contig.
        version: Version identifier embedded in sequence IDs (e.g. "1.0").
        reference_genome: Reference genome to use. Defaults to "GRCh38".
        tmp_dir: Temporary directory for Hail checkpoints and (when not retained) per-contig
            intermediate TSVs.
        polars_chunk_size: Maximum number of rows per polars read batch. Bounds the in-process
            DataFrame size when streaming each per-contig TSV into DuckDB.
        retain_per_contig_tsvs: If True, write per-contig TSVs alongside the duckdb output rather
            than into `tmp_dir`.
        force: If True, overwrite an existing duckdb output. Otherwise raise FileExistsError.
        spark_driver_memory_gb: Memory in GB to allocate to the Spark driver.
        spark_executor_memory_gb: Memory in GB to allocate to the Spark executor.
    """
    assert_path_is_readable(in_table_pairs_tsv)
    assert_path_is_readable(reference_fasta)
    assert_path_is_readable(reference_fasta.with_suffix(".fai"))
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

    # determine per-contig TSV paths and fail fast on writability when retained
    per_contig_tsv_dir: Path = output_base.parent if retain_per_contig_tsvs else tmp_dir
    per_contig_tsvs: dict[str, Path] = {
        tp.contig: per_contig_tsv_dir
        / f"{output_base.name}.haplotypes_gnomad_merge.{tp.contig}.tsv.bgz"
        for tp in table_pairs
    }
    if retain_per_contig_tsvs:
        for tsv_path in per_contig_tsvs.values():
            assert_path_is_writable(tsv_path)

    os.environ["PYSPARK_SUBMIT_ARGS"] = (
        f"--driver-memory {spark_driver_memory_gb}g "
        f"--executor-memory {spark_executor_memory_gb}g "
        "pyspark-shell"
    )
    hl.init(tmp_dir=str(tmp_dir))

    # Bootstrap pop legends from the first available source on each side. The haplotype side may
    # be absent on some rows (or on every row, e.g. a sites-only run).
    first_with_hap: TablePair | None = next(
        (tp for tp in table_pairs if tp.haplotype_table_path is not None), None
    )
    hgdp_pops_legend: list[str] = (
        hl.read_table(str(first_with_hap.haplotype_table_path)).pops.collect()[0]
        if first_with_hap is not None
        else []
    )
    gnomad_pops_legend: list[str] = hl.read_table(
        str(table_pairs[0].sites_table_path)
    ).pops.collect()[0]
    # All pairs must share the same pops legends so a single remap into the joint legend is valid
    # for every contig; otherwise the exported gnomAD_AF_* columns would be misaligned. Rows
    # without a haplotype table are skipped on the haplotype-side check.
    for tp in table_pairs[1:]:
        tp_gnomad_pops: list[str] = hl.read_table(str(tp.sites_table_path)).pops.collect()[0]
        if tp_gnomad_pops != gnomad_pops_legend:
            raise ValueError(
                f"gnomAD pops legend mismatch for contig {tp.contig}: "
                f"{tp_gnomad_pops} vs {gnomad_pops_legend}."
            )
    for tp in table_pairs:
        if tp is first_with_hap or tp.haplotype_table_path is None:
            continue
        tp_hgdp_pops: list[str] = hl.read_table(str(tp.haplotype_table_path)).pops.collect()[0]
        if tp_hgdp_pops != hgdp_pops_legend:
            raise ValueError(
                f"HGDP haplotype pops legend mismatch for contig {tp.contig}: "
                f"{tp_hgdp_pops} vs {hgdp_pops_legend}."
            )
    # Joint legend: gnomAD pops in their original order, then any HGDP-only pops appended.
    joint_pops_legend: list[str] = list(gnomad_pops_legend) + [
        p for p in hgdp_pops_legend if p not in gnomad_pops_legend
    ]
    hgdp_to_joint: list[int] = [joint_pops_legend.index(p) for p in hgdp_pops_legend]
    gnomad_to_joint: list[int] = [joint_pops_legend.index(p) for p in gnomad_pops_legend]
    # Inverse remaps for reshuffling each source's `gnomad_freqs` inner array into joint order.
    hgdp_at_joint: list[int] = [
        hgdp_pops_legend.index(p) if p in hgdp_pops_legend else -1 for p in joint_pops_legend
    ]
    gnomad_at_joint: list[int] = [
        gnomad_pops_legend.index(p) if p in gnomad_pops_legend else -1 for p in joint_pops_legend
    ]

    hl.get_reference(reference_genome).add_sequence(str(reference_fasta))

    with duckdb.connect(str(out_duckdb_file)) as conn:
        sequence_id_offset: int = 0
        created_table: bool = False
        for table_pair in table_pairs:
            contig_seq_ht = build_contig_sequences_table(
                table_pair=table_pair,
                window_size=window_size,
                version=version,
                sequence_id_offset=sequence_id_offset,
                hgdp_to_joint=hgdp_to_joint,
                gnomad_to_joint=gnomad_to_joint,
                hgdp_at_joint=hgdp_at_joint,
                gnomad_at_joint=gnomad_at_joint,
            )
            contig_tsv: Path = per_contig_tsvs[table_pair.contig]
            export_sequences_table_to_tsv(
                ht=contig_seq_ht,
                out_file=contig_tsv,
                joint_pops_legend=joint_pops_legend,
            )

            contig_rows: int = 0
            for df in iter_dataframe_chunks(
                tsv=contig_tsv,
                joint_pops_legend=joint_pops_legend,
                chunk_size=polars_chunk_size,
            ):
                if not created_table:
                    conn.execute("CREATE TABLE sequences AS SELECT * FROM df")
                    created_table = True
                else:
                    conn.execute("INSERT INTO sequences SELECT * FROM df")
                contig_rows += df.height
                sequence_id_offset += df.height

            if not retain_per_contig_tsvs and contig_tsv.exists():
                contig_tsv.unlink()

            logger.info(
                f"Appended {contig_rows} rows for contig {table_pair.contig} "
                f"(running total: {sequence_id_offset})"
            )

        if not created_table:
            raise ValueError("No sequence rows were produced; cannot create DuckDB index.")

        conn.execute("CREATE INDEX idx_sequence_id ON sequences(sequence_id)")
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
        Tuple of (checkpointed Hail table, population legend list).
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


def iter_dataframe_chunks(
    *,
    tsv: Path,
    joint_pops_legend: list[str],
    chunk_size: int,
) -> Iterator[polars.DataFrame]:
    """
    Yield polars DataFrames of up to `chunk_size` rows from a sequences TSV.

    The `sequence_id` and `gnomAD_AF_*` columns are explicitly typed as strings so that
    schema inference cannot misread comma-delimited per-variant AFs as floats.

    Args:
        tsv: Path to the sequences TSV (bgz-compressed).
        joint_pops_legend: Ordered list of population codes used to name `gnomAD_AF_{pop}`
            columns; must match what `export_sequences_table_to_tsv` wrote.
        chunk_size: Maximum rows per yielded DataFrame.

    Yields:
        Polars DataFrame batches read from `tsv`.
    """
    schema_overrides: dict[str, type[polars.DataType]] = {
        "sequence_id": polars.String,
        **{f"gnomAD_AF_{pop}": polars.String for pop in joint_pops_legend},
    }
    # Hail's TSV export emits "NA" for missing scalar fields; "null" is included for
    # robustness against other writers.
    lf = polars.scan_csv(
        tsv,
        separator="\t",
        schema_overrides=schema_overrides,
        null_values=["NA", "null"],
    )
    for df in lf.collect_batches(chunk_size=chunk_size):
        if df.height > 0:
            yield df
