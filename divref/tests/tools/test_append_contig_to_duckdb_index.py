"""Tests for the append_contig_to_duckdb_index tool."""

import shutil
from pathlib import Path

import duckdb
import hail as hl
import polars
import pytest

from divref.tools.append_contig_to_duckdb_index import _add_compatibility_flag
from divref.tools.append_contig_to_duckdb_index import _stream_tsv_into_sequences
from divref.tools.append_contig_to_duckdb_index import append_contig_to_duckdb_index
from divref.tools.append_contig_to_duckdb_index import export_sequences_table_to_tsv
from divref.tools.append_contig_to_duckdb_index import iter_dataframe_chunks
from divref.tools.init_duckdb_index import init_duckdb_index


def _write_table_pairs_tsv(path: Path, rows: list[tuple[str, str, str]]) -> Path:
    """Write a table-pairs TSV with a (contig, haplotype_table_path, sites_table_path) per row."""
    lines = ["contig\thaplotype_table_path\tsites_table_path"]
    lines += [f"{contig}\t{hap}\t{sites}" for contig, hap, sites in rows]
    path.write_text("\n".join(lines) + "\n")
    return path


def _table_pairs(datadir: Path, tmp_path: Path) -> Path:
    """Build the two-contig table-pairs TSV (chr1 = full pair, chrX = sites-only)."""
    return _write_table_pairs_tsv(
        tmp_path / "table_pairs.tsv",
        rows=[
            (
                "chr1",
                str(datadir / "chr1_100001_200000_haplotypes.ht"),
                str(datadir / "chr1_100001_200000.gnomad_afs.ht"),
            ),
            (
                "chrX",
                "",
                str(datadir / "chrX_50000000_50025000.gnomad_afs.ht"),
            ),
        ],
    )


def _reference_fasta(datadir: Path) -> Path:
    """Path to the committed masked mini-reference FASTA."""
    return datadir / "test_reference.chr1_chrX.fa.gz"


def _db_path(output_base: Path) -> Path:
    """The DuckDB index path for a given output base."""
    return Path(f"{output_base}.haplotypes_gnomad_merge.index.duckdb")


def test_append_single_contig_creates_sequences(
    hail_context: None,  # noqa: ARG001
    datadir: Path,
    tmp_path: Path,
) -> None:
    """Appending one contig creates `sequences` with rows whose ids start at DR-9.9-0."""
    table_pairs_tsv = _table_pairs(datadir, tmp_path)
    output_base = tmp_path / "idx"
    init_duckdb_index(
        in_table_pairs_tsv=table_pairs_tsv,
        output_base=output_base,
        version="9.9",
        window_size=25,
        force=True,
    )
    append_contig_to_duckdb_index(
        in_table_pairs_tsv=table_pairs_tsv,
        contig="chr1",
        output_base=output_base,
        reference_fasta=_reference_fasta(datadir),
        window_size=25,
        version="9.9",
    )

    with duckdb.connect(str(_db_path(output_base))) as conn:
        count_row = conn.execute("SELECT COUNT(*) FROM sequences").fetchone()
        assert count_row is not None
        assert count_row[0] > 0
        first_row = conn.execute(
            "SELECT sequence_id FROM sequences ORDER BY sequence_id LIMIT 1"
        ).fetchone()
        assert first_row is not None
        assert first_row[0] == "DR-9.9-0"


def test_append_two_contigs_contiguous_ids(
    hail_context: None,  # noqa: ARG001
    datadir: Path,
    tmp_path: Path,
) -> None:
    """Appending chr1 then chrX yields contiguous DR-9.9-0..N-1 ids continuing the offset."""
    table_pairs_tsv = _table_pairs(datadir, tmp_path)
    output_base = tmp_path / "idx"
    reference_fasta = _reference_fasta(datadir)
    init_duckdb_index(
        in_table_pairs_tsv=table_pairs_tsv,
        output_base=output_base,
        version="9.9",
        window_size=25,
        force=True,
    )

    append_contig_to_duckdb_index(
        in_table_pairs_tsv=table_pairs_tsv,
        contig="chr1",
        output_base=output_base,
        reference_fasta=reference_fasta,
        window_size=25,
        version="9.9",
    )
    with duckdb.connect(str(_db_path(output_base))) as conn:
        chr1_count_row = conn.execute("SELECT COUNT(*) FROM sequences").fetchone()
        assert chr1_count_row is not None
        chr1_count: int = chr1_count_row[0]

    append_contig_to_duckdb_index(
        in_table_pairs_tsv=table_pairs_tsv,
        contig="chrX",
        output_base=output_base,
        reference_fasta=reference_fasta,
        window_size=25,
        version="9.9",
    )

    with duckdb.connect(str(_db_path(output_base))) as conn:
        total_row = conn.execute("SELECT COUNT(*) FROM sequences").fetchone()
        assert total_row is not None
        total: int = total_row[0]
        assert total > chr1_count

        ids = [
            row[0]
            for row in conn.execute(
                "SELECT sequence_id FROM sequences ORDER BY CAST(SPLIT_PART("
                "sequence_id, '-', 3) AS INTEGER)"
            ).fetchall()
        ]
        # The ids form the contiguous set DR-9.9-0 .. DR-9.9-(total-1).
        assert ids == [f"DR-9.9-{i}" for i in range(total)]

        # chrX's rows are exactly those beyond chr1's count; its first id continues the offset.
        chrx_first_row = conn.execute(
            "SELECT sequence_id FROM sequences WHERE contig = 'chrX' "
            "ORDER BY CAST(SPLIT_PART(sequence_id, '-', 3) AS INTEGER) LIMIT 1"
        ).fetchone()
        assert chrx_first_row is not None
        assert chrx_first_row[0] == f"DR-9.9-{chr1_count}"


def test_append_legend_mismatch_raises(
    hail_context: None,  # noqa: ARG001
    datadir: Path,
    tmp_path: Path,
) -> None:
    """A stored legend that disagrees with the contig's source legend raises ValueError."""
    table_pairs_tsv = _table_pairs(datadir, tmp_path)
    output_base = tmp_path / "idx"
    init_duckdb_index(
        in_table_pairs_tsv=table_pairs_tsv,
        output_base=output_base,
        version="9.9",
        window_size=25,
        force=True,
    )

    # Corrupt a stored legend in a copy of the DB so the chr1 source legends no longer match.
    db = _db_path(output_base)
    corrupt_output_base = tmp_path / "corrupt"
    corrupt_db = _db_path(corrupt_output_base)
    shutil.copy(db, corrupt_db)
    with duckdb.connect(str(corrupt_db)) as conn:
        conn.execute("DROP TABLE gnomad_variant_pops_legend")
        conn.execute(
            "CREATE TABLE gnomad_variant_pops_legend AS SELECT ? AS pops_legend",
            ['["totally", "wrong"]'],
        )

    with pytest.raises(ValueError):
        append_contig_to_duckdb_index(
            in_table_pairs_tsv=table_pairs_tsv,
            contig="chr1",
            output_base=corrupt_output_base,
            reference_fasta=_reference_fasta(datadir),
            window_size=25,
            version="9.9",
        )


def _make_sequences_ht() -> hl.Table:
    """
    Build a 2-row, 2-pop synthetic sequences HT mirroring the schema fed to the exporter.

    Row 0 has entries for both joint pops in `all_pop_freqs`; row 1 has only pop 0, simulating
    a source that doesn't cover pop 1 — exercises the `dict.get(i, missing_struct)` lookup.
    """
    freq_struct = hl.tstruct(AF=hl.tfloat64, AC=hl.tint32)
    pop_freq_struct = hl.tstruct(
        pop=hl.tint32,
        empirical_AC=hl.tint64,
        empirical_AF=hl.tfloat64,
        fraction_phased=hl.tfloat64,
        estimated_gnomad_AF=hl.tfloat64,
    )
    schema = hl.tstruct(
        sequence=hl.tstr,
        sequence_length=hl.tint32,
        sequence_id=hl.tstr,
        n_variants=hl.tint32,
        contig=hl.tstr,
        start=hl.tint32,
        end=hl.tint32,
        popmax_empirical_AF=hl.tfloat64,
        popmax_empirical_AC=hl.tint64,
        source=hl.tstr,
        estimated_gnomad_AF=hl.tfloat64,
        fraction_phased=hl.tfloat64,
        max_pop=hl.tint32,
        variant_strs=hl.tarray(hl.tstr),
        gnomad_freqs=hl.tarray(hl.tarray(freq_struct)),
        all_pop_freqs=hl.tarray(pop_freq_struct),
    )
    rows = [
        {
            "sequence": "ACGT",
            "sequence_length": 4,
            "sequence_id": "DR-1.0-0",
            "n_variants": 2,
            "contig": "chr1",
            "start": 99,
            "end": 130,
            "popmax_empirical_AF": 0.42,
            "popmax_empirical_AC": 42,
            "source": "HGDP_haplotype",
            "estimated_gnomad_AF": 0.07,
            "fraction_phased": 1.40,
            "max_pop": 0,
            "variant_strs": ["chr1:100:A:T", "chr1:120:C:G"],
            # Multi-variant row where the amr AF is present on variant 0 but missing on
            # variant 1. Exercises the per-element substitution + `hl.delimit` join: the cell
            # must stay a comma-delimited "0.04000,NA" and not collapse to a single "NA".
            "gnomad_freqs": [
                [{"AF": 0.05, "AC": 5}, {"AF": 0.04, "AC": 4}],
                [{"AF": 0.07, "AC": 7}, {"AF": None, "AC": 0}],
            ],
            "all_pop_freqs": [
                {
                    "pop": 0,
                    "empirical_AC": 42,
                    "empirical_AF": 0.42,
                    "fraction_phased": 1.40,
                    "estimated_gnomad_AF": 0.07,
                },
                {
                    "pop": 1,
                    "empirical_AC": 19,
                    "empirical_AF": 0.38,
                    "fraction_phased": 1.52,
                    "estimated_gnomad_AF": 0.05,
                },
            ],
        },
        {
            "sequence": "GGCC",
            "sequence_length": 4,
            "sequence_id": "DR-1.0-1",
            "n_variants": 1,
            "contig": "chr1",
            "start": 199,
            "end": 230,
            "popmax_empirical_AF": 0.10,
            "popmax_empirical_AC": 5,
            "source": "gnomAD_variant",
            "estimated_gnomad_AF": 0.10,
            "fraction_phased": 1.0,
            "max_pop": 0,
            "variant_strs": ["chr1:200:G:C"],
            # Single-variant row with the amr per-pop AF deliberately missing. Exercises the
            # writer's `hl.if_else(hl.is_defined(...), ..., hl.literal("NA"))` branch in
            # export_sequences_table_to_tsv: the resulting `gnomAD_AF_amr` cell degenerates
            # to a bare "NA" token, which polars' `null_values` would otherwise convert to
            # None before iter_dataframe_chunks's `fill_null("NA")` restores it.
            "gnomad_freqs": [
                [{"AF": 0.10, "AC": 10}, {"AF": None, "AC": 0}],
            ],
            "all_pop_freqs": [
                {
                    "pop": 0,
                    "empirical_AC": 5,
                    "empirical_AF": 0.10,
                    "fraction_phased": 1.0,
                    "estimated_gnomad_AF": 0.10,
                },
            ],
        },
    ]
    return hl.Table.parallelize(rows, schema=schema)


def test_export_sequences_table_to_tsv_per_pop_columns(
    hail_context: None,  # noqa: ARG001
    tmp_path: Path,
) -> None:
    """Per-pop flat columns are emitted for every joint pop, with nulls where source absent."""
    ht = _make_sequences_ht()
    out_file = tmp_path / "sequences.tsv.bgz"
    joint_pops_legend = ["afr", "amr"]

    export_sequences_table_to_tsv(ht=ht, out_file=out_file, joint_pops_legend=joint_pops_legend)

    # Read via the production loader so the writer + reader pipeline is end-to-end validated
    # (in particular, the `fill_null("NA")` step that restores bare-"NA" cells after polars'
    # `null_values` coercion).
    chunks = list(
        iter_dataframe_chunks(tsv=out_file, joint_pops_legend=joint_pops_legend, chunk_size=10)
    )
    assert len(chunks) == 1
    df = chunks[0]

    # Renamed scalar columns are present; pre-rename names are gone.
    assert "popmax_estimated_gnomad_AF" in df.columns
    assert "popmax_fraction_phased" in df.columns
    assert "estimated_gnomad_AF" not in df.columns
    assert "fraction_phased" not in df.columns

    # `max_pop` integer index resolved to label via the legend.
    assert df["max_pop"].to_list() == ["afr", "afr"]

    # Per-pop flat columns present for every joint pop.
    for pop in joint_pops_legend:
        for field in (
            "empirical_AC",
            "empirical_AF",
            "fraction_phased",
            "estimated_gnomAD_haplotype_AF",
        ):
            assert f"{field}_{pop}" in df.columns

    # gnomAD_AF_{pop} columns still emitted (comma-delimited per-variant strings).
    assert df["gnomAD_AF_afr"][0] == "0.05000,0.07000"
    # Variant 0 has an amr AF, variant 1 does not: the cell preserves both elements as a
    # comma-delimited string rather than collapsing the whole array to a single "NA".
    assert df["gnomAD_AF_amr"][0] == "0.04000,NA"
    # Row 1 has a missing amr AF on its single variant: the writer emits a bare "NA" cell,
    # polars' null_values turns that into None, and iter_dataframe_chunks restores "NA".
    assert df["gnomAD_AF_afr"][1] == "0.10000"
    assert df["gnomAD_AF_amr"][1] == "NA"

    # Row 0: both pops have data.
    assert df["empirical_AC_afr"][0] == 42
    assert df["empirical_AF_afr"][0] == pytest.approx(0.42)
    assert df["fraction_phased_afr"][0] == pytest.approx(1.40)
    assert df["estimated_gnomAD_haplotype_AF_afr"][0] == pytest.approx(0.07)
    assert df["empirical_AC_amr"][0] == 19
    assert df["empirical_AF_amr"][0] == pytest.approx(0.38)
    assert df["fraction_phased_amr"][0] == pytest.approx(1.52)
    assert df["estimated_gnomAD_haplotype_AF_amr"][0] == pytest.approx(0.05)

    # Row 1: only pop 0 present in all_pop_freqs; pop 1's flat columns must be null.
    assert df["empirical_AC_afr"][1] == 5
    assert df["empirical_AF_afr"][1] == pytest.approx(0.10)
    assert df["fraction_phased_afr"][1] == pytest.approx(1.0)
    assert df["estimated_gnomAD_haplotype_AF_afr"][1] == pytest.approx(0.10)
    assert df["empirical_AC_amr"][1] is None
    assert df["empirical_AF_amr"][1] is None
    assert df["fraction_phased_amr"][1] is None
    assert df["estimated_gnomAD_haplotype_AF_amr"][1] is None

    # Scalar `popmax_*` columns match the max_pop entry (pop 0 in both rows).
    assert df["popmax_empirical_AF"][0] == pytest.approx(0.42)
    assert df["popmax_fraction_phased"][0] == pytest.approx(1.40)
    assert df["popmax_estimated_gnomad_AF"][0] == pytest.approx(0.07)


def test_stream_empty_tsv_creates_sequences_table(tmp_path: Path) -> None:
    """
    A contig that yields zero rows still leaves a valid (empty) sequences table behind.

    Guards the finalize step: it must always find a `sequences` table, even if the first appended
    contig produced no rows. No Hail needed — this exercises the polars/DuckDB streaming path
    directly with a header-only TSV.
    """
    joint_pops_legend = ["afr"]
    # Mirror the real export header (export_sequences_table_to_tsv), including the `variants`,
    # `source`, and `n_variants` columns the compatibility-flag step reads.
    columns = [
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
        "popmax_estimated_gnomad_AF",
        "popmax_fraction_phased",
        "max_pop",
        "variants",
        "gnomAD_AF_afr",
        "empirical_AC_afr",
        "empirical_AF_afr",
        "fraction_phased_afr",
        "estimated_gnomAD_haplotype_AF_afr",
    ]
    header_only_tsv = tmp_path / "empty.sequences.tsv"
    header_only_tsv.write_text("\t".join(columns) + "\n")

    db = tmp_path / "idx.duckdb"
    with duckdb.connect(str(db)) as conn:
        appended = _stream_tsv_into_sequences(
            conn,
            tsv=header_only_tsv,
            joint_pops_legend=joint_pops_legend,
            chunk_size=100,
        )
        assert appended == 0
        exists = conn.execute(
            "SELECT 1 FROM information_schema.tables WHERE table_name = 'sequences'"
        ).fetchone()
        assert exists is not None
        count = conn.execute("SELECT COUNT(*) FROM sequences").fetchone()
        assert count is not None
        assert count[0] == 0


def test_stream_tsv_flags_overlapping_haplotype(tmp_path: Path) -> None:
    """
    Streaming a TSV with an overlapping haplotype persists its incompatibility reason.

    The committed e2e test data has no overlapping haplotypes, so this injects the case at the
    TSV -> DuckDB boundary where the flag is applied: an overlapping HGDP haplotype (a SNP at a
    deleted base) is flagged `snp_in_deletion`, while a clean haplotype and a `gnomAD_variant`
    row stay `PASS`.
    """
    columns = [
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
        "popmax_estimated_gnomad_AF",
        "popmax_fraction_phased",
        "max_pop",
        "variants",
        "gnomAD_AF_afr",
        "empirical_AC_afr",
        "empirical_AF_afr",
        "fraction_phased_afr",
        "estimated_gnomAD_haplotype_AF_afr",
    ]
    rows = [
        # overlapping: a SNP at a base the deletion removes -> snp_in_deletion
        [
            "ACGT",
            "4",
            "DR-9.9-0",
            "2",
            "chr1",
            "274",
            "326",
            "0.05",
            "7",
            "HGDP_haplotype",
            "0.004",
            "0.5",
            "afr",
            "chr1:300:AT:A,chr1:301:T:A",
            "0.1,0.2",
            "7",
            "0.05",
            "0.5",
            "0.004",
        ],
        # clean two-variant haplotype -> PASS
        [
            "ACGT",
            "4",
            "DR-9.9-1",
            "2",
            "chr1",
            "174",
            "235",
            "0.05",
            "5",
            "HGDP_haplotype",
            "0.004",
            "0.5",
            "afr",
            "chr1:200:A:T,chr1:210:C:G",
            "0.1,0.2",
            "5",
            "0.05",
            "0.5",
            "0.004",
        ],
        # single gnomAD variant -> PASS
        [
            "A",
            "1",
            "DR-9.9-2",
            "1",
            "chr1",
            "24",
            "76",
            "0.1",
            "100",
            "gnomAD_variant",
            "0.1",
            "1.0",
            "afr",
            "chr1:50:A:T",
            "0.1",
            "100",
            "0.1",
            "1.0",
            "0.1",
        ],
    ]
    tsv = tmp_path / "with_rows.sequences.tsv"
    tsv.write_text("\t".join(columns) + "\n" + "\n".join("\t".join(r) for r in rows) + "\n")

    db = tmp_path / "idx.duckdb"
    with duckdb.connect(str(db)) as conn:
        appended = _stream_tsv_into_sequences(
            conn, tsv=tsv, joint_pops_legend=["afr"], chunk_size=100
        )
        assert appended == 3
        flags = dict(conn.execute("SELECT sequence_id, haplotype_filter FROM sequences").fetchall())
    assert flags == {
        "DR-9.9-0": "snp_in_deletion",
        "DR-9.9-1": "PASS",
        "DR-9.9-2": "PASS",
    }


def test_add_compatibility_flag_values() -> None:
    """PASS for gnomAD/single/clean rows; the incompatibility reason for an overlapping pair."""
    df = polars.DataFrame({
        "variants": [
            "chr1:300:AT:A,chr1:301:T:A",  # SNP at a deleted base
            "chr1:200:A:T,chr1:210:C:G",  # clean haplotype
            "chr1:50:A:T",  # gnomAD single variant
        ],
        "source": ["HGDP_haplotype", "HGDP_haplotype", "gnomAD_variant"],
        "n_variants": [2, 2, 1],
    })
    out = _add_compatibility_flag(df)
    assert out["haplotype_filter"].to_list() == ["snp_in_deletion", "PASS", "PASS"]


def test_add_compatibility_flag_empty_frame() -> None:
    """An empty frame (used to create the table) still gains a String haplotype_filter column."""
    df = polars.DataFrame(
        {"variants": [], "source": [], "n_variants": []},
        schema={
            "variants": polars.String,
            "source": polars.String,
            "n_variants": polars.Int64,
        },
    )
    out = _add_compatibility_flag(df)
    assert out.height == 0
    assert out.schema["haplotype_filter"] == polars.String
