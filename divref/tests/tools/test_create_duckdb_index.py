"""Tests for the create_duckdb_index tool."""

from pathlib import Path

import hail as hl
import pytest

from divref.tools.create_duckdb_index import export_sequences_table_to_tsv
from divref.tools.create_duckdb_index import iter_dataframe_chunks


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
