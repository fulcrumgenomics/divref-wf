from pathlib import Path

import polars as pl
import pysam
import pytest
from pydantic import ValidationError

from divref.duckdb_index import sequences_tsv_columns
from divref.tools.create_duckdb_from_tsv import SourceMetadata
from divref.tools.create_duckdb_from_tsv import _format_raw_number
from divref.tools.create_duckdb_from_tsv import _hail_argmax
from divref.tools.create_duckdb_from_tsv import build_sequences_frame
from divref.tools.create_duckdb_from_tsv import read_and_validate_variants
from divref.tools.create_duckdb_from_tsv import read_source_metadata
from divref.tools.create_duckdb_from_tsv import validate_variants_header

# NB the fixture is `…fa.gz`, so its index is `…fa.fai` (samtools/pysam append `.fai` to the full
# name). Derive it from the `.fa.gz` path -- `("…fa.gz").with_suffix(".fai")` replaces `.gz` and
# yields `…fa.fai`. Do NOT apply `.with_suffix(".fai")` to a `…fa` name: that replaces `.fa` and
# yields `…fai`, which does not exist. Simplest is the literal `…fa.fai`.
_FIXTURE_FA = "test_reference.chr1_chrX.fa.gz"
_FIXTURE_FAI = "test_reference.chr1_chrX.fa.fai"


def test_read_source_metadata_happy(datadir: Path) -> None:
    meta = read_source_metadata(datadir / "tsv_source" / "source_meta.yml")
    assert meta.source_name == "test_cohort"
    assert meta.version == "9.9"
    assert meta.reference_genome == "GRCh38"
    assert meta.populations == ["afr", "eas"]


@pytest.mark.parametrize(
    "overrides, match",
    [
        pytest.param({"reference_genome": "GRCh37"}, "GRCh38", id="non-grch38-rejected"),
        pytest.param({"populations": []}, "at least one", id="empty-populations-rejected"),
        pytest.param(
            {"populations": ["afr", "afr"]}, "unique", id="duplicate-populations-rejected"
        ),
        pytest.param({"version": 1.10}, "version", id="non-string-version-rejected"),
        pytest.param(
            {"source_name": "my cohort"},
            "identifier|source_name",
            id="non-identifier-source-name-rejected",
        ),
        pytest.param(
            {"populations": ["af-r"]},
            "identifier|population",
            id="non-identifier-population-rejected",
        ),
    ],
)
def test_source_metadata_refusals(overrides: dict[str, object], match: str) -> None:
    base = {
        "source_name": "s",
        "version": "1.0",
        "reference_genome": "GRCh38",
        "populations": ["afr"],
    }
    with pytest.raises(ValidationError, match=match):
        SourceMetadata.model_validate({**base, **overrides})


def test_header_ok() -> None:
    validate_variants_header(
        ["contig", "pos", "ref", "alt", "AC_afr", "AF_afr", "AC_eas", "AF_eas"], ["afr", "eas"]
    )


def test_header_unknown_pop_rejected() -> None:
    with pytest.raises(ValueError, match="unexpected|unknown"):
        validate_variants_header(
            ["contig", "pos", "ref", "alt", "AC_afr", "AF_afr", "AC_sas", "AF_sas"], ["afr"]
        )


def test_header_missing_pop_rejected() -> None:
    with pytest.raises(ValueError, match="missing"):
        validate_variants_header(
            ["contig", "pos", "ref", "alt", "AC_afr", "AF_afr"], ["afr", "eas"]
        )


def test_read_variants_happy(datadir: Path) -> None:
    df = read_and_validate_variants(datadir / "tsv_source" / "variants.tsv", ["afr", "eas"])
    assert df.height == 4
    assert df["AC_afr"].dtype == pl.Int64
    assert df["AF_afr"].dtype == pl.Float64


@pytest.mark.parametrize(
    "row, match",
    [
        pytest.param("chr1\t0\tA\tG\t1\t0.1\t1\t0.1", "pos", id="pos-lt-1"),
        pytest.param("chr1\t\tA\tG\t1\t0.1\t1\t0.1", "pos|missing|null", id="null-pos"),
        pytest.param("\t10\tA\tG\t1\t0.1\t1\t0.1", "contig|missing|null", id="null-contig"),
        pytest.param("chr1\t10\tA\tX\t1\t0.1\t1\t0.1", "ref/alt|base", id="bad-base"),
        pytest.param("chr1\t10\tA\tG\t1\t1.5\t1\t0.1", "0, 1|range", id="af-out-of-range"),
        pytest.param("chr1\t10\tA\tG\t-1\t0.1\t1\t0.1", ">= 0|negative", id="negative-ac"),
        pytest.param("chr1\t10\tA\tG\t\t\t\t", "at least one", id="all-missing-af"),
        # both-or-neither per pop: a lone AC or lone AF is a refusal (no Hail oracle for the
        # asymmetric rendering).
        pytest.param("chr1\t10\tA\tG\t1\t0.1\t\t0.1", "both|AC.*AF", id="lone-af-eas"),
        pytest.param("chr1\t10\tA\tG\t1\t0.1\t1\t", "both|AC.*AF", id="lone-ac-eas"),
    ],
)
def test_variant_row_refusals(tmp_path: Path, row: str, match: str) -> None:
    p = tmp_path / "v.tsv"
    p.write_text("contig\tpos\tref\talt\tAC_afr\tAF_afr\tAC_eas\tAF_eas\n" + row + "\n")
    with pytest.raises(ValueError, match=match):
        read_and_validate_variants(p, ["afr", "eas"])


@pytest.mark.parametrize(
    "value, expected",
    [
        # An integer-valued AF keeps its decimal point (Hail's Double export, not Python %g).
        pytest.param(1.0, "1.0", id="integer-valued-af-keeps-decimal-point"),
        # 5 sig figs of 3.0312e-05 is itself; magnitude < 1e-4 renders in scientific notation.
        # Cross-checked against duckdb_index_golden/sequences.chr1_chrX.tsv's empirical_AF_afr.
        pytest.param(3.0312e-05, "3.0312e-05", id="small-af-renders-scientific-notation"),
        # 0.123456789 -> first 5 significant digits 12345, 6th digit 6 rounds the 5th up: 0.12346.
        pytest.param(0.123456789, "0.12346", id="value-needs-5-sig-fig-rounding"),
        pytest.param(None, None, id="none-passes-through"),
    ],
)
def test_format_raw_number(value: float | None, expected: str | None) -> None:
    assert _format_raw_number(value) == expected


@pytest.mark.parametrize(
    "values, expected_index",
    [
        pytest.param([None, 0.2, 0.9, 0.3], 2, id="skip-none"),
        pytest.param([0.5, 0.5, 0.1], 0, id="tie-first-index-wins"),
        pytest.param([0.1, 0.4, 0.2], 1, id="normal-max"),
        pytest.param([None, None], None, id="all-none-raises"),
    ],
)
def test_hail_argmax(values: list[float | None], expected_index: int | None) -> None:
    if expected_index is None:
        with pytest.raises(ValueError, match="popmax"):
            _hail_argmax(values)
    else:
        assert _hail_argmax(values) == expected_index


def _one_variant_wide(contig: str, pos: int, ref: str, alt: str) -> pl.DataFrame:
    """Build a 1-row wide variants frame (population `afr` only)."""
    return pl.DataFrame({
        "contig": [contig],
        "pos": [pos],
        "ref": [ref],
        "alt": [alt],
        "AC_afr": [10],
        "AF_afr": [0.1],
    })


@pytest.mark.parametrize(
    "kind",
    [
        pytest.param("snp", id="snp"),
        pytest.param("ins", id="ins"),
        pytest.param("del", id="del"),
    ],
)
def test_builder_snp_ins_del(datadir: Path, kind: str) -> None:
    fa = pysam.FastaFile(str(datadir / _FIXTURE_FA), filepath_index=str(datadir / _FIXTURE_FAI))
    w, pos, contig = 25, 100_100, "chr1"  # interior position within the fixture's covered range
    # ref must match the reference bases at pos (build_sequences_frame now validates this), so
    # derive ref -- and a consistent alt -- from the reference rather than hardcoding bases.
    if kind == "snp":
        ref = fa.fetch(contig, pos - 1, pos)
        alt = "A" if ref != "A" else "C"  # any base different from ref
    elif kind == "ins":
        ref = fa.fetch(contig, pos - 1, pos)
        alt = ref + "A"
    else:
        ref = fa.fetch(contig, pos - 1, pos + 1)
        alt = ref[0]
    left = fa.fetch(contig, pos - 1 - w, pos - 1)
    right = fa.fetch(contig, pos - 1 + len(ref), pos - 1 + len(ref) + w)
    expected_seq = left + alt + right  # bases as-is; NO .upper()
    df = _one_variant_wide(contig, pos, ref, alt)
    row = build_sequences_frame(
        df=df,
        populations=["afr"],
        reference=datadir / _FIXTURE_FA,
        window_size=w,
        version="9.9",
        source="s",
        sequence_id_offset=0,
    ).row(0, named=True)
    assert row["sequence"] == expected_seq
    assert row["start"] == pos - 1 - w
    assert row["end"] == pos - 1 + len(ref) + w
    assert row["sequence_length"] == 2 * w + len(alt)


def test_builder_rejects_edge_of_contig_variant(datadir: Path) -> None:
    df = _one_variant_wide("chr1", 5, "A", "G")  # pos-1-w = -21 < 0 for w=25; fails before any
    # ref/reference comparison, so the (possibly wrong) hardcoded ref here doesn't matter.
    with pytest.raises(ValueError, match="contig bounds|out of"):
        build_sequences_frame(
            df=df,
            populations=["afr"],
            reference=datadir / _FIXTURE_FA,
            window_size=25,
            version="9.9",
            source="s",
            sequence_id_offset=0,
        )


def test_builder_rejects_ref_reference_mismatch(datadir: Path) -> None:
    fa = pysam.FastaFile(str(datadir / _FIXTURE_FA), filepath_index=str(datadir / _FIXTURE_FAI))
    pos, contig = 100_100, "chr1"
    true_ref = fa.fetch(contig, pos - 1, pos)
    wrong_ref = "A" if true_ref != "A" else "C"  # deliberately does not match the reference
    df = _one_variant_wide(contig, pos, wrong_ref, "G")
    with pytest.raises(ValueError, match="reference bases"):
        build_sequences_frame(
            df=df,
            populations=["afr"],
            reference=datadir / _FIXTURE_FA,
            window_size=25,
            version="9.9",
            source="s",
            sequence_id_offset=0,
        )


# Popmax tie-break and null-skip are NOT covered by the golden (it has 0 AF ties; null-AF pops
# appear in only 10 rows), so pin `hl.argmax`'s two edge behaviours here with hand-built rows.
def test_builder_popmax_tie_breaks_on_first_index(datadir: Path) -> None:
    # Two pops share the max AF. hl.argmax(unique=False) returns the LOWEST index, so max_pop
    # must be the first such pop in legend order (afr), with its AC.
    df = pl.DataFrame({
        "contig": ["chr1"],
        "pos": [100_100],
        "ref": ["T"],  # true reference base at chr1:100100
        "alt": ["G"],
        "AC_afr": [10],
        "AF_afr": [0.20],
        "AC_eas": [7],
        "AF_eas": [0.20],
    })
    row = build_sequences_frame(
        df=df,
        populations=["afr", "eas"],
        reference=datadir / _FIXTURE_FA,
        window_size=25,
        version="9.9",
        source="s",
        sequence_id_offset=0,
    ).row(0, named=True)
    assert row["max_pop"] == "afr"
    assert row["popmax_empirical_AC"] == 10


def test_builder_popmax_skips_null_af_pops(datadir: Path) -> None:
    # afr has no data (both AC/AF empty); argmax must skip the missing element and pick eas,
    # matching hl.argmax's missing-element semantics. afr's gnomAD_AF cell renders "NA".
    df = pl.DataFrame(
        {
            "contig": ["chr1"],
            "pos": [100_100],
            "ref": ["T"],  # true reference base at chr1:100100
            "alt": ["G"],
            "AC_afr": [None],
            "AF_afr": [None],
            "AC_eas": [3],
            "AF_eas": [0.05],
        },
        schema_overrides={"AC_afr": pl.Int64, "AF_afr": pl.Float64},
    )
    row = build_sequences_frame(
        df=df,
        populations=["afr", "eas"],
        reference=datadir / _FIXTURE_FA,
        window_size=25,
        version="9.9",
        source="s",
        sequence_id_offset=0,
    ).row(0, named=True)
    assert row["max_pop"] == "eas"
    assert row["s_AF_afr"] == "NA"


def test_builder_empty_variants_frame_returns_typed_empty_frame(datadir: Path) -> None:
    df = pl.DataFrame(
        schema={"contig": pl.String, "pos": pl.Int64, "ref": pl.String, "alt": pl.String}
    )
    built = build_sequences_frame(
        df=df,
        populations=["afr"],
        reference=datadir / _FIXTURE_FA,
        window_size=25,
        version="9.9",
        source="s",
        sequence_id_offset=0,
    )
    assert built.height == 0
    assert built.columns == sequences_tsv_columns(
        ["afr"], af_prefix="s", popmax_estimated_col="popmax_estimated_s_AF"
    )


def test_builder_rejects_all_null_popmax(datadir: Path) -> None:
    df = pl.DataFrame(
        # ref "T" is the true reference base at chr1:100100 (this test targets the all-null-AF
        # ValueError, not the ref/reference-mismatch ValueError, so ref must be valid here).
        {"contig": ["chr1"], "pos": [100_100], "ref": ["T"], "alt": ["G"], "AF_afr": [None]},
        schema_overrides={"AF_afr": pl.Float64},
    )
    with pytest.raises(ValueError, match="popmax"):
        build_sequences_frame(
            df=df,
            populations=["afr"],
            reference=datadir / _FIXTURE_FA,
            window_size=25,
            version="9.9",
            source="s",
            sequence_id_offset=0,
        )
