from pathlib import Path

import polars as pl
import pytest
from pydantic import ValidationError

from divref.tools.create_duckdb_from_tsv import SourceMetadata
from divref.tools.create_duckdb_from_tsv import read_and_validate_variants
from divref.tools.create_duckdb_from_tsv import read_source_metadata
from divref.tools.create_duckdb_from_tsv import validate_variants_header


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
