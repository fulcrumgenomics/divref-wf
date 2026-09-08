"""Build a DivRef DuckDB index from a wide single-variant TSV (pure Python; no Hail)."""

import logging
import re
from collections.abc import Sequence
from pathlib import Path
from typing import Literal

import polars
import yaml
from pydantic import BaseModel
from pydantic import ConfigDict
from pydantic import StrictStr
from pydantic import field_validator

logger = logging.getLogger(__name__)

# source_name and each population are interpolated into column names (`<source>_AF_<pop>`,
# `empirical_AC_<pop>`), so both must be valid identifiers.
_IDENTIFIER = re.compile(r"[A-Za-z][A-Za-z0-9_]*")


class SourceMetadata(BaseModel):
    """Sidecar metadata for one TSV variation source."""

    model_config = ConfigDict(frozen=True)

    source_name: StrictStr
    # StrictStr: an unquoted YAML `version: 1.10` (float 1.1) would corrupt DR-{version}-N.
    version: StrictStr
    reference_genome: Literal["GRCh38"]
    populations: list[str]

    @field_validator("source_name")
    @classmethod
    def _valid_identifier(cls, value: str) -> str:
        if _IDENTIFIER.fullmatch(value) is None:
            raise ValueError("source_name must be a valid identifier ^[A-Za-z][A-Za-z0-9_]*$.")
        return value

    @field_validator("populations")
    @classmethod
    def _non_empty_unique(cls, value: list[str]) -> list[str]:
        if len(value) < 1:
            raise ValueError("populations must contain at least one population.")
        if len(set(value)) != len(value):
            raise ValueError("populations must be unique.")
        for population in value:
            if _IDENTIFIER.fullmatch(population) is None:
                raise ValueError(
                    f"population {population!r} must be a valid identifier ^[A-Za-z][A-Za-z0-9_]*$."
                )
        return value


def read_source_metadata(path: Path) -> SourceMetadata:
    """
    Read and validate a source_meta.yml sidecar.

    Args:
        path: Path to the source_meta.yml file.

    Returns:
        The validated `SourceMetadata`.

    Raises:
        ValueError: If the YAML document does not parse to a mapping.
    """
    with path.open() as handle:
        data: object = yaml.safe_load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"{path} is not a YAML mapping.")
    return SourceMetadata.model_validate(data)


_FIXED_COLUMNS = frozenset({"contig", "pos", "ref", "alt"})


def validate_variants_header(columns: list[str], populations: list[str]) -> None:
    """
    Assert a variants-TSV header matches the fixed columns plus the population legend.

    The header must contain `contig, pos, ref, alt` plus `AC_<pop>` and `AF_<pop>` for each
    population in `populations`, in any order. Extra columns that are not `AC_`/`AF_` (e.g. `rsid`)
    are permitted and ignored; an `AC_`/`AF_` column naming an unknown population is rejected.

    Args:
        columns: The TSV header, as parsed column names.
        populations: The population codes from the `source_meta.yml` legend.

    Raises:
        ValueError: If a required column is missing, or an `AC_`/`AF_` column names a
            population outside `populations`.
    """
    expected = set(_FIXED_COLUMNS)
    for pop in populations:
        expected |= {f"AC_{pop}", f"AF_{pop}"}
    present = set(columns)
    missing = expected - present
    if missing:
        raise ValueError(f"variants TSV is missing required columns: {sorted(missing)}")
    unexpected = {c for c in present - expected if c.startswith(("AC_", "AF_"))}
    if unexpected:
        raise ValueError(
            f"variants TSV has unexpected/unknown population columns: {sorted(unexpected)}"
        )


_CONTIG_PATTERN = r"^chr([1-9]|1[0-9]|2[0-2]|X|Y)$"
_BASE_PATTERN = r"^[ACGTN]+$"


def read_and_validate_variants(path: Path, populations: list[str]) -> polars.DataFrame:
    """
    Read the wide variants TSV with a legend-driven schema and validate it on the fly.

    Validates the header against `populations`, then checks: the fixed columns are non-null;
    `pos` is >= 1; `contig` is a GRCh38 main-contig token; `ref`/`alt` match `^[ACGTN]+$`; each
    population's `AC_`/`AF_` pair is both-present-or-both-null with `AC_ >= 0` and
    `AF_ in [0, 1]`; and each row has at least one defined AF.

    Args:
        path: Path to the wide variants TSV.
        populations: The population codes from the `source_meta.yml` legend.

    Returns:
        A `polars.DataFrame` with columns `contig, pos, ref, alt` and, for each population,
        `AC_<pop>: Int64` / `AF_<pop>: Float64`.

    Raises:
        ValueError: If the header does not match `populations`, or any row fails validation.
    """
    overrides: dict[str, type[polars.DataType]] = {
        "contig": polars.String,
        "pos": polars.Int64,
        "ref": polars.String,
        "alt": polars.String,
    }
    for pop in populations:
        overrides[f"AC_{pop}"] = polars.Int64
        overrides[f"AF_{pop}"] = polars.Float64
    header = polars.read_csv(path, separator="\t", n_rows=0).columns
    validate_variants_header(header, populations)
    df = polars.read_csv(path, separator="\t", schema_overrides=overrides, null_values=["", "NA"])

    # 1. The four fixed columns are required and non-null (empty pos/contig must not slip through
    #    the range/regex checks below, which are null-tolerant).
    for col in ("contig", "pos", "ref", "alt"):
        _validate(df, df[col].is_not_null(), f"{col} is required and must be non-empty")
    # 2. Fixed-field ranges/shape. str.contains is vectorized and null-safe (null -> caught above).
    _validate(df, df["pos"] >= 1, "pos must be >= 1")
    _validate(
        df, df["contig"].str.contains(_CONTIG_PATTERN), "contig must be a GRCh38 main-contig token"
    )
    for col in ("ref", "alt"):
        _validate(
            df, df[col].str.contains(_BASE_PATTERN), f"{col} must be A/C/G/T/N bases (^[ACGTN]+$)"
        )
    # 3. Per-pop: AC/AF must both be present or neither for each population, with range guards.
    for pop in populations:
        af, ac = df[f"AF_{pop}"], df[f"AC_{pop}"]
        _validate(
            df,
            af.is_null() == ac.is_null(),
            f"AC_{pop} and AF_{pop} must both be present or both empty",
        )
        _validate(df, af.is_null() | ((af >= 0.0) & (af <= 1.0)), f"AF_{pop} must be in [0, 1]")
        _validate(df, ac.is_null() | (ac >= 0), f"AC_{pop} must be >= 0")
    # 4. At least one defined AF per row (computed once).
    any_af = df.select(
        polars.any_horizontal(polars.col(f"AF_{p}").is_not_null() for p in populations).alias("m")
    )["m"]
    _validate(df, any_af, "each variant must have at least one defined AF")
    return df


def _validate(df: polars.DataFrame, mask: polars.Series, message: str) -> None:
    """
    Raise a `ValueError` naming the first offending variant if `mask` is not all-True.

    Args:
        df: The variants frame being validated; must have `contig, pos, ref, alt` columns.
        mask: A boolean, row-aligned series. A null counts as a failure (`fill_null(False)`).
        message: The failure reason to report.

    Raises:
        ValueError: If any element of `mask` is false or null.
    """
    bad = df.filter(~mask.fill_null(False))
    if bad.height > 0:
        first = bad.row(0, named=True)
        raise ValueError(
            f"{message} ({bad.height} row(s) failed). First offending variant: "
            f"{first['contig']}:{first['pos']}:{first['ref']}:{first['alt']}"
        )


def _format_raw_number(value: float | None) -> str | None:
    """
    Render a "raw" (unformatted) AF or fraction-phased value as Hail's TSV exporter does.

    Verified against every gnomAD single-variant row's true (Hail-table) AF in the committed
    golden: Hail's `Double` -> text export for these columns rounds to 5 SIGNIFICANT figures
    (not decimal places), then renders that rounded value with the shortest round-trip decimal
    text -- switching to scientific notation for small magnitudes and always keeping a decimal
    point (`1.0`, not `1`). `%.5g` alone does not reproduce this: Python's `%g` strips the
    decimal point off whole numbers (`1`, not `1.0`), so the 5-sig-fig rounding and the
    decimal-text rendering must be two separate steps. This differs from
    `polars.Series.cast(pl.String)` on a plain `Float64` column, which neither rounds to 5
    significant figures nor ever emits scientific notation. Unlike the `%.5f`-formatted
    `<source>_AF_<pop>` column, `popmax_empirical_AF`, `empirical_AF_<pop>`,
    `fraction_phased_*`, and `estimated_<source>_haplotype_AF_<pop>` are passed through this way
    rather than left as a native `Float64` column.

    Args:
        value: The value, or `None` if this population/row has no data.

    Returns:
        The 5-significant-figure decimal text of `value`, or `None` if `value` is `None`.
    """
    if value is None:
        return None
    return str(float(f"{value:.5g}"))


def _hail_argmax(values: Sequence[float | None]) -> int:
    """
    Replicate `hl.argmax`'s index-of-maximum semantics over per-population AFs.

    Null AFs are skipped. On a tie, the first (lowest) index wins, matching Hail's
    `hl.argmax(unique=False)` default used by `build_gnomad_variant_table_entries`.

    Args:
        values: Per-population AF values, in population-legend order; a population with no
            data is `None`.

    Returns:
        The index of the maximum non-null value.

    Raises:
        ValueError: If every value is `None`.
    """
    best_index: int | None = None
    best_value: float | None = None
    for index, value in enumerate(values):
        if value is None:
            continue
        if best_value is None or value > best_value:
            best_value = value
            best_index = index
    if best_index is None:
        raise ValueError("Cannot compute popmax: every population AF is null for this variant.")
    return best_index
