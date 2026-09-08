"""Build a DivRef DuckDB index from a wide single-variant TSV (pure Python; no Hail)."""

import logging
import re
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import polars
import pysam
import yaml
from pydantic import BaseModel
from pydantic import ConfigDict
from pydantic import StrictStr
from pydantic import field_validator

from divref.duckdb_index import sequences_tsv_columns

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
    Render a "raw" AF or fraction-phased value as Hail's TSV exporter does.

    Round to 5 significant figures, keep a decimal point (`1.0`, not `1`), and use scientific
    notation for small magnitudes. Round-tripping through `float()` restores the decimal point
    that `%.5g` strips off whole numbers.

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


def _as_int(value: object) -> int:
    """Narrow a polars row cell (known to be a Python `int`) to `int` without a `cast()`."""
    if not isinstance(value, int):
        raise TypeError(f"Expected int, got {type(value).__name__}: {value!r}")
    return value


def _as_str(value: object) -> str:
    """Narrow a polars row cell (known to be a Python `str`) to `str` without a `cast()`."""
    if not isinstance(value, str):
        raise TypeError(f"Expected str, got {type(value).__name__}: {value!r}")
    return value


def _as_optional_float(value: object) -> float | None:
    """Narrow a polars row cell (a `float` or `None`) without a `cast()`."""
    if value is None:
        return None
    if not isinstance(value, float):
        raise TypeError(f"Expected float or None, got {type(value).__name__}: {value!r}")
    return value


def _as_optional_int(value: object) -> int | None:
    """Narrow a polars row cell (an `int` or `None`) without a `cast()`."""
    if value is None:
        return None
    if not isinstance(value, int):
        raise TypeError(f"Expected int or None, got {type(value).__name__}: {value!r}")
    return value


# eq=False so the frozen dataclass keeps a (default, identity-based) __hash__ despite its `dict`
# field; without it, frozen+eq would synthesize a __hash__ that raises TypeError on the
# unhashable `values` dict. Never compared or hashed by value -- only field-accessed.
@dataclass(frozen=True, kw_only=True, eq=False)
class _SequenceRow:
    """
    One built sequences-row, plus the `(start, variants)` sort key.

    Attributes:
        values: The row's column values, keyed by column name. Does not yet contain
            `sequence_id`; `build_sequences_frame` fills it in after sorting.
        start: The row's window start (0-based), used as the primary sort key.
        variants: The row's `variants` string, used as the sort tie-break.
    """

    values: dict[str, object]
    start: int
    variants: str


def _build_sequence_row(
    *,
    record: dict[str, object],
    populations: list[str],
    fasta: pysam.FastaFile,
    window_size: int,
    source: str,
) -> _SequenceRow:
    """
    Build one variant's sequences row, matching `build_gnomad_variant_table_entries`.

    Args:
        record: One row of the wide variants frame (`contig, pos, ref, alt` plus `AC_<pop>` /
            `AF_<pop>` per population), as returned by `polars.DataFrame.iter_rows(named=True)`.
        populations: Ordered population codes; drives per-pop column names and popmax order.
        fasta: Open reference FASTA reader.
        window_size: Flanking reference context size on each side of the variant.
        source: Source name; prefixes the annotation columns (Addendum A).

    Returns:
        The built `_SequenceRow`.

    Raises:
        ValueError: If the variant's flanking window falls outside the contig's bounds, or if
            the supplied `ref` does not match the reference bases at `contig:pos`.
    """
    contig = _as_str(record["contig"])
    pos = _as_int(record["pos"])
    ref = _as_str(record["ref"])
    alt = _as_str(record["alt"])
    variants = f"{contig}:{pos}:{ref}:{alt}"

    start = pos - 1 - window_size
    end = pos - 1 + len(ref) + window_size
    contig_length = fasta.get_reference_length(contig)
    if start < 0 or end > contig_length:
        raise ValueError(
            f"Variant {variants} window [{start}, {end}) is out of contig bounds "
            f"[0, {contig_length}) for {contig}."
        )
    # Bases as-is (no .upper()) throughout: matches Hail's case-preserving reference reader.
    # A mismatch here means the input TSV's ref (or its pos) is wrong relative to this
    # reference -- fail loud rather than silently building a sequence around the wrong context.
    reference_bases = fasta.fetch(contig, pos - 1, pos - 1 + len(ref))
    if reference_bases != ref:
        raise ValueError(
            f"Variant {variants}: reference bases at {contig}:{pos} are {reference_bases!r}, "
            f"not the supplied ref {ref!r}."
        )
    left = fasta.fetch(contig, start, pos - 1)
    right = fasta.fetch(contig, pos - 1 + len(ref), end)
    sequence = left + alt + right

    afs = [_as_optional_float(record[f"AF_{pop}"]) for pop in populations]
    popmax_index = _hail_argmax(afs)
    max_pop = populations[popmax_index]
    popmax_af = afs[popmax_index]
    popmax_ac = _as_optional_int(record[f"AC_{max_pop}"])

    values: dict[str, object] = {
        "sequence": sequence,
        "sequence_length": len(sequence),
        "n_variants": 1,
        "contig": contig,
        "start": start,
        "end": end,
        "popmax_empirical_AF": _format_raw_number(popmax_af),
        "popmax_empirical_AC": popmax_ac,
        "source": source,
        f"popmax_estimated_{source}_AF": _format_raw_number(popmax_af),
        "popmax_fraction_phased": _format_raw_number(1.0),
        "max_pop": max_pop,
        "variants": variants,
    }
    for pop, af in zip(populations, afs, strict=True):
        ac = _as_optional_int(record[f"AC_{pop}"])
        values[f"{source}_AF_{pop}"] = "NA" if af is None else f"{af:.5f}"
        values[f"empirical_AC_{pop}"] = ac
        values[f"empirical_AF_{pop}"] = _format_raw_number(af)
        values[f"fraction_phased_{pop}"] = _format_raw_number(1.0)
        values[f"estimated_{source}_haplotype_AF_{pop}"] = _format_raw_number(af)

    return _SequenceRow(values=values, start=start, variants=variants)


def build_sequences_frame(
    *,
    df: polars.DataFrame,
    populations: list[str],
    reference: Path,
    window_size: int,
    version: str,
    source: str,
    sequence_id_offset: int,
) -> polars.DataFrame:
    """
    Build the per-variant "sequences" rows for one wide single-variant TSV.

    Reproduces `build_gnomad_variant_table_entries`'s single-variant semantics in pure Python:
    for each variant, fetches flanking reference context with `pysam`, computes popmax by
    replicating `hl.argmax` (skip null AFs, first index wins on ties), and formats per-population
    columns to match the Hail TSV exporter byte-for-byte. The annotation columns are prefixed
    with `source` (`<source>_AF_<pop>`, `estimated_<source>_haplotype_AF_<pop>`,
    `popmax_estimated_<source>_AF`); `empirical_*` and `fraction_phased_*` are not prefixed.

    Args:
        df: Wide variants frame with `contig, pos, ref, alt` and, per population, `AC_<pop>` /
            `AF_<pop>` (each nullable).
        populations: Ordered population codes; drives per-pop column names and popmax order.
        reference: Path to the (optionally bgzipped) reference FASTA. Its index is resolved as
            `reference.with_suffix(".fai")`.
        window_size: Flanking reference context size on each side of the variant.
        version: Version identifier baked into `sequence_id` (`DR-{version}-{n}`).
        source: Source name; prefixes the annotation columns (Addendum A) and fills the
            `source` column.
        sequence_id_offset: Added to each row's post-sort index to form `sequence_id`.

    Returns:
        A `polars.DataFrame` with columns exactly `sequences_tsv_columns(populations,
        af_prefix=source, popmax_estimated_col=f"popmax_estimated_{source}_AF")`, sorted by
        `(start, variants)` with `sequence_id` assigned in that order.

    Raises:
        ValueError: If any variant's flanking window falls outside the contig's bounds, or its
            `ref` does not match the reference bases at its `contig:pos`.
    """
    columns = sequences_tsv_columns(
        populations, af_prefix=source, popmax_estimated_col=f"popmax_estimated_{source}_AF"
    )
    if df.height == 0:
        return polars.DataFrame({column: [] for column in columns})

    with pysam.FastaFile(
        str(reference), filepath_index=str(reference.with_suffix(".fai"))
    ) as fasta:
        rows = [
            _build_sequence_row(
                record=record,
                populations=populations,
                fasta=fasta,
                window_size=window_size,
                source=source,
            )
            for record in df.iter_rows(named=True)
        ]

    # (start, variants) tie-break equals Hail's (min_pos, source, variants) within a single
    # source (append_contig_to_duckdb_index.py:265-271).
    rows.sort(key=lambda row: (row.start, row.variants))
    for index, row in enumerate(rows):
        row.values["sequence_id"] = f"DR-{version}-{sequence_id_offset + index}"

    return polars.DataFrame([row.values for row in rows]).select(columns)
