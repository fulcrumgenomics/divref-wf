from pathlib import Path

import pytest
from pydantic import ValidationError

from divref.tools.create_duckdb_from_tsv import SourceMetadata
from divref.tools.create_duckdb_from_tsv import read_source_metadata


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
        # A YAML `version: 1.10` (unquoted) parses as float 1.1 and would silently corrupt the
        # resource identity DR-{version}-N; StrictStr must reject the non-string scalar.
        pytest.param({"version": 1.10}, "version", id="non-string-version-rejected"),
        # Addendum A: source_name becomes a column-name prefix + DuckDB base name.
        pytest.param(
            {"source_name": "my cohort"},
            "identifier|source_name",
            id="non-identifier-source-name-rejected",
        ),
        # Addendum A: each population becomes a column-name suffix (`empirical_AC_<pop>`).
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
