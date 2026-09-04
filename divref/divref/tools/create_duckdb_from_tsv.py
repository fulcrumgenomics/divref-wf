"""Build a DivRef DuckDB index from a wide single-variant TSV (pure Python; no Hail)."""

import logging
import re
from pathlib import Path
from typing import Literal

import yaml
from pydantic import BaseModel
from pydantic import ConfigDict
from pydantic import StrictStr
from pydantic import field_validator

logger = logging.getLogger(__name__)

# Addendum A: source_name is a column-name prefix (`<source>_AF_<pop>`) and the DuckDB base name;
# each population is a column-name suffix (`empirical_AC_<pop>`). Both must be bare identifiers.
_IDENTIFIER = re.compile(r"[A-Za-z][A-Za-z0-9_]*")


class SourceMetadata(BaseModel):
    """Sidecar metadata for one TSV variation source."""

    model_config = ConfigDict(frozen=True)

    source_name: StrictStr
    # StrictStr so an unquoted YAML `version: 1.10` (float 1.1) is rejected, not coerced to "1.1"
    # and baked into sequence_id DR-{version}-N.
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
