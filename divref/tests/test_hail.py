"""Tests for divref.hail.hail_init input validation (no Hail context required)."""

from pathlib import Path

import pytest

from divref.hail import hail_init


def test_hail_init_missing_credentials_file_raises(tmp_path: Path) -> None:
    """GCS mode with a provided-but-missing credentials file fails loudly before starting Hail."""
    with pytest.raises(FileNotFoundError, match="credentials file not found"):
        hail_init(gcs_credentials_path=tmp_path / "does_not_exist.json", use_s3=False)


def test_hail_init_requires_credentials_for_gcs() -> None:
    """GCS mode with no credentials path raises ValueError."""
    with pytest.raises(ValueError, match="gcs_credentials_path is required"):
        hail_init(gcs_credentials_path=None, use_s3=False)
