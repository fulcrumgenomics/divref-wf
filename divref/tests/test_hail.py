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


@pytest.mark.parametrize(("driver_gb", "executor_gb"), [(0, 1), (1, 0)])
def test_hail_init_rejects_sub_1gb_memory(driver_gb: int, executor_gb: int) -> None:
    """
    Spark driver/executor memory below 1GB fails loudly before starting Hail.

    The memory guards run before any credentials/JAR resolution, so ``use_s3=True`` isolates them.
    """
    with pytest.raises(ValueError, match="at least 1GB"):
        hail_init(
            spark_driver_memory_gb=driver_gb,
            spark_executor_memory_gb=executor_gb,
            use_s3=True,
        )
