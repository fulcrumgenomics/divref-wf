"""Dry-run integration check for `workflows/generate_divref_from_tsv.smk`."""

import shutil
import subprocess
from pathlib import Path

import pytest

# tests/ -> divref/ -> repo root, where `workflows/` lives.
_REPO_ROOT = Path(__file__).resolve().parents[2]


@pytest.mark.skipif(
    shutil.which("pixi") is None,
    reason="needs the pixi-managed snakemake environment (run under `pixi run pytest`)",
)
def test_workflow_dag_dry_run() -> None:
    """Dry-run the TSV-source workflow against the committed test config and check its DAG."""
    result = subprocess.run(
        [
            "pixi",
            "run",
            "snakemake",
            "-s",
            "workflows/generate_divref_from_tsv.smk",
            "--configfile",
            "workflows/config/config_tsv_test.yml",
            "-n",
        ],
        cwd=_REPO_ROOT,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    # Snakemake may emit the dry-run job listing on either stream across versions; check both.
    dag = result.stdout + result.stderr
    assert "create_tsv_index" in dag
    assert "create_tsv_fasta" in dag
