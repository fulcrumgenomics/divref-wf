"""Tests for the extract_sample_metadata tool."""

from pathlib import Path
from unittest.mock import patch

import hail as hl

from divref.tools.extract_sample_metadata import extract_sample_metadata


def test_extract_sample_metadata(
    hail_context: None,  # noqa: ARG001
    datadir: Path,
    tmp_path: Path,
) -> None:
    """Happy-path: extract sample metadata from local test data."""
    in_samples = str(datadir / "hgdp_1kg_sample_metadata.ht")
    out_sa = tmp_path / "sa.ht"

    with patch("divref.tools.extract_sample_metadata.hail_init"):
        extract_sample_metadata(
            in_gnomad_hgdp_sample_data=in_samples,
            out_sample_metadata=out_sa,
        )

    sa = hl.read_table(str(out_sa))
    sa_count = sa.count()
    assert sa_count == 4151

    # Each sample should have a pop and sex_karyotype field
    first_sample = sa.head(1).collect()[0]
    assert hasattr(first_sample, "pop")
    assert isinstance(first_sample.pop, str)
    assert hasattr(first_sample, "sex_karyotype")
    assert isinstance(first_sample.sex_karyotype, str)

    # Sex karyotype values include the expected typed and aneuploid labels.
    sex_counts = dict(sa.aggregate(hl.agg.counter(sa.sex_karyotype)))
    assert sex_counts["XX"] > 0
    assert sex_counts["XY"] > 0
    # The downstream aneuploidy filter in `compute_haplotypes` only has meaningful regression
    # coverage if the fixture exposes at least one non-XX/XY sample. Guard the fixture itself
    # so a future regeneration that silently drops aneuploidies doesn't quietly weaken tests.
    non_canonical_karyotypes = set(sex_counts) - {"XX", "XY"}
    assert non_canonical_karyotypes, (
        "fixture has no aneuploid/ambiguous samples — the compute_haplotypes "
        "aneuploidy filter test loses meaning"
    )
