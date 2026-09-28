"""Tests for the extract_gnomad_afs tool."""

from pathlib import Path
from unittest.mock import patch

import hail as hl

from divref import defaults
from divref.tools.extract_gnomad_afs import extract_gnomad_afs


def test_extract_gnomad_afs(
    hail_context: None,  # noqa: ARG001
    datadir: Path,
    tmp_path: Path,
) -> None:
    """Happy-path: extract AFs from local test variant and sample data."""
    in_sites = str(datadir / "chr1_100001_200000.ht")
    out_va = tmp_path / "va.ht"

    with patch("divref.tools.extract_gnomad_afs.hail_init"):
        extract_gnomad_afs(
            in_gnomad_sites_table=in_sites,
            out_variant_annotation_table=out_va,
            contig="chr1",
            freq_threshold=0.001,
            populations=defaults.POPULATIONS,
            reference_genome=defaults.REFERENCE_GENOME,
        )

    va = hl.read_table(str(out_va))
    va_count = va.count()
    assert va_count == 966

    # Globals should contain the requested populations
    pops = hl.eval(va.globals.pops)
    assert pops == defaults.POPULATIONS

    # Each row should have pop_freqs with one entry per population
    first_row = va.head(1).collect()[0]
    assert len(first_row.pop_freqs) == len(defaults.POPULATIONS)
    assert all(hasattr(f, "AF") for f in first_row.pop_freqs)

    # All retained variants must exceed freq_threshold in at least one population
    min_max_af = va.aggregate(hl.agg.min(hl.max(va.pop_freqs.map(lambda x: x.AF))))
    assert min_max_af >= 0.001


def test_extract_gnomad_afs_records_prefilter_threshold(
    hail_context: None,  # noqa: ARG001
    datadir: Path,
    tmp_path: Path,
) -> None:
    """The pre-filter threshold is recorded for `compute_haplotypes` to pass on."""
    out_va = tmp_path / "va.ht"
    with patch("divref.tools.extract_gnomad_afs.hail_init"):
        extract_gnomad_afs(
            in_gnomad_sites_table=str(datadir / "chr1_100001_200000.ht"),
            out_variant_annotation_table=out_va,
            contig="chr1",
            freq_threshold=0.01,
            populations=defaults.POPULATIONS,
            reference_genome=defaults.REFERENCE_GENOME,
        )
    recorded = hl.eval(hl.read_table(str(out_va)).index_globals().prefilter_parameters)
    assert recorded == hl.Struct(freq_threshold=0.01)
