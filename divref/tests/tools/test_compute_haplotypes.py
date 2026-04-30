"""Tests for the compute_haplotypes tool."""

from dataclasses import dataclass
from pathlib import Path
from typing import NamedTuple
from unittest.mock import patch

import hail as hl
import pytest

from divref.tools.compute_haplotypes import _enumerate_subfragments
from divref.tools.compute_haplotypes import _form_parent_blocks
from divref.tools.compute_haplotypes import compute_haplotypes


class CarrierPos(NamedTuple):
    """Helper class to define a carrier position."""

    position: int
    row_index: int
    ref_len: int


@dataclass
class CarrierInfo:
    """
    Helper class to define carrier info from a sample.

    `left` and `right` can be any order; the function under test sorts internally.
    """

    col_index: int
    pop_int: int
    left: list[CarrierPos]
    right: list[CarrierPos]


def _make_cols_ht(
    samples: list[CarrierInfo],
) -> hl.Table:
    """
    Construct a synthetic cols Table for testing `_form_parent_blocks`.

    Args:
        samples: list of `CarrierInfo`

    Returns:
        Hail Table with the schema expected by `_form_parent_blocks`.
    """
    carrier_type = hl.tstruct(
        locus=hl.tstruct(contig=hl.tstr, position=hl.tint32),
        row_idx=hl.tint64,
        ref_len=hl.tint32,
    )
    row_type = hl.tstruct(
        col_idx=hl.tint32,
        pop_int=hl.tint32,
        left_carriers=hl.tarray(carrier_type),
        right_carriers=hl.tarray(carrier_type),
    )
    rows = [
        {
            "col_idx": s.col_index,
            "pop_int": s.pop_int,
            "left_carriers": [
                {
                    "locus": {"contig": "chr1", "position": c.position},
                    "row_idx": c.row_index,
                    "ref_len": c.ref_len,
                }
                for c in s.left
            ],
            "right_carriers": [
                {
                    "locus": {"contig": "chr1", "position": c.position},
                    "row_idx": c.row_index,
                    "ref_len": c.ref_len,
                }
                for c in s.right
            ],
        }
        for s in samples
    ]
    return hl.Table.parallelize(rows, schema=row_type)


def test_form_parent_blocks_single_block(hail_context: None) -> None:  # noqa: ARG001
    """Three carriers within window_size=25: one parent block emitted, sorted by position."""
    cols_ht = _make_cols_ht([
        CarrierInfo(
            col_index=0,
            pop_int=7,
            left=[CarrierPos(130, 2, 1), CarrierPos(100, 0, 1), CarrierPos(110, 1, 1)],
            right=[],
        )
    ])
    result = sorted(_form_parent_blocks(cols_ht, window_size=25).collect(), key=lambda r: r.strand)
    assert len(result) == 1
    row = result[0]
    assert row.col_idx == 0
    assert row.pop_int == 7
    assert row.strand == 0
    assert [v.row_idx for v in row.parent_block] == [0, 1, 2]
    assert [v.locus.position for v in row.parent_block] == [100, 110, 130]


def test_form_parent_blocks_singleton_dropped(hail_context: None) -> None:  # noqa: ARG001
    """A single carrier on a strand produces no parent block (length-< 2 filter)."""
    cols_ht = _make_cols_ht([
        CarrierInfo(col_index=0, pop_int=0, left=[CarrierPos(100, 0, 1)], right=[])
    ])
    result = _form_parent_blocks(cols_ht, window_size=25).collect()
    assert result == []


def test_form_parent_blocks_splits_at_large_gap(hail_context: None) -> None:  # noqa: ARG001
    """Gap of exactly window_size triggers a split; both halves preserved when length ≥ 2."""
    # gaps 9, 25, 3
    cols_ht = _make_cols_ht([
        CarrierInfo(
            col_index=0,
            pop_int=0,
            left=[
                CarrierPos(100, 0, 1),
                CarrierPos(110, 1, 1),
                CarrierPos(136, 2, 1),
                CarrierPos(140, 3, 1),
            ],
            right=[],
        )
    ])
    rows = sorted(
        _form_parent_blocks(cols_ht, window_size=25).collect(),
        key=lambda r: r.parent_block[0].locus.position,
    )
    assert len(rows) == 2
    assert [v.row_idx for v in rows[0].parent_block] == [0, 1]
    assert [v.row_idx for v in rows[1].parent_block] == [2, 3]


def test_form_parent_blocks_isolates_singleton(hail_context: None) -> None:  # noqa: ARG001
    """A carrier that breaks adjacency on both sides becomes a singleton and is dropped."""
    # gaps 99, 9
    cols_ht = _make_cols_ht([
        CarrierInfo(
            col_index=0,
            pop_int=0,
            left=[CarrierPos(100, 0, 1), CarrierPos(200, 1, 1), CarrierPos(210, 2, 1)],
            right=[],
        )
    ])
    rows = _form_parent_blocks(cols_ht, window_size=25).collect()
    assert len(rows) == 1
    assert [v.row_idx for v in rows[0].parent_block] == [1, 2]


def test_form_parent_blocks_ref_len_closes_gap(hail_context: None) -> None:  # noqa: ARG001
    """Ref-allele length is subtracted from gap (matches `variant_distance` semantics)."""
    # 50 bp deletion at position 100 closes the 50 bp positional gap to V at position 150.
    # gap = 150 - 100 - 50 = 0
    cols_ht = _make_cols_ht([
        CarrierInfo(
            col_index=0, pop_int=0, left=[CarrierPos(100, 0, 50), CarrierPos(150, 1, 1)], right=[]
        )
    ])
    rows = _form_parent_blocks(cols_ht, window_size=25).collect()
    assert len(rows) == 1
    assert [v.row_idx for v in rows[0].parent_block] == [0, 1]


def test_form_parent_blocks_left_and_right_independent(
    hail_context: None,  # noqa: ARG001
) -> None:
    """Left- and right-strand carriers are not interleaved when forming blocks."""
    cols_ht = _make_cols_ht([
        CarrierInfo(
            col_index=5,
            pop_int=2,
            left=[CarrierPos(100, 0, 1), CarrierPos(110, 1, 1)],
            right=[CarrierPos(200, 2, 1), CarrierPos(210, 3, 1)],
        )
    ])
    rows = sorted(_form_parent_blocks(cols_ht, window_size=25).collect(), key=lambda r: r.strand)
    assert len(rows) == 2
    assert rows[0].strand == 0
    assert [v.row_idx for v in rows[0].parent_block] == [0, 1]
    assert rows[1].strand == 1
    assert [v.row_idx for v in rows[1].parent_block] == [2, 3]


@dataclass
class ParentInfo:
    """Helper to define a single parent block for testing `_enumerate_subfragments`."""

    col_index: int
    pop_int: int
    strand: int
    block: list[CarrierPos]


def _make_parents_ht(parents: list[ParentInfo]) -> hl.Table:
    """
    Construct a synthetic parents Hail Table for testing `_enumerate_subfragments`.

    Args:
        parents: list of `ParentInfo` rows; each becomes one (sample, strand, parent block) row.

    Returns:
        Hail Table with the schema produced by `_form_parent_blocks`.
    """
    carrier_type = hl.tstruct(
        locus=hl.tstruct(contig=hl.tstr, position=hl.tint32),
        row_idx=hl.tint64,
        ref_len=hl.tint32,
    )
    row_type = hl.tstruct(
        col_idx=hl.tint32,
        pop_int=hl.tint32,
        strand=hl.tint32,
        parent_block=hl.tarray(carrier_type),
    )
    rows = [
        {
            "col_idx": p.col_index,
            "pop_int": p.pop_int,
            "strand": p.strand,
            "parent_block": [
                {
                    "locus": {"contig": "chr1", "position": c.position},
                    "row_idx": c.row_index,
                    "ref_len": c.ref_len,
                }
                for c in p.block
            ],
        }
        for p in parents
    ]
    return hl.Table.parallelize(rows, schema=row_type)


def test_enumerate_subfragments_length_two(hail_context: None) -> None:  # noqa: ARG001
    """A parent of length 2 emits one sub-fragment (the full block itself)."""
    parents = _make_parents_ht([
        ParentInfo(
            col_index=0,
            pop_int=0,
            strand=0,
            block=[CarrierPos(100, 0, 1), CarrierPos(110, 1, 1)],
        )
    ])
    rows = _enumerate_subfragments(parents).collect()
    assert len(rows) == 1
    assert [v.row_idx for v in rows[0].sub_fragment] == [0, 1]


def test_enumerate_subfragments_length_three(hail_context: None) -> None:  # noqa: ARG001
    """A parent of length 3 emits 3 sub-fragments: [0,1], [0,1,2], [1,2]."""
    parents = _make_parents_ht([
        ParentInfo(
            col_index=0,
            pop_int=0,
            strand=0,
            block=[CarrierPos(100, 0, 1), CarrierPos(110, 1, 1), CarrierPos(120, 2, 1)],
        )
    ])
    rows = _enumerate_subfragments(parents).collect()
    seen = sorted([tuple(v.row_idx for v in r.sub_fragment) for r in rows])
    assert seen == [(0, 1), (0, 1, 2), (1, 2)]


def test_enumerate_subfragments_length_four(hail_context: None) -> None:  # noqa: ARG001
    """A parent of length 4 emits 6 sub-fragments — every contiguous slice of length ≥ 2."""
    parents = _make_parents_ht([
        ParentInfo(
            col_index=0,
            pop_int=0,
            strand=0,
            block=[
                CarrierPos(100, 0, 1),
                CarrierPos(110, 1, 1),
                CarrierPos(120, 2, 1),
                CarrierPos(130, 3, 1),
            ],
        )
    ])
    rows = _enumerate_subfragments(parents).collect()
    seen = sorted([tuple(v.row_idx for v in r.sub_fragment) for r in rows])
    assert seen == [
        (0, 1),
        (0, 1, 2),
        (0, 1, 2, 3),
        (1, 2),
        (1, 2, 3),
        (2, 3),
    ]


def test_enumerate_subfragments_preserves_metadata(hail_context: None) -> None:  # noqa: ARG001
    """Each sub-fragment row inherits col_idx, pop_int, strand from its parent."""
    parents = _make_parents_ht([
        ParentInfo(
            col_index=42,
            pop_int=3,
            strand=1,
            block=[CarrierPos(100, 0, 1), CarrierPos(110, 1, 1), CarrierPos(120, 2, 1)],
        )
    ])
    rows = _enumerate_subfragments(parents).collect()
    assert len(rows) == 3
    assert all(r.col_idx == 42 for r in rows)
    assert all(r.pop_int == 3 for r in rows)
    assert all(r.strand == 1 for r in rows)


def test_enumerate_subfragments_multiple_parents(hail_context: None) -> None:  # noqa: ARG001
    """Sub-fragments from different parents are emitted independently."""
    parents = _make_parents_ht([
        ParentInfo(
            col_index=0,
            pop_int=0,
            strand=0,
            block=[CarrierPos(100, 0, 1), CarrierPos(110, 1, 1)],
        ),
        ParentInfo(
            col_index=1,
            pop_int=0,
            strand=0,
            block=[CarrierPos(200, 5, 1), CarrierPos(210, 6, 1), CarrierPos(220, 7, 1)],
        ),
    ])
    rows = _enumerate_subfragments(parents).collect()
    by_sample: dict[int, list[tuple[int, ...]]] = {}
    for r in rows:
        by_sample.setdefault(r.col_idx, []).append(tuple(v.row_idx for v in r.sub_fragment))
    assert sorted(by_sample[0]) == [(0, 1)]
    assert sorted(by_sample[1]) == [(5, 6), (5, 6, 7), (6, 7)]


def test_form_parent_blocks_multiple_samples(hail_context: None) -> None:  # noqa: ARG001
    """Each sample produces its own parent blocks; shared positions don't cross sample lines."""
    cols_ht = _make_cols_ht([
        CarrierInfo(
            col_index=0,
            pop_int=0,
            left=[CarrierPos(100, 0, 1), CarrierPos(110, 1, 1), CarrierPos(120, 2, 1)],
            right=[],
        ),
        CarrierInfo(
            col_index=1,
            pop_int=1,
            left=[CarrierPos(78, 99, 1), CarrierPos(100, 0, 1), CarrierPos(110, 1, 1)],
            right=[],
        ),
    ])
    rows = sorted(_form_parent_blocks(cols_ht, window_size=25).collect(), key=lambda r: r.col_idx)
    assert len(rows) == 2
    assert rows[0].col_idx == 0
    assert [v.row_idx for v in rows[0].parent_block] == [0, 1, 2]
    assert rows[1].col_idx == 1
    assert [v.row_idx for v in rows[1].parent_block] == [99, 0, 1]


@pytest.mark.parametrize(
    "variant_freq_threshold,haplotype_freq_threshold,expected_count",
    [
        (0.0, 0.0, 517),
        (0.005, 0.0, 295),
        (0.0, 0.005, 30),
        (0.005, 0.005, 33),
    ],
)
def test_compute_haplotypes(
    hail_context: None,  # noqa: ARG001
    datadir: Path,
    tmp_path: Path,
    variant_freq_threshold: float,
    haplotype_freq_threshold: float,
    expected_count: int,
) -> None:
    """Happy-path: compute haplotypes from test VCF with gnomAD annotations."""
    # --- act ---
    in_sites = datadir / "chr1_100001_200000.gnomad_afs.ht"
    in_samples = datadir / "hgdp_1kg_sample_metadata.extract.ht"
    vcf_path = datadir / "chr1_100001_200000.vcf.gz"
    output_base = tmp_path / "haplos"

    with patch("divref.tools.compute_haplotypes.hl.init"):
        compute_haplotypes(
            vcfs_path=vcf_path,
            gnomad_va_file=in_sites,
            gnomad_sa_file=in_samples,
            window_size=5000,
            variant_freq_threshold=variant_freq_threshold,
            haplotype_freq_threshold=haplotype_freq_threshold,
            output_base=output_base,
            temp_dir=tmp_path / "hail_tmp",
        )

    # --- assert ---
    result = hl.read_table(f"{output_base}.ht")
    result_count = result.count()
    assert result_count == expected_count

    results: list[hl.Struct] = result.collect()

    # Each haplotype has >= 2 variants (single-variant haplotypes are filtered)
    assert all(len(r.haplotype) >= 2 for r in results)
    assert all(len(r.variants) == len(r.haplotype) for r in results)
    assert all(len(r.gnomad_freqs) == len(r.haplotype) for r in results)

    # Population frequency summary fields are present
    assert all(hasattr(r, "max_pop") for r in results)
    assert all(hasattr(r, "max_empirical_AF") for r in results)
    assert all(hasattr(r, "max_empirical_AC") for r in results)
    assert all(hasattr(r, "all_pop_freqs") for r in results)
    assert all(len(r.all_pop_freqs) > 0 for r in results)


def test_compute_haplotypes_no_variants(
    hail_context: None,  # noqa: ARG001
    datadir: Path,
    tmp_path: Path,
) -> None:
    """All variants are filtered out."""
    # --- act ---
    in_sites = datadir / "chr1_100001_200000.gnomad_afs.ht"
    in_samples = datadir / "hgdp_1kg_sample_metadata.extract.ht"
    vcf_path = datadir / "chr1_100001_200000.vcf.gz"
    output_base = tmp_path / "haplos"

    with (
        patch("divref.tools.compute_haplotypes.hl.init"),
        pytest.raises(ValueError, match="No variants found with minimum population AF"),
    ):
        compute_haplotypes(
            vcfs_path=vcf_path,
            gnomad_va_file=in_sites,
            gnomad_sa_file=in_samples,
            window_size=5000,
            variant_freq_threshold=1,
            haplotype_freq_threshold=0,
            output_base=output_base,
            temp_dir=tmp_path / "hail_tmp",
        )
