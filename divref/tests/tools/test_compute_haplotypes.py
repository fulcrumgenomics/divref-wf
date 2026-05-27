"""Tests for the compute_haplotypes tool."""

from dataclasses import dataclass
from pathlib import Path
from typing import NamedTuple
from unittest.mock import patch

import hail as hl
import pytest

from divref.tools.compute_haplotypes import _aggregate_containment_ac
from divref.tools.compute_haplotypes import _apply_containment_dedup
from divref.tools.compute_haplotypes import _attach_component_info
from divref.tools.compute_haplotypes import _compute_metrics
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
        carriers=hl.tarray(carrier_type),
        strand=hl.tint32,
    )

    def _carrier_dict(c: CarrierPos) -> dict[str, object]:
        return {
            "locus": {"contig": "chr1", "position": c.position},
            "row_idx": c.row_index,
            "ref_len": c.ref_len,
        }

    rows: list[dict[str, object]] = []
    for s in samples:
        for strand, group in ((0, s.left), (1, s.right)):
            rows.append({
                "col_idx": s.col_index,
                "pop_int": s.pop_int,
                "carriers": [_carrier_dict(c) for c in group],
                "strand": strand,
            })
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


@dataclass
class SubFragmentInfo:
    """Helper to define one sub-fragment row for testing `_aggregate_containment_ac`."""

    col_index: int
    pop_int: int
    strand: int
    sub_fragment: list[CarrierPos]


def _make_subfragments_ht(rows: list[SubFragmentInfo]) -> hl.Table:
    """Construct a synthetic sub-fragments Table with the schema of `_enumerate_subfragments`."""
    carrier_type = hl.tstruct(
        locus=hl.tstruct(contig=hl.tstr, position=hl.tint32),
        row_idx=hl.tint64,
        ref_len=hl.tint32,
    )
    row_type = hl.tstruct(
        col_idx=hl.tint32,
        pop_int=hl.tint32,
        strand=hl.tint32,
        sub_fragment=hl.tarray(carrier_type),
    )
    return hl.Table.parallelize(
        [
            {
                "col_idx": r.col_index,
                "pop_int": r.pop_int,
                "strand": r.strand,
                "sub_fragment": [
                    {
                        "locus": {"contig": "chr1", "position": c.position},
                        "row_idx": c.row_index,
                        "ref_len": c.ref_len,
                    }
                    for c in r.sub_fragment
                ],
            }
            for r in rows
        ],
        schema=row_type,
    )


def test_aggregate_containment_ac_single_row(hail_context: None) -> None:  # noqa: ARG001
    """One sub-fragment from one sample in pop 0 → AC vector with 1 at index 0, 0 elsewhere."""
    sf = _make_subfragments_ht([
        SubFragmentInfo(
            col_index=0,
            pop_int=0,
            strand=0,
            sub_fragment=[CarrierPos(100, 0, 1), CarrierPos(110, 1, 1)],
        )
    ])
    rows = _aggregate_containment_ac(sf, n_pops=3).collect()
    assert len(rows) == 1
    assert list(rows[0].haplotype) == [0, 1]
    assert list(rows[0].per_pop_AC) == [1, 0, 0]


def test_aggregate_containment_ac_canonical_shared_subfragment(
    hail_context: None,  # noqa: ARG001
) -> None:
    """
    Canonical regression for AC summation across parents that share a sub-fragment.

    A's parent [V0_id=0,V1_id=1,V2_id=2] and B's parent [V0_id=99,V1_id=0,V2_id=1] both in
    pop 0; the shared sub-fragment [V1_id=0,V2_id=1] appears in both parents and should
    aggregate to AC=2 in pop 0.
    """
    # All sub-fragments enumerable from each parent (per `_enumerate_subfragments`):
    a_subs = [[0, 1], [0, 1, 2], [1, 2]]
    b_subs = [[99, 0], [99, 0, 1], [0, 1]]

    rows_in: list[SubFragmentInfo] = []
    for h in a_subs:
        rows_in.append(
            SubFragmentInfo(
                col_index=0,
                pop_int=0,
                strand=0,
                sub_fragment=[CarrierPos(position=100 + r, row_index=r, ref_len=1) for r in h],
            )
        )
    for h in b_subs:
        rows_in.append(
            SubFragmentInfo(
                col_index=1,
                pop_int=0,
                strand=0,
                sub_fragment=[CarrierPos(position=78 + r, row_index=r, ref_len=1) for r in h],
            )
        )

    result = _aggregate_containment_ac(_make_subfragments_ht(rows_in), n_pops=2).collect()
    by_hap: dict[tuple[int, ...], list[int]] = {
        tuple(r.haplotype): list(r.per_pop_AC) for r in result
    }

    # The shared [V1, V2] (row_idx tuple (0, 1)) must have AC=2 in pop 0.
    assert by_hap[(0, 1)] == [2, 0]
    # Each parent-unique full block has AC=1.
    assert by_hap[(0, 1, 2)] == [1, 0]
    assert by_hap[(99, 0, 1)] == [1, 0]
    # The other proper sub-fragments [V2_id=1, V2_id=2] and [V0_id=99, V1_id=0] each have AC=1.
    assert by_hap[(1, 2)] == [1, 0]
    assert by_hap[(99, 0)] == [1, 0]


def test_aggregate_containment_ac_multi_population(hail_context: None) -> None:  # noqa: ARG001
    """Two samples in different populations contribute to different slots of the AC vector."""
    rows_in = [
        SubFragmentInfo(
            col_index=0,
            pop_int=0,
            strand=0,
            sub_fragment=[CarrierPos(100, 0, 1), CarrierPos(110, 1, 1)],
        ),
        SubFragmentInfo(
            col_index=1,
            pop_int=2,
            strand=1,
            sub_fragment=[CarrierPos(100, 0, 1), CarrierPos(110, 1, 1)],
        ),
    ]
    result = _aggregate_containment_ac(_make_subfragments_ht(rows_in), n_pops=4).collect()
    assert len(result) == 1
    assert list(result[0].haplotype) == [0, 1]
    assert list(result[0].per_pop_AC) == [1, 0, 1, 0]


def _make_hap_table(rows: list[tuple[list[int], list[int]]]) -> hl.Table:
    """
    Construct a synthetic haplotype Table for testing `_attach_component_info`.

    Args:
        rows: list of `(haplotype, per_pop_AC)` tuples.

    Returns:
        Hail Table keyed by `haplotype` with fields `haplotype` (array<int64>) and
        `per_pop_AC` (array<int64>).
    """
    row_type = hl.tstruct(haplotype=hl.tarray(hl.tint64), per_pop_AC=hl.tarray(hl.tint64))
    return hl.Table.parallelize(
        [{"haplotype": h, "per_pop_AC": ac} for (h, ac) in rows], schema=row_type
    ).key_by("haplotype")


def _make_variants_ht(
    variants: list[tuple[int, int, str, str, list[float], dict[int, int]]],
) -> hl.Table:
    """
    Construct a synthetic variants Table mirroring the schema from `mt.rows().select(...)`.

    Args:
        variants: list of `(row_idx, position, ref, alt, freq_per_pop_AF, ANs_by_pop)`. The
            simplified `freq` field stores just AF per population (not the full struct);
            `frequencies_by_pop` is a dict of pop_int → struct(AN=...). These are sufficient
            for the lookup-mechanic tests; `_compute_metrics` is tested separately with
            richer fields.

    Returns:
        Hail Table keyed by (locus, alleles) with fields `row_idx`, `locus`, `alleles`, `freq`,
        and `frequencies_by_pop`.
    """
    freq_type = hl.tstruct(AF=hl.tfloat64)
    fbp_value_type = hl.tstruct(AN=hl.tint32)
    row_type = hl.tstruct(
        locus=hl.tstruct(contig=hl.tstr, position=hl.tint32),
        alleles=hl.tarray(hl.tstr),
        row_idx=hl.tint64,
        freq=hl.tarray(freq_type),
        frequencies_by_pop=hl.tdict(hl.tint32, fbp_value_type),
    )
    rows = [
        {
            "locus": {"contig": "chr1", "position": pos},
            "alleles": [ref, alt],
            "row_idx": row_idx,
            "freq": [{"AF": af} for af in afs],
            "frequencies_by_pop": {p: {"AN": an} for p, an in ans.items()},
        }
        for (row_idx, pos, ref, alt, afs, ans) in variants
    ]
    return hl.Table.parallelize(rows, schema=row_type).key_by("locus", "alleles")


def test_attach_component_info_basic_lookup(hail_context: None) -> None:  # noqa: ARG001
    """Each haplotype row gets variants/gnomad_freqs/frequencies_by_pop in haplotype order."""
    hap_table = _make_hap_table([([0, 2], [1, 0]), ([1], [3, 0])])
    variants_ht = _make_variants_ht([
        (0, 100, "A", "T", [0.1, 0.2], {0: 100, 1: 200}),
        (1, 110, "C", "G", [0.3, 0.4], {0: 102, 1: 198}),
        (2, 120, "G", "A", [0.5, 0.6], {0: 99, 1: 201}),
    ])
    rows = sorted(
        _attach_component_info(hap_table, variants_ht).collect(),
        key=lambda r: list(r.haplotype),
    )

    assert len(rows) == 2
    # haplotype [0, 2] picks up variants at positions 100 and 120
    assert [v.locus.position for v in rows[0].variants] == [100, 120]
    assert [list(v.alleles) for v in rows[0].variants] == [["A", "T"], ["G", "A"]]
    assert [[s.AF for s in r] for r in rows[0].gnomad_freqs] == [[0.1, 0.2], [0.5, 0.6]]
    assert [{p: s.AN for p, s in fbp.items()} for fbp in rows[0].frequencies_by_pop] == [
        {0: 100, 1: 200},
        {0: 99, 1: 201},
    ]
    # haplotype [1] picks up variant at position 110
    assert [v.locus.position for v in rows[1].variants] == [110]


def test_attach_component_info_preserves_per_pop_ac(hail_context: None) -> None:  # noqa: ARG001
    """`per_pop_AC` is left untouched by component lookup."""
    hap_table = _make_hap_table([([0, 1], [3, 5])])
    variants_ht = _make_variants_ht([
        (0, 100, "A", "T", [0.1], {0: 10}),
        (1, 110, "C", "G", [0.2], {0: 12}),
    ])
    rows = _attach_component_info(hap_table, variants_ht).collect()
    assert len(rows) == 1
    assert list(rows[0].per_pop_AC) == [3, 5]


def test_attach_component_info_array_lengths_match_haplotype(
    hail_context: None,  # noqa: ARG001
) -> None:
    """variants/gnomad_freqs/frequencies_by_pop are parallel-length to haplotype."""
    hap_table = _make_hap_table([([0, 1, 2, 3], [1, 0])])
    variants_ht = _make_variants_ht([
        (0, 100, "A", "T", [0.1], {0: 10}),
        (1, 110, "C", "G", [0.2], {0: 11}),
        (2, 120, "G", "A", [0.3], {0: 12}),
        (3, 130, "T", "C", [0.4], {0: 13}),
    ])
    rows = _attach_component_info(hap_table, variants_ht).collect()
    assert len(rows[0].variants) == 4
    assert len(rows[0].gnomad_freqs) == 4
    assert len(rows[0].frequencies_by_pop) == 4


@dataclass
class MetricsRow:
    """Helper to construct one row for testing `_compute_metrics`."""

    haplotype: list[int]
    per_pop_ac: list[int]
    # One inner list per variant, each entry is the gnomAD AF for that variant in that pop.
    variant_pop_af: list[list[float]]
    # One inner dict per variant, mapping pop_int -> AN.
    variant_pop_an: list[dict[int, int]]


def _make_metrics_input_ht(rows: list[MetricsRow]) -> hl.Table:
    """Construct a synthetic hap_table with the schema expected by `_compute_metrics`."""
    freq_type = hl.tstruct(AF=hl.tfloat64)
    fbp_value_type = hl.tstruct(AN=hl.tint32)
    row_type = hl.tstruct(
        haplotype=hl.tarray(hl.tint64),
        per_pop_AC=hl.tarray(hl.tint64),
        variants=hl.tarray(
            hl.tstruct(
                locus=hl.tstruct(contig=hl.tstr, position=hl.tint32),
                alleles=hl.tarray(hl.tstr),
            )
        ),
        gnomad_freqs=hl.tarray(hl.tarray(freq_type)),
        frequencies_by_pop=hl.tarray(hl.tdict(hl.tint32, fbp_value_type)),
    )
    table_rows = [
        {
            "haplotype": r.haplotype,
            "per_pop_AC": r.per_pop_ac,
            "variants": [
                {
                    "locus": {"contig": "chr1", "position": 100 + i},
                    "alleles": ["A", "T"],
                }
                for i in range(len(r.haplotype))
            ],
            "gnomad_freqs": [[{"AF": af} for af in afs] for afs in r.variant_pop_af],
            "frequencies_by_pop": [
                {p: {"AN": an} for p, an in ans.items()} for ans in r.variant_pop_an
            ],
        }
        for r in rows
    ]
    return hl.Table.parallelize(table_rows, schema=row_type)


def test_compute_metrics_single_population(hail_context: None) -> None:  # noqa: ARG001
    """One pop, two variants: AF=AC/min_AN, max_pop=0, fraction_phased computed correctly."""
    ht = _make_metrics_input_ht([
        MetricsRow(
            haplotype=[0, 1],
            per_pop_ac=[3],
            variant_pop_af=[[0.5], [0.4]],
            variant_pop_an=[{0: 100}, {0: 200}],
        )
    ])
    rows = _compute_metrics(ht, n_pops=1).collect()
    assert len(rows) == 1
    row = rows[0]
    assert row.max_pop == 0
    # min_AN over variants for pop 0 = min(100, 200) = 100. AF = 3/100 = 0.03.
    assert row.max_empirical_AF == pytest.approx(0.03)
    assert row.max_empirical_AC == 3
    # min_variant_frequency = min(0.5, 0.4) = 0.4
    assert row.min_variant_frequency == pytest.approx(0.4)
    assert row.fraction_phased == pytest.approx(0.03 / 0.4)
    # estimated_gnomad_AF = min(0.5 * fp, 0.4 * fp) = 0.4 * fp
    assert row.estimated_gnomad_AF == pytest.approx(0.4 * 0.03 / 0.4)


def test_compute_metrics_picks_max_population(hail_context: None) -> None:  # noqa: ARG001
    """max_pop is the population with the highest empirical AF."""
    ht = _make_metrics_input_ht([
        MetricsRow(
            haplotype=[0],
            # pop 0: AC=2, AN=200 -> AF=0.01.  pop 1: AC=5, AN=100 -> AF=0.05.  pop 2: AC=0.
            per_pop_ac=[2, 5, 0],
            variant_pop_af=[[0.3, 0.2, 0.1]],
            variant_pop_an=[{0: 200, 1: 100, 2: 50}],
        )
    ])
    rows = _compute_metrics(ht, n_pops=3).collect()
    assert rows[0].max_pop == 1
    assert rows[0].max_empirical_AC == 5
    assert rows[0].max_empirical_AF == pytest.approx(0.05)
    # all_pop_freqs sorted by AF desc: pop 1 (0.05), pop 0 (0.01), pop 2 (0.0)
    assert [s.pop for s in rows[0].all_pop_freqs] == [1, 0, 2]


def test_compute_metrics_zero_an_yields_missing(hail_context: None) -> None:  # noqa: ARG001
    """Population with AN=0 across all variants gets missing empirical_AF and is not chosen."""
    ht = _make_metrics_input_ht([
        MetricsRow(
            haplotype=[0],
            # pop 0: AC=1, AN=0 -> AF missing.  pop 1: AC=1, AN=10 -> AF=0.1.
            per_pop_ac=[1, 1],
            variant_pop_af=[[0.5, 0.4]],
            variant_pop_an=[{0: 0, 1: 10}],
        )
    ])
    rows = _compute_metrics(ht, n_pops=2).collect()
    assert rows[0].max_pop == 1
    assert rows[0].max_empirical_AF == pytest.approx(0.1)
    # all_pop_freqs: pop 1 first (defined AF), pop 0 last (missing AF sorts to end)
    assert rows[0].all_pop_freqs[0].pop == 1
    assert rows[0].all_pop_freqs[0].empirical_AF == pytest.approx(0.1)
    assert rows[0].all_pop_freqs[1].pop == 0
    assert rows[0].all_pop_freqs[1].empirical_AF is None


def test_compute_metrics_tie_breaks_to_smallest_index(
    hail_context: None,  # noqa: ARG001
) -> None:
    """When two pops tie for max AF, argmax returns the smallest index."""
    ht = _make_metrics_input_ht([
        MetricsRow(
            haplotype=[0],
            per_pop_ac=[1, 1],  # both AC=1 over AN=100 -> AF=0.01 each
            variant_pop_af=[[0.2, 0.2]],
            variant_pop_an=[{0: 100, 1: 100}],
        )
    ])
    rows = _compute_metrics(ht, n_pops=2).collect()
    assert rows[0].max_pop == 0


def test_apply_containment_dedup_canonical_three_rows(
    hail_context: None,  # noqa: ARG001
) -> None:
    """
    Canonical example: A's parent [V0,V1,V2] and B's parent [V99,V0,V1] both pop 0.

    All sub-fragments enumerated with containment AC:
      [V0, V1, V2]    AC=(1, ...) — A only
      [V99, V0, V1]   AC=(1, ...) — B only
      [V0, V1]        AC=(2, ...) — both contain it
      [V1, V2]        AC=(1, ...) — A only, contained in [V0,V1,V2] AC=(1) → dropped
      [V99, V0]       AC=(1, ...) — B only, contained in [V99,V0,V1] AC=(1) → dropped

    Expect 3 rows out: the two full blocks (AC=1 each) and the shared 2-element overlap (AC=2).
    """
    hap_table = _make_hap_table([
        ([0, 1, 2], [1, 0]),  # A's full block
        ([99, 0, 1], [1, 0]),  # B's full block
        ([0, 1], [2, 0]),  # shared overlap (different AC, NOT subsumed)
        ([1, 2], [1, 0]),  # sub of A only — should drop (same AC as longer [0,1,2])
        ([99, 0], [1, 0]),  # sub of B only — should drop (same AC as longer [99,0,1])
    ])
    result = sorted(_apply_containment_dedup(hap_table).collect(), key=lambda r: list(r.haplotype))
    haps = [tuple(r.haplotype) for r in result]
    assert haps == [(0, 1), (0, 1, 2), (99, 0, 1)]
    by_hap: dict[tuple[int, ...], list[int]] = {
        tuple(r.haplotype): list(r.per_pop_AC) for r in result
    }
    assert by_hap[(0, 1, 2)] == [1, 0]
    assert by_hap[(99, 0, 1)] == [1, 0]
    assert by_hap[(0, 1)] == [2, 0]


def test_apply_containment_dedup_identical_full_block(
    hail_context: None,  # noqa: ARG001
) -> None:
    """
    Two samples emit identical full block [V0,V1,V2] (AC=2 in pop 0).

    Sub-fragments [V0,V1] and [V1,V2] also have AC=2 (each is contained in both parents).
    All three have AC=2; the proper sub-fragments are subsumed by the full block.
    Expect a single row out: the full block.
    """
    hap_table = _make_hap_table([
        ([0, 1, 2], [2]),
        ([0, 1], [2]),
        ([1, 2], [2]),
    ])
    result = _apply_containment_dedup(hap_table).collect()
    assert len(result) == 1
    assert list(result[0].haplotype) == [0, 1, 2]


def test_apply_containment_dedup_mixed_ac(hail_context: None) -> None:  # noqa: ARG001
    """
    A's parent [V0,V1,V2,V3] (AC=1) and B,C both with [V0,V1,V2] (AC=2 → 3 with A).

    Containment AC for each unique haplotype:
      [V0,V1,V2,V3]  AC=1 (A only)
      [V0,V1,V2]     AC=3 (A as proper sub + B + C exact)
      Other proper sub-fragments [V0,V1], [V1,V2], [V2,V3], [V1,V2,V3] each have some AC value
        but those subsumed by [V0,V1,V2] (AC=3) are dropped if they share AC=3.

    Expect [V0,V1,V2,V3] AC=1 and [V0,V1,V2] AC=3 to survive.
    """
    hap_table = _make_hap_table([
        ([0, 1, 2, 3], [1]),  # A's full block, unique
        ([0, 1, 2], [3]),  # appears in A (sub) + B (full) + C (full)
        ([0, 1], [3]),  # in A (sub), B (sub), C (sub) — same AC as longer [0,1,2] → drop
        ([1, 2], [3]),  # in A (sub), B (sub), C (sub) — same AC → drop
        ([1, 2, 3], [1]),  # in A only as proper sub. Same AC as [0,1,2,3] → drop.
        ([2, 3], [1]),  # in A only. Same AC as [0,1,2,3] → drop. (Also as longer [1,2,3].)
    ])
    result = sorted(_apply_containment_dedup(hap_table).collect(), key=lambda r: len(r.haplotype))
    haps = [tuple(r.haplotype) for r in result]
    assert haps == [(0, 1, 2), (0, 1, 2, 3)]


def test_apply_containment_dedup_no_op_when_unique_acs(
    hail_context: None,  # noqa: ARG001
) -> None:
    """If every row has a distinct AC vector, nothing is dropped (no AC-equal pairs)."""
    hap_table = _make_hap_table([
        ([0, 1, 2], [1, 0]),
        ([0, 1], [2, 0]),
        ([1, 2], [3, 0]),
    ])
    result = _apply_containment_dedup(hap_table).collect()
    assert len(result) == 3


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
        (0.0, 0.0, 457),
        (0.005, 0.0, 283),
        (0.0, 0.005, 43),
        (0.005, 0.005, 35),
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


@pytest.mark.parametrize(
    "variant_freq_threshold,haplotype_freq_threshold,expected_count",
    [
        (0.0, 0.0, 1195),
        (0.005, 0.0, 884),
        (0.0, 0.005, 796),
        (0.005, 0.005, 670),
    ],
)
def test_compute_haplotypes_chrx_nonpar(
    hail_context: None,  # noqa: ARG001
    datadir: Path,
    tmp_path: Path,
    variant_freq_threshold: float,
    haplotype_freq_threshold: float,
    expected_count: int,
) -> None:
    """Compute haplotypes on a chrX non-PAR window with sex-aware ploidy correction."""
    in_sites = datadir / "chrX_50000000_50100000.gnomad_afs.ht"
    in_samples = datadir / "hgdp_1kg_sample_metadata.extract.ht"
    vcf_path = datadir / "chrX_50000000_50100000.vcf.gz"
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

    result = hl.read_table(f"{output_base}.ht")
    assert result.count() == expected_count

    results: list[hl.Struct] = result.collect()
    assert all(len(r.haplotype) >= 2 for r in results)
    assert all(len(r.variants) == len(r.haplotype) for r in results)


def test_compute_haplotypes_chrx_aneuploid_filter_and_haploid_male_an(
    hail_context: None,  # noqa: ARG001
    datadir: Path,
    tmp_path: Path,
) -> None:
    """
    Aneuploidy/ambiguous samples are excluded; chrX non-PAR males contribute haploid AN.

    With the committed sample metadata fixture (1905 XX + 2216 XY + 30 aneuploid/ambiguous)
    and the population restriction (5 populations: afr, amr, eas, sas, nfe), a fully-
    genotyped chrX non-PAR variant yields a per-variant AN of exactly 5532 summed across
    populations: 2 alleles per pop-assigned XX sample + 1 per pop-assigned XY sample. The
    diploid-everywhere bound (no haploid correction) would be 2 * 4121 = 8242; the all-
    samples-included bound (no aneuploidy filter) would also exceed 5532. A regression in
    either fix would shift this value, so we pin it exactly.
    """
    in_sites = datadir / "chrX_50000000_50100000.gnomad_afs.ht"
    in_samples = datadir / "hgdp_1kg_sample_metadata.extract.ht"
    vcf_path = datadir / "chrX_50000000_50100000.vcf.gz"
    output_base = tmp_path / "haplos"

    with patch("divref.tools.compute_haplotypes.hl.init"):
        compute_haplotypes(
            vcfs_path=vcf_path,
            gnomad_va_file=in_sites,
            gnomad_sa_file=in_samples,
            window_size=5000,
            variant_freq_threshold=0.0,
            haplotype_freq_threshold=0.0,
            output_base=output_base,
            temp_dir=tmp_path / "hail_tmp",
        )

    variants = hl.read_table(f"{output_base}.variants.ht")
    # max-over-variants of (sum-over-pops of AN). With the aneuploidy filter + haploid-male
    # correction, every fully-genotyped variant has exactly this value.
    summed_an_max: int = variants.aggregate(
        hl.agg.max(hl.sum(variants.frequencies_by_pop.values().map(lambda v: v.AN)))
    )
    assert summed_an_max == 5532, (
        f"sum(AN) per variant did not match the expected haploid-male non-PAR total: "
        f"{summed_an_max}"
    )
