"""Tests for divref.haplotype_compat."""

import pytest

from divref.haplotype_compat import Variant
from divref.haplotype_compat import classify_haplotype
from divref.haplotype_compat import classify_pair
from divref.haplotype_compat import compatibility_flag
from divref.haplotype_compat import count_bypass_resolutions
from divref.haplotype_compat import end_coordinate_shortfall
from divref.haplotype_compat import parse_variants_string
from divref.haplotype_compat import start_coordinate_shortfall
from divref.haplotype_compat import variant_distance
from divref.haplotype_compat import variants_overlap


def _v(token: str) -> Variant:
    contig, pos, ref, alt = token.split(":")
    return (contig, int(pos), ref, alt)


def test_parse_variants_string_sorts_and_roundtrips() -> None:
    """Tokens are parsed and ordered by (position, reference length)."""
    parsed = parse_variants_string("chr1:1744031:AAC:A,chr1:1744033:C:A")
    assert parsed == [("chr1", 1744031, "AAC", "A"), ("chr1", 1744033, "C", "A")]


def test_parse_variants_string_sorts_same_position_by_ref_length() -> None:
    """Same-position variants order by ascending reference length."""
    parsed = parse_variants_string("chr1:600:AACAC:A,chr1:600:AAC:A")
    assert [len(v[2]) for v in parsed] == [3, 5]


@pytest.mark.parametrize(
    ("a", "b", "expected"),
    [
        ("chr1:1744031:AAC:A", "chr1:1744033:C:A", -1),
        ("chr1:1:AA:T", "chr1:3:A:T", 0),
        ("chr1:1:A:T", "chr1:3:A:T", 1),
    ],
)
def test_variant_distance(a: str, b: str, expected: int) -> None:
    """Distance is negative on overlap, zero when a deletion closes the gap, positive otherwise."""
    assert variant_distance(_v(a), _v(b)) == expected


@pytest.mark.parametrize(
    ("a", "b", "expected"),
    [
        ("chr1:300:AT:A", "chr1:301:T:A", "snp_in_deletion"),
        ("chr1:400:TGG:T", "chr1:402:G:GTTTT", "indel_in_deletion"),
        ("chr1:500:AAAG:A", "chr1:501:AAG:A", "overlapping_deletions"),
        ("chr1:600:AAC:A", "chr1:600:AACAC:A", "same_position"),
        ("chr1:700:AT:ATGG", "chr1:701:C:G", "insertion_anchor_conflict"),
        ("chr1:200:A:T", "chr1:210:C:G", None),  # distance >= 0: compatible
        ("chr1:1:AA:T", "chr1:3:A:T", None),  # distance == 0: deletion closes gap, keep
    ],
)
def test_classify_pair(a: str, b: str, expected: str | None) -> None:
    """Each overlap shape maps to its reason; compatible pairs return None."""
    assert classify_pair(_v(a), _v(b)) == expected


def test_classify_haplotype_clean() -> None:
    """A haplotype with no overlapping pairs has no reasons."""
    variants = parse_variants_string("chr2:1000:A:T,chr2:1005:C:G,chr2:1010:G:A")
    assert classify_haplotype(variants) == []


def test_classify_haplotype_single_reason() -> None:
    """A single overlapping pair yields one reason."""
    variants = parse_variants_string("chr1:300:AT:A,chr1:301:T:A")
    assert classify_haplotype(variants) == ["snp_in_deletion"]


def test_classify_haplotype_catches_nonadjacent_overlap() -> None:
    """A long deletion swallowing a non-adjacent variant is still flagged (adjacency suffices)."""
    variants = parse_variants_string("chr1:100:AAAAAAAAAAAAAAAAAAAA:A,chr1:105:C:T,chr1:110:G:A")
    assert classify_haplotype(variants)


def test_classify_haplotype_multiple_reasons() -> None:
    """Distinct overlap shapes in one haplotype are all reported."""
    variants = parse_variants_string("chr1:300:AT:A,chr1:301:T:A,chr1:600:AAC:A,chr1:600:AACAC:A")
    assert set(classify_haplotype(variants)) == {"snp_in_deletion", "same_position"}


@pytest.mark.parametrize(
    ("variants_str", "expected"),
    [
        ("chr1:100:A:T,chr1:200:C:G", "PASS"),  # clean
        ("chr1:50:A:T", "PASS"),  # single variant
        ("chr1:300:AT:A,chr1:301:T:A", "snp_in_deletion"),
        (
            "chr1:300:AT:A,chr1:301:T:A,chr1:600:AAC:A,chr1:600:AACAC:A",
            "same_position;snp_in_deletion",
        ),
    ],
)
def test_compatibility_flag(variants_str: str, expected: str) -> None:
    """PASS for compatible/single rows; sorted ';'-joined reasons otherwise."""
    assert compatibility_flag(variants_str) == expected


@pytest.mark.parametrize(
    ("a", "b", "expected"),
    [
        ("chr1:300:AT:A", "chr1:301:T:A", True),
        ("chr1:600:AAC:A", "chr1:600:AACAC:A", True),
        ("chr1:200:A:T", "chr1:210:C:G", False),
        ("chr1:1:AA:T", "chr1:3:A:T", False),
    ],
)
def test_variants_overlap(a: str, b: str, expected: bool) -> None:
    """Overlap detection is order-independent and excludes distance-0 touching pairs."""
    assert variants_overlap(_v(a), _v(b)) is expected
    assert variants_overlap(_v(b), _v(a)) is expected


def test_count_bypass_resolutions_length2_conflict_recovers_nothing() -> None:
    """A length-2 conflict has no >= 2-variant resolution."""
    assert count_bypass_resolutions(parse_variants_string("chr1:300:AT:A,chr1:301:T:A")) == 0


def test_count_bypass_resolutions_clean_is_one() -> None:
    """A clean haplotype is its own single resolution."""
    assert count_bypass_resolutions(parse_variants_string("chr1:100:A:T,chr1:200:C:G")) == 1


def test_count_bypass_resolutions_single_internal_conflict_is_two() -> None:
    """A single conflicting pair between clean flanks splits two ways."""
    variants = parse_variants_string("chr1:100:A:T,chr1:300:AT:A,chr1:301:T:A,chr1:400:G:C")
    assert count_bypass_resolutions(variants) == 2


def test_count_bypass_resolutions_three_mutually_exclusive_alleles() -> None:
    """Three mutually overlapping deletions plus a clean flank give three resolutions."""
    variants = parse_variants_string(
        "chr1:100:A:T,chr1:500:AAAAAAAA:A,chr1:502:AAAA:A,chr1:504:AA:A"
    )
    assert count_bypass_resolutions(variants) == 3


def test_end_coordinate_shortfall_undershoots_on_early_deletion() -> None:
    """An early deletion reaching past the last variant makes the stored end too short."""
    variants = parse_variants_string("chr1:100:AAAAA:A,chr1:102:C:T")
    assert end_coordinate_shortfall(variants, window_size=25, stored_end=127) == 2


def test_end_coordinate_shortfall_zero_when_last_variant_reaches_furthest() -> None:
    """No shortfall when the last-by-position variant reaches furthest right."""
    variants = parse_variants_string("chr10:112867606:AT:A,chr10:112867607:T:A")
    stored_end = 112867607 + 25
    assert end_coordinate_shortfall(variants, window_size=25, stored_end=stored_end) == 0


def test_start_coordinate_shortfall_zero_for_correct_left_edge() -> None:
    """The left edge is structurally correct, so the start shortfall is zero."""
    variants = parse_variants_string("chr1:100:AAAAA:A,chr1:102:C:T")
    assert start_coordinate_shortfall(variants, window_size=25, stored_start=74) == 0
