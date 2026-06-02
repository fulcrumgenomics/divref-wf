"""Tests for divref.haplotype_compat."""

import pytest

from divref.haplotype_compat import Variant
from divref.haplotype_compat import classify_haplotype
from divref.haplotype_compat import classify_pair
from divref.haplotype_compat import compatibility_flag
from divref.haplotype_compat import end_extends_past_rightmost_variant
from divref.haplotype_compat import parse_variants_string
from divref.haplotype_compat import variant_distance
from divref.haplotype_compat import variant_kind


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
        # Genuine conflicts (cannot co-occur on one chromosome).
        ("chr1:300:AT:A", "chr1:301:T:A", "snp_in_deletion"),  # SNP on a deleted base
        ("chr1:500:AAAG:A", "chr1:501:AAG:A", "overlapping_deletions"),
        ("chr1:700:AT:ATGG", "chr1:701:C:G", "insertion_anchor_conflict"),
        ("chr1:600:AAC:A", "chr1:600:AACAC:A", "same_position_deletion"),  # deletion + deletion
        ("chr1:100:C:CA", "chr1:100:C:CGG", "same_position_insertion"),  # insertion + insertion
        ("chr1:100:C:A", "chr1:100:C:G", "same_position_snp"),  # snp + snp
        # deletion + insertion at one site: reciprocal (insert == delete) vs distinct alleles.
        (
            "chr1:100:CAA:C",
            "chr1:100:C:CAA",
            "same_position_reciprocal_insertion_deletion",
        ),  # deletes "AA", inserts "AA" -> nets to reference
        ("chr1:100:CA:C", "chr1:100:C:CGG", "same_position_insertion_deletion"),  # "A" != "GG"
        # Insertion anchored at a deleted base: contested anchor -> conflict.
        ("chr1:400:TGG:T", "chr1:402:G:GTTTT", "insertion_in_deletion"),
        # Composable overlaps -> None.
        ("chr1:100:C:A", "chr1:100:C:CAA", None),  # snp + insertion at one site
        ("chr1:100:C:A", "chr1:100:CA:C", None),  # snp + deletion at one site
        # Non-overlapping.
        ("chr1:200:A:T", "chr1:210:C:G", None),  # distance >= 0
        ("chr1:1:AA:T", "chr1:3:A:T", None),  # distance == 0: deletion closes gap
    ],
)
def test_classify_pair(a: str, b: str, expected: str | None) -> None:
    """Genuine conflicts map to a reason; compatible/composable pairs return None."""
    assert classify_pair(_v(a), _v(b)) == expected


@pytest.mark.parametrize(
    ("ref", "alt", "expected"),
    [
        ("C", "T", "snp"),
        ("AT", "A", "deletion"),
        ("A", "ATG", "insertion"),
        ("AT", "GC", "mnv"),
    ],
)
def test_variant_kind(ref: str, alt: str, expected: str) -> None:
    """Allele lengths classify a variant as snp / deletion / insertion / mnv."""
    assert variant_kind(ref, alt) == expected


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
    assert set(classify_haplotype(variants)) == {"snp_in_deletion", "same_position_deletion"}


@pytest.mark.parametrize(
    ("variants_str", "expected"),
    [
        ("chr1:100:A:T,chr1:200:C:G", "PASS"),  # clean
        ("chr1:50:A:T", "PASS"),  # single variant
        ("chr1:400:TGG:T,chr1:402:G:GTTTT", "insertion_in_deletion"),  # insertion at a deleted base
        ("chr1:100:C:A,chr1:100:C:CAA", "PASS"),  # composable snp + insertion at one site
        ("chr1:300:AT:A,chr1:301:T:A", "snp_in_deletion"),
        (
            "chr1:300:AT:A,chr1:301:T:A,chr1:600:AAC:A,chr1:600:AACAC:A",
            "same_position_deletion;snp_in_deletion",
        ),
        # An early deletion reaches past the last variant: the end-extension flag joins the reason.
        (
            "chr1:100:AAAAA:A,chr1:102:C:T",
            "end_extends_past_rightmost_variant;snp_in_deletion",
        ),
    ],
)
def test_compatibility_flag(variants_str: str, expected: str) -> None:
    """PASS for compatible/single rows; sorted ';'-joined flags otherwise."""
    assert compatibility_flag(variants_str) == expected


@pytest.mark.parametrize(
    ("variants_str", "expected"),
    [
        # An early deletion (100, ref AAAAA -> end 104) reaches past the rightmost variant (102).
        ("chr1:100:AAAAA:A,chr1:102:C:T", True),
        # The rightmost-by-position variant reaches furthest: no extension.
        ("chr1:100:A:T,chr1:200:GG:G", False),
        # Two same-position deletions: the longer one is the rightmost (ties by ref length), so the
        # end already comes from it -- not an extension past the rightmost variant.
        ("chr1:600:AAC:A,chr1:600:AACAC:A", False),
        ("chr1:100:A:T,chr1:200:C:G", False),  # clean
    ],
)
def test_end_extends_past_rightmost_variant(variants_str: str, expected: bool) -> None:
    """True only when an earlier, longer-reference variant reaches past the rightmost variant."""
    assert end_extends_past_rightmost_variant(parse_variants_string(variants_str)) is expected
