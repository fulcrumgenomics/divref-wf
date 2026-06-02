"""Tests for shared Hail utilities in haplotype.py."""

from pathlib import Path
from typing import Any
from unittest.mock import patch

import hail as hl
import pytest

from divref.haplotype import get_haplo_sequence
from divref.haplotype import haplo_coordinates
from divref.haplotype import to_hashable_items
from divref.haplotype import variant_distance

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------


def _make_variant(position: int, ref: str, alt: str, contig: str = "chr1") -> hl.Struct:
    return hl.Struct(locus=hl.Struct(contig=contig, position=position), alleles=[ref, alt])


def _make_haplotype_table(variant_positions: list[tuple[str, int, str, str]]) -> hl.Table:
    variant_type = hl.tstruct(
        locus=hl.tstruct(contig=hl.tstr, position=hl.tint32), alleles=hl.tarray(hl.tstr)
    )
    row_type = hl.tstruct(
        variants=hl.tarray(variant_type),
        haplotype=hl.tarray(hl.tstr),
        gnomad_freqs=hl.tarray(hl.tfloat64),
    )
    variants = [
        {"locus": {"contig": contig, "position": pos}, "alleles": [ref, alt]}
        for contig, pos, ref, alt in variant_positions
    ]
    return hl.Table.parallelize(
        [
            {
                "variants": variants,
                "haplotype": [str(i) for i in range(len(variants))],
                "gnomad_freqs": [0.1] * len(variants),
            }
        ],
        schema=row_type,
    )


# ---------------------------------------------------------------------------
# get_haplo_sequence
# ---------------------------------------------------------------------------


def _create_reference_mock(reference_sequence: str) -> Any:
    """
    Create a mock for hl.get_sequence that returns substrings of a fixed reference.

    The mock accepts the same arguments as hl.get_sequence and returns the
    appropriate substring of the provided reference string.

    Args:
        reference_sequence: The reference string to use for subsequence extraction.

    Returns:
        A callable that mimics hl.get_sequence using the provided reference.
    """

    def mock_get_sequence(
        _contig: str,
        position: int,
        before: int = 0,
        after: int = 0,
        reference_genome: Any = None,  # noqa: ARG001
    ) -> Any:
        return hl.str(reference_sequence)[position - before : position + after + 1]

    return mock_get_sequence


def test_get_haplo_sequence_edge_cases(hail_context: None) -> None:  # noqa: ARG001
    """Test get_haplo_sequence with SNPs, insertions, and deletions."""
    reference = "01234567891"

    two_snps = [_make_variant(4, "A", "T"), _make_variant(6, "G", "C")]
    two_insertions = [_make_variant(4, "A", "AT"), _make_variant(6, "G", "GC")]
    two_deletions = [_make_variant(4, "AT", "A"), _make_variant(7, "GC", "G")]

    mock_get_sequence = _create_reference_mock(reference)

    with patch("hail.get_sequence", side_effect=mock_get_sequence):
        assert hl.eval(get_haplo_sequence(context_size=2, variants=two_snps)) == "23T5C78"
        assert hl.eval(get_haplo_sequence(context_size=2, variants=two_insertions)) == "23AT5GC78"
        assert hl.eval(get_haplo_sequence(context_size=2, variants=two_deletions)) == "23A6G91"


def test_get_haplo_sequence_deletion_consumes_interior_variant(
    hail_context: None,  # noqa: ARG001
) -> None:
    """
    A variant inside a deletion is consumed, and trailing context reaches the max reference end.

    A deletion at 4 (ref len 5, spans 4-8, alt `D`) plus a SNP at 6 inside it: the SNP's single alt
    base falls entirely within the consumed deletion, so it contributes nothing; the cursor reaches
    the deletion's end (9) and the trailing context runs to 9 + context. Result `23` + `D` + ref
    `9A`. (This pair stays flagged `snp_in_deletion`; the sequence is merely well-defined.)
    """
    reference = "0123456789A"
    variants = [_make_variant(4, "AAAAA", "D"), _make_variant(6, "C", "T")]
    with patch("hail.get_sequence", side_effect=_create_reference_mock(reference)):
        assert hl.eval(get_haplo_sequence(context_size=2, variants=variants)) == "23D9A"


@pytest.mark.parametrize(
    ("specs", "expected", "note"),
    [
        # Composable overlaps: the alt is composed onto the reference, not concatenated.
        ([(4, "AAA", "D"), (6, "G", "GTT")], "23DTT78", "insertion composes inside a deletion"),
        ([(4, "X", "Y"), (4, "X", "XZZ")], "23YZZ56", "snp + insertion at one site"),
        ([(4, "X", "Y"), (4, "XBC", "X")], "23Y78", "snp + deletion at one site"),
        # Genuinely-incompatible overlaps still resolve to a defined (flagged) sequence.
        ([(4, "AAA", "A"), (4, "AAAA", "A")], "23A89", "two deletions: the longer one wins"),
        ([(4, "A", "ABB"), (4, "A", "ACC")], "23ABBCC56", "two insertions: both are applied"),
        ([(4, "A", "T"), (4, "A", "G")], "23T56", "two snps: the first one wins"),
        ([(4, "A", "ATT"), (4, "AB", "A")], "23ATT67", "insertion + deletion: insert then skip"),
    ],
)
def test_get_haplo_sequence_overlap_composition(
    specs: list[tuple[int, str, str]],
    expected: str,
    note: str,
    hail_context: None,  # noqa: ARG001
) -> None:
    """Overlapping variants are composed onto the reference with a cursor, not concatenated."""
    reference = "0123456789A"
    variants = [_make_variant(pos, ref, alt) for pos, ref, alt in specs]
    with patch("hail.get_sequence", side_effect=_create_reference_mock(reference)):
        assert hl.eval(get_haplo_sequence(context_size=2, variants=variants)) == expected, note


def test_get_haplo_sequence_single(
    datadir: Path,
    hail_reference_genome: hl.ReferenceGenome,
    hail_context: None,  # noqa: ARG001
) -> None:
    """get_haplo_sequence should return the correct haplotype sequence."""
    test_fasta: Path = datadir / "test.fa"
    test_fai: Path = datadir / "test.fa.fai"

    hail_reference_genome.add_sequence(str(test_fasta), str(test_fai))

    variant: hl.Struct = _make_variant(position=100, ref="A", alt="C")
    haplo_seq = get_haplo_sequence(
        context_size=2, variants=[variant], reference_genome=hail_reference_genome.name
    )
    assert hl.eval(haplo_seq) == "CCCTC"


def test_get_haplo_sequence_invalid_reference_genome_raises(
    hail_context: None,  # noqa: ARG001
) -> None:
    """get_haplo_sequence should raise KeyError when given an unregistered reference genome."""
    variant = _make_variant(position=100, ref="A", alt="C")
    with pytest.raises(KeyError, match="nonexistent_genome"):
        get_haplo_sequence(
            context_size=2, variants=[variant], reference_genome="nonexistent_genome"
        )


def test_get_haplo_sequence_empty_list_raises() -> None:
    """get_haplo_sequence should raise ValueError when given an empty list."""
    with pytest.raises(ValueError, match="at least one variant"):
        get_haplo_sequence(context_size=2, variants=[])


def test_get_haplo_sequence_empty_tuple_raises() -> None:
    """get_haplo_sequence should raise ValueError when given an empty tuple."""
    with pytest.raises(ValueError, match="at least one variant"):
        get_haplo_sequence(context_size=2, variants=())


# ---------------------------------------------------------------------------
# to_hashable_items
# ---------------------------------------------------------------------------


def test_to_hashable_items_empty() -> None:
    """to_hashable_items should return an empty tuple for an empty dict."""
    assert to_hashable_items({}) == ()


def test_to_hashable_items_single_entry() -> None:
    """to_hashable_items should return a one-element tuple for a single-entry dict."""
    assert to_hashable_items({"key": "value"}) == (("key", "value"),)


def test_to_hashable_items_sorted_by_key() -> None:
    """to_hashable_items should return items sorted by key regardless of insertion order."""
    assert to_hashable_items({"b": 2, "a": 1, "c": 3}) == (("a", 1), ("b", 2), ("c", 3))


# ---------------------------------------------------------------------------
# variant_distance
# ---------------------------------------------------------------------------


def test_variant_distance_adjacent_snps(hail_context: None) -> None:  # noqa: ARG001
    """SNP at 100, next SNP at 101: distance = 101 - 100 - len("A") = 0."""
    assert (
        hl.eval(variant_distance(_make_variant(100, "A", "T"), _make_variant(101, "C", "G"))) == 0
    )


def test_variant_distance_snps_with_gap(hail_context: None) -> None:  # noqa: ARG001
    """SNP at 100, next SNP at 103: 2 reference bases separate them."""
    assert (
        hl.eval(variant_distance(_make_variant(100, "A", "T"), _make_variant(103, "C", "G"))) == 2
    )


def test_variant_distance_deletion_closes_gap(hail_context: None) -> None:  # noqa: ARG001
    """Deletion AT→A at 100 (consumes 2 ref bases), next variant at 102: distance = 0."""
    assert (
        hl.eval(variant_distance(_make_variant(100, "AT", "A"), _make_variant(102, "C", "G"))) == 0
    )


# ---------------------------------------------------------------------------
# haplo_coordinates
# ---------------------------------------------------------------------------


def test_haplo_coordinates_snp(hail_context: None) -> None:  # noqa: ARG001
    """Single SNP at position P with window W: start = P - W, end = P + 1 + W."""
    variants = hl.array([_make_variant(100, "A", "T")])
    coords = hl.eval(haplo_coordinates(10, variants))
    assert coords.start == 89
    assert coords.end == 110


def test_haplo_coordinates_insertion(hail_context: None) -> None:  # noqa: ARG001
    """Insertion (ref len 1) has the same start/end as a SNP at the same position."""
    variants = hl.array([_make_variant(100, "A", "ACGT")])
    coords = hl.eval(haplo_coordinates(10, variants))
    assert coords.start == 89
    assert coords.end == 110


def test_haplo_coordinates_deletion(hail_context: None) -> None:  # noqa: ARG001
    """Deletion with ref len 4 at position 100 with window 10: end = 100 - 1 + 4 + 10 = 113."""
    variants = hl.array([_make_variant(100, "ACGT", "A")])
    coords = hl.eval(haplo_coordinates(10, variants))
    assert coords.start == 89
    assert coords.end == 113


def test_haplo_coordinates_multi_variant(hail_context: None) -> None:  # noqa: ARG001
    """Start uses first variant; end uses last variant + ref allele length + window."""
    variants = hl.array([
        _make_variant(100, "A", "T"),
        _make_variant(200, "GG", "G"),
    ])
    coords = hl.eval(haplo_coordinates(10, variants))
    assert coords.start == 89  # 100 - 1 - 10
    assert coords.end == 211  # 200 -1 + 2 + 10


def test_haplo_coordinates_early_deletion_extends_past_last_variant(
    hail_context: None,  # noqa: ARG001
) -> None:
    """
    `end` uses the rightmost reference end over all variants, not the last-by-position one.

    A deletion at 100 spans 100-104 (ref len 5); the SNP at 102 is last by position but reaches
    only to 102, so `end` must use the deletion's reach (ref end 105).
    """
    variants = hl.array([
        _make_variant(100, "AAAAA", "A"),
        _make_variant(102, "C", "T"),
    ])
    coords = hl.eval(haplo_coordinates(10, variants))
    assert coords.start == 89  # 100 - 1 - 10
    assert coords.end == 114  # max_ref_end (100 + 5) - 1 + 10


def test_haplo_coordinates_matches_sequence_length(hail_context: None) -> None:  # noqa: ARG001
    """For a SNP-only haplotype, end - start must equal len(sequence)."""
    reference = "0" * 300

    mock_get_sequence = _create_reference_mock(reference)

    with patch("hail.get_sequence", side_effect=mock_get_sequence):
        window = 10
        variants = hl.array([_make_variant(150, "A", "T")])
        seq = hl.eval(
            get_haplo_sequence(context_size=window, variants=[_make_variant(150, "A", "T")])
        )
        coords = hl.eval(haplo_coordinates(window, variants))
        assert coords.end - coords.start == len(seq)
