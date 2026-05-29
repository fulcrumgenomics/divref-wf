"""Tests for remap_divref models and helper functions."""

from typing import Any

import pytest

from divref.tools.remap_divref import Haplotype
from divref.tools.remap_divref import ReferenceMapping
from divref.tools.remap_divref import Variant
from divref.tools.remap_divref import _parse_pop_freqs

_DEFAULT_GNOMAD_AFS: dict[str, str] = {
    "afr": "0.1,0.2,0.3",
    "amr": "0.15,0.25,0.35",
    "eas": "0.12,0.22,0.32",
    "nfe": "0.11,0.21,0.31",
    "sas": "0.13,0.23,0.33",
}

_DEFAULT_ESTIMATED_GNOMAD_AF_PER_POP: dict[str, float] = {
    "afr": 0.05,
    "amr": 0.15,
    "eas": 0.08,
    "nfe": 0.10,
    "sas": 0.07,
}


def create_haplotype(
    sequence_id: str = "test_hap",
    sequence: str = "ACGT",
    sequence_length: int = 4,
    n_variants: int = 3,
    variants: str = "1:500:A:T,1:505:C:G,1:510:T:A",
    gnomad_afs: dict[str, str] | None = None,
    estimated_gnomad_af_per_pop: dict[str, float | None] | None = None,
    **kwargs: Any,
) -> Haplotype:
    """
    Create a Haplotype instance with sensible defaults for testing.

    Args:
        sequence_id: Haplotype sequence identifier.
        sequence: Haplotype sequence string.
        sequence_length: Length of the sequence.
        n_variants: Number of variants in the haplotype.
        variants: Comma-separated variant strings (chrom:pos:ref:alt).
        gnomad_afs: Mapping from pop label to comma-delimited per-variant AF string. Defaults to
            a five-population (afr/amr/eas/nfe/sas) dict if not provided.
        estimated_gnomad_af_per_pop: Mapping from pop label to per-pop scalar estimated gnomAD
            AF. Defaults to a five-population dict if not provided.
        **kwargs: Additional fields passed to Haplotype.

    Returns:
        A Haplotype instance for use in tests.
    """
    defaults: dict[str, Any] = {
        "popmax_fraction_phased": 1.0,
        "popmax_empirical_AF": 0.25,
        "popmax_empirical_AC": 1000,
        "popmax_estimated_gnomad_AF": 0.15,
        "max_pop": "amr",
        "source": "test_source",
    }
    defaults.update(kwargs)
    return Haplotype(
        sequence_id=sequence_id,
        sequence=sequence,
        sequence_length=sequence_length,
        n_variants=n_variants,
        variants=variants,
        gnomad_afs=dict(_DEFAULT_GNOMAD_AFS) if gnomad_afs is None else gnomad_afs,
        estimated_gnomad_af_per_pop=(
            dict(_DEFAULT_ESTIMATED_GNOMAD_AF_PER_POP)
            if estimated_gnomad_af_per_pop is None
            else estimated_gnomad_af_per_pop
        ),
        **defaults,
    )


def test_basic_snp_variants_involved() -> None:
    """
    With context_size=10, haplotype pos 0 = reference pos 490.

    Variants at ref positions 500, 505, 510 map to haplotype positions 10, 15, 20.
    Query window [12, 17) overlaps only the second variant (pos 15..16).
    """
    haplotype = create_haplotype(sequence_id="hap1", variants="1:500:A:T,1:505:C:G,1:510:T:A")
    rm: ReferenceMapping = haplotype.reference_mapping(12, 17, 10)

    assert rm.first_variant_index == 1
    assert rm.last_variant_index == 1
    assert rm.population_frequencies == {
        "afr": [0.1, 0.2, 0.3],
        "amr": [0.15, 0.25, 0.35],
        "eas": [0.12, 0.22, 0.32],
        "nfe": [0.11, 0.21, 0.31],
        "sas": [0.13, 0.23, 0.33],
    }


def test_deletion_shifts_coordinates() -> None:
    """
    Second variant is a deletion (CC→G), which shifts subsequent haplotype positions.

    Query [14, 20) spans the deletion and the following SNP.
    """
    haplotype = create_haplotype(sequence_id="hap2", variants="1:500:A:T,1:505:CC:G,1:510:T:A")
    rm = haplotype.reference_mapping(14, 20, 10)

    assert rm.first_variant_index == 1
    assert rm.last_variant_index == 2


def test_insertion_shifts_coordinates() -> None:
    """
    First variant is an insertion (T→TTT), which shifts subsequent haplotype positions.

    Query [15, 18) lands in the expanded haplotype space after the insertion.
    """
    haplotype = create_haplotype(sequence_id="hap3", variants="1:500:T:TTT,1:505:C:G,1:510:T:A")
    rm = haplotype.reference_mapping(15, 18, 10)

    assert rm.first_variant_index == 1
    assert rm.last_variant_index == 1


def test_no_variants_in_range() -> None:
    """Query [0, 5) is entirely in the flanking context before any variant."""
    haplotype = create_haplotype(sequence_id="hap4", variants="1:500:A:T,1:510:C:G")
    rm = haplotype.reference_mapping(0, 5, 10)

    assert rm.first_variant_index is None
    assert rm.last_variant_index is None


def test_complex_multi_indel_mapping() -> None:
    """Three variants: deletion, insertion, deletion. Query [9, 22) spans all three."""
    haplotype = create_haplotype(
        sequence_id="hap5",
        variants="1:500:AT:A,1:505:C:CTT,1:510:GGG:T",
        gnomad_afs={
            "afr": "0.05,0.15,0.25",
            "amr": "0.06,0.16,0.26",
            "eas": "0.07,0.17,0.27",
            "nfe": "0.04,0.14,0.24",
            "sas": "0.03,0.13,0.23",
        },
    )
    rm = haplotype.reference_mapping(9, 22, 10)

    assert rm.first_variant_index == 0
    assert rm.last_variant_index == 2


def test_large_insertion_with_null_frequencies() -> None:
    """
    Single large insertion; null gnomAD frequencies should be parsed as 0.0.

    context_size=25, variant at ref pos 90349349.
    Query [22, 42) starts before the insertion (ref start = 90349349 - 3 = 90349346).
    """
    alt = "TATGCAAGTGTCATCAGATGAATTGATGACATTTTTGTCAAGTTTAAGCACTGAAAGAACAAACCTCTAAATC"
    haplotype = create_haplotype(
        sequence_id="hap6",
        sequence="ACTACTATCTATATCATCTACTACTACTATCATCATCATCAT",
        sequence_length=42,
        n_variants=1,
        variants=f"chr12:90349349:T:{alt}",
        gnomad_afs={
            "afr": "null",
            "amr": "null",
            "eas": "null",
            "nfe": "null",
            "sas": "null",
        },
    )
    rm = haplotype.reference_mapping(22, 42, 25)

    assert rm.first_variant_index == 0
    assert rm.last_variant_index == 0
    assert rm.start == 90349346
    assert rm.end == 90349350
    assert rm.population_frequencies == {
        "afr": [0.0],
        "amr": [0.0],
        "eas": [0.0],
        "nfe": [0.0],
        "sas": [0.0],
    }


# ---------------------------------------------------------------------------
# Variant.render
# ---------------------------------------------------------------------------


def test_variant_render_snp() -> None:
    """Variant.render should format a SNP as 'chr:pos:ref:alt'."""
    assert Variant(chromosome="chr1", position=100, reference="A", alternate="T").render() == (
        "chr1:100:A:T"
    )


def test_variant_render_insertion() -> None:
    """Variant.render should format an insertion as 'chr:pos:ref:alt'."""
    assert Variant(chromosome="chr2", position=200, reference="A", alternate="ATG").render() == (
        "chr2:200:A:ATG"
    )


def test_variant_render_deletion() -> None:
    """Variant.render should format a deletion as 'chr:pos:ref:alt'."""
    assert Variant(chromosome="chrX", position=300, reference="ATG", alternate="A").render() == (
        "chrX:300:ATG:A"
    )


# ---------------------------------------------------------------------------
# Haplotype.parsed_variants and .contig
# ---------------------------------------------------------------------------


def test_parsed_variants_single() -> None:
    """parsed_variants should parse a single variant string into one Variant."""
    hap = create_haplotype(variants="chr1:100:A:T", n_variants=1)
    vs = hap.parsed_variants()
    assert len(vs) == 1
    assert vs[0].chromosome == "chr1"
    assert vs[0].position == 100
    assert vs[0].reference == "A"
    assert vs[0].alternate == "T"


def test_parsed_variants_multiple() -> None:
    """parsed_variants should parse multiple comma-separated variants."""
    hap = create_haplotype(variants="chr1:100:A:T,chr1:200:CC:G", n_variants=2)
    vs = hap.parsed_variants()
    assert len(vs) == 2
    assert vs[1].position == 200
    assert vs[1].reference == "CC"


def test_parsed_variants_cached() -> None:
    """Second call to parsed_variants returns the exact same list object (no re-parsing)."""
    hap = create_haplotype(variants="chr1:100:A:T,chr1:200:C:G", n_variants=2)
    assert hap.parsed_variants() is hap.parsed_variants()


def test_contig_returns_first_variant_chromosome() -> None:
    """Contig should return the chromosome of the first parsed variant."""
    hap = create_haplotype(variants="chr5:100:A:T,chr5:200:C:G", n_variants=2)
    assert hap.contig() == "chr5"


# ---------------------------------------------------------------------------
# ReferenceMapping.variants_involved_str
# ---------------------------------------------------------------------------


def test_variants_involved_str_empty() -> None:
    """variants_involved_str should return an empty string when no variants are involved."""
    rm = ReferenceMapping(
        chromosome="chr1",
        start=100,
        end=200,
        variants_involved=[],
        first_variant_index=None,
        last_variant_index=None,
        population_frequencies={},
    )
    assert rm.variants_involved_str() == ""


def test_variants_involved_str_single() -> None:
    """variants_involved_str should return a single rendered variant string."""
    rm = ReferenceMapping(
        chromosome="chr1",
        start=100,
        end=200,
        variants_involved=[Variant(chromosome="chr1", position=150, reference="A", alternate="T")],
        first_variant_index=0,
        last_variant_index=0,
        population_frequencies={},
    )
    assert rm.variants_involved_str() == "chr1:150:A:T"


def test_variants_involved_str_multiple() -> None:
    """variants_involved_str should return comma-separated rendered variant strings."""
    rm = ReferenceMapping(
        chromosome="chr1",
        start=100,
        end=300,
        variants_involved=[
            Variant(chromosome="chr1", position=150, reference="A", alternate="T"),
            Variant(chromosome="chr1", position=200, reference="CC", alternate="C"),
        ],
        first_variant_index=0,
        last_variant_index=1,
        population_frequencies={},
    )
    assert rm.variants_involved_str() == "chr1:150:A:T,chr1:200:CC:C"


# ---------------------------------------------------------------------------
# _parse_pop_freqs
# ---------------------------------------------------------------------------


def test_parse_pop_freqs_all_floats() -> None:
    """_parse_pop_freqs should parse a comma-separated string of floats."""
    assert _parse_pop_freqs("0.1,0.2,0.3") == [0.1, 0.2, 0.3]


def test_parse_pop_freqs_nulls_become_zero() -> None:
    """_parse_pop_freqs should replace 'null' entries with 0.0."""
    assert _parse_pop_freqs("0.1,null,0.3") == [0.1, 0.0, 0.3]


def test_parse_pop_freqs_single_null() -> None:
    """_parse_pop_freqs should return [0.0] for a single 'null' entry."""
    assert _parse_pop_freqs("null") == [0.0]


def test_parse_pop_freqs_na_treated_as_missing() -> None:
    """'NA' (Hail's default TSV missing token) should parse to 0.0 alongside 'null'."""
    assert _parse_pop_freqs("0.1,NA,null,0.4") == [0.1, 0.0, 0.0, 0.4]


def test_from_row_splits_per_pop_columns() -> None:
    """Haplotype.from_row should pull per-pop columns into the dict fields by pop label."""
    pops_legend = ["afr", "amr", "eas"]
    row: dict[str, Any] = {
        "sequence_id": "row_hap",
        "sequence": "ACGT",
        "sequence_length": 4,
        "n_variants": 1,
        "popmax_fraction_phased": 1.0,
        "popmax_empirical_AF": 0.5,
        "popmax_empirical_AC": 10,
        "popmax_estimated_gnomad_AF": 0.5,
        "max_pop": "afr",
        "variants": "chr1:100:A:T",
        "source": "test",
        "gnomAD_AF_afr": "0.5",
        "gnomAD_AF_amr": "0.3",
        "gnomAD_AF_eas": "NA",
        "estimated_gnomAD_haplotype_AF_afr": 0.5,
        "estimated_gnomAD_haplotype_AF_amr": 0.3,
        "estimated_gnomAD_haplotype_AF_eas": None,
    }
    hap = Haplotype.from_row(row, pops_legend)
    assert hap.gnomad_afs == {"afr": "0.5", "amr": "0.3", "eas": "NA"}
    assert list(hap.gnomad_afs.keys()) == pops_legend
    assert hap.estimated_gnomad_af_per_pop == {"afr": 0.5, "amr": 0.3, "eas": None}
    assert list(hap.estimated_gnomad_af_per_pop.keys()) == pops_legend


# ---------------------------------------------------------------------------
# Error paths: malformed variant strings
# ---------------------------------------------------------------------------


def test_parsed_variants_empty_string_raises() -> None:
    """variants="" yields one empty token which cannot be split into 4 fields."""
    hap = create_haplotype(variants="", n_variants=0)
    with pytest.raises(ValueError):
        hap.parsed_variants()


def test_parsed_variants_too_few_fields_raises() -> None:
    """'chr1:100:A' has only 3 colon-delimited fields; unpacking into 4 raises ValueError."""
    hap = create_haplotype(variants="chr1:100:A", n_variants=1)
    with pytest.raises(ValueError):
        hap.parsed_variants()


def test_contig_malformed_variants_raises() -> None:
    """contig() delegates to parsed_variants(), so a malformed string propagates the error."""
    hap = create_haplotype(variants="", n_variants=0)
    with pytest.raises(ValueError):
        hap.contig()
