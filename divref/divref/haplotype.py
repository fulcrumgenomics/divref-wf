"""Shared utilities for Hail-based DivRef pipeline tools."""

from typing import Hashable
from typing import TypeVar

import hail as hl

from divref import defaults

_V = TypeVar("_V", bound=Hashable)
"""Type variable for hashable dictionary values used in to_hashable_items."""


def to_hashable_items(d: dict[str, _V]) -> tuple[tuple[str, _V], ...]:
    """
    Convert a dictionary to a sorted tuple of items for use as a hashable key.

    Args:
        d: Dictionary with hashable values to convert.

    Returns:
        Sorted tuple of (key, value) pairs.
    """
    return tuple(sorted(d.items()))


def get_haplo_sequence(
    context_size: int,
    variants: hl.Expression,
    reference_genome: str = defaults.REFERENCE_GENOME,
) -> hl.Expression:
    """
    Construct a haplotype sequence string with flanking genomic context.

    Builds a sequence by combining alternate alleles from each variant with
    intervening reference sequence, bounded by context_size flanking bases on
    each side.

    Args:
        context_size: Number of reference bases to include flanking each end.
        variants: Hail array expression of variant structs with locus and alleles fields.
            Must contain at least one variant.
        reference_genome: Name of the reference genome. Defaults to "GRCh38".

    Returns:
        Hail string expression representing the full haplotype sequence.

    Raises:
        ValueError: If variants is a Python sequence with no elements.
    """
    if isinstance(variants, (list, tuple)) and len(variants) == 0:
        raise ValueError(
            "get_haplo_sequence requires at least one variant; received an empty sequence"
        )
    sorted_variants = hl.sorted(variants, key=lambda x: x.locus.position)
    min_variant = sorted_variants[0]
    max_variant = sorted_variants[-1]
    min_pos = min_variant.locus.position
    max_pos = max_variant.locus.position
    max_variant_size = hl.len(max_variant.alleles[0])
    full_context = hl.get_sequence(
        min_variant.locus.contig,
        min_pos,
        before=context_size,
        after=(max_pos - min_pos + max_variant_size + context_size - 1),
        reference_genome=reference_genome,
    )

    # (min_pos - index_translation) equals context_size, mapping locus positions to string indices
    index_translation = min_pos - context_size

    def get_chunk_until_next_variant(i: hl.Expression) -> hl.Expression:
        """
        Return the alternate allele plus intervening reference bases up to the next variant.

        Args:
            i: Hail integer expression indexing into sorted_variants.

        Returns:
            Hail string expression for the alternate allele concatenated with the
            reference bases between this variant and the next (or the trailing context).
        """
        v = sorted_variants[i]
        variant_size = hl.len(v.alleles[0])
        reference_buffer_size = hl.if_else(
            i == hl.len(sorted_variants) - 1,
            context_size,
            sorted_variants[i + 1].locus.position - (v.locus.position + variant_size),
        )
        start = v.locus.position - index_translation + variant_size
        return v.alleles[1] + full_context[start : start + reference_buffer_size]

    return full_context[:context_size] + hl.delimit(
        hl.range(hl.len(sorted_variants)).map(get_chunk_until_next_variant),
        "",
    )


def variant_distance(v1: hl.Expression, v2: hl.Expression) -> hl.Expression:
    """
    Calculate the number of reference bases between two variants.

    For example: 1:1:A:T and 1:3:A:T have distance 1 (one base between them).
    1:1:AA:T and 1:3:A:T have distance 0 (deletion closes the gap).

    Args:
        v1: First variant Hail struct with locus and alleles fields.
        v2: Second variant Hail struct with locus and alleles fields.

    Returns:
        Hail int32 expression for the number of reference bases between v1 and v2.
    """
    return v2.locus.position - v1.locus.position - hl.len(v1.alleles[0])


def haplo_coordinates(
    window_size: int,
    variants: hl.Expression,
) -> hl.Expression:
    """
    Compute the 0-based half-open reference genome coordinates of a haplotype sequence window.

    The variant position coordinate is 1-based.

    The window spans from `window_size` bases before the first variant to `window_size` bases
    after the rightmost reference base touched by any variant. `end` uses the maximum reference
    end over all variants, not just the last-by-position variant: an earlier deletion can have a
    reference allele that extends past the last variant (which happens only for haplotypes whose
    variants overlap, i.e. those flagged incompatible). For a haplotype with no overlapping
    variants the maximum reference end is the last variant's, so `end` is unchanged.

    Args:
        window_size: Number of flanking reference bases on each side (same value passed to
            get_haplo_sequence).
        variants: Hail array expression of variant structs with locus and alleles fields.

    Returns:
        Hail struct expression with int32 fields `start` (inclusive) and `end` (exclusive).
    """
    sorted_variants = hl.sorted(variants, key=lambda x: x.locus.position)
    min_variant = sorted_variants[0]
    max_ref_end = hl.max(sorted_variants.map(lambda v: v.locus.position + hl.len(v.alleles[0])))
    return hl.struct(
        start=min_variant.locus.position - 1 - window_size,
        end=max_ref_end - 1 + window_size,
    )
