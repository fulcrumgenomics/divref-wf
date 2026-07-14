"""Shared utilities for Hail-based DivRef pipeline tools."""

import hail as hl

from divref import defaults


def _max_reference_end(variants: hl.Expression) -> hl.Expression:
    """
    Return one past the rightmost reference base any variant in the haplotype touches.

    Computed as the maximum over variants of `locus.position + len(ref_allele)`. For a haplotype
    whose variants do not overlap this equals the last-by-position variant's reference end, but an
    earlier deletion can have a reference allele that reaches further right. Shared by
    `haplo_coordinates` and `get_haplo_sequence` so the stored `end` coordinate and the emitted
    sequence span always agree.

    Args:
        variants: Hail array expression of variant structs with locus and alleles fields.

    Returns:
        Hail int expression: the maximum `locus.position + len(ref_allele)` over the variants.
    """
    return hl.max(variants.map(lambda v: v.locus.position + hl.len(v.alleles[0])))


def get_haplo_sequence(
    context_size: int,
    variants: hl.Expression,
    reference_genome: str = defaults.REFERENCE_GENOME,
) -> hl.Expression:
    """
    Construct a haplotype sequence string with flanking genomic context.

    Composes the variants onto the reference left to right with a running reference cursor, rather
    than concatenating alternate alleles: each variant contributes the reference between the cursor
    and its position, then its alternate allele, and advances the cursor past its reference allele.
    A variant that starts inside already-consumed reference (its position is before the cursor,
    i.e. it overlaps an earlier variant) contributes only the portion of its alternate allele past
    the cursor. This resolves every overlap to a defined sequence -- e.g. a deletion plus an
    insertion anchored on a deleted base inserts the new bases without re-adding the deleted anchor.
    Most overlaps are incompatibilities flagged elsewhere (a SNP on a deleted base, an insertion in
    a deletion, two deletions, two alleles at one position); the cursor still resolves them to a
    defined (if not meaningful) sequence via that tie-break. A SNP co-located with an indel is the
    one genuinely compatible overlap.

    The window spans `context_size` reference bases before the first variant through `context_size`
    bases after the rightmost reference end across all variants (see `_max_reference_end`), matching
    `haplo_coordinates`.

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
    # Tertiary tiebreak on alt length so a substitution (len(alt) == 1) composes before a pure
    # insertion (len(alt) > 1) at the same position. A remaining tie (same position, ref length,
    # and alt length) falls back to input order, but that requires two alleles at one site, which
    # cannot occur on a single phased haplotype.
    sorted_variants = hl.sorted(
        variants, key=lambda x: (x.locus.position, hl.len(x.alleles[0]), hl.len(x.alleles[1]))
    )
    min_pos = sorted_variants[0].locus.position
    max_ref_end = _max_reference_end(sorted_variants)
    full_context = hl.get_sequence(
        sorted_variants[0].locus.contig,
        min_pos,
        before=context_size,
        after=(max_ref_end - min_pos + context_size - 1),
        reference_genome=reference_genome,
    )
    # A reference position `p` maps to full_context index `p - translation`.
    translation = min_pos - context_size

    def compose(acc: hl.Expression, v: hl.Expression) -> hl.Expression:
        """
        Fold one variant into the (seq, cursor) accumulator using the reference cursor.

        Args:
            acc: Struct with `seq` (sequence so far) and `cursor` (next unconsumed reference pos).
            v: The variant struct to apply.

        Returns:
            The updated `(seq, cursor)` struct.
        """
        cursor = acc.cursor
        v_pos = v.locus.position
        v_ref_len = hl.len(v.alleles[0])
        ref_gap = full_context[cursor - translation : hl.max(cursor, v_pos) - translation]
        overlap = hl.max(0, cursor - v_pos)
        alt_contribution = v.alleles[1][hl.min(overlap, hl.len(v.alleles[1])) :]
        return hl.struct(
            seq=acc.seq + ref_gap + alt_contribution,
            cursor=hl.max(cursor, v_pos + v_ref_len),
        )

    composed = hl.fold(compose, hl.struct(seq="", cursor=min_pos), sorted_variants)
    trailing = full_context[
        composed.cursor - translation : composed.cursor - translation + context_size
    ]
    return full_context[:context_size] + composed.seq + trailing


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
    # Same sort key as get_haplo_sequence so the two stay in lockstep; only the minimum position
    # (`[0].locus.position`) is read here, and `_max_reference_end` is order-independent.
    sorted_variants = hl.sorted(
        variants, key=lambda x: (x.locus.position, hl.len(x.alleles[0]), hl.len(x.alleles[1]))
    )
    min_variant = sorted_variants[0]
    return hl.struct(
        start=min_variant.locus.position - 1 - window_size,
        end=_max_reference_end(sorted_variants) - 1 + window_size,
    )
