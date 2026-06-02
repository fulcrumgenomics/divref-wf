"""
Variant-compatibility classification for DivRef haplotypes.

A haplotype is *incompatible* when two of its component variants have overlapping reference
spans, so they cannot co-occur on a single chromosome (e.g. a SNP at a base a deletion removes,
or two overlapping repeat-contraction deletions). These come from upstream phasing errors at
tandem repeats. This module provides pure functions to detect and classify such cases.

Used by `append_contig_to_duckdb_index` to compute the per-row `haplotype_filter` flag.
"""

# A component variant: (contig, pos, ref, alt). pos is 1-based.
Variant = tuple[str, int, str, str]

# Reason labels in precedence order (first match wins in classify_pair) and stable output order.
# A SNP co-located with an indel is NOT listed: the SNP substitutes the shared base and the indel
# acts on/after it, so the overlap composes to a well-defined haplotype (see classify_pair). An
# insertion anchored inside a deletion IS a conflict (`insertion_in_deletion`): its anchor base is
# one the deletion removes, so the two make conflicting claims about that base.
REASONS: tuple[str, ...] = (
    "same_position_snp",
    "same_position_deletion",
    "same_position_insertion",
    "same_position_reciprocal_insertion_deletion",
    "same_position_insertion_deletion",
    "same_position_other",
    "snp_in_deletion",
    "overlapping_deletions",
    "insertion_in_deletion",
    "insertion_anchor_conflict",
    "other_overlap",
)

# An extra `haplotype_filter` flag that is NOT a classify_pair reason: set when the haplotype's
# reference end is determined by an earlier, longer-reference variant rather than the rightmost
# one, so the window `end` differs from a rightmost-variant-only calculation (see
# `end_extends_past_rightmost_variant`). compatibility_flag appends it to the reason set; it always
# co-occurs with an overlap reason, since only an overlap lets an earlier variant reach furthest.
END_EXTENDS_FLAG = "end_extends_past_rightmost_variant"


def variant_kind(ref: str, alt: str) -> str:
    """Classify a variant as `snp`, `deletion`, `insertion`, or `mnv` from its alleles."""
    if len(ref) == 1 and len(alt) == 1:
        return "snp"
    if len(ref) > len(alt):
        return "deletion"
    if len(alt) > len(ref):
        return "insertion"
    return "mnv"


def parse_variants_string(s: str) -> list[Variant]:
    """
    Parse a comma-separated `chr:pos:ref:alt` list into variants sorted by (pos, len(ref)).

    Args:
        s: Comma-separated variant tokens, each `contig:position:ref:alt`.

    Returns:
        Variants sorted by ascending position, then by reference-allele length.
    """
    out: list[Variant] = []
    for token in s.split(","):
        contig, pos, ref, alt = token.split(":")
        out.append((contig, int(pos), ref, alt))
    out.sort(key=lambda v: (v[1], len(v[2])))
    return out


def variant_distance(v1: Variant, v2: Variant) -> int:
    """
    Reference bases between two variants; negative when their reference spans overlap.

    Args:
        v1: The earlier (or equal-position) variant.
        v2: The later variant.

    Returns:
        `v2.pos - v1.pos - len(v1.ref)`. Zero means adjacent (a deletion closes the gap).
    """
    return v2[1] - v1[1] - len(v1[2])


def _same_position_reason(v1: Variant, v2: Variant, k1: str, k2: str) -> str:
    """
    Name a genuine same-position conflict by its variant-type pair.

    Args:
        v1: The first variant at the shared position.
        v2: The second variant at the shared position.
        k1: `variant_kind` of `v1`.
        k2: `variant_kind` of `v2`.

    Returns:
        A `same_position_*` reason label. A deletion+insertion pair is `reciprocal` when the
        insertion inserts exactly the bases the deletion removes (applying both nets back to the
        reference); otherwise it is a distinct allele pair. Pairs involving an MNV fall to `other`.
    """
    kinds = {k1, k2}
    if kinds == {"snp"}:
        return "same_position_snp"
    if kinds == {"deletion"}:
        return "same_position_deletion"
    if kinds == {"insertion"}:
        return "same_position_insertion"
    if kinds == {"deletion", "insertion"}:
        ins, dele = (v1, v2) if k1 == "insertion" else (v2, v1)
        inserted = ins[3][len(ins[2]) :]  # insertion alt minus its anchor prefix
        deleted = dele[2][len(dele[3]) :]  # deletion ref minus its anchor prefix
        if inserted == deleted:
            return "same_position_reciprocal_insertion_deletion"
        return "same_position_insertion_deletion"
    return "same_position_other"


def classify_pair(v1: Variant, v2: Variant) -> str | None:
    """
    Classify an adjacent variant pair into an incompatibility reason.

    Returns `None` when the pair is compatible: either it does not overlap (`variant_distance >= 0`)
    or it is a SNP co-located with an indel, which composes to a well-defined haplotype (the SNP
    substitutes the shared base and the indel acts on/after it, so the junction survives). A reason
    is returned for every other overlap, including an insertion anchored inside a deletion, whose
    anchor base the deletion removes (a contested base, resolvable only via a tie-break).

    Args:
        v1: The earlier variant by position.
        v2: The later variant by position.

    Returns:
        A reason label from `REASONS`, or `None` if the pair is compatible/composable.
    """
    if variant_distance(v1, v2) >= 0:
        return None
    pos1, pos2 = v1[1], v2[1]
    k1 = variant_kind(v1[2], v1[3])
    k2 = variant_kind(v2[2], v2[3])
    if pos1 == pos2:
        # A SNP changes the shared first base while a co-located indel acts on/after it, so the
        # two compose to a well-defined haplotype. Any other same-position pair is two alleles at
        # one site and cannot co-occur on one chromosome.
        if "snp" in (k1, k2) and ("insertion" in (k1, k2) or "deletion" in (k1, k2)):
            return None
        return _same_position_reason(v1, v2, k1, k2)
    # pos1 < pos2 and the spans overlap: v2 starts inside v1's reference allele.
    if k1 == "deletion":
        if k2 == "snp":
            return "snp_in_deletion"  # SNP changes a base the deletion removes
        if k2 == "deletion":
            return "overlapping_deletions"
        if k2 == "insertion":
            # The insertion's anchor base is one the deletion removes: conflicting claims, and the
            # insert junction is interior to the deleted span (composed only via a tie-break).
            return "insertion_in_deletion"
        return "other_overlap"  # MNV inside a deletion (not observed in practice)
    if k1 == "insertion":
        return "insertion_anchor_conflict"
    return "other_overlap"


def classify_haplotype(variants: list[Variant]) -> list[str]:
    """
    Return the incompatibility reasons from every overlapping consecutive pair.

    Adjacent pairs are sufficient: positions are non-decreasing and each reference allele is
    >= 1 bp, so if no consecutive pair overlaps then no pair overlaps at all. A long deletion
    swallowing a non-adjacent downstream variant is still flagged at its own adjacent boundary.

    Args:
        variants: Position-sorted component variants (e.g. from `parse_variants_string`).

    Returns:
        One reason per incompatible adjacent pair (may repeat reasons or be empty).
    """
    reasons: list[str] = []
    for v1, v2 in zip(variants, variants[1:], strict=False):
        reason = classify_pair(v1, v2)
        if reason is not None:
            reasons.append(reason)
    return reasons


def end_extends_past_rightmost_variant(variants: list[Variant]) -> bool:
    """
    Whether the haplotype's reference end is set by a variant other than the rightmost one.

    The rightmost variant is the one with the largest position, ties broken by the longest
    reference allele -- the variant a naive `end` runs to. An earlier long deletion can reach
    further right than that, making the true end (the maximum reference end over all variants)
    larger. This is True exactly for the haplotypes whose `end` a rightmost-variant-only window
    would truncate; it implies an overlap, so it never occurs on an otherwise-compatible haplotype.

    Args:
        variants: Component variants of the haplotype.

    Returns:
        True iff some variant's reference reaches past the rightmost (largest-position, then
        longest-reference) variant's reference end.
    """
    rightmost = max(variants, key=lambda v: (v[1], len(v[2])))
    rightmost_ref_end = rightmost[1] - 1 + len(rightmost[2])
    max_ref_end = max(v[1] - 1 + len(v[2]) for v in variants)
    return max_ref_end > rightmost_ref_end


def compatibility_flag(variants_str: str) -> str:
    """
    Compute the VCF-style `haplotype_filter` value for a haplotype's `variants` string.

    Args:
        variants_str: Comma-separated `chr:pos:ref:alt` component variants.

    Returns:
        `"PASS"` if the haplotype has no flags, else the `;`-joined sorted flags: the distinct
        incompatibility reasons plus `end_extends_past_rightmost_variant` when an earlier variant's
        reference reaches past the rightmost one (so the `end` differs from a rightmost-only
        window). Single-variant rows are always `"PASS"`.
    """
    variants = parse_variants_string(variants_str)
    flags = set(classify_haplotype(variants))
    if end_extends_past_rightmost_variant(variants):
        flags.add(END_EXTENDS_FLAG)
    return ";".join(sorted(flags)) if flags else "PASS"
