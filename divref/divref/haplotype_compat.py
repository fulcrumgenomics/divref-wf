"""
Variant-compatibility classification for DivRef haplotypes.

A haplotype is *incompatible* when two of its component variants have overlapping reference
spans, so they cannot co-occur on a single chromosome (e.g. a SNP at a base a deletion removes,
or two overlapping repeat-contraction deletions). These come from upstream phasing errors at
tandem repeats. This module provides pure functions to detect and classify such cases.

Shared by `append_contig_to_duckdb_index` (to compute the per-row `haplotype_filter` flag) and by
`scripts/evaluate_haplotype_incompatibility.py` (to audit a built index).
"""

# A component variant: (contig, pos, ref, alt). pos is 1-based.
Variant = tuple[str, int, str, str]

# Reason labels in precedence order (first match wins in classify_pair) and stable output order.
REASONS: tuple[str, ...] = (
    "same_position",
    "snp_in_deletion",
    "overlapping_deletions",
    "indel_in_deletion",
    "insertion_anchor_conflict",
    "other_overlap",
)


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


def classify_pair(v1: Variant, v2: Variant) -> str | None:
    """
    Classify an adjacent variant pair into an incompatibility reason.

    Args:
        v1: The earlier variant by position.
        v2: The later variant by position.

    Returns:
        A reason label from `REASONS`, or `None` if the pair is compatible (distance >= 0).
    """
    if variant_distance(v1, v2) >= 0:
        return None
    pos1, ref1, alt1 = v1[1], v1[2], v1[3]
    pos2, ref2, alt2 = v2[1], v2[2], v2[3]
    v1_del = len(ref1) > len(alt1)
    v2_del = len(ref2) > len(alt2)
    v2_snp = len(ref2) == 1 and len(alt2) == 1
    v1_ins = len(alt1) > len(ref1)
    if pos1 == pos2:
        return "same_position"
    if v1_del and v2_snp and pos1 < pos2 < pos1 + len(ref1):
        return "snp_in_deletion"
    if v1_del and v2_del:
        return "overlapping_deletions"
    if v1_del and pos1 < pos2 < pos1 + len(ref1):
        return "indel_in_deletion"
    if v1_ins:
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


def compatibility_flag(variants_str: str) -> str:
    """
    Compute the VCF-style `haplotype_filter` value for a haplotype's `variants` string.

    Args:
        variants_str: Comma-separated `chr:pos:ref:alt` component variants.

    Returns:
        `"PASS"` if no adjacent component variants overlap (including single-variant rows), else
        the `;`-joined sorted distinct incompatibility reasons.
    """
    reasons = classify_haplotype(parse_variants_string(variants_str))
    return ";".join(sorted(set(reasons))) if reasons else "PASS"


def variants_overlap(v1: Variant, v2: Variant) -> bool:
    """
    Return whether two variants' reference spans overlap (order-independent).

    Args:
        v1: A variant.
        v2: Another variant.

    Returns:
        True iff the later variant starts within the earlier variant's reference span.
    """
    earlier, later = (v1, v2) if v1[1] <= v2[1] else (v2, v1)
    return later[1] < earlier[1] + len(earlier[2])


def count_bypass_resolutions(variants: list[Variant]) -> int:
    """
    Count distinct maximal pairwise-compatible sub-haplotypes (>= 2 variants).

    Models the "explode the conflict" alternative to dropping: a resolution keeps a maximal set
    of variants with no two overlapping (a maximal independent set in the overlap graph), valid
    only if it retains >= 2 variants. For an incompatible haplotype, 0 means exploding recovers
    nothing (e.g. a length-2 conflict yields only singletons), and larger values quantify the
    combinatorial blow-up at multi-conflict repeat loci.

    Args:
        variants: Component variants of the haplotype.

    Returns:
        The number of distinct maximal compatible resolutions of size >= 2.
    """
    # Bron-Kerbosch maximal-clique enumeration is worst-case exponential, but DivRef haplotypes are
    # short (a handful of variants, never more than a few dozen even at dense repeat loci), so the
    # graph is tiny and enumeration is cheap. No size guard is needed for this bounded input.
    n = len(variants)
    if n < 2:
        return 0
    # Complement adjacency: i ~ j iff the two variants do NOT overlap (i.e. are compatible).
    compatible: list[set[int]] = [set() for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if not variants_overlap(variants[i], variants[j]):
                compatible[i].add(j)
                compatible[j].add(i)

    count = 0

    def bron_kerbosch(r: set[int], p: set[int], x: set[int]) -> None:
        # Maximal cliques in the complement graph == maximal independent sets in the overlap graph.
        nonlocal count
        if not p and not x:
            if len(r) >= 2:
                count += 1
            return
        pivot = max(p | x, key=lambda u: len(compatible[u] & p))
        for v in list(p - compatible[pivot]):
            bron_kerbosch(r | {v}, p & compatible[v], x & compatible[v])
            p = p - {v}
            x = x | {v}

    bron_kerbosch(set(), set(range(n)), set())
    return count


def end_coordinate_shortfall(variants: list[Variant], window_size: int, stored_end: int) -> int:
    """
    Compute bp by which a stored 0-based-exclusive `end` fails to cover all variants' references.

    Args:
        variants: Component variants of the haplotype.
        window_size: Flanking reference-context size.
        stored_end: The 0-based-exclusive `end` recorded for the haplotype.

    Returns:
        `max_v(pos - 1 + len(ref)) + window_size - stored_end`; > 0 means deleted reference is
        truncated out of the window.
    """
    rightmost_ref_end = max(v[1] - 1 + len(v[2]) for v in variants)  # 0-based exclusive
    return rightmost_ref_end + window_size - stored_end


def start_coordinate_shortfall(variants: list[Variant], window_size: int, stored_start: int) -> int:
    """
    Compute bp by which a stored 0-based-inclusive `start` is too far right (should be 0).

    Reference alleles only extend rightward, so the leftmost touched base is the min-position
    variant and `start` is structurally correct. This is a defensive check.

    Args:
        variants: Component variants of the haplotype.
        window_size: Flanking reference-context size.
        stored_start: The 0-based-inclusive `start` recorded for the haplotype.

    Returns:
        `stored_start - (min_v(pos - 1) - window_size)`; should always be 0.
    """
    leftmost_ref_start = min(v[1] - 1 for v in variants)  # 0-based inclusive
    required_start = leftmost_ref_start - window_size
    return stored_start - required_start
