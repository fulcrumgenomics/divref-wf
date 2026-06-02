"""Compare old vs new compute_haplotypes haplotype catalogs (DuckDB vs DuckDB).

Both the original (DivRef 1.1) and the new pipeline expose haplotypes as `sequences`
rows with `source = 'HGDP_haplotype'`, carrying a `variants` string
(`chr:pos:ref:alt,...`) and a `popmax_empirical_AC` column. This compares the two
catalogs over the autosomes (the contigs both algorithms compute), partitions them
into shared / old-only / new-only, analyses old-only sub-fragments and new-only
variant provenance, and compares AC for shared haplotypes.

Emits a stdout summary plus a machine-readable summary TSV
(`data/analysis/compute_haplotypes/algo_comparison.summary.tsv`); downstream plotting
(`compare_haplotypes_venn.R`) and the blog read counts from the TSV, not from
hardcoded values. Pure Python + DuckDB, no Spark.

    pixi run python scripts/compare_haplotypes.py                 # autosomes chr1-chr22
    pixi run python scripts/compare_haplotypes.py --contigs chr22 # chr22 deep-dive
"""

import argparse
import sys
from collections import Counter
from collections import defaultdict
from pathlib import Path

import duckdb

DEFAULT_OLD_DUCKDB = "data/analysis/input/DivRef-v1.1.haplotypes_gnomad_merge.index.duckdb"
DEFAULT_NEW_DUCKDB = "data/work/output/hgdp_1kg.haplotypes_gnomad_merge.index.duckdb"
DEFAULT_CONTIGS = [f"chr{i}" for i in range(1, 23)]
SUMMARY_TSV = Path("data/analysis/compute_haplotypes/algo_comparison.summary.tsv")

# A haplotype key: an ordered tuple of (contig, position, (ref, alt)) variants. Variant
# tuples embed the contig, so keys from different contigs never collide.
Variant = tuple[str, int, tuple[str, str]]
VariantTuple = tuple[Variant, ...]


def parse_variants_string(s: str) -> VariantTuple:
    """Parse a comma-separated `chr:pos:ref:alt` list into a variant tuple."""
    out: list[Variant] = []
    for token in s.split(","):
        contig, pos, ref, alt = token.split(":")
        out.append((contig, int(pos), (ref, alt)))
    return tuple(out)


def key_to_variants_string(key: VariantTuple) -> str:
    """Reconstruct the `chr:pos:ref:alt,...` string for a haplotype key (stored order)."""
    return ",".join(f"{contig}:{pos}:{ref}:{alt}" for contig, pos, (ref, alt) in key)


def has_overlapping_pair(key: VariantTuple) -> bool:
    """Return whether two component variants' reference spans overlap.

    Sorts by (position, reference length) and checks each adjacent pair's distance
    (`pos2 - pos1 - len(ref1)`); a negative value means the reference spans overlap. Both builders
    construct overlapping variants by composition (the new cursor builder differently from the
    original concatenation), so only non-overlapping haplotypes are guaranteed an identical
    sequence -- a co-located SNP+indel is compatible (PASS) yet built differently by each.
    """
    ordered = sorted(key, key=lambda v: (v[1], len(v[2][0])))
    return any(
        later[1] - earlier[1] - len(earlier[2][0]) < 0
        for earlier, later in zip(ordered, ordered[1:], strict=False)
    )


def is_contiguous_sub(short: VariantTuple, long_: VariantTuple) -> bool:
    """True iff `short` is a strict contiguous sub-array of `long_`."""
    n, m = len(short), len(long_)
    if n >= m:
        return False
    for i in range(m - n + 1):
        if long_[i : i + n] == short:
            return True
    return False


def build_map(
    duckdb_path: str, contigs: list[str], side: str
) -> dict[VariantTuple, tuple[int, str]]:
    """Read HGDP_haplotype rows for `contigs` into a {key: (popmax_empirical_AC, sequence)} map.

    On a repeated variant tuple, keeps the row with the larger `popmax_empirical_AC` and logs a
    line to stderr naming the tuple and both AC values, so a silent merge is never missed. Both
    index builders already deduplicate haplotypes, so any collision flags a regression.
    """
    placeholders = ",".join(["?"] * len(contigs))
    con = duckdb.connect(str(Path(duckdb_path).resolve()), read_only=True)
    query_rows = con.execute(
        "SELECT variants, popmax_empirical_AC, sequence FROM sequences "
        f"WHERE source = 'HGDP_haplotype' AND contig IN ({placeholders})",  # noqa: S608
        contigs,
    ).fetchall()
    con.close()

    result: dict[VariantTuple, tuple[int, str]] = {}
    collisions = 0
    for variants_str, ac, sequence in query_rows:
        key = parse_variants_string(variants_str)
        if key in result:
            collisions += 1
            print(
                f"DUPLICATE TUPLE [{side}]: {variants_str} popmax_AC {result[key][0]} vs {ac}",
                file=sys.stderr,
            )
            if ac > result[key][0]:
                result[key] = (ac, sequence)
        else:
            result[key] = (ac, sequence)
    print(
        f"{side}: {len(query_rows)} rows -> {len(result)} unique haplotypes "
        f"({collisions} duplicate tuples)"
    )
    return result


def compare_shared_sequences(
    shared: set[VariantTuple],
    old_map: dict[VariantTuple, tuple[int, str]],
    new_map: dict[VariantTuple, tuple[int, str]],
) -> tuple[int, int, int, int, list[str]]:
    """Compare old vs new stored sequences for shared haplotypes, split by whether they overlap.

    A shared haplotype whose component variants do not overlap must build to the identical sequence
    under both algorithms; a mismatch there is a regression. A haplotype with an overlapping pair
    (a composable SNP co-located with an indel, or a flagged incompatibility) is composed by the
    new cursor builder differently from the original concatenation, so it may legitimately differ.

    Args:
        shared: Haplotype keys present in both maps.
        old_map: Original-pipeline `{key: (popmax_empirical_AC, sequence)}`.
        new_map: New-pipeline `{key: (popmax_empirical_AC, sequence)}`.

    Returns:
        `(independent_match, independent_mismatch, overlapping_match, overlapping_differ,
        independent_mismatch_examples)`; a non-empty examples list (capped at 10) is a regression
        to investigate, since non-overlapping variants must build the same sequence.
    """
    independent_match = 0
    independent_mismatch = 0
    overlapping_match = 0
    overlapping_differ = 0
    independent_mismatch_examples: list[str] = []
    for k in shared:
        sequences_match = old_map[k][1] == new_map[k][1]
        if has_overlapping_pair(k):
            if sequences_match:
                overlapping_match += 1
            else:
                overlapping_differ += 1
        elif sequences_match:
            independent_match += 1
        else:
            independent_mismatch += 1
            if len(independent_mismatch_examples) < 10:
                independent_mismatch_examples.append(key_to_variants_string(k))
    return (
        independent_match,
        independent_mismatch,
        overlapping_match,
        overlapping_differ,
        independent_mismatch_examples,
    )


def count_contig_haplotypes(duckdb_path: str, contig: str) -> int:
    """Count HGDP_haplotype rows for a single contig (used for the chrX coverage note)."""
    con = duckdb.connect(str(Path(duckdb_path).resolve()), read_only=True)
    n = con.execute(
        "SELECT COUNT(*) FROM sequences WHERE source = 'HGDP_haplotype' AND contig = ?",
        [contig],
    ).fetchone()[0]
    con.close()
    return int(n)


def main() -> None:
    """Run the DuckDB-to-DuckDB haplotype comparison and write the summary TSV."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--old-duckdb", default=DEFAULT_OLD_DUCKDB, help="Original-pipeline (DivRef 1.1) DuckDB."
    )
    parser.add_argument(
        "--new-duckdb", default=DEFAULT_NEW_DUCKDB, help="New whole-genome DuckDB index."
    )
    parser.add_argument(
        "--contigs",
        nargs="+",
        default=DEFAULT_CONTIGS,
        help="Contigs to compare (default: autosomes chr1-chr22).",
    )
    args = parser.parse_args()

    old_map = build_map(args.old_duckdb, args.contigs, "old")
    new_map = build_map(args.new_duckdb, args.contigs, "new")

    old_keys = set(old_map)
    new_keys = set(new_map)
    shared = old_keys & new_keys
    old_only = old_keys - new_keys
    new_only = new_keys - old_keys

    print(f"Old count: {len(old_keys)}")
    print(f"New count: {len(new_keys)}")
    print(f"Shared:    {len(shared)}")
    print(f"Old only:  {len(old_only)}")
    print(f"New only:  {len(new_only)}")

    # --- Old-only sub-fragment analysis (inverted-index accelerated) --------------
    #
    # For each old-only haplotype, find any new haplotype that strictly contains it as a
    # contiguous sub-array. An enclosing new haplotype must contain every variant of the
    # old-only key, so it lies in the intersection of the per-variant new-key sets; the
    # rarest variant's set is the cheapest superset to scan. This bounds the work by how
    # many new haplotypes actually share a variant, rather than O(old_only x new_keys).
    #
    # - both_have_longer:    enclosing new ALSO emitted by old -> both algorithms find the
    #                        longer fragment; only new recognises the short one as redundant.
    # - only_new_has_longer: enclosing new NOT in old -> old emitted the short fragment but
    #                        never accumulated enough AC for the longer one, while new's
    #                        containment counting credited it with carriers from every parent.
    new_keys_by_variant: dict[Variant, set[VariantTuple]] = defaultdict(set)
    for nk in new_keys:
        for v in nk:
            new_keys_by_variant[v].add(nk)

    old_only_subfragment_of_new = 0
    old_only_subfragment_both_have_longer = 0
    old_only_subfragment_only_new_has_longer = 0
    for k in old_only:
        rarest = min(k, key=lambda v: len(new_keys_by_variant.get(v, ())))
        candidates = new_keys_by_variant.get(rarest, set())
        enclosing_news = [
            nk for nk in candidates if len(nk) > len(k) and is_contiguous_sub(k, nk)
        ]
        if not enclosing_news:
            continue
        old_only_subfragment_of_new += 1
        if any(nk in old_keys for nk in enclosing_news):
            old_only_subfragment_both_have_longer += 1
        else:
            old_only_subfragment_only_new_has_longer += 1
    old_only_not_subfragment = len(old_only) - old_only_subfragment_of_new

    print()
    print("=== Old-only sub-fragment analysis ===")
    print(
        f"Old-only that are proper contiguous sub-fragments of some new haplotype: "
        f"{old_only_subfragment_of_new} / {len(old_only)}"
    )
    print(
        f"  enclosing new also in old (both algorithms find longer): "
        f"{old_only_subfragment_both_have_longer}"
    )
    print(
        f"  enclosing new NOT in old (new's containment counting recovered AC):  "
        f"{old_only_subfragment_only_new_has_longer}"
    )
    print(f"  old-only NOT a sub-fragment of any new (residual): {old_only_not_subfragment}")

    # --- New-only decomposition ---------------------------------------------------
    old_all_variants: set[Variant] = set()
    for k in old_keys:
        old_all_variants.update(k)

    new_only_all_variants_in_old = 0
    new_only_mixed_variants = 0
    new_only_all_novel = 0
    all_novel_lengths: list[int] = []
    for k in new_only:
        in_old = [v in old_all_variants for v in k]
        if all(in_old):
            new_only_all_variants_in_old += 1
        elif not any(in_old):
            new_only_all_novel += 1
            all_novel_lengths.append(len(k))
        else:
            new_only_mixed_variants += 1

    all_novel_length_counts = Counter(all_novel_lengths)
    new_only_all_novel_length_2 = all_novel_length_counts.get(2, 0)
    new_only_all_novel_max_length = max(all_novel_lengths) if all_novel_lengths else 0

    print()
    print("=== New-only decomposition ===")
    print(f"All variants present somewhere in old:        {new_only_all_variants_in_old}")
    print(f"Mix of shared and novel variants:             {new_only_mixed_variants}")
    print(f"All variants novel (none in any old haplo):   {new_only_all_novel}")
    print(f"  length distribution: {dict(sorted(all_novel_length_counts.items()))}")
    print(f"  length-2 count:      {new_only_all_novel_length_2}")
    print(f"  max length:          {new_only_all_novel_max_length}")

    # --- AC comparison for shared haplotypes --------------------------------------
    shared_same_ac = 0
    shared_new_higher_ac = 0
    shared_old_higher_ac = 0
    for k in shared:
        old_ac = old_map[k][0]
        new_ac = new_map[k][0]
        if old_ac == new_ac:
            shared_same_ac += 1
        elif new_ac > old_ac:
            shared_new_higher_ac += 1
        else:
            shared_old_higher_ac += 1

    print()
    print("=== AC comparison for shared haplotypes ===")
    print(f"Shared with same AC:        {shared_same_ac}")
    print(f"Shared where new > old AC:  {shared_new_higher_ac}")
    print(f"Shared where old > new AC:  {shared_old_higher_ac}")

    # --- Sequence equality for shared haplotypes ----------------------------------
    # A non-overlapping shared haplotype that differs is a regression, printed loudly. Overlapping
    # ones (composable co-located variants or flagged incompatibilities) are composed differently
    # by the new cursor builder than by the original concatenation, so differences are expected.
    (
        independent_seq_match,
        independent_seq_mismatch,
        overlapping_seq_match,
        overlapping_seq_differ,
        independent_mismatch_examples,
    ) = compare_shared_sequences(shared, old_map, new_map)

    print()
    print("=== Sequence equality for shared haplotypes ===")
    print(f"Non-overlapping shared, sequence matches:  {independent_seq_match}")
    print(f"Non-overlapping shared, sequence MISMATCH: {independent_seq_mismatch}")
    print(f"Overlapping shared, sequence matches:      {overlapping_seq_match}")
    print(f"Overlapping shared, sequence differs:      {overlapping_seq_differ}")
    if independent_mismatch_examples:
        print("  !! non-overlapping sequence mismatches (unexpected -- investigate):")
        for variants_str in independent_mismatch_examples:
            print(f"     {variants_str}")

    # --- chrX coverage note (separate from the autosome algorithm comparison) -----
    # Every chrX haplotype is new-only by construction: the original workflow computed no
    # chrX haplotypes, so this is a coverage difference, not an algorithm difference.
    new_chrx_haplotypes = count_contig_haplotypes(args.new_duckdb, "chrX")
    print()
    print("=== chrX coverage (separate from the algorithm comparison) ===")
    print(f"New chrX haplotypes (old workflow computed none): {new_chrx_haplotypes}")

    # --- Summary TSV --------------------------------------------------------------
    SUMMARY_TSV.parent.mkdir(parents=True, exist_ok=True)
    summary_rows: list[tuple[str, int]] = [
        ("shared", len(shared)),
        ("old_only", len(old_only)),
        ("new_only", len(new_only)),
        ("old_only_subfragment_of_new", old_only_subfragment_of_new),
        ("old_only_subfragment_both_have_longer", old_only_subfragment_both_have_longer),
        ("old_only_subfragment_only_new_has_longer", old_only_subfragment_only_new_has_longer),
        ("old_only_not_subfragment", old_only_not_subfragment),
        ("new_only_all_variants_in_old", new_only_all_variants_in_old),
        ("new_only_mixed_variants", new_only_mixed_variants),
        ("new_only_all_novel", new_only_all_novel),
        ("new_only_all_novel_length_2", new_only_all_novel_length_2),
        ("new_only_all_novel_max_length", new_only_all_novel_max_length),
        ("shared_same_ac", shared_same_ac),
        ("shared_new_higher_ac", shared_new_higher_ac),
        ("shared_old_higher_ac", shared_old_higher_ac),
        ("shared_independent_seq_match", independent_seq_match),
        ("shared_independent_seq_mismatch", independent_seq_mismatch),
        ("shared_overlapping_seq_match", overlapping_seq_match),
        ("shared_overlapping_seq_differ", overlapping_seq_differ),
        ("new_chrX_haplotypes", new_chrx_haplotypes),
    ]
    with SUMMARY_TSV.open("w") as f:
        f.write("metric\tvalue\n")
        for name, value in summary_rows:
            f.write(f"{name}\t{value}\n")
    print()
    print(f"wrote summary TSV: {SUMMARY_TSV}")


if __name__ == "__main__":
    main()
