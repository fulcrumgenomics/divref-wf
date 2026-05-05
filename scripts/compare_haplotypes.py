"""Compare old vs new compute_haplotypes outputs by variant tuples.

Old haplotypes come from the DuckDB index produced by the old pipeline
(`sequences` table, `source = 'HGDP_haplotype'`). New haplotypes come directly
from the new `compute_haplotypes` Hail table.
"""

from collections import Counter
from pathlib import Path

import duckdb
import hail as hl

OLD_DUCKDB = (
    "data/analysis/compute_haplotypes/test_data_old/hgdp_1kg.haplotypes_gnomad_merge.index.duckdb"
)
NEW_HT = (
    "data/analysis/compute_haplotypes/test_data_new/hgdp_1kg.haplotypes.chr22.ht"
)
CONTIG = "chr22"


def parse_variants_string(s: str) -> tuple:
    """Parse a comma-separated `chr:pos:ref:alt` list into a variant tuple."""
    out = []
    for token in s.split(","):
        contig, pos, ref, alt = token.split(":")
        out.append((contig, int(pos), (ref, alt)))
    return tuple(out)


def variant_tuple_from_ht_row(row: hl.Struct) -> tuple:
    return tuple(
        (v.locus.contig, v.locus.position, tuple(v.alleles)) for v in row.variants
    )


# Load old (DuckDB) ----
con = duckdb.connect(str(Path(OLD_DUCKDB).resolve()), read_only=True)
old_query = con.execute(
    "SELECT variants, popmax_empirical_AC FROM sequences "
    "WHERE source = 'HGDP_haplotype' AND contig = ?",
    [CONTIG],
).fetchall()
con.close()

old_map: dict[tuple, int] = {}
for variants_str, ac in old_query:
    key = parse_variants_string(variants_str)
    # If the same haplotype appears more than once (shouldn't, but defensively
    # handle the old two-pass union), keep the max AC.
    if key in old_map:
        old_map[key] = max(old_map[key], ac)
    else:
        old_map[key] = ac

# Load new (Hail table) ----
hl.init(quiet=True)
new = hl.read_table(NEW_HT)
new_rows = new.collect()
new_map = {variant_tuple_from_ht_row(r): r.max_empirical_AC for r in new_rows}

old_keys = set(old_map.keys())
new_keys = set(new_map.keys())

shared = old_keys & new_keys
old_only = old_keys - new_keys
new_only = new_keys - old_keys

print(f"Old count: {len(old_keys)}")
print(f"New count: {len(new_keys)}")
print(f"Shared: {len(shared)}")
print(f"Old only: {len(old_only)}")
print(f"New only: {len(new_only)}")


def length_dist(keys: set) -> dict:
    return dict(sorted(Counter(len(k) for k in keys).items()))


print()
print("=== Length distribution ===")
print(f"Old all   : {length_dist(old_keys)}")
print(f"New all   : {length_dist(new_keys)}")
print(f"Shared    : {length_dist(shared)}")
print(f"Old only  : {length_dist(old_only)}")
print(f"New only  : {length_dist(new_only)}")


def is_contiguous_sub(short: tuple, long_: tuple) -> bool:
    n, m = len(short), len(long_)
    if n >= m:
        return False
    for i in range(m - n + 1):
        if long_[i : i + n] == short:
            return True
    return False


# Are old-only haplotypes proper sub-fragments of any NEW haplotype? -> would mean
# the new algorithm's containment dedup dropped them as redundant.
new_keys_by_len = sorted(new_keys, key=len, reverse=True)
sub_of_new_count = 0
for k in old_only:
    for nk in new_keys_by_len:
        if len(nk) <= len(k):
            break
        if is_contiguous_sub(k, nk):
            sub_of_new_count += 1
            break
print()
print(
    f"Old-only that are proper sub-fragments of some new haplotype: "
    f"{sub_of_new_count} / {len(old_only)}"
)

# Are old-only haplotypes proper sub-fragments of any OTHER OLD haplotype?
# -> would mean the new algorithm could have dropped them via containment dedup
# even if they didn't appear in the new output.
old_keys_by_len = sorted(old_keys, key=len, reverse=True)
sub_of_old_count = 0
for k in old_only:
    for ok in old_keys_by_len:
        if len(ok) <= len(k):
            break
        if is_contiguous_sub(k, ok):
            sub_of_old_count += 1
            break
print(
    f"Old-only that are proper sub-fragments of some other old haplotype: "
    f"{sub_of_old_count} / {len(old_only)}"
)

# Are old-only haplotypes' variant SETS contained in any new haplotype?
new_var_sets = [set(k) for k in new_keys]
old_subset_of_some_new = 0
old_disjoint_from_all_new = 0
for k in old_only:
    s = set(k)
    if any(s <= ns for ns in new_var_sets):
        old_subset_of_some_new += 1
    elif all(s.isdisjoint(ns) for ns in new_var_sets):
        old_disjoint_from_all_new += 1
print()
print(
    f"Old-only with variants ⊆ some new haplotype's variants: "
    f"{old_subset_of_some_new} / {len(old_only)}"
)
print(
    f"Old-only completely variant-disjoint from every new haplotype: "
    f"{old_disjoint_from_all_new} / {len(old_only)}"
)


# How many old-only haplotypes have any single variant that is NEVER in any new haplotype?
new_all_variants = set()
for k in new_keys:
    new_all_variants.update(k)
old_with_variant_not_in_new = 0
for k in old_only:
    if not all(v in new_all_variants for v in k):
        old_with_variant_not_in_new += 1
print(
    f"Old-only containing at least one variant absent from every new haplotype: "
    f"{old_with_variant_not_in_new} / {len(old_only)}"
)


# Mutually exclusive partition of old-only haplotypes, ordered from most-missing
# (variants entirely absent from new) to least-missing (would-be present in new
# except that containment dedup dropped them).
#
# Tiers:
#   1. fully variant-disjoint: no variant in this haplotype appears in any new haplotype.
#   2. partially missing variants: at least one variant absent from every new haplotype,
#      but not all (some variants do appear in some new haplotype somewhere).
#   3. variants split across new haplotypes: every variant appears somewhere in new,
#      but no single new haplotype contains all of them as a subset.
#   4. variant-set subsumed by some new haplotype, but not contiguous: all variants
#      fit inside one new haplotype's variant set, yet this old key is not a
#      contiguous sub-array of any new haplotype (different ordering or non-
#      adjacency-contiguous selection).
#   5. proper contiguous sub-fragment of some new haplotype: would have been emitted
#      by new but was dropped by containment dedup (same per-pop AC as the enclosing
#      fragment).
new_var_sets_by_len = sorted(((nk, set(nk)) for nk in new_keys), key=lambda t: len(t[0]), reverse=True)
new_all_variants: set = set()
for nk in new_keys:
    new_all_variants.update(nk)


def categorize_old_only(k: tuple) -> str:
    if all(v not in new_all_variants for v in k):
        return "1_fully_variant_disjoint"
    if any(v not in new_all_variants for v in k):
        return "2_partially_missing_variants"
    s = set(k)
    if not any(s <= ns for _, ns in new_var_sets_by_len):
        return "3_variants_split_across_new_haplotypes"
    for nk, _ in new_var_sets_by_len:
        if len(nk) <= len(k):
            break
        if is_contiguous_sub(k, nk):
            return "5_proper_contiguous_subfragment_of_new"
    return "4_subset_of_some_new_but_not_contiguous"


tier_counts = Counter(categorize_old_only(k) for k in old_only)
print()
print("=== Old-only mutually exclusive breakdown (most missing → least missing) ===")
total = 0
for tier in (
    "1_fully_variant_disjoint",
    "2_partially_missing_variants",
    "3_variants_split_across_new_haplotypes",
    "4_subset_of_some_new_but_not_contiguous",
    "5_proper_contiguous_subfragment_of_new",
):
    c = tier_counts.get(tier, 0)
    total += c
    print(f"  {tier}: {c}")
print(f"  total: {total} (old-only count: {len(old_only)})")


# Inspect a few old-only haplotypes — print the variant positions
print()
print("=== First 5 old-only haplotypes (positions) ===")
for k in list(old_only)[:5]:
    positions = [v[1] for v in k]
    gaps = [positions[i] - positions[i - 1] for i in range(1, len(positions))]
    print(f"  n={len(k)}, positions={positions}, gaps={gaps}")

print()
print("=== First 5 new-only haplotypes (positions) ===")
for k in list(new_only)[:5]:
    positions = [v[1] for v in k]
    gaps = [positions[i] - positions[i - 1] for i in range(1, len(positions))]
    print(f"  n={len(k)}, positions={positions}, gaps={gaps}")

# Are the AC values different for shared haplotypes?
# Old: popmax_empirical_AC from the DuckDB sequences table.
# New: max_empirical_AC from the Hail table.
print()
print("=== AC comparison for shared haplotypes ===")
same_count = 0
new_higher = 0
old_higher = 0
for k in shared:
    old_ac = old_map[k]
    new_ac = new_map[k]
    if old_ac == new_ac:
        same_count += 1
    elif new_ac > old_ac:
        new_higher += 1
    else:
        old_higher += 1
print(f"Shared with same AC: {same_count}")
print(f"Shared where new > old: {new_higher}")
print(f"Shared where old > new: {old_higher}")
