"""Compare old vs new compute_haplotypes outputs by variant tuples."""

from collections import Counter

import hail as hl

hl.init(quiet=True)

old = hl.read_table(
    "/Users/ameynert/fulcrum/repos/divref-wf/data/analysis/compute_haplotypes/test_data_old/chr1_100001_200000_haplotypes.ht"
)
new = hl.read_table(
    "/Users/ameynert/fulcrum/repos/divref-wf/data/analysis/compute_haplotypes/test_data_new/chr1_100001_200000_haplotypes.ht"
)


def variant_tuple(row: hl.Struct) -> tuple:
    return tuple(
        (v.locus.contig, v.locus.position, tuple(v.alleles)) for v in row.variants
    )


old_rows = old.collect()
new_rows = new.collect()

old_map = {variant_tuple(r): r for r in old_rows}
new_map = {variant_tuple(r): r for r in new_rows}

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
# Old has max_empirical_AC; new also has max_empirical_AC (int64 vs int32).
print()
print("=== AC comparison for shared haplotypes ===")
ac_diffs = []
for k in shared:
    old_ac = old_map[k].max_empirical_AC
    new_ac = new_map[k].max_empirical_AC
    ac_diffs.append((old_ac, new_ac))
diff_counts = Counter((o, n) for o, n in ac_diffs if o != n)
same_count = sum(1 for o, n in ac_diffs if o == n)
new_higher = sum(1 for o, n in ac_diffs if n > o)
old_higher = sum(1 for o, n in ac_diffs if o > n)
print(f"Shared with same max_empirical_AC: {same_count}")
print(f"Shared where new > old: {new_higher}")
print(f"Shared where old > new: {old_higher}")
