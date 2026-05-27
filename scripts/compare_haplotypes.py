"""Compare old vs new compute_haplotypes outputs by variant tuples.

Old haplotypes come from the DuckDB index produced by the old pipeline
(`sequences` table, `source = 'HGDP_haplotype'`). New haplotypes come directly
from the new `compute_haplotypes` Hail table.
"""

import random
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
SITES_HT = "data/work/inputs/hgdp_1kg.sites.chr22.ht"
CONTIG = "chr22"
HAPLOTYPE_FREQ_THRESHOLD = 0.005
ADJACENCY_WINDOW = 25


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
    "SELECT variants, popmax_empirical_AC, estimated_gnomad_AF "
    "FROM sequences WHERE source = 'HGDP_haplotype' AND contig = ?",
    [CONTIG],
).fetchall()
con.close()

old_map: dict[tuple, int] = {}
old_est_af: dict[tuple, float] = {}
for variants_str, ac, est_af in old_query:
    key = parse_variants_string(variants_str)
    # If the same haplotype appears more than once (shouldn't, but defensively
    # handle the old two-pass union), keep the max AC.
    if key in old_map:
        old_map[key] = max(old_map[key], ac)
        old_est_af[key] = max(old_est_af[key], est_af)
    else:
        old_map[key] = ac
        old_est_af[key] = est_af

# Load new (Hail table) ----
hl.init(quiet=True)
new = hl.read_table(NEW_HT)
new_rows = new.collect()
new_map = {variant_tuple_from_ht_row(r): r.max_empirical_AC for r in new_rows}
new_est_af: dict[tuple, float] = {
    variant_tuple_from_ht_row(r): r.estimated_gnomad_AF for r in new_rows
}

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

# Investigate tiers 1 & 2 (variants absent from any new haplotype): is exclusion
# from new driven by (a) the new algorithm's tighter 25 bp adjacency rule cutting
# pairs that old's 100 bp bins kept (structural), or (b) old's two-pass-union AC
# inflation pushing `estimated_gnomad_AF` over the 0.005 threshold for haplotypes
# that, correctly counted, would fail it (filter)?
#
# Discriminator: the ref-aware gap used by new's adjacency rule
#   gap_i = pos[i+1] - pos[i] - len(ref[i])
# For a length-2 haplotype with gap < window, new MUST consider it (parent block
# of length 2 forms in any sample carrying both alts). Its absence from new
# haplotypes can only be explained by the threshold/AN filter — i.e., a count
# difference. For gap >= window, new's adjacency rule cuts the pair, so the
# absence is structural — bin-vs-adjacency, not AC.
#
# We then look up component variants in the upstream sites Hail table to fold in
# the smallest gnomAD per-pop AF, and report what fraction of "filter-rejected"
# tier-1+2 haplotypes have an old `estimated_gnomad_AF` that would fall below
# threshold if AC were halved (a proxy for two-pass-union double counting being
# halved away in the new single-pass count).


def max_ref_aware_gap(k: tuple) -> int:
    if len(k) < 2:
        return 0
    gaps: list[int] = []
    for i in range(1, len(k)):
        prev = k[i - 1]
        curr = k[i]
        ref_prev_len = len(prev[2][0])
        gaps.append(curr[1] - prev[1] - ref_prev_len)
    return max(gaps)


tier12_keys = [
    k
    for k in old_only
    if categorize_old_only(k)
    in {"1_fully_variant_disjoint", "2_partially_missing_variants"}
]
tier12_structural = [k for k in tier12_keys if max_ref_aware_gap(k) >= ADJACENCY_WINDOW]
tier12_filterable = [k for k in tier12_keys if max_ref_aware_gap(k) < ADJACENCY_WINDOW]

print()
print("=== Tier 1+2 root-cause split ===")
print(f"Total tier 1+2: {len(tier12_keys)}")
print(
    f"  structural (max ref-aware gap >= {ADJACENCY_WINDOW} bp; new's adjacency cuts): "
    f"{len(tier12_structural)}"
)
print(
    f"  filter-only (max ref-aware gap < {ADJACENCY_WINDOW} bp; new sees them, must have"
    f" failed threshold): {len(tier12_filterable)}"
)


def quantiles(values: list[float]) -> str:
    if not values:
        return "n=0"
    s = sorted(values)
    n = len(s)
    return (
        f"n={n} min={s[0]:.4g} p10={s[n // 10]:.4g} p50={s[n // 2]:.4g} "
        f"p90={s[min(n - 1, n * 9 // 10)]:.4g} max={s[-1]:.4g}"
    )


print()
print(
    f"=== Tier 1+2 filter-only subgroup: would AC halving drop them below "
    f"{HAPLOTYPE_FREQ_THRESHOLD}? ==="
)
filterable_afs = [
    old_est_af[k] for k in tier12_filterable if k in old_est_af and old_est_af[k] is not None
]
print(f"old estimated_gnomad_AF distribution: {quantiles(filterable_afs)}")
halved = sum(1 for af in filterable_afs if af / 2 < HAPLOTYPE_FREQ_THRESHOLD)
quartered = sum(1 for af in filterable_afs if af / 4 < HAPLOTYPE_FREQ_THRESHOLD)
already_below = sum(1 for af in filterable_afs if af < HAPLOTYPE_FREQ_THRESHOLD)
print(
    f"  already < {HAPLOTYPE_FREQ_THRESHOLD} (shouldn't be — sanity check): {already_below}"
)
print(
    f"  AC/2 would drop est_af below {HAPLOTYPE_FREQ_THRESHOLD}: "
    f"{halved} / {len(filterable_afs)}"
)
print(
    f"  AC/4 would drop est_af below {HAPLOTYPE_FREQ_THRESHOLD}: "
    f"{quartered} / {len(filterable_afs)}"
)

# For comparison: same distribution for shared haplotypes that ARE in new
shared_afs = [old_est_af[k] for k in shared if k in old_est_af and old_est_af[k] is not None]
print(f"\nshared (passed both): {quantiles(shared_afs)}")

# Look up component variants in the upstream sites Hail table to compute, per
# tier-1+2 filter-only haplotype, the smallest gnomAD AF across its components in
# any single population — a lower bound on what new's `estimated_gnomad_AF` could
# have been before applying fraction_phased. If even this lower bound is well
# below threshold, the haplotype was already on a knife-edge.
print()
print("=== Looking up component variant pop_freqs from sites HT ===")
sites = hl.read_table(SITES_HT)
pops_legend: list[str] = sites.globals.pops.collect()[0]
n_pops = len(pops_legend)

# Build a driver-side dict: (pos, ref, alt) -> per-pop AFs
sites_collected = sites.aggregate(
    hl.agg.collect(
        hl.struct(
            pos=sites.locus.position,
            ref=sites.alleles[0],
            alt=sites.alleles[1],
            af=sites.pop_freqs.map(lambda x: x.AF),
        )
    )
)
sites_map: dict[tuple, list[float]] = {(s.pos, s.ref, s.alt): list(s.af) for s in sites_collected}

# For each filter-only tier-1+2 haplotype, find the per-population min component AF and
# track the largest such min (analogous to choosing max_pop). If this max-of-pop-mins
# is also small, then old's estimated_gnomad_AF = fraction_phased * (max-of-pop-mins),
# so a high old `estimated_gnomad_AF` implies a high `fraction_phased` — the
# AC-driven multiplier — and halving fraction_phased halves estimated_gnomad_AF.
not_in_sites = 0
component_min_pop_afs: list[float] = []
implied_fraction_phased: list[float] = []
for k in tier12_filterable:
    component_afs: list[list[float]] = []
    missing = False
    for v in k:
        key = (v[1], v[2][0], v[2][1])
        if key not in sites_map:
            missing = True
            break
        component_afs.append(sites_map[key])
    if missing:
        not_in_sites += 1
        continue
    # min across components per pop, then max across pops
    per_pop_mins = [min(afs[p] for afs in component_afs) for p in range(n_pops)]
    max_pop_min = max(per_pop_mins)
    component_min_pop_afs.append(max_pop_min)
    if max_pop_min > 0 and k in old_est_af:
        implied_fraction_phased.append(old_est_af[k] / max_pop_min)

print(f"filter-only haplotypes with all components in sites HT: "
      f"{len(tier12_filterable) - not_in_sites} / {len(tier12_filterable)}")
print(
    f"max-of-pop-mins of component gnomAD AF: {quantiles(component_min_pop_afs)}"
)
print(
    f"implied fraction_phased = old_est_af / (max-of-pop-mins): "
    f"{quantiles(implied_fraction_phased)}"
)
print(
    f"  if fraction_phased halved, est_af / 2 < {HAPLOTYPE_FREQ_THRESHOLD}: "
    f"{sum(1 for fp, af in zip(implied_fraction_phased, component_min_pop_afs, strict=False) if (fp / 2) * af < HAPLOTYPE_FREQ_THRESHOLD)}"
    f" / {len(implied_fraction_phased)}"
)


# ---------------------------------------------------------------------------
# Per-tier breakdown of the halving stats (blog TODO 3).
#
# The combined "Tier 1+2 filter-only" stat above doesn't distinguish whether
# the threshold-edge behavior is uniform across tiers, or whether one tier is
# driving the signal. Split and report separately.
# ---------------------------------------------------------------------------
tier1_keys = {k for k in old_only if categorize_old_only(k) == "1_fully_variant_disjoint"}
tier2_keys = {k for k in old_only if categorize_old_only(k) == "2_partially_missing_variants"}
tier1_filterable = [k for k in tier12_filterable if k in tier1_keys]
tier2_filterable = [k for k in tier12_filterable if k in tier2_keys]


def halving_breakdown(label: str, keys: list[tuple]) -> None:
    afs = [old_est_af[k] for k in keys if k in old_est_af and old_est_af[k] is not None]
    halved_n = sum(1 for af in afs if af / 2 < HAPLOTYPE_FREQ_THRESHOLD)
    quartered_n = sum(1 for af in afs if af / 4 < HAPLOTYPE_FREQ_THRESHOLD)
    print(f"\n{label} (n={len(afs)}):")
    print(f"  est_af distribution: {quantiles(afs)}")
    print(f"  AC/2 drops below {HAPLOTYPE_FREQ_THRESHOLD}: {halved_n} / {len(afs)}")
    print(f"  AC/4 drops below {HAPLOTYPE_FREQ_THRESHOLD}: {quartered_n} / {len(afs)}")


print()
print("=== Filter-only halving stats split by tier (TODO 3) ===")
halving_breakdown("Tier 1 (fully variant-disjoint) filter-only", tier1_filterable)
halving_breakdown("Tier 2 (partially missing variants) filter-only", tier2_filterable)


# ---------------------------------------------------------------------------
# Tier 2 "exclude variants that passed" analysis (blog TODO 4).
#
# Tier 2 haplotypes have a mix of variants that appear somewhere in the new
# output and variants that don't. Reduce each tier-2 haplotype to just its
# "absent" variants (those NOT in any new haplotype), recompute the lower
# bound on est_af using only the remaining variants' per-pop AFs, and ask
# whether halving AC would still drop the reduced haplotype below threshold.
#
# If yes for most: the haplotype was rescued by the inflated AC; the absent
# variants alone wouldn't have cleared 0.005 even at full empirical AC.
# ---------------------------------------------------------------------------
print()
print("=== Tier 2 'exclude variants that passed' analysis (TODO 4) ===")
tier2_reduced_results: list[tuple[int, float, float, float]] = []
tier2_zero_variants_remaining = 0
tier2_missing_sites = 0
for k in tier2_filterable:
    if k not in old_est_af or old_est_af[k] is None:
        continue
    absent_variants = tuple(v for v in k if v not in new_all_variants)
    if not absent_variants:
        # Definitionally impossible for tier 2, but guard.
        tier2_zero_variants_remaining += 1
        continue
    # Look up component AFs for the absent variants only.
    component_afs: list[list[float]] = []
    missing = False
    for v in absent_variants:
        key = (v[1], v[2][0], v[2][1])
        if key not in sites_map:
            missing = True
            break
        component_afs.append(sites_map[key])
    if missing:
        tier2_missing_sites += 1
        continue
    per_pop_mins_absent = [min(afs[p] for afs in component_afs) for p in range(n_pops)]
    max_pop_min_absent = max(per_pop_mins_absent)
    # Use the same implied fraction_phased that explained the full-haplotype
    # est_af, and apply it to the reduced-variant-set lower bound. This is an
    # upper-bound estimate on what the reduced haplotype's est_af could be —
    # the true value for the reduced haplotype, computed end-to-end, can only
    # be lower because the empirical_AC over a subset of variants is ≤ the
    # original.
    fp = old_est_af[k] / max_pop_min_absent if max_pop_min_absent > 0 else 0.0
    reduced_full_est_af = fp * max_pop_min_absent
    reduced_halved_est_af = reduced_full_est_af / 2
    tier2_reduced_results.append((
        len(absent_variants),
        max_pop_min_absent,
        reduced_full_est_af,
        reduced_halved_est_af,
    ))

print(f"tier 2 filter-only with components in sites HT: {len(tier2_reduced_results)}")
print(f"  zero absent variants (sanity): {tier2_zero_variants_remaining}")
print(f"  missing component AFs: {tier2_missing_sites}")
if tier2_reduced_results:
    reduced_full = [r[2] for r in tier2_reduced_results]
    reduced_halved = [r[3] for r in tier2_reduced_results]
    fail_full = sum(1 for af in reduced_full if af < HAPLOTYPE_FREQ_THRESHOLD)
    fail_halved = sum(1 for af in reduced_halved if af < HAPLOTYPE_FREQ_THRESHOLD)
    print(f"  reduced est_af (upper bound) distribution: {quantiles(reduced_full)}")
    print(
        f"  reduced est_af < {HAPLOTYPE_FREQ_THRESHOLD} at full AC: "
        f"{fail_full} / {len(reduced_full)}"
    )
    print(
        f"  reduced est_af < {HAPLOTYPE_FREQ_THRESHOLD} at half AC: "
        f"{fail_halved} / {len(reduced_halved)}"
    )


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


# ---------------------------------------------------------------------------
# Positional distribution of absent variants within Tier 2 haplotypes.
#
# For each tier-2 haplotype, classify the absent-variant index pattern:
#  - "prefix": absent variants form a contiguous run starting at index 0.
#  - "suffix": absent variants form a contiguous run ending at index n-1.
#  - "interior": all absent variants strictly between present variants (not at
#    either end).
#  - "split": absent variants in two or more disjoint runs.
#
# Also report the count of absent variants and the normalized index of the
# leftmost / rightmost absent variant.
# ---------------------------------------------------------------------------
print()
print("=== Tier 2 positional distribution of absent variants ===")


def classify_positions(absent_indices: list[int], n: int) -> str:
    if not absent_indices:
        return "empty"
    sorted_idx = sorted(absent_indices)
    # Detect runs of consecutive indices.
    runs: list[list[int]] = [[sorted_idx[0]]]
    for i in sorted_idx[1:]:
        if i == runs[-1][-1] + 1:
            runs[-1].append(i)
        else:
            runs.append([i])
    if len(runs) > 1:
        return "split"
    run = runs[0]
    starts_at_0 = run[0] == 0
    ends_at_n_minus_1 = run[-1] == n - 1
    if starts_at_0 and ends_at_n_minus_1:
        # All variants absent — wouldn't be tier 2, but guard.
        return "full"
    if starts_at_0:
        return "prefix"
    if ends_at_n_minus_1:
        return "suffix"
    return "interior"


tier2_positional = Counter()
tier2_n_absent = []
tier2_first_absent_frac = []
tier2_last_absent_frac = []
for k in tier2_filterable:
    n = len(k)
    absent_idx = [i for i, v in enumerate(k) if v not in new_all_variants]
    tier2_positional[classify_positions(absent_idx, n)] += 1
    tier2_n_absent.append(len(absent_idx))
    if absent_idx:
        # Use (n-1) as denominator so the value is in [0, 1] for length >= 2.
        denom = max(n - 1, 1)
        tier2_first_absent_frac.append(absent_idx[0] / denom)
        tier2_last_absent_frac.append(absent_idx[-1] / denom)

print(f"  n tier 2 filter-only: {len(tier2_filterable)}")
print("  pattern:")
for label in ("prefix", "suffix", "interior", "split", "full", "empty"):
    if tier2_positional[label]:
        pct = 100 * tier2_positional[label] / len(tier2_filterable)
        print(f"    {label}: {tier2_positional[label]} ({pct:.1f}%)")
print(f"  haplotype length distribution: {dict(sorted(Counter(len(k) for k in tier2_filterable).items()))}")
print(f"  absent count per haplotype: {dict(sorted(Counter(tier2_n_absent).items()))}")
print(f"  normalized first-absent index: {quantiles(tier2_first_absent_frac)}")
print(f"  normalized last-absent index: {quantiles(tier2_last_absent_frac)}")


# ---------------------------------------------------------------------------
# Haplotype length distribution by tier (within original-only) and the new-only set.
# ---------------------------------------------------------------------------
print()
print("=== Length distribution by tier ===")
by_tier: dict[str, list[tuple]] = {
    "1_fully_variant_disjoint": [],
    "2_partially_missing_variants": [],
    "3_variants_split_across_new_haplotypes": [],
    "4_subset_of_some_new_but_not_contiguous": [],
    "5_proper_contiguous_subfragment_of_new": [],
}
for k in old_only:
    by_tier[categorize_old_only(k)].append(k)


def fmt_len_dist(keys: list[tuple]) -> str:
    if not keys:
        return "n=0"
    lens = Counter(len(k) for k in keys)
    total = sum(lens.values())
    parts = [f"{n}:{c}" for n, c in sorted(lens.items())]
    return f"n={total} " + " ".join(parts)


print(f"  Tier 1 (fully variant-disjoint):           {fmt_len_dist(by_tier['1_fully_variant_disjoint'])}")
print(f"  Tier 2 (partially missing variants):       {fmt_len_dist(by_tier['2_partially_missing_variants'])}")
print(f"  Tier 3 (split across new haplotypes):      {fmt_len_dist(by_tier['3_variants_split_across_new_haplotypes'])}")
print(f"  Tier 4 (subset, not contiguous):           {fmt_len_dist(by_tier['4_subset_of_some_new_but_not_contiguous'])}")
print(f"  Tier 5 (contiguous sub-fragment of new):   {fmt_len_dist(by_tier['5_proper_contiguous_subfragment_of_new'])}")
print(f"  New-only (for reference):                  {fmt_len_dist(list(new_only))}")
print(f"  Shared (for reference):                    {fmt_len_dist(list(shared))}")


# ---------------------------------------------------------------------------
# Emit per-haplotype estimated_gnomad_AF rows to a TSV for plotting (one row
# per haplotype). Categories follow the same partition as the length-dist
# plot, plus Tiers 3 and 4 included so the reader can see where the small
# tiers fall. For original-only tiers, est_af is from the original algorithm
# (the value that made the haplotype pass the original's filter). For Shared
# and New-only, est_af is from the new algorithm.
# ---------------------------------------------------------------------------
af_tsv_path = Path("data/analysis/compute_haplotypes/algo_comparison.est_af.tsv")
af_tsv_path.parent.mkdir(parents=True, exist_ok=True)


def tier_label(tier_key: str) -> str:
    return {
        "1_fully_variant_disjoint": "Tier 1 (fully variant-disjoint)",
        "2_partially_missing_variants": "Tier 2 (partially missing variants)",
        "3_variants_split_across_new_haplotypes": "Tier 3 (split across new)",
        "4_subset_of_some_new_but_not_contiguous": "Tier 4 (subset, not contiguous)",
        "5_proper_contiguous_subfragment_of_new": "Tier 5 (sub-fragment of new)",
    }[tier_key]


with af_tsv_path.open("w") as f:
    f.write("category\test_af\n")
    for k in shared:
        af = new_est_af.get(k)
        if af is not None:
            f.write(f"Shared (both algorithms)\t{af}\n")
    for k in new_only:
        af = new_est_af.get(k)
        if af is not None:
            f.write(f"New-only (found only by new)\t{af}\n")
    for k in old_only:
        af = old_est_af.get(k)
        if af is None:
            continue
        f.write(f"{tier_label(categorize_old_only(k))}\t{af}\n")
print(f"\nwrote {af_tsv_path}")


# ---------------------------------------------------------------------------
# Show concrete Tier 4 examples: for each Tier 4 haplotype, find the new
# haplotype(s) whose variant set subsumes it, and print position tuples so
# the reader can see the (v1, v3)-vs-(v1, v2, v3) pattern in real data.
# ---------------------------------------------------------------------------
def fmt_variant(v: tuple) -> str:
    """Format a variant tuple `(contig, pos, (ref, alt))` as `pos:ref>alt`."""
    return f"{v[1]}:{v[2][0]}>{v[2][1]}"


print()
print("=== Tier 4 examples: original-only haplotype with its subsuming new haplotype(s) ===")
tier4_keys = [k for k in old_only if categorize_old_only(k) == "4_subset_of_some_new_but_not_contiguous"]
for i, k in enumerate(tier4_keys):
    set_k = set(k)
    print(f"\nTier 4 example #{i + 1}: original-only haplotype")
    print(f"  variants: [{', '.join(fmt_variant(v) for v in k)}]")
    print(f"  old AC={old_map.get(k)}, old est_af={old_est_af.get(k):.4g}")
    # Find new haplotypes whose variant set ⊇ set_k.
    subsuming = [nk for nk in new_keys if set_k <= set(nk)]
    subsuming.sort(key=lambda nk: (len(nk), nk))
    for j, nk in enumerate(subsuming[:3]):
        extras = [v for v in nk if v not in set_k]
        extras_str = ", ".join(fmt_variant(v) for v in extras)
        print(f"  subsuming new haplotype #{j + 1}: variants [{', '.join(fmt_variant(v) for v in nk)}]")
        print(f"    extra variants in new not in original: [{extras_str}]")
        print(f"    new AC={new_map.get(nk)}, new est_af={new_est_af.get(nk):.4g}")
    if len(subsuming) > 3:
        print(f"  ... and {len(subsuming) - 3} more subsuming new haplotypes")


# ---------------------------------------------------------------------------
# Tier 3 examples: for each Tier 3 haplotype, show where its component
# variants end up in the new output. Every variant is in some new haplotype,
# but they don't co-occur in a single new haplotype.
# ---------------------------------------------------------------------------
print()
print("=== Tier 3 examples: original-only haplotype with the new haplotypes its variants land in ===")
tier3_keys = [k for k in old_only if categorize_old_only(k) == "3_variants_split_across_new_haplotypes"]
# Build a quick lookup: variant -> list of new haplotypes containing it.
new_by_variant: dict[tuple, list[tuple]] = {}
for nk in new_keys:
    for v in nk:
        new_by_variant.setdefault(v, []).append(nk)

for i, k in enumerate(tier3_keys):
    print(f"\nTier 3 example #{i + 1}: original-only haplotype")
    print(f"  variants: [{', '.join(fmt_variant(v) for v in k)}]")
    print(f"  old AC={old_map.get(k)}, old est_af={old_est_af.get(k):.4g}")
    for v in k:
        new_homes = sorted(new_by_variant.get(v, []), key=len)[:2]
        print(f"  variant {fmt_variant(v)}:")
        for nk in new_homes:
            print(f"    appears in new haplotype [{', '.join(fmt_variant(vv) for vv in nk)}], AC={new_map.get(nk)}, est_af={new_est_af.get(nk):.4g}")
        if not new_homes:
            print("    (no new haplotype contains this variant — would be Tier 1/2, not 3)")


# ---------------------------------------------------------------------------
# Tier 1, 2, and 5 random examples. Use a fixed seed so re-runs are stable;
# sample without replacement up to N per tier (clipped to tier size).
# ---------------------------------------------------------------------------
TIER_SAMPLE_SIZE = 10
TIER_SAMPLE_SEED = 0

tier1_all = [k for k in old_only if categorize_old_only(k) == "1_fully_variant_disjoint"]
tier2_all = [k for k in old_only if categorize_old_only(k) == "2_partially_missing_variants"]
tier5_all = [k for k in old_only if categorize_old_only(k) == "5_proper_contiguous_subfragment_of_new"]

rng = random.Random(TIER_SAMPLE_SEED)
tier1_sample = rng.sample(tier1_all, min(TIER_SAMPLE_SIZE, len(tier1_all)))
tier2_sample = rng.sample(tier2_all, min(TIER_SAMPLE_SIZE, len(tier2_all)))
tier5_sample = rng.sample(tier5_all, min(TIER_SAMPLE_SIZE, len(tier5_all)))


def print_tier_header(name: str, n_total: int, n_sampled: int) -> None:
    print()
    print(f"=== {name} examples ({n_sampled} random of {n_total}, seed={TIER_SAMPLE_SEED}) ===")


print_tier_header(
    "Tier 1 (fully variant-disjoint from every new haplotype)",
    len(tier1_all), len(tier1_sample),
)
for i, k in enumerate(tier1_sample):
    print(f"\nTier 1 example #{i + 1}: original-only haplotype")
    print(f"  variants: [{', '.join(fmt_variant(v) for v in k)}]")
    print(f"  old AC={old_map.get(k)}, old est_af={old_est_af.get(k):.4g}")
    print("  (no variant of this haplotype appears in any new haplotype)")


print_tier_header(
    "Tier 2 (partially missing variants)",
    len(tier2_all), len(tier2_sample),
)
for i, k in enumerate(tier2_sample):
    print(f"\nTier 2 example #{i + 1}: original-only haplotype")
    print(f"  variants: [{', '.join(fmt_variant(v) for v in k)}]")
    print(f"  old AC={old_map.get(k)}, old est_af={old_est_af.get(k):.4g}")
    present = [v for v in k if v in new_all_variants]
    absent = [v for v in k if v not in new_all_variants]
    print(f"  variants present in new ({len(present)}): [{', '.join(fmt_variant(v) for v in present)}]")
    print(f"  variants absent from every new haplotype ({len(absent)}): [{', '.join(fmt_variant(v) for v in absent)}]")
    # For each present variant, show the shortest new haplotype containing it.
    for v in present:
        new_homes = sorted(new_by_variant.get(v, []), key=len)[:1]
        if new_homes:
            nk = new_homes[0]
            print(f"    {fmt_variant(v)} appears in new haplotype [{', '.join(fmt_variant(vv) for vv in nk)}], AC={new_map.get(nk)}, est_af={new_est_af.get(nk):.4g}")


print_tier_header(
    "Tier 5 (proper contiguous sub-fragment of some new haplotype)",
    len(tier5_all), len(tier5_sample),
)
for i, k in enumerate(tier5_sample):
    print(f"\nTier 5 example #{i + 1}: original-only haplotype")
    print(f"  variants: [{', '.join(fmt_variant(v) for v in k)}]")
    print(f"  old AC={old_map.get(k)}, old est_af={old_est_af.get(k):.4g}")
    # Find the shortest new haplotype that contains k as a contiguous sub-array.
    enclosing = []
    for nk in new_keys:
        if len(nk) > len(k) and is_contiguous_sub(k, nk):
            enclosing.append(nk)
    enclosing.sort(key=lambda nk: (len(nk), nk))
    if enclosing:
        # Show the shortest one (the immediate enclosing fragment).
        nk = enclosing[0]
        extras = [v for v in nk if v not in set(k)]
        print(f"  enclosing new haplotype: [{', '.join(fmt_variant(v) for v in nk)}]")
        print(f"    extra variants in new: [{', '.join(fmt_variant(v) for v in extras)}]")
        print(f"    new AC={new_map.get(nk)}, new est_af={new_est_af.get(nk):.4g}")
        if len(enclosing) > 1:
            print(f"    (+ {len(enclosing) - 1} more enclosing new haplotypes)")
