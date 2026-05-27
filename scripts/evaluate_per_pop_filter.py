"""Evaluate the impact of changing the haplotype inclusion filter from a single
scalar `estimated_gnomad_AF >= threshold` to a per-population threshold
(`max over pops of estimated_gnomad_AF_by_pop >= threshold`).

Expected inputs:
  NEW_HT          - compute_haplotypes output Hail table, generated with the
                    `estimated_gnomad_AF_by_pop` column populated and the
                    haplotype-AF threshold lowered (e.g. 0.0) so we can see
                    haplotypes the current scalar filter would have dropped.
  OLD_DUCKDB      - original-pipeline DuckDB index (same dataset as the Hail
                    table; used to cross-tabulate against the Tier 1-5 buckets).
  HAPLOTYPE_FREQ_THRESHOLD - the inclusion threshold (default 0.005).

Outputs (stdout):
  - Counts of haplotypes that pass under each filter rule.
  - Cross-tabulation of rescued haplotypes against the old/new partition
    (Shared, New-only, Original-only Tiers 1-5).
  - Per-population breakdown of which pop is rescuing each rescued haplotype.
  - A handful of example rescued haplotypes with full per-pop est_af vectors.
"""

from collections import Counter
from pathlib import Path

import duckdb
import hail as hl

OLD_DUCKDB = (
    "data/analysis/compute_haplotypes/test_data_old/hgdp_1kg.haplotypes_gnomad_merge.index.duckdb"
)
# Re-generated with the per-pop est_af patch + threshold lowered. Path can be
# overridden by the caller; the default mirrors the file that
# `compare_haplotypes.py` uses.
NEW_HT = (
    "data/analysis/compute_haplotypes/test_data_new/hgdp_1kg.haplotypes.chr22.ht"
)
CONTIG = "chr22"
HAPLOTYPE_FREQ_THRESHOLD = 0.005
N_EXAMPLES = 10
# Populations the production DivRef bundle filters on. The variant annotation table used
# to generate the test fixture may include additional gnomAD populations (e.g. `mid`,
# `fin`, `oth`); restricting the per-pop rescue to these 5 keeps the analysis aligned
# with what the production bundle actually emits.
DIVREF_POPS = {"afr", "amr", "eas", "sas", "nfe"}


def parse_variants_string(s: str) -> tuple:
    out = []
    for token in s.split(","):
        contig, pos, ref, alt = token.split(":")
        out.append((contig, int(pos), (ref, alt)))
    return tuple(out)


def variant_tuple_from_ht_row(row: hl.Struct) -> tuple:
    return tuple(
        (v.locus.contig, v.locus.position, tuple(v.alleles)) for v in row.variants
    )


def fmt_variant(v: tuple) -> str:
    return f"{v[1]}:{v[2][0]}>{v[2][1]}"


# Load original-pipeline DuckDB rows so we can cross-tab rescued haplotypes
# against the existing Tier 1-5 categorization.
con = duckdb.connect(str(Path(OLD_DUCKDB).resolve()), read_only=True)
old_rows = con.execute(
    "SELECT variants, popmax_empirical_AC, estimated_gnomad_AF "
    "FROM sequences WHERE source = 'HGDP_haplotype' AND contig = ?",
    [CONTIG],
).fetchall()
con.close()

old_keys: set[tuple] = set()
for variants_str, _, _ in old_rows:
    old_keys.add(parse_variants_string(variants_str))


# Load new HT.
hl.init(quiet=True)
new = hl.read_table(NEW_HT)
new_rows = new.collect()

# Build a dict: haplotype-key -> (scalar_est_af, est_af_by_pop, pops_legend, length).
# Also collect a globals.pops so we can label the rescuing pop.
pops_legend: list[str] = new.globals.pops.collect()[0]

# Read the per-pop quantities out of `all_pop_freqs`, which after the option-(A) patch
# bundles (empirical_AC, empirical_AF, fraction_phased, estimated_gnomad_AF) per pop.
# `all_pop_freqs` is sorted by `empirical_AF` desc; for the rescue analysis we index by
# `pop` so we can align across haplotypes and against the legend.
new_data: dict[tuple, dict] = {}
for r in new_rows:
    key = variant_tuple_from_ht_row(r)
    by_pop_est_af: list[float | None] = [None] * len(pops_legend)
    for entry in r.all_pop_freqs:
        by_pop_est_af[entry.pop] = entry.estimated_gnomad_AF
    new_data[key] = {
        "scalar_est_af": r.estimated_gnomad_AF,
        "by_pop_est_af": by_pop_est_af,
        "max_empirical_AC": r.max_empirical_AC,
        "length": len(key),
    }

print(f"Loaded {len(new_data)} haplotypes from new HT")
print(f"Pops legend: {pops_legend}")


# Apply both filter rules.
def passes_scalar(d: dict) -> bool:
    return d["scalar_est_af"] is not None and d["scalar_est_af"] >= HAPLOTYPE_FREQ_THRESHOLD


def passes_per_pop(d: dict, allowed_pop_indices: set[int]) -> bool:
    by_pop = [
        v for p, v in enumerate(d["by_pop_est_af"]) if v is not None and p in allowed_pop_indices
    ]
    if not by_pop:
        return False
    return max(by_pop) >= HAPLOTYPE_FREQ_THRESHOLD


keys = set(new_data.keys())
scalar_pass = {k for k in keys if passes_scalar(new_data[k])}

divref_pop_indices = {i for i, p in enumerate(pops_legend) if p in DIVREF_POPS}
all_pop_indices = set(range(len(pops_legend)))
per_pop_divref_pass = {k for k in keys if passes_per_pop(new_data[k], divref_pop_indices)}
per_pop_all_pass = {k for k in keys if passes_per_pop(new_data[k], all_pop_indices)}

rescued_divref = per_pop_divref_pass - scalar_pass
rescued_all = per_pop_all_pass - scalar_pass
dropped_divref = scalar_pass - per_pop_divref_pass

print()
print(f"=== Threshold filter comparison (threshold = {HAPLOTYPE_FREQ_THRESHOLD}) ===")
print(f"  DivRef pops considered for per-pop filter: {sorted(DIVREF_POPS)}")
print(f"  Extra gnomAD pops in HT (not used by DivRef): "
      f"{sorted(set(pops_legend) - DIVREF_POPS)}")
print(f"  passes scalar (current):                       {len(scalar_pass)}")
print(f"  passes per-pop, DivRef pops only (proposed):   {len(per_pop_divref_pass)}")
print(f"  passes per-pop, all pops in HT:                {len(per_pop_all_pass)}")
print(f"  rescued by per-pop (DivRef pops only):         {len(rescued_divref)}")
print(f"  rescued by per-pop (all pops in HT):           {len(rescued_all)}")
print(f"  dropped by per-pop (DivRef pops):              {len(dropped_divref)}  (should be 0)")

# Continue downstream analyses using the DivRef-pops-only rescue set, since that's
# what would actually change the production bundle.
rescued = rescued_divref


# Cross-tab the rescued set against old DuckDB membership.
rescued_in_old = rescued & old_keys
rescued_not_in_old = rescued - old_keys

print()
print("=== Rescued haplotypes vs original pipeline output ===")
print(f"  rescued AND in original DuckDB: {len(rescued_in_old)}  (i.e., currently 'Tier *' / Shared)")
print(f"  rescued but NOT in original DuckDB: {len(rescued_not_in_old)}  (genuinely new)")


# Per-pop attribution: for each rescued haplotype, identify the DivRef pop with the
# highest per-pop est_af (the "rescuing pop") and report the count breakdown.
rescuing_pop_counts: Counter = Counter()
for k in rescued:
    by_pop = new_data[k]["by_pop_est_af"]
    best = -1.0
    best_p = None
    for p, v in enumerate(by_pop):
        if v is None or p not in divref_pop_indices:
            continue
        if v > best:
            best = v
            best_p = p
    if best_p is not None:
        rescuing_pop_counts[pops_legend[best_p]] += 1

print()
print("=== Rescuing population breakdown (DivRef pops only) ===")
for pop, count in sorted(rescuing_pop_counts.items(), key=lambda t: -t[1]):
    print(f"  {pop}: {count}")


# Length distribution of rescued haplotypes.
print()
print("=== Length distribution of rescued haplotypes ===")
length_counts = Counter(new_data[k]["length"] for k in rescued)
for n in sorted(length_counts):
    print(f"  length {n}: {length_counts[n]}")


# Show a few example rescued haplotypes with their per-pop est_af vectors (DivRef pops only).
def best_divref_est_af(k: tuple) -> float:
    by_pop = new_data[k]["by_pop_est_af"]
    vals = [v for p, v in enumerate(by_pop) if v is not None and p in divref_pop_indices]
    return max(vals) if vals else 0.0


print()
print(f"=== First {N_EXAMPLES} rescued haplotypes (sorted by max DivRef per-pop est_af desc) ===")
rescued_sorted = sorted(rescued, key=lambda k: -best_divref_est_af(k))
for i, k in enumerate(rescued_sorted[:N_EXAMPLES]):
    d = new_data[k]
    by_pop_str = ", ".join(
        f"{pops_legend[p]}={v:.4g}" if v is not None else f"{pops_legend[p]}=NA"
        for p, v in enumerate(d["by_pop_est_af"])
    )
    in_old = "in original DuckDB" if k in old_keys else "genuinely new"
    print(f"\nrescued example #{i + 1} ({in_old}):")
    print(f"  variants: [{', '.join(fmt_variant(v) for v in k)}]")
    print(f"  scalar est_af: {d['scalar_est_af']:.4g}")
    print(f"  per-pop est_af: {by_pop_str}")
    print(f"  AC (max pop): {d['max_empirical_AC']}")
