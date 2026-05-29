"""Compare old vs new compute_haplotypes outputs by variant tuples.

Old haplotypes come from the DuckDB index produced by the original pipeline
(`sequences` table, `source = 'HGDP_haplotype'`).
New haplotypes come directly from the new `compute_haplotypes` Hail table.

Emits two things:

1. Stdout summary covering every count the blog cites.
2. A machine-readable summary TSV
   (`data/analysis/compute_haplotypes/algo_comparison.summary.tsv`) with one
   `metric\tvalue` row per cited count. Downstream plotting and the blog
   reproduction document read from this TSV; nothing downstream hardcodes any
   of these counts.
"""

import argparse
from collections import Counter
from pathlib import Path

import duckdb
import hail as hl

DEFAULT_OLD_DUCKDB = (
    "data/analysis/input/DivRef-v1.1.haplotypes_gnomad_merge.index.duckdb"
)
DEFAULT_NEW_HT = "data/work/haplotypes/hgdp_1kg.haplotypes.chr22.ht"
CONTIG = "chr22"
SUMMARY_TSV = Path("data/analysis/compute_haplotypes/algo_comparison.summary.tsv")

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--old-duckdb",
    default=DEFAULT_OLD_DUCKDB,
    help=(
        "Path to the original-pipeline DuckDB index (defaults to the file "
        "produced by `workflows/compare_divref_gnomad.smk`'s "
        "`download_divref_index` rule)."
    ),
)
parser.add_argument(
    "--new-ht",
    default=DEFAULT_NEW_HT,
    help="Path to the new algorithm's compute_haplotypes Hail table for the contig.",
)
args = parser.parse_args()
OLD_DUCKDB = args.old_duckdb
NEW_HT = args.new_ht


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


def is_contiguous_sub(short: tuple, long_: tuple) -> bool:
    """True iff `short` is a strict contiguous sub-array of `long_`."""
    n, m = len(short), len(long_)
    if n >= m:
        return False
    for i in range(m - n + 1):
        if long_[i : i + n] == short:
            return True
    return False


# Load old (DuckDB) -------------------------------------------------------------
con = duckdb.connect(str(Path(OLD_DUCKDB).resolve()), read_only=True)
old_query = con.execute(
    "SELECT variants, popmax_empirical_AC "
    "FROM sequences WHERE source = 'HGDP_haplotype' AND contig = ?",
    [CONTIG],
).fetchall()
con.close()

old_map: dict[tuple, int] = {}
for variants_str, ac in old_query:
    key = parse_variants_string(variants_str)
    # Defensive: keep the larger AC if a haplotype shows up twice.
    if key in old_map:
        old_map[key] = max(old_map[key], ac)
    else:
        old_map[key] = ac

# Load new (Hail table) ---------------------------------------------------------
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
print(f"Shared:    {len(shared)}")
print(f"Old only:  {len(old_only)}")
print(f"New only:  {len(new_only)}")


# Old-only sub-fragment analysis -----------------------------------------------
#
# For each old-only haplotype, find any new haplotype that strictly contains it
# as a contiguous sub-array.  When such an enclosing new exists, ask whether
# any of those enclosing new haplotypes is itself present in the old catalog.
#
# - both_have_longer:    enclosing new ALSO emitted by old -> both algorithms
#                        find the longer fragment; only new recognises the
#                        short one as redundant.
# - only_new_has_longer: enclosing new NOT in old -> old emitted the short
#                        fragment but never accumulated enough AC for the
#                        longer one to pass its AF filter, while new's
#                        containment counting credited the longer fragment
#                        with carriers from every parent block.
new_keys_by_len = sorted(new_keys, key=len, reverse=True)
old_only_subfragment_of_new = 0
old_only_subfragment_both_have_longer = 0
old_only_subfragment_only_new_has_longer = 0
for k in old_only:
    enclosing_news = [
        nk for nk in new_keys_by_len if len(nk) > len(k) and is_contiguous_sub(k, nk)
    ]
    if not enclosing_news:
        continue
    old_only_subfragment_of_new += 1
    if any(nk in old_keys for nk in enclosing_news):
        old_only_subfragment_both_have_longer += 1
    else:
        old_only_subfragment_only_new_has_longer += 1

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


# New-only decomposition --------------------------------------------------------
#
# Partition the new-only haplotypes by how their variants relate to the old
# catalog's variant universe.  Also report the length distribution of the
# all-novel bucket so the blog's "21 of 36 are 2-variant, max 8" is derivable.
old_all_variants: set = set()
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
print(
    f"All variants present somewhere in old:        {new_only_all_variants_in_old}"
)
print(
    f"Mix of shared and novel variants:             {new_only_mixed_variants}"
)
print(
    f"All variants novel (none in any old haplo):   {new_only_all_novel}"
)
print(f"  length distribution: {dict(sorted(all_novel_length_counts.items()))}")
print(f"  length-2 count:      {new_only_all_novel_length_2}")
print(f"  max length:          {new_only_all_novel_max_length}")


# AC comparison for shared haplotypes ------------------------------------------
shared_same_ac = 0
shared_new_higher_ac = 0
shared_old_higher_ac = 0
for k in shared:
    old_ac = old_map[k]
    new_ac = new_map[k]
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


# Summary TSV ------------------------------------------------------------------
SUMMARY_TSV.parent.mkdir(parents=True, exist_ok=True)
rows: list[tuple[str, int]] = [
    ("shared", len(shared)),
    ("old_only", len(old_only)),
    ("new_only", len(new_only)),
    ("old_only_subfragment_of_new", old_only_subfragment_of_new),
    ("old_only_subfragment_both_have_longer", old_only_subfragment_both_have_longer),
    (
        "old_only_subfragment_only_new_has_longer",
        old_only_subfragment_only_new_has_longer,
    ),
    ("new_only_all_variants_in_old", new_only_all_variants_in_old),
    ("new_only_mixed_variants", new_only_mixed_variants),
    ("new_only_all_novel", new_only_all_novel),
    ("new_only_all_novel_length_2", new_only_all_novel_length_2),
    ("new_only_all_novel_max_length", new_only_all_novel_max_length),
    ("shared_same_ac", shared_same_ac),
    ("shared_new_higher_ac", shared_new_higher_ac),
    ("shared_old_higher_ac", shared_old_higher_ac),
]
with SUMMARY_TSV.open("w") as f:
    f.write("metric\tvalue\n")
    for name, value in rows:
        f.write(f"{name}\t{value}\n")
print()
print(f"wrote summary TSV: {SUMMARY_TSV}")
