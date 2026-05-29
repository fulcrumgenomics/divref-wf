"""Count per-chromosome (v1, v2, v3) genotype combinations at the Tier 4 example
positions on chr22, aggregated across all samples that pass the workflow's
pop/sex filters.

Positions (1-based):
  v1 = chr22:24626868 C>T
  v2 = chr22:24626884 T>C
  v3 = chr22:24626892 A>G

Output:
  - Table of all 8 (v1, v2, v3) ∈ {ref, alt}^3 combinations with per-population
    counts, summed over (sample, strand) pairs (one chromosome per strand).
  - Sanity check: AC for each individual variant.
"""

from pathlib import Path

import hail as hl

from divref.hail import hail_init

CONTIG = "chr22"
POSITIONS: dict[str, int] = {
    "v1": 24626868,
    "v2": 24626884,
    "v3": 24626892,
}
PHASED_BCF = (
    "gs://gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2/"
    f"hgdp1kgp_{CONTIG}.filtered.SNV_INDEL.phased.shapeit5.bcf"
)
SAMPLE_HT = "data/work/inputs/hgdp_1kg.sample_meta.extract.ht"
DIVREF_POPS = {"afr", "amr", "eas", "sas", "nfe"}

hail_init(
    gcs_credentials_path=Path("~/.config/gcloud/application_default_credentials.json").expanduser()
)

# Read sample meta to filter to (pop ∈ DIVREF_POPS, sex_karyotype ∈ {XX, XY})
sa = hl.read_table(SAMPLE_HT)
sa = sa.filter(hl.literal(DIVREF_POPS).contains(sa.pop))

# Import the BCF for just the three positions.
intervals = [
    hl.parse_locus_interval(
        f"{CONTIG}:{min(POSITIONS.values())}-{max(POSITIONS.values()) + 1}",
        reference_genome="GRCh38",
    )
]
mt = hl.import_vcf(
    PHASED_BCF,
    reference_genome="GRCh38",
    force_bgz=True,
    array_elements_required=False,
)
mt = hl.filter_intervals(mt, intervals)
mt = mt.filter_rows(hl.literal(set(POSITIONS.values())).contains(mt.locus.position))

# Attach pop and filter to the kept samples.
mt = mt.annotate_cols(pop=sa[mt.col_key].pop)
mt = mt.filter_cols(hl.is_defined(mt.pop))

# Collect each variant's (left, right) genotypes per sample.
mt = mt.annotate_entries(
    left=mt.GT[0],
    right=mt.GT[1],
)

# Build a per-(sample, strand) row indexed by position.
rows = mt.entries().select("locus", "pop", "left", "right").collect()

# Aggregate: for each (sample, strand), one chromosome.
# combo[v_name] = 0 (ref) or 1 (alt)
import itertools
from collections import Counter, defaultdict

per_chrom_combos: dict[tuple[str, int], dict[str, int]] = defaultdict(dict)
sample_pops: dict[str, str] = {}
pos_to_vname = {pos: name for name, pos in POSITIONS.items()}

for r in rows:
    s = r.s
    pop = r.pop
    sample_pops[s] = pop
    vname = pos_to_vname[r.locus.position]
    per_chrom_combos[(s, 0)][vname] = int(r.left)
    per_chrom_combos[(s, 1)][vname] = int(r.right)

combo_counts_overall: Counter = Counter()
combo_counts_per_pop: dict[str, Counter] = defaultdict(Counter)

for (s, strand), combo in per_chrom_combos.items():
    # Skip chromosomes lacking a genotype at any of the three positions.
    if any(combo.get(v) is None for v in POSITIONS):
        continue
    key = tuple(combo[v] for v in ("v1", "v2", "v3"))
    combo_counts_overall[key] += 1
    combo_counts_per_pop[sample_pops[s]][key] += 1


def fmt_key(k: tuple[int, int, int]) -> str:
    return " ".join("alt" if x else "ref" for x in k)


print(f"=== {CONTIG} (v1, v2, v3) chromosome counts ===")
print(f"  v1={POSITIONS['v1']}, v2={POSITIONS['v2']}, v3={POSITIONS['v3']}")
print(f"  total chromosomes after filters: {sum(combo_counts_overall.values())}")
print()
print(f"{'v1':>5} {'v2':>5} {'v3':>5}   total  " + "  ".join(
    f"{p:>5}" for p in sorted(DIVREF_POPS)
))
all_combos = list(itertools.product([0, 1], repeat=3))
for k in sorted(all_combos, key=lambda c: combo_counts_overall.get(c, 0), reverse=True):
    cells = "  ".join(
        f"{combo_counts_per_pop[p].get(k, 0):>5}" for p in sorted(DIVREF_POPS)
    )
    print(
        f"{('alt' if k[0] else 'ref'):>5} "
        f"{('alt' if k[1] else 'ref'):>5} "
        f"{('alt' if k[2] else 'ref'):>5}  "
        f"{combo_counts_overall.get(k, 0):>6}  {cells}"
    )

print()
print("=== Single-variant AC (alt count, all chromosomes) ===")
for vname in ("v1", "v2", "v3"):
    idx = ("v1", "v2", "v3").index(vname)
    ac = sum(c for k, c in combo_counts_overall.items() if k[idx] == 1)
    print(f"  {vname} ({POSITIONS[vname]}): AC = {ac}")
