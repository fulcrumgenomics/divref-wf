"""Evaluate HGDP_haplotype rows in a DivRef DuckDB index for variant problems.

Two dimensions, reported with per-category counts and sampled examples so each can be addressed:

  A. Incompatibility taxonomy: adjacent component variants whose reference spans overlap
     (`variant_distance < 0`) cannot co-occur on one chromosome. Each incompatible adjacent
     pair is classified into a reason bucket (e.g. a SNP inside a deletion).
  B. End-coordinate adequacy: haplotypes whose stored `end` window is too short to span the
     full reference allele of an earlier deletion, so the deleted bases (and their trailing
     context) are truncated out of the FASTA record and the `end` column.

Dimension B is a strict subset of Dimension A (an end-undershoot implies an overlap), but it
isolates the cases where the emitted output is actually corrupted. See
docs/superpowers/specs/2026-06-01-haplotype-variant-incompatibility-eval-design.md.

Pure Python + DuckDB, read-only, no Spark.

    pixi run python scripts/evaluate_haplotype_incompatibility.py                 # whole genome
    pixi run python scripts/evaluate_haplotype_incompatibility.py --contigs chr22 # one contig
"""

import argparse
import sys
from collections import Counter
from collections import defaultdict
from dataclasses import dataclass
from dataclasses import field
from pathlib import Path

import duckdb

DEFAULT_DUCKDB = "data/work/output/hgdp_1kg.haplotypes_gnomad_merge.index.duckdb"
DEFAULT_OUT_DIR = Path("data/analysis/haplotype_incompatibility")

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
    """Parse a comma-separated `chr:pos:ref:alt` list, sorted by (pos, len(ref))."""
    out: list[Variant] = []
    for token in s.split(","):
        contig, pos, ref, alt = token.split(":")
        out.append((contig, int(pos), ref, alt))
    out.sort(key=lambda v: (v[1], len(v[2])))
    return out


def variant_distance(v1: Variant, v2: Variant) -> int:
    """Reference bases between v1 and v2; < 0 means their reference spans overlap."""
    return v2[1] - v1[1] - len(v1[2])


def classify_pair(v1: Variant, v2: Variant) -> str | None:
    """Reason label for an incompatible adjacent pair, or None if compatible (distance >= 0)."""
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
    """Reasons from every incompatible consecutive pair (position-sorted); may repeat / be empty.

    Adjacent pairs are sufficient: positions are non-decreasing and each reference allele is
    >= 1 bp, so if no consecutive pair overlaps then no pair overlaps at all. A long deletion
    swallowing a non-adjacent downstream variant is still flagged at its own adjacent boundary.
    """
    reasons: list[str] = []
    for v1, v2 in zip(variants, variants[1:]):
        reason = classify_pair(v1, v2)
        if reason is not None:
            reasons.append(reason)
    return reasons


def variants_overlap(v1: Variant, v2: Variant) -> bool:
    """True iff two variants' reference spans overlap (order-independent)."""
    earlier, later = (v1, v2) if v1[1] <= v2[1] else (v2, v1)
    return later[1] < earlier[1] + len(earlier[2])


def count_bypass_resolutions(variants: list[Variant]) -> int:
    """Distinct maximal pairwise-compatible sub-haplotypes (>= 2 variants) for this haplotype.

    Models the "explode the conflict" alternative to dropping: a resolution keeps a maximal set
    of variants with no two overlapping (a maximal independent set in the overlap graph), and is
    a valid haplotype only if it retains >= 2 variants. Returns the count of such resolutions.

    Interpretation for an INCOMPATIBLE haplotype:
      - 0  -> exploding recovers nothing (e.g. a length-2 conflict: both resolutions are
              singletons, which are not haplotypes and are already in the gnomAD_variant track).
      - 1  -> a single unambiguous resolution.
      - >1 -> exploding would spawn this many speculative candidates from one (mis-phased) block;
              larger values quantify the combinatorial blow-up at multi-conflict repeat loci.
    """
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
    """bp by which the stored 0-based-exclusive `end` fails to cover every variant's reference
    span plus window_size of trailing context. > 0 means deleted reference is truncated.

    The stored `end` is set from the last-by-position variant, but an earlier deletion can reach
    further right. window_size cancels in the comparison, but is included so the result is the
    literal bp shortfall against the stored value.
    """
    rightmost_ref_end = max(v[1] - 1 + len(v[2]) for v in variants)  # 0-based exclusive
    return rightmost_ref_end + window_size - stored_end


def start_coordinate_shortfall(variants: list[Variant], window_size: int, stored_start: int) -> int:
    """bp by which the stored 0-based-inclusive `start` is too far right (should always be 0).

    Reference alleles only extend rightward, so the leftmost touched base is the min-position
    variant and `start` is structurally correct. This is a defensive check.
    """
    leftmost_ref_start = min(v[1] - 1 for v in variants)  # 0-based inclusive
    required_start = leftmost_ref_start - window_size
    return stored_start - required_start


@dataclass
class Summary:
    """Aggregated evaluation results across all scanned haplotypes."""

    window_size: int
    total_haplotypes: int = 0
    haplotypes_with_any_incompatibility: int = 0
    reason_haplotype_counts: Counter = field(default_factory=Counter)
    reason_pair_counts: Counter = field(default_factory=Counter)
    reason_contig_counts: dict[str, Counter] = field(default_factory=lambda: defaultdict(Counter))
    end_undershoot_count: int = 0
    end_undershoot_by_reason: Counter = field(default_factory=Counter)
    end_undershoot_max_bp: int = 0
    end_undershoot_hist: Counter = field(default_factory=Counter)
    start_undershoot_anomalies: int = 0
    # Length distribution of incompatible haplotypes, and the "explode instead of drop" sizing.
    incompatible_length_hist: Counter = field(default_factory=Counter)
    bypass_resolutions_hist: Counter = field(default_factory=Counter)
    bypass_resolutions_total: int = 0
    bypass_resolutions_max: int = 0
    recoverable_haplotypes: int = 0  # incompatible haplotypes with >= 1 bypass resolution
    examples: dict[str, list[tuple]] = field(default_factory=lambda: defaultdict(list))


def _read_window_size(con: duckdb.DuckDBPyConnection) -> int:
    return int(con.execute("SELECT window_size FROM window_size").fetchone()[0])


def load_and_classify(
    con: duckdb.DuckDBPyConnection,
    examples_per: int = 20,
    contigs: list[str] | None = None,
) -> Summary:
    """Scan HGDP_haplotype rows and aggregate both evaluation dimensions."""
    summary = Summary(window_size=_read_window_size(con))

    query = (
        'SELECT contig, n_variants, variants, popmax_empirical_AC, start, "end" '
        "FROM sequences WHERE source = 'HGDP_haplotype' AND n_variants >= 2"
    )
    params: list[str] = []
    if contigs:
        placeholders = ",".join(["?"] * len(contigs))
        query += f" AND contig IN ({placeholders})"  # noqa: S608
        params = contigs

    def _add_example(bucket: str, row: tuple, shortfall: int = 0) -> None:
        if len(summary.examples[bucket]) < examples_per:
            contig, n_variants, variants_str, ac, start, end = row
            summary.examples[bucket].append(
                (bucket, contig, n_variants, ac, start, end, shortfall, variants_str)
            )

    cur = con.execute(query, params)
    while True:
        batch = cur.fetchmany(50_000)
        if not batch:
            break
        for row in batch:
            contig, _n_variants, variants_str, _ac, start, end = row
            variants = parse_variants_string(variants_str)
            summary.total_haplotypes += 1

            reasons = classify_haplotype(variants)
            distinct_reasons = set(reasons)
            if reasons:
                summary.haplotypes_with_any_incompatibility += 1
                summary.incompatible_length_hist[len(variants)] += 1
                resolutions = count_bypass_resolutions(variants)
                summary.bypass_resolutions_hist[resolutions] += 1
                summary.bypass_resolutions_total += resolutions
                summary.bypass_resolutions_max = max(summary.bypass_resolutions_max, resolutions)
                if resolutions >= 1:
                    summary.recoverable_haplotypes += 1
                for reason in reasons:
                    summary.reason_pair_counts[reason] += 1
                for reason in distinct_reasons:
                    summary.reason_haplotype_counts[reason] += 1
                    summary.reason_contig_counts[reason][contig] += 1
                    _add_example(reason, row)

            end_short = end_coordinate_shortfall(variants, summary.window_size, end)
            if end_short > 0:
                summary.end_undershoot_count += 1
                summary.end_undershoot_max_bp = max(summary.end_undershoot_max_bp, end_short)
                summary.end_undershoot_hist[end_short] += 1
                for reason in distinct_reasons or {"none"}:
                    summary.end_undershoot_by_reason[reason] += 1
                _add_example("end_undershoot", row, shortfall=end_short)

            start_short = start_coordinate_shortfall(variants, summary.window_size, start)
            if start_short > 0:
                summary.start_undershoot_anomalies += 1
                _add_example("start_anomaly", row, shortfall=start_short)

    return summary


def write_summary_tsv(summary: Summary, path: Path) -> None:
    """Write the machine-readable `metric<TAB>value` summary."""
    rows: list[tuple[str, int]] = [
        ("window_size", summary.window_size),
        ("total_haplotypes", summary.total_haplotypes),
        ("haplotypes_with_any_incompatibility", summary.haplotypes_with_any_incompatibility),
    ]
    for reason in REASONS:
        rows.append((f"haplotype_count.{reason}", summary.reason_haplotype_counts.get(reason, 0)))
    for reason in REASONS:
        rows.append((f"pair_count.{reason}", summary.reason_pair_counts.get(reason, 0)))
    rows.append(("haplotypes_with_end_undershoot", summary.end_undershoot_count))
    rows.append(("end_undershoot_max_bp", summary.end_undershoot_max_bp))
    rows.append(("start_undershoot_anomalies", summary.start_undershoot_anomalies))
    for reason in REASONS:
        rows.append(
            (f"end_undershoot.{reason}", summary.end_undershoot_by_reason.get(reason, 0))
        )

    # "Explode instead of drop" sizing.
    n_len2 = summary.incompatible_length_hist.get(2, 0)
    rows.append(("incompatible_length2", n_len2))
    rows.append(("incompatible_length_ge3", summary.haplotypes_with_any_incompatibility - n_len2))
    rows.append(("recoverable_by_explode", summary.recoverable_haplotypes))
    rows.append(("bypass_resolutions_total", summary.bypass_resolutions_total))
    rows.append(("bypass_resolutions_max", summary.bypass_resolutions_max))
    for length in sorted(summary.incompatible_length_hist):
        rows.append((f"incompatible_length.{length}", summary.incompatible_length_hist[length]))

    with path.open("w") as f:
        f.write("metric\tvalue\n")
        for name, value in rows:
            f.write(f"{name}\t{value}\n")


def write_examples_tsv(summary: Summary, path: Path) -> None:
    """Write sampled example rows per bucket for manual review."""
    with path.open("w") as f:
        f.write("dimension\tbucket\tcontig\tn_variants\tpopmax_empirical_AC\tstart\tend\t")
        f.write("shortfall_bp\tvariants\n")
        for bucket, rows in summary.examples.items():
            dimension = "B" if bucket in ("end_undershoot", "start_anomaly") else "A"
            for ex in rows:
                _bucket, contig, n_variants, ac, start, end, shortfall, variants_str = ex
                f.write(
                    f"{dimension}\t{bucket}\t{contig}\t{n_variants}\t{ac}\t{start}\t{end}\t"
                    f"{shortfall}\t{variants_str}\n"
                )


def print_summary(summary: Summary) -> None:
    """Print a human-readable summary, surfacing catch-all / anomaly buckets loudly."""
    print(f"window_size: {summary.window_size}")
    print(f"total HGDP_haplotype rows (n_variants >= 2): {summary.total_haplotypes}")
    print(f"haplotypes with any incompatibility:         {summary.haplotypes_with_any_incompatibility}")
    print("\n=== Dimension A: incompatibility reasons (haplotype counts / pair counts) ===")
    for reason in REASONS:
        print(
            f"  {reason:28s} {summary.reason_haplotype_counts.get(reason, 0):>10d}"
            f" / {summary.reason_pair_counts.get(reason, 0):>10d}"
        )
    print("\n=== Dimension B: end-coordinate undershoot ===")
    print(f"  haplotypes with end undershoot: {summary.end_undershoot_count}")
    print(f"  max shortfall (bp):             {summary.end_undershoot_max_bp}")
    print(f"  by incompatibility reason:      {dict(summary.end_undershoot_by_reason)}")
    if summary.start_undershoot_anomalies:
        print(f"  !! start-coordinate anomalies (expected 0): {summary.start_undershoot_anomalies}")

    n_len2 = summary.incompatible_length_hist.get(2, 0)
    n_ge3 = summary.haplotypes_with_any_incompatibility - n_len2
    print("\n=== Explode-vs-drop sizing (incompatible haplotypes) ===")
    print(f"  length distribution: {dict(sorted(summary.incompatible_length_hist.items()))}")
    print(f"  length-2 (explode recovers nothing):     {n_len2}")
    print(f"  length>=3 (explode could yield bypasses): {n_ge3}")
    print(f"  recoverable by explode (>=1 resolution):  {summary.recoverable_haplotypes}")
    print(f"  total bypass candidates exploding spawns: {summary.bypass_resolutions_total}")
    print(f"  max bypass candidates from one haplotype: {summary.bypass_resolutions_max}")
    print(f"  bypass-count distribution: {dict(sorted(summary.bypass_resolutions_hist.items()))}")

    if summary.examples.get("other_overlap"):
        print("\n!! 'other_overlap' examples (unclassified — refine the taxonomy):")
        for ex in summary.examples["other_overlap"]:
            print(f"    {ex[1]}  {ex[-1]}")
    if summary.examples.get("start_anomaly"):
        print("\n!! 'start_anomaly' examples (unexpected left-edge truncation):")
        for ex in summary.examples["start_anomaly"]:
            print(f"    {ex[1]}  {ex[-1]}")


def main() -> None:
    """Run the evaluation against a DivRef DuckDB index and write TSV + stdout summaries."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--duckdb", default=DEFAULT_DUCKDB, help="DivRef DuckDB index path.")
    parser.add_argument("--out-dir", default=DEFAULT_OUT_DIR, type=Path, help="Output directory.")
    parser.add_argument(
        "--examples-per-reason", default=20, type=int, help="Sampled examples per bucket."
    )
    parser.add_argument(
        "--contigs", nargs="+", default=None, help="Restrict to these contigs (default: all)."
    )
    args = parser.parse_args()

    try:
        con = duckdb.connect(str(Path(args.duckdb).resolve()), read_only=True)
    except duckdb.IOException as e:
        print(
            f"ERROR: could not open {args.duckdb} read-only (is it locked by another "
            f"process?): {e}",
            file=sys.stderr,
        )
        sys.exit(1)

    try:
        summary = load_and_classify(con, args.examples_per_reason, contigs=args.contigs)
    finally:
        con.close()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_summary_tsv(summary, args.out_dir / "summary.tsv")
    write_examples_tsv(summary, args.out_dir / "examples.tsv")
    print_summary(summary)
    print(f"\nwrote {args.out_dir / 'summary.tsv'} and {args.out_dir / 'examples.tsv'}")


if __name__ == "__main__":
    main()
