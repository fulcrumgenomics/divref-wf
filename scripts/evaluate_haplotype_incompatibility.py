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

from divref.haplotype_compat import REASONS
from divref.haplotype_compat import classify_haplotype
from divref.haplotype_compat import count_bypass_resolutions
from divref.haplotype_compat import end_coordinate_shortfall
from divref.haplotype_compat import parse_variants_string
from divref.haplotype_compat import start_coordinate_shortfall

DEFAULT_DUCKDB = "data/work/output/hgdp_1kg.haplotypes_gnomad_merge.index.duckdb"
DEFAULT_OUT_DIR = Path("data/analysis/haplotype_incompatibility")

# One HGDP_haplotype query row: (contig, n_variants, variants, popmax_empirical_AC, start, end).
QueryRow = tuple[str, int, str, int, int, int]


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
    """
    Scan the HGDP_haplotype rows of a DivRef index and aggregate both evaluation dimensions.

    Reads every `sequences` row with `source = 'HGDP_haplotype'` and `n_variants >= 2`, classifies
    each haplotype's adjacent-variant incompatibilities (Dimension A) and its end-coordinate
    undershoot (Dimension B), and accumulates per-reason counts, the length / bypass-resolution
    distributions, and up to `examples_per` sampled rows per bucket.

    Args:
        con: Open DuckDB connection to a DivRef index (read-only is sufficient).
        examples_per: Maximum number of sampled example rows retained per bucket.
        contigs: Restrict the scan to these contigs; None scans all contigs.

    Returns:
        A populated `Summary`.

    Raises:
        duckdb.Error: If the `window_size` or `sequences` query fails (e.g. a missing table).
        ValueError: If a row's `variants` string is malformed (from `parse_variants_string`).
    """
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

    def _add_example(bucket: str, row: QueryRow, shortfall: int = 0) -> None:
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
    """
    Write the machine-readable `metric<TAB>value` summary TSV.

    Args:
        summary: The aggregated results to serialise.
        path: Output path for the summary TSV; its parent directory must already exist.

    Returns:
        None.

    Raises:
        OSError: If `path` cannot be opened for writing.
    """
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
    """
    Write the sampled per-bucket example rows as a TSV for manual review.

    Emits one row per retained example with a `dimension` column (`A` for incompatibility-reason
    buckets, `B` for the coordinate buckets) plus the example's contig, variant, and coordinate
    fields.

    Args:
        summary: The aggregated results whose `examples` are written.
        path: Output path for the examples TSV; its parent directory must already exist.

    Returns:
        None.

    Raises:
        OSError: If `path` cannot be opened for writing.
    """
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
    """
    Print a human-readable summary to stdout, surfacing catch-all / anomaly buckets loudly.

    Args:
        summary: The aggregated results to print.

    Returns:
        None.
    """
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
    """
    CLI entry point: evaluate a DivRef DuckDB index and write the summary + examples TSVs.

    Parses `--duckdb`, `--out-dir`, `--examples-per-reason`, and `--contigs`, opens the index
    read-only, runs `load_and_classify`, writes both TSVs into the output directory, and prints
    the summary to stdout.

    Returns:
        None.

    Raises:
        SystemExit: Exits with status 1 if the index cannot be opened read-only (e.g. it is
            locked by another process).
    """
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
