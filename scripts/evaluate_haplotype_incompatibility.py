"""Audit HGDP_haplotype rows in a DivRef DuckDB index for variant incompatibilities.

Re-derives each haplotype's incompatibility reasons from its `variants` string -- the same
classification the build persists into `haplotype_filter` -- and reports the per-reason
distribution with sampled examples. This is an independent cross-check of the persisted column
(read straight from the DB) and surfaces any `other_overlap` rows the taxonomy doesn't yet name.

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
from divref.haplotype_compat import parse_variants_string

DEFAULT_DUCKDB = "data/work/output/hgdp_1kg.haplotypes_gnomad_merge.index.duckdb"
DEFAULT_OUT_DIR = Path("data/analysis/haplotype_incompatibility")


@dataclass
class Summary:
    """Aggregated incompatibility-reason counts across all scanned haplotypes."""

    total_haplotypes: int = 0
    haplotypes_with_any_incompatibility: int = 0
    # Distinct reasons per haplotype (a multi-reason haplotype counts under each of its reasons).
    reason_haplotype_counts: Counter = field(default_factory=Counter)
    reason_pair_counts: Counter = field(default_factory=Counter)  # one per overlapping pair
    examples: dict[str, list[str]] = field(default_factory=lambda: defaultdict(list))


def load_and_classify(
    con: duckdb.DuckDBPyConnection,
    examples_per: int = 20,
    contigs: list[str] | None = None,
) -> Summary:
    """
    Scan the HGDP_haplotype rows of a DivRef index and tally incompatibility reasons.

    Reads every `sequences` row with `source = 'HGDP_haplotype'` and `n_variants >= 2`, classifies
    each haplotype's adjacent-variant incompatibilities, and accumulates per-reason haplotype and
    pair counts plus up to `examples_per` sampled `variants` strings per reason.

    Args:
        con: Open DuckDB connection to a DivRef index (read-only is sufficient).
        examples_per: Maximum number of sampled example rows retained per reason.
        contigs: Restrict the scan to these contigs; None scans all contigs.

    Returns:
        A populated `Summary`.

    Raises:
        duckdb.Error: If the `sequences` query fails (e.g. a missing table).
        ValueError: If a row's `variants` string is malformed (from `parse_variants_string`).
    """
    summary = Summary()
    query = "SELECT variants FROM sequences WHERE source = 'HGDP_haplotype' AND n_variants >= 2"
    params: list[str] = []
    if contigs:
        placeholders = ",".join(["?"] * len(contigs))
        query += f" AND contig IN ({placeholders})"  # noqa: S608
        params = contigs

    cur = con.execute(query, params)
    while True:
        batch = cur.fetchmany(50_000)
        if not batch:
            break
        for (variants_str,) in batch:
            summary.total_haplotypes += 1
            reasons = classify_haplotype(parse_variants_string(variants_str))
            if not reasons:
                continue
            summary.haplotypes_with_any_incompatibility += 1
            for reason in reasons:
                summary.reason_pair_counts[reason] += 1
            for reason in set(reasons):
                summary.reason_haplotype_counts[reason] += 1
                if len(summary.examples[reason]) < examples_per:
                    summary.examples[reason].append(variants_str)
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
        ("total_haplotypes", summary.total_haplotypes),
        ("haplotypes_with_any_incompatibility", summary.haplotypes_with_any_incompatibility),
    ]
    for reason in REASONS:
        rows.append((f"haplotype_count.{reason}", summary.reason_haplotype_counts.get(reason, 0)))
    for reason in REASONS:
        rows.append((f"pair_count.{reason}", summary.reason_pair_counts.get(reason, 0)))

    with path.open("w") as f:
        f.write("metric\tvalue\n")
        for name, value in rows:
            f.write(f"{name}\t{value}\n")


def write_examples_tsv(summary: Summary, path: Path) -> None:
    """
    Write the sampled per-reason example `variants` strings as a TSV for manual review.

    Args:
        summary: The aggregated results whose `examples` are written.
        path: Output path for the examples TSV; its parent directory must already exist.

    Returns:
        None.

    Raises:
        OSError: If `path` cannot be opened for writing.
    """
    with path.open("w") as f:
        f.write("reason\tvariants\n")
        for reason, variant_strs in summary.examples.items():
            for variants_str in variant_strs:
                f.write(f"{reason}\t{variants_str}\n")


def print_summary(summary: Summary) -> None:
    """
    Print a human-readable summary to stdout, surfacing the `other_overlap` catch-all loudly.

    Args:
        summary: The aggregated results to print.

    Returns:
        None.
    """
    print(f"total HGDP_haplotype rows (n_variants >= 2): {summary.total_haplotypes}")
    print(f"haplotypes with any incompatibility: {summary.haplotypes_with_any_incompatibility}")
    print("\n=== incompatibility reasons (haplotype counts / pair counts) ===")
    for reason in REASONS:
        print(
            f"  {reason:44s} {summary.reason_haplotype_counts.get(reason, 0):>8d}"
            f" / {summary.reason_pair_counts.get(reason, 0):>8d}"
        )
    if summary.examples.get("other_overlap"):
        print("\n!! 'other_overlap' examples (unclassified -- refine the taxonomy):")
        for variants_str in summary.examples["other_overlap"]:
            print(f"    {variants_str}")


def main() -> None:
    """
    CLI entry point: audit a DivRef DuckDB index and write the summary + examples TSVs.

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
        "--examples-per-reason", default=20, type=int, help="Sampled examples per reason."
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
