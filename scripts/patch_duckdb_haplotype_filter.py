"""Migrate a DivRef DuckDB built before the haplotype_filter / end-coordinate fix into a copy.

This is a one-off, Hail-free migration of an index produced by the OLD `append_contig_to_duckdb_index`
(no `haplotype_filter` column, last-variant-only `end`). It writes a fixed COPY (leaving the source
untouched) so the original and fixed indexes can be diffed as an integration test, without re-running
the workflow. The permanent fix lives in `append_contig_to_duckdb_index`; this reuses the same
`divref.haplotype_compat` logic, so the patched copy is identical to what a rebuild would produce:

  - adds `haplotype_filter` (VCF-style: PASS, else the ';'-joined incompatibility reason), and
  - corrects `end` for haplotypes whose deletion span was truncated
    (`new_end = old_end + end_coordinate_shortfall`, equal to the fixed `haplo_coordinates.end`).

    pixi run python scripts/patch_duckdb_haplotype_filter.py \
        --in-duckdb data/work/output/hgdp_1kg.haplotypes_gnomad_merge.index.duckdb \
        --out-duckdb data/analysis/hgdp_1kg.index.fixed.duckdb
"""

import argparse
import shutil
from pathlib import Path

import duckdb

from divref.haplotype_compat import compatibility_flag
from divref.haplotype_compat import end_coordinate_shortfall
from divref.haplotype_compat import parse_variants_string

DEFAULT_IN = "data/work/output/hgdp_1kg.haplotypes_gnomad_merge.index.duckdb"
DEFAULT_OUT = "data/analysis/hgdp_1kg.index.fixed.duckdb"


def patch_duckdb(in_duckdb: Path, out_duckdb: Path) -> tuple[int, int]:
    """
    Write a fixed copy of `in_duckdb` with `haplotype_filter` added and `end` corrected.

    Args:
        in_duckdb: Source DivRef DuckDB (left untouched).
        out_duckdb: Destination path for the fixed copy (overwritten if present).

    Returns:
        `(n_flagged, n_end_fixed)`: how many haplotypes were flagged non-PASS and how many had
        their `end` corrected.
    """
    out_duckdb.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(in_duckdb, out_duckdb)

    con = duckdb.connect(str(out_duckdb))
    try:
        window_size = int(con.execute("SELECT window_size FROM window_size").fetchone()[0])
        # Every row defaults to PASS (gnomAD_variant and compatible HGDP rows stay PASS); only
        # incompatible haplotypes and end-undershoots need an explicit UPDATE.
        con.execute("ALTER TABLE sequences ADD COLUMN IF NOT EXISTS haplotype_filter VARCHAR")
        con.execute("UPDATE sequences SET haplotype_filter = 'PASS'")

        filter_updates: list[tuple[str, str]] = []
        end_updates: list[tuple[int, str]] = []
        cur = con.execute(
            'SELECT sequence_id, variants, "end" FROM sequences '
            "WHERE source = 'HGDP_haplotype' AND n_variants >= 2"
        )
        while True:
            batch = cur.fetchmany(100_000)
            if not batch:
                break
            for sequence_id, variants_str, end in batch:
                flag = compatibility_flag(variants_str)
                if flag != "PASS":
                    filter_updates.append((flag, sequence_id))
                shortfall = end_coordinate_shortfall(
                    parse_variants_string(variants_str), window_size, end
                )
                if shortfall > 0:
                    end_updates.append((end + shortfall, sequence_id))

        con.executemany(
            "UPDATE sequences SET haplotype_filter = ? WHERE sequence_id = ?", filter_updates
        )
        con.executemany('UPDATE sequences SET "end" = ? WHERE sequence_id = ?', end_updates)
    finally:
        con.close()

    return len(filter_updates), len(end_updates)


def main() -> None:
    """Patch a DivRef DuckDB into a fixed copy and report what changed."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--in-duckdb", default=DEFAULT_IN, type=Path, help="Source index.")
    parser.add_argument("--out-duckdb", default=DEFAULT_OUT, type=Path, help="Fixed-copy output.")
    args = parser.parse_args()

    n_flagged, n_end_fixed = patch_duckdb(args.in_duckdb, args.out_duckdb)
    print(f"wrote fixed copy: {args.out_duckdb}")
    print(f"haplotypes flagged non-PASS: {n_flagged}")
    print(f"end coordinates corrected:   {n_end_fixed}")


if __name__ == "__main__":
    main()
