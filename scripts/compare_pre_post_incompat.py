"""Compare a pre-fix DivRef build against the re-generated (post-fix) build.

Validates the incompatibility-flag + end-coordinate fix end-to-end by diffing the old index/FASTAs
(saved under data/work/output/pre_incompat) against the freshly rebuilt ones in data/work/output.

Expected differences (and nothing else):
  - the new index has a `haplotype_filter` column (PASS, else the incompatibility reason),
  - `end` and `sequence` differ only for the end-undershoot haplotypes (the sequence is longer by
    the deletion's truncated span), and those rows are all flagged non-PASS,
  - per-contig FASTAs differ only for contigs that contain such haplotypes.

    pixi run python scripts/compare_pre_post_incompat.py
"""

import argparse
import hashlib
from pathlib import Path

import duckdb

DEFAULT_OLD_DUCKDB = "data/work/output/pre_incompat/hgdp_1kg.haplotypes_gnomad_merge.index.duckdb"
DEFAULT_NEW_DUCKDB = "data/work/output/hgdp_1kg.haplotypes_gnomad_merge.index.duckdb"
DEFAULT_OLD_FASTA_DIR = "data/work/output/pre_incompat"
DEFAULT_NEW_FASTA_DIR = "data/work/output"
DEFAULT_FASTA_BASE = "hgdp_1kg.haplotypes_gnomad_merge"
DEFAULT_CONTIGS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


def compare_duckdbs(old_path: Path, new_path: Path) -> dict[str, object]:
    """Diff the `sequences` tables of the old and new indexes by `sequence_id`."""
    con = duckdb.connect()
    con.execute(f"ATTACH '{old_path.resolve()}' AS old (READ_ONLY)")
    con.execute(f"ATTACH '{new_path.resolve()}' AS new (READ_ONLY)")

    def scalar(sql: str) -> int:
        return int(con.execute(sql).fetchone()[0])

    metrics: dict[str, object] = {
        "old_row_count": scalar("SELECT COUNT(*) FROM old.sequences"),
        "new_row_count": scalar("SELECT COUNT(*) FROM new.sequences"),
        "sequence_ids_only_in_old": scalar(
            "SELECT COUNT(*) FROM (SELECT sequence_id FROM old.sequences "
            "EXCEPT SELECT sequence_id FROM new.sequences)"
        ),
        "sequence_ids_only_in_new": scalar(
            "SELECT COUNT(*) FROM (SELECT sequence_id FROM new.sequences "
            "EXCEPT SELECT sequence_id FROM old.sequences)"
        ),
    }
    # Rows that changed value (joined on the shared sequence_id).
    join = (
        "FROM new.sequences n JOIN old.sequences o ON n.sequence_id = o.sequence_id "  # noqa: S608
    )
    metrics["rows_end_changed"] = scalar(f'SELECT COUNT(*) {join} WHERE n."end" != o."end"')
    metrics["rows_sequence_changed"] = scalar(
        f"SELECT COUNT(*) {join} WHERE n.sequence != o.sequence"
    )
    metrics["rows_changed"] = scalar(
        f'SELECT COUNT(*) {join} WHERE n."end" != o."end" OR n.sequence != o.sequence'
    )
    # Every changed row must be flagged non-PASS (the fix only touches incompatible haplotypes).
    metrics["changed_rows_that_are_pass"] = scalar(
        f'SELECT COUNT(*) {join} WHERE (n."end" != o."end" OR n.sequence != o.sequence) '
        "AND n.haplotype_filter = 'PASS'"
    )
    metrics["new_flagged_non_pass"] = scalar(
        "SELECT COUNT(*) FROM new.sequences WHERE haplotype_filter != 'PASS'"
    )
    filter_rows = con.execute(
        "SELECT haplotype_filter, COUNT(*) FROM new.sequences GROUP BY 1 ORDER BY 2 DESC"
    ).fetchall()
    metrics["new_filter_distribution"] = {row[0]: int(row[1]) for row in filter_rows}
    con.close()
    return metrics


def _fasta_md5(path: Path) -> str | None:
    """MD5 of a FASTA file, or None if it does not exist."""
    if not path.exists():
        return None
    h = hashlib.md5()  # noqa: S324 - non-cryptographic file-equality check
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def compare_fastas(
    old_dir: Path, new_dir: Path, contigs: list[str], base: str
) -> dict[str, list[str]]:
    """Per-contig FASTA equality by MD5; report which contigs changed / are missing."""
    identical: list[str] = []
    changed: list[str] = []
    missing: list[str] = []
    for contig in contigs:
        old_md5 = _fasta_md5(old_dir / f"{base}.{contig}.fasta")
        new_md5 = _fasta_md5(new_dir / f"{base}.{contig}.fasta")
        if old_md5 is None or new_md5 is None:
            missing.append(contig)
        elif old_md5 == new_md5:
            identical.append(contig)
        else:
            changed.append(contig)
    return {"identical": identical, "changed": changed, "missing": missing}


def main() -> None:
    """Run the old-vs-new comparison and print a pass/fail-oriented report."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--old-duckdb", default=DEFAULT_OLD_DUCKDB, type=Path)
    parser.add_argument("--new-duckdb", default=DEFAULT_NEW_DUCKDB, type=Path)
    parser.add_argument("--old-fasta-dir", default=DEFAULT_OLD_FASTA_DIR, type=Path)
    parser.add_argument("--new-fasta-dir", default=DEFAULT_NEW_FASTA_DIR, type=Path)
    parser.add_argument("--fasta-base", default=DEFAULT_FASTA_BASE)
    parser.add_argument("--contigs", nargs="+", default=DEFAULT_CONTIGS)
    args = parser.parse_args()

    db = compare_duckdbs(args.old_duckdb, args.new_duckdb)
    print("=== DuckDB comparison (old vs new) ===")
    for key in (
        "old_row_count",
        "new_row_count",
        "sequence_ids_only_in_old",
        "sequence_ids_only_in_new",
        "rows_end_changed",
        "rows_sequence_changed",
        "rows_changed",
        "changed_rows_that_are_pass",
        "new_flagged_non_pass",
    ):
        print(f"  {key:28s} {db[key]}")
    print(f"  new_filter_distribution: {db['new_filter_distribution']}")

    ok = (
        db["old_row_count"] == db["new_row_count"]
        and db["sequence_ids_only_in_old"] == 0
        and db["sequence_ids_only_in_new"] == 0
        and db["changed_rows_that_are_pass"] == 0
    )
    print(f"  => row sets identical and all changed rows flagged: {'PASS' if ok else 'FAIL'}")

    fasta = compare_fastas(args.old_fasta_dir, args.new_fasta_dir, args.contigs, args.fasta_base)
    print("\n=== FASTA comparison (per contig) ===")
    print(f"  identical contigs: {len(fasta['identical'])}")
    print(f"  changed contigs:   {fasta['changed']}")
    if fasta["missing"]:
        print(f"  !! missing FASTAs: {fasta['missing']}")


if __name__ == "__main__":
    main()
