"""Tests for scripts/patch_duckdb_haplotype_filter.py.

Run with the uv-managed divref env (has pytest + duckdb):

    uv run --directory divref pytest --no-cov \\
        ../scripts/tests/test_patch_duckdb_haplotype_filter.py
"""

import sys
from pathlib import Path

import duckdb

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import patch_duckdb_haplotype_filter as patch  # noqa: E402


def _make_source_db(path: Path) -> None:
    con = duckdb.connect(str(path))
    con.execute("CREATE TABLE window_size (window_size INTEGER)")
    con.execute("INSERT INTO window_size VALUES (25)")
    con.execute(
        "CREATE TABLE sequences "
        '(sequence_id VARCHAR, variants VARCHAR, source VARCHAR, n_variants INTEGER, "end" BIGINT)'
    )
    rows = [
        ("s1", "chr1:200:A:T,chr1:210:C:G", "HGDP_haplotype", 2, 235),  # clean, end ok
        ("s2", "chr1:300:AT:A,chr1:301:T:A", "HGDP_haplotype", 2, 326),  # flagged, end ok
        ("s3", "chr1:100:AAAAA:A,chr1:102:C:T", "HGDP_haplotype", 2, 127),  # flagged + undershoot
        ("s4", "chr1:50:A:T", "gnomAD_variant", 1, 76),  # gnomAD -> PASS
    ]
    con.executemany("INSERT INTO sequences VALUES (?, ?, ?, ?, ?)", rows)
    con.close()


def test_patch_duckdb_flags_and_fixes_end(tmp_path: Path) -> None:
    src = tmp_path / "src.duckdb"
    dst = tmp_path / "fixed.duckdb"
    _make_source_db(src)

    n_flagged, n_end_fixed = patch.patch_duckdb(src, dst)
    assert n_flagged == 2  # s2, s3
    assert n_end_fixed == 1  # s3 only

    con = duckdb.connect(str(dst), read_only=True)
    try:
        flags = dict(con.execute("SELECT sequence_id, haplotype_filter FROM sequences").fetchall())
        ends = dict(con.execute('SELECT sequence_id, "end" FROM sequences').fetchall())
    finally:
        con.close()

    assert flags == {
        "s1": "PASS",
        "s2": "snp_in_deletion",
        "s3": "snp_in_deletion",
        "s4": "PASS",
    }
    assert ends["s3"] == 129  # 127 + shortfall(2)
    assert ends["s1"] == 235  # unchanged

    # The source is left untouched (no haplotype_filter column).
    src_con = duckdb.connect(str(src), read_only=True)
    try:
        cols = [d[0] for d in src_con.execute("SELECT * FROM sequences LIMIT 0").description]
    finally:
        src_con.close()
    assert "haplotype_filter" not in cols
