"""Tests for scripts/evaluate_haplotype_incompatibility.py (DuckDB load + aggregation).

The pure classification functions now live in `divref.haplotype_compat` and are tested in
`divref/tests/test_haplotype_compat.py`; this file covers the DuckDB load and Summary aggregation.

Run with the uv-managed divref env, which has both pytest and duckdb:

    uv run --directory divref pytest --no-cov \\
        ../scripts/tests/test_evaluate_haplotype_incompatibility.py
"""

import sys
from pathlib import Path

import duckdb
import pytest

# The script lives in scripts/ (standalone, not a package); make it importable.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import evaluate_haplotype_incompatibility as ev  # noqa: E402


@pytest.fixture()
def fixture_con() -> duckdb.DuckDBPyConnection:
    con = duckdb.connect()
    con.execute("CREATE TABLE window_size (window_size INTEGER)")
    con.execute("INSERT INTO window_size VALUES (25)")
    con.execute(
        "CREATE TABLE sequences (contig VARCHAR, n_variants INTEGER, variants VARCHAR, "
        'popmax_empirical_AC INTEGER, start BIGINT, "end" BIGINT, source VARCHAR)'
    )
    rows = [
        # contig, n_variants, variants, AC, start, end, source
        ("chr1", 2, "chr1:200:A:T,chr1:210:C:G", 5, 174, 235, "HGDP_haplotype"),  # clean
        ("chr1", 2, "chr1:300:AT:A,chr1:301:T:A", 7, 274, 326, "HGDP_haplotype"),  # snp_in_del
        ("chr1", 2, "chr1:400:TGG:T,chr1:402:G:GTTTT", 4, 374, 427, "HGDP_haplotype"),  # indel
        ("chr1", 2, "chr1:500:AAAG:A,chr1:501:AAG:A", 9, 474, 528, "HGDP_haplotype"),  # ovl del
        ("chr1", 2, "chr1:600:AAC:A,chr1:600:AACAC:A", 3, 574, 629, "HGDP_haplotype"),  # same_pos
        ("chr1", 2, "chr1:100:AAAAA:A,chr1:102:C:T", 6, 74, 127, "HGDP_haplotype"),  # undershoot
        ("chr2", 3, "chr2:1000:A:T,chr2:1005:C:G,chr2:1010:G:A", 8, 974, 1035, "HGDP_haplotype"),
        ("chr1", 1, "chr1:50:A:T", 11, 24, 76, "gnomAD_variant"),  # ignored
        ("chr1", 1, "chr1:60:A:T", 12, 34, 86, "gnomAD_variant"),  # ignored
    ]
    con.executemany("INSERT INTO sequences VALUES (?, ?, ?, ?, ?, ?, ?)", rows)
    return con


def test_load_and_classify_counts(fixture_con: duckdb.DuckDBPyConnection) -> None:
    s = ev.load_and_classify(fixture_con, examples_per=20)
    assert s.window_size == 25
    assert s.total_haplotypes == 7  # only HGDP_haplotype rows with n_variants >= 2
    assert s.haplotypes_with_any_incompatibility == 5
    assert dict(s.reason_haplotype_counts) == {
        "snp_in_deletion": 2,
        "indel_in_deletion": 1,
        "overlapping_deletions": 1,
        "same_position": 1,
    }
    assert dict(s.reason_pair_counts) == {
        "snp_in_deletion": 2,
        "indel_in_deletion": 1,
        "overlapping_deletions": 1,
        "same_position": 1,
    }
    assert s.end_undershoot_count == 1
    assert s.end_undershoot_max_bp == 2
    assert dict(s.end_undershoot_by_reason) == {"snp_in_deletion": 1}
    assert s.start_undershoot_anomalies == 0
    # All five incompatible fixture haplotypes are length-2, so exploding recovers nothing.
    assert dict(s.incompatible_length_hist) == {2: 5}
    assert s.recoverable_haplotypes == 0
    assert s.bypass_resolutions_total == 0
    assert s.bypass_resolutions_max == 0


def test_write_summary_tsv(tmp_path: Path, fixture_con: duckdb.DuckDBPyConnection) -> None:
    s = ev.load_and_classify(fixture_con, examples_per=20)
    out = tmp_path / "summary.tsv"
    ev.write_summary_tsv(s, out)
    text = out.read_text()
    assert "metric\tvalue\n" in text
    assert "total_haplotypes\t7\n" in text
    assert "haplotype_count.snp_in_deletion\t2\n" in text
    assert "haplotypes_with_end_undershoot\t1\n" in text
    assert "end_undershoot_max_bp\t2\n" in text
