"""Tests for scripts/evaluate_haplotype_incompatibility.py.

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


def _v(token: str) -> ev.Variant:
    contig, pos, ref, alt = token.split(":")
    return (contig, int(pos), ref, alt)


# --- Task 1: parser + distance ----------------------------------------------------------------


def test_parse_variants_string_sorts_and_roundtrips() -> None:
    parsed = ev.parse_variants_string("chr1:1744031:AAC:A,chr1:1744033:C:A")
    assert parsed == [("chr1", 1744031, "AAC", "A"), ("chr1", 1744033, "C", "A")]


def test_parse_variants_string_sorts_same_position_by_ref_length() -> None:
    parsed = ev.parse_variants_string("chr1:600:AACAC:A,chr1:600:AAC:A")
    assert [len(v[2]) for v in parsed] == [3, 5]


@pytest.mark.parametrize(
    ("a", "b", "expected"),
    [
        ("chr1:1744031:AAC:A", "chr1:1744033:C:A", -1),
        ("chr1:1:AA:T", "chr1:3:A:T", 0),
        ("chr1:1:A:T", "chr1:3:A:T", 1),
    ],
)
def test_variant_distance(a: str, b: str, expected: int) -> None:
    assert ev.variant_distance(_v(a), _v(b)) == expected


# --- Task 2: classify_pair --------------------------------------------------------------------


@pytest.mark.parametrize(
    ("a", "b", "expected"),
    [
        ("chr1:300:AT:A", "chr1:301:T:A", "snp_in_deletion"),
        ("chr1:400:TGG:T", "chr1:402:G:GTTTT", "indel_in_deletion"),
        ("chr1:500:AAAG:A", "chr1:501:AAG:A", "overlapping_deletions"),
        ("chr1:600:AAC:A", "chr1:600:AACAC:A", "same_position"),
        ("chr1:700:AT:ATGG", "chr1:701:C:G", "insertion_anchor_conflict"),
        ("chr1:200:A:T", "chr1:210:C:G", None),  # distance >= 0: compatible
        ("chr1:1:AA:T", "chr1:3:A:T", None),  # distance == 0: deletion closes gap, keep
    ],
)
def test_classify_pair(a: str, b: str, expected: str | None) -> None:
    assert ev.classify_pair(_v(a), _v(b)) == expected


# --- Task 3: classify_haplotype ---------------------------------------------------------------


def test_classify_haplotype_clean() -> None:
    variants = ev.parse_variants_string("chr2:1000:A:T,chr2:1005:C:G,chr2:1010:G:A")
    assert ev.classify_haplotype(variants) == []


def test_classify_haplotype_single_reason() -> None:
    variants = ev.parse_variants_string("chr1:300:AT:A,chr1:301:T:A")
    assert ev.classify_haplotype(variants) == ["snp_in_deletion"]


def test_classify_haplotype_catches_nonadjacent_overlap() -> None:
    # A long deletion at 100 (spans 100-119) swallows a variant at 110 that is NOT adjacent
    # to it in the array (a SNP at 105 sits between). The overlap is still flagged, because the
    # deletion overlaps the 105 SNP too (adjacent-pair sufficiency).
    variants = ev.parse_variants_string(
        "chr1:100:AAAAAAAAAAAAAAAAAAAA:A,chr1:105:C:T,chr1:110:G:A"
    )
    assert ev.classify_haplotype(variants)  # non-empty -> flagged


def test_classify_haplotype_multiple_reasons() -> None:
    # del+snp overlap, then a clean gap, then a same-position pair.
    variants = ev.parse_variants_string(
        "chr1:300:AT:A,chr1:301:T:A,chr1:600:AAC:A,chr1:600:AACAC:A"
    )
    assert set(ev.classify_haplotype(variants)) == {"snp_in_deletion", "same_position"}


# --- explode-vs-drop sizing: overlap + bypass resolutions -------------------------------------


@pytest.mark.parametrize(
    ("a", "b", "expected"),
    [
        ("chr1:300:AT:A", "chr1:301:T:A", True),       # snp inside deletion
        ("chr1:600:AAC:A", "chr1:600:AACAC:A", True),  # same position
        ("chr1:200:A:T", "chr1:210:C:G", False),       # clear gap
        ("chr1:1:AA:T", "chr1:3:A:T", False),          # distance 0: touching, not overlapping
    ],
)
def test_variants_overlap(a: str, b: str, expected: bool) -> None:
    assert ev.variants_overlap(_v(a), _v(b)) is expected
    assert ev.variants_overlap(_v(b), _v(a)) is expected  # order-independent


def test_count_bypass_resolutions_length2_conflict_recovers_nothing() -> None:
    # Exploding a length-2 conflict yields only singletons -> no >=2-variant resolution.
    variants = ev.parse_variants_string("chr1:300:AT:A,chr1:301:T:A")
    assert ev.count_bypass_resolutions(variants) == 0


def test_count_bypass_resolutions_clean_is_one() -> None:
    variants = ev.parse_variants_string("chr1:100:A:T,chr1:200:C:G")
    assert ev.count_bypass_resolutions(variants) == 1


def test_count_bypass_resolutions_single_internal_conflict_is_two() -> None:
    # clean A, conflicting B/C pair, clean D -> resolutions {A,B,D} and {A,C,D}.
    variants = ev.parse_variants_string("chr1:100:A:T,chr1:300:AT:A,chr1:301:T:A,chr1:400:G:C")
    assert ev.count_bypass_resolutions(variants) == 2


def test_count_bypass_resolutions_three_mutually_exclusive_alleles() -> None:
    # Three deletions all overlapping each other (a repeat ladder): at most one kept, plus the
    # clean flank -> three single-allele resolutions.
    variants = ev.parse_variants_string(
        "chr1:100:A:T,chr1:500:AAAAAAAA:A,chr1:502:AAAA:A,chr1:504:AA:A"
    )
    assert ev.count_bypass_resolutions(variants) == 3


# --- Task 3b: coordinate-adequacy functions ---------------------------------------------------


def test_end_coordinate_shortfall_undershoots_on_early_deletion() -> None:
    # del at 100 spans 100-104 (rightmost ref end 0-based excl = 100-1+5 = 104); last-by-position
    # variant is the SNP at 102, so the stored end = 102-1+1+25 = 127. Required = 104+25 = 129.
    variants = ev.parse_variants_string("chr1:100:AAAAA:A,chr1:102:C:T")
    assert ev.end_coordinate_shortfall(variants, window_size=25, stored_end=127) == 2


def test_end_coordinate_shortfall_zero_when_last_variant_reaches_furthest() -> None:
    variants = ev.parse_variants_string("chr10:112867606:AT:A,chr10:112867607:T:A")
    # rightmost ref end 0-based excl = 112867607-1+1 = 112867607; stored end = that + 25.
    stored_end = 112867607 + 25
    assert ev.end_coordinate_shortfall(variants, window_size=25, stored_end=stored_end) == 0


def test_start_coordinate_shortfall_zero_for_correct_left_edge() -> None:
    variants = ev.parse_variants_string("chr1:100:AAAAA:A,chr1:102:C:T")
    # leftmost 0-based base = 100-1 = 99; required start = 99 - 25 = 74.
    assert ev.start_coordinate_shortfall(variants, window_size=25, stored_start=74) == 0


# --- Task 4: aggregation over an in-memory DuckDB fixture -------------------------------------


@pytest.fixture()
def fixture_con() -> duckdb.DuckDBPyConnection:
    con = duckdb.connect()
    con.execute("CREATE TABLE window_size (window_size INTEGER)")
    con.execute("INSERT INTO window_size VALUES (25)")
    con.execute(
        'CREATE TABLE sequences (contig VARCHAR, n_variants INTEGER, variants VARCHAR, '
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
    con.executemany(
        'INSERT INTO sequences VALUES (?, ?, ?, ?, ?, ?, ?)',
        rows,
    )
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
