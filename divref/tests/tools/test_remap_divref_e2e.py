"""End-to-end tests for the remap_divref() orchestration (CSV in, DuckDB lookup, CSV out)."""

import json
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import duckdb
import pandas as pd
import pytest

from divref.tools.remap_divref import _find_index_in_cwd
from divref.tools.remap_divref import _get_index_connection
from divref.tools.remap_divref import remap_divref

_POPS: tuple[str, ...] = ("afr", "amr")
_VARIANTS: str = "chr1:500:A:T,chr1:505:C:G,chr1:510:T:A"


def _seq_row(
    sequence_id: str,
    variants: str = _VARIANTS,
    *,
    haplotype_filter: str = "PASS",
    pops: Iterable[str] = _POPS,
) -> dict[str, Any]:
    """Build one `sequences` row dict with the columns remap_divref reads."""
    row: dict[str, Any] = {
        "sequence_id": sequence_id,
        "sequence": "ACGT",
        "sequence_length": 4,
        "n_variants": len(variants.split(",")),
        "popmax_fraction_phased": 1.0,
        "popmax_empirical_AF": 0.25,
        "popmax_empirical_AC": 1000,
        "popmax_estimated_gnomad_AF": 0.15,
        "max_pop": "amr",
        "variants": variants,
        "source": "HGDP_haplotype",
        "haplotype_filter": haplotype_filter,
    }
    for pop in pops:
        row[f"gnomAD_AF_{pop}"] = "0.1,0.2,0.3"
        row[f"estimated_gnomAD_haplotype_AF_{pop}"] = 0.05
    return row


def _build_index(
    db_path: Path,
    *,
    seq_rows: list[dict[str, Any]],
    version: str = "9.9",
    window_size: int = 10,
    pops: Iterable[str] = _POPS,
    skip: Iterable[str] = (),
) -> None:
    """
    Build a minimal DivRef DuckDB index (metadata tables + a `sequences` table).

    Args:
        db_path: Path to the DuckDB file to create.
        seq_rows: `sequences` rows (e.g. from `_seq_row`).
        version: Value for the VERSION table.
        window_size: Value for the window_size table.
        pops: Population labels for the joint_pops_legend table.
        skip: Metadata table names to NOT create (to exercise error paths).
    """
    skip = set(skip)
    with duckdb.connect(str(db_path)) as conn:
        if "VERSION" not in skip:
            conn.execute("CREATE TABLE VERSION AS SELECT ? AS version", [version])
        if "window_size" not in skip:
            conn.execute("CREATE TABLE window_size AS SELECT ? AS window_size", [window_size])
        if "joint_pops_legend" not in skip:
            conn.execute(
                "CREATE TABLE joint_pops_legend AS SELECT ? AS pops_legend",
                [json.dumps(list(pops))],
            )
        seq_df = pd.DataFrame(seq_rows)
        conn.register("seq_df", seq_df)
        conn.execute("CREATE TABLE sequences AS SELECT * FROM seq_df")


def _write_calitas_tsv(path: Path, rows: list[dict[str, Any]]) -> Path:
    """Write a CALITAS-style input TSV (the `chromosome` column holds the DivRef sequence id)."""
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    return path


def _calitas_row(
    sequence_id: str,
    *,
    coordinate_start: int,
    coordinate_end: int,
    strand: str = "+",
    padded_target: str = "ACGTA",
    unpadded_target_sequence: str = "ACGTA",
) -> dict[str, Any]:
    """Build one CALITAS input row."""
    return {
        "chromosome": sequence_id,
        "coordinate_start": coordinate_start,
        "coordinate_end": coordinate_end,
        "strand": strand,
        "padded_target": padded_target,
        "unpadded_target_sequence": unpadded_target_sequence,
    }


def test_remap_plus_strand_appends_reference_coordinates(tmp_path: Path) -> None:
    """A '+' strand hit is remapped to reference coordinates with the full metadata columns."""
    db_path = tmp_path / "index.duckdb"
    _build_index(db_path, seq_rows=[_seq_row("hapA", haplotype_filter="snp_in_deletion")])
    input_tsv = _write_calitas_tsv(
        tmp_path / "calitas.tsv",
        [_calitas_row("hapA", coordinate_start=12, coordinate_end=17)],
    )
    output_tsv = tmp_path / "out.tsv"

    remap_divref(input_path=input_tsv, output_path=output_tsv, index_path=db_path)

    out = pd.read_csv(output_tsv, sep="\t")
    row = out.iloc[0]
    # Original DivRef-space coordinates preserved under divref_* columns.
    assert row["divref_sequence_id"] == "hapA"
    assert row["divref_start"] == 12
    assert row["divref_end"] == 17
    # Remapped reference coordinates replace the original chromosome/start/end columns.
    assert row["chromosome"] == "chr1"
    assert row["coordinate_start"] == 502
    assert row["coordinate_end"] == 507
    assert row["genome_build"] == "DivRef-v9.9"
    assert row["all_variants"] == _VARIANTS
    assert row["variants_involved"] == "chr1:505:C:G"
    assert row["n_variants_involved"] == 1
    assert row["max_pop"] == "amr"
    assert row["variant_source"] == "HGDP_haplotype"
    # haplotype_filter is carried through verbatim.
    assert row["haplotype_filter"] == "snp_in_deletion"
    # One gnomAD_AF_/estimated column per joint-legend population.
    assert row["gnomAD_AF_afr"] == "0.1,0.2,0.3"
    assert row["estimated_gnomAD_haplotype_AF_amr"] == 0.05


def test_remap_minus_strand_adjusts_start_by_padding(tmp_path: Path) -> None:
    """On the '-' strand the padding-length adjustment is subtracted from the start coordinate."""
    db_path = tmp_path / "index.duckdb"
    _build_index(db_path, seq_rows=[_seq_row("hapA")])
    # padded_target has 2 more (gap-stripped) bases than the unpadded target, so padded_len_adj=2;
    # on the '-' strand start 14 -> 12, mapping the same [12, 19) window as a hand-computed case.
    input_tsv = _write_calitas_tsv(
        tmp_path / "calitas.tsv",
        [
            _calitas_row(
                "hapA",
                coordinate_start=14,
                coordinate_end=19,
                strand="-",
                padded_target="ACGTAA",
                unpadded_target_sequence="ACGT",
            )
        ],
    )
    output_tsv = tmp_path / "out.tsv"

    remap_divref(input_path=input_tsv, output_path=output_tsv, index_path=db_path)

    row = pd.read_csv(output_tsv, sep="\t").iloc[0]
    assert row["divref_start"] == 14  # original, unadjusted
    assert row["coordinate_start"] == 502
    assert row["coordinate_end"] == 509


def test_remap_missing_required_field_raises(tmp_path: Path) -> None:
    """An input TSV missing a required column raises ValueError before touching the index."""
    db_path = tmp_path / "index.duckdb"
    _build_index(db_path, seq_rows=[_seq_row("hapA")])
    # Drop the 'strand' column.
    bad = _calitas_row("hapA", coordinate_start=12, coordinate_end=17)
    del bad["strand"]
    input_tsv = _write_calitas_tsv(tmp_path / "calitas.tsv", [bad])

    with pytest.raises(ValueError, match="Required fields"):
        remap_divref(input_path=input_tsv, output_path=tmp_path / "out.tsv", index_path=db_path)


def test_remap_unknown_haplotype_raises(tmp_path: Path) -> None:
    """A CALITAS row whose sequence id is absent from the index raises RuntimeError."""
    db_path = tmp_path / "index.duckdb"
    _build_index(db_path, seq_rows=[_seq_row("hapA")])
    input_tsv = _write_calitas_tsv(
        tmp_path / "calitas.tsv",
        [_calitas_row("hapMISSING", coordinate_start=12, coordinate_end=17)],
    )

    with pytest.raises(RuntimeError, match="Unable to find haplotype"):
        remap_divref(input_path=input_tsv, output_path=tmp_path / "out.tsv", index_path=db_path)


@pytest.mark.parametrize("missing_table", ["VERSION", "window_size", "joint_pops_legend"])
def test_remap_missing_metadata_table_raises(tmp_path: Path, missing_table: str) -> None:
    """A missing metadata table raises a RuntimeError naming that table."""
    db_path = tmp_path / "index.duckdb"
    _build_index(db_path, seq_rows=[_seq_row("hapA")], skip=[missing_table])
    input_tsv = _write_calitas_tsv(
        tmp_path / "calitas.tsv",
        [_calitas_row("hapA", coordinate_start=12, coordinate_end=17)],
    )

    with pytest.raises(RuntimeError, match=missing_table):
        remap_divref(input_path=input_tsv, output_path=tmp_path / "out.tsv", index_path=db_path)


def test_remap_closes_the_index_connection(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """remap_divref closes the DuckDB connection it opens (no leaked handle)."""
    db_path = tmp_path / "index.duckdb"
    _build_index(db_path, seq_rows=[_seq_row("hapA")])
    input_tsv = _write_calitas_tsv(
        tmp_path / "calitas.tsv",
        [_calitas_row("hapA", coordinate_start=12, coordinate_end=17)],
    )

    opened: list[duckdb.DuckDBPyConnection] = []
    real_connect = duckdb.connect

    def tracking_connect(*args: Any, **kwargs: Any) -> duckdb.DuckDBPyConnection:
        conn = real_connect(*args, **kwargs)
        opened.append(conn)
        return conn

    monkeypatch.setattr("divref.tools.remap_divref.duckdb.connect", tracking_connect)
    remap_divref(input_path=input_tsv, output_path=tmp_path / "out.tsv", index_path=db_path)

    assert opened, "expected remap_divref to open a DuckDB connection"
    # Executing on a closed connection raises; this fails if the connection was left open.
    with pytest.raises(duckdb.ConnectionException):
        opened[-1].execute("SELECT 1")


def test_get_index_connection_discovers_index_in_cwd(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """With no index_path, the connection is opened from a .duckdb discovered under cwd."""
    db_path = tmp_path / "nested" / "index.duckdb"
    db_path.parent.mkdir()
    _build_index(db_path, seq_rows=[_seq_row("hapA")])
    monkeypatch.chdir(tmp_path)

    conn = _get_index_connection(None)
    try:
        assert conn.execute("SELECT * FROM VERSION").fetchone() is not None
    finally:
        conn.close()


def test_find_index_in_cwd_returns_none_when_absent(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """_find_index_in_cwd returns None (and _get_index_connection raises) when no index exists."""
    monkeypatch.chdir(tmp_path)
    assert _find_index_in_cwd() is None
    with pytest.raises(RuntimeError, match="Unable to find a DuckDB index"):
        _get_index_connection(None)
