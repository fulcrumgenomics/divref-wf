"""Inspect the 6 chr22 haplotypes that the original algorithm emits but the new
algorithm does not, against the new algorithm's intermediate Hail tables.

Sources:
  variants.ht  — per-variant call_stats (frequencies_by_pop) and gnomAD freqs
  hap_ac.ht    — per-unique-haplotype containment AC, pre-filter, pre-dedup
  parents.ht   — per-(sample, strand) parent block compositions

For each case, computes the new algorithm's empirical_AF, fraction_phased,
and estimated_gnomad_AF by hand (same formula as `_compute_metrics`), so the
analysis works even if the final .ht hasn't been written yet.
"""

import argparse
from collections import Counter
from pathlib import Path

import duckdb
import hail as hl

DEFAULT_BASE = "data/work/haplotypes/hgdp_1kg.haplotypes.chr22"
GNOMAD_AFS_HT_FALLBACK = "data/work/inputs/hgdp_1kg.sites.chr22.ht"
DEFAULT_OLD_DUCKDB = (
    "data/analysis/input/DivRef-v1.1.haplotypes_gnomad_merge.index.duckdb"
)

CASES = [
    ("Adjacent #1", [(19951135, ("G", "A")), (19951136, ("A", "G"))]),
    (
        "VNTR",
        [
            (40457473, ("AGAAAGAAAGAAAGAAG", "A")),
            (40457477, ("AGAAAGAAAGAAG", "A")),
        ],
    ),
    ("Short-gap", [(32402562, ("G", "A")), (32402565, ("TCAG", "T"))]),
    ("Intermediate", [(22627350, ("T", "C")), (22627366, ("G", "A"))]),
    (
        "Triple+skip",
        [
            (24600987, ("A", "G")),
            (24601005, ("T", "C")),
            (24601022, ("A", "C")),
        ],
    ),
    ("Non-contig", [(24626868, ("C", "T")), (24626892, ("A", "G"))]),
]


def main() -> None:
    import os

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--base",
        default=DEFAULT_BASE,
        help=(
            "Base path of the compute_haplotypes output bundle (without the "
            ".variants.ht / .hap_ac.ht / .parents.ht / .ht suffix). "
            f"Default: {DEFAULT_BASE}"
        ),
    )
    parser.add_argument(
        "--old-duckdb",
        default=DEFAULT_OLD_DUCKDB,
        help=(
            "Path to the original-pipeline DuckDB index (defaults to the file "
            "produced by `workflows/compare_divref_gnomad.smk`'s "
            "`download_divref_index` rule)."
        ),
    )
    args = parser.parse_args()
    base = args.base
    old_duckdb = args.old_duckdb

    os.environ["PYSPARK_SUBMIT_ARGS"] = "--driver-memory 8g --executor-memory 8g pyspark-shell"
    hl.init(quiet=True)

    ht_variants = hl.read_table(f"{base}.variants.ht")
    ht_hap_ac = hl.read_table(f"{base}.hap_ac.ht")
    ht_parents = hl.read_table(f"{base}.parents.ht")
    if hl.hadoop_exists(f"{base}.ht"):
        ht_final = hl.read_table(f"{base}.ht")
    else:
        ht_final = None

    if hl.hadoop_exists(GNOMAD_AFS_HT_FALLBACK):
        pops_legend = hl.read_table(GNOMAD_AFS_HT_FALLBACK).globals.pops.collect()[0]
    else:
        pops_legend = ht_variants.globals.pops.collect()[0]
    print(f"pops legend: {pops_legend}")
    n_pops = len(pops_legend)

    # Filter variants.ht Hail-side to a window around each case (±100bp) so that
    # parent-block neighbors are also resolvable when we print them later.
    case_positions = sorted({pos for _label, expected in CASES for pos, _ in expected})
    case_locus_alleles = {
        (pos, alleles): True for _label, expected in CASES for pos, alleles in expected
    }
    window = 100
    intervals = [
        hl.parse_locus_interval(
            f"chr22:{max(1, pos - window)}-{pos + window}", reference_genome="GRCh38"
        )
        for pos in case_positions
    ]
    print(
        f"filtering variants.ht to ±{window}bp windows around {len(case_positions)} case positions ..."
    )
    ht_variants_filtered = hl.filter_intervals(ht_variants, intervals)
    variant_rows = ht_variants_filtered.collect()
    pa_to_v = {}
    rid_to_pa = {}
    for r in variant_rows:
        pa = (r.locus.position, tuple(r.alleles))
        rid_to_pa[r.row_idx] = pa
        if pa in case_locus_alleles:
            pa_to_v[pa] = r
    print(
        f"  variants in window: {len(rid_to_pa)};  "
        f"resolved {len(pa_to_v)} of {len(case_locus_alleles)} case variants"
    )

    # Resolve each case's variants to row_idx tuples.
    case_rids_per_case = []
    for _label, expected in CASES:
        rids = [pa_to_v[pa].row_idx if pa in pa_to_v else None for pa in expected]
        case_rids_per_case.append(rids)

    # Build the set of case haplotype keys (each is a tuple of int64 row_idxs).
    case_keys = []
    for _label, expected in CASES:
        rids = [pa_to_v[pa].row_idx if pa in pa_to_v else None for pa in expected]
        if None not in rids:
            case_keys.append(tuple(rids))
    case_keys_lit = hl.literal(
        {",".join(str(r) for r in k) for k in case_keys}, dtype=hl.tset(hl.tstr)
    )
    print(f"filtering hap_ac.ht to {len(case_keys)} case haplotype keys ...")
    ht_hap_ac_filtered = ht_hap_ac.annotate(
        _hap_str=hl.delimit(ht_hap_ac.haplotype.map(hl.str), ",")
    )
    ht_hap_ac_filtered = ht_hap_ac_filtered.filter(case_keys_lit.contains(ht_hap_ac_filtered._hap_str))
    hap_ac_by_hap = {tuple(r.haplotype): r for r in ht_hap_ac_filtered.collect()}
    print(f"  hap_ac rows matched: {len(hap_ac_by_hap)}")

    # Same filter on the final HT (if it exists) so we can cross-check the manual est_af
    # computation against the value Hail actually wrote.
    final_by_hap: dict = {}
    if ht_final is not None:
        print(f"filtering final HT to {len(case_keys)} case haplotype keys ...")
        ht_final_filtered = ht_final.annotate(
            _hap_str=hl.delimit(ht_final.haplotype.map(hl.str), ",")
        )
        ht_final_filtered = ht_final_filtered.filter(
            case_keys_lit.contains(ht_final_filtered._hap_str)
        )
        final_by_hap = {tuple(r.haplotype): r for r in ht_final_filtered.collect()}
        print(f"  final HT rows matched: {len(final_by_hap)}")
    else:
        print("  final HT not present yet; skipping cross-check")

    all_case_rids = set()
    for rids in case_rids_per_case:
        for r in rids:
            if r is not None:
                all_case_rids.add(r)
    print(f"filtering parents.ht to blocks touching any of {len(all_case_rids)} case row_idxs ...")
    rid_set_lit = hl.literal(all_case_rids, dtype=hl.tset(hl.tint64))
    parents_touched = ht_parents.filter(
        ht_parents.parent_block.any(lambda c: rid_set_lit.contains(c.row_idx))
    ).collect()
    print(f"parent blocks touching any case variant: {len(parents_touched)}\n")

    def fmt_v(pa: tuple) -> str:
        pos, (ref, alt) = pa
        return f"{pos}:{ref}>{alt}"

    def fmt_block(rids: tuple) -> str:
        return ", ".join(fmt_v(rid_to_pa[r]) if r in rid_to_pa else f"rid={r}" for r in rids)

    # Look up each case's row in the OLD DuckDB so we can compare max_pop labels.
    print(f"querying OLD DuckDB at {old_duckdb} ...")
    con = duckdb.connect(str(Path(old_duckdb).resolve()), read_only=True)
    old_lookup: dict = {}
    for label, expected in CASES:
        variants_str = ",".join(
            f"chr22:{pos}:{ref}:{alt}" for pos, (ref, alt) in expected
        )
        rows = con.execute(
            "SELECT max_pop, popmax_empirical_AC, popmax_empirical_AF, "
            "fraction_phased, estimated_gnomad_AF "
            "FROM sequences "
            "WHERE source = 'HGDP_haplotype' AND variants = ?",
            [variants_str],
        ).fetchall()
        old_lookup[label] = rows[0] if rows else None
    con.close()

    for (label, expected), case_rids in zip(CASES, case_rids_per_case, strict=True):
        print(f"==== {label} ====")
        print(f"  variants: {', '.join(fmt_v(e) for e in expected)}")
        if None in case_rids:
            missing = [expected[i] for i, r in enumerate(case_rids) if r is None]
            print(f"  not all variants are in variants.ht: missing {missing}")
            print()
            continue

        old_row = old_lookup.get(label)
        if old_row:
            mp_old, ac_old, af_old, fp_old, est_old = old_row
            print(
                f"  OLD DUCKDB: max_pop={mp_old} popmax_AC={ac_old} "
                f"popmax_AF={af_old:.5f} fp={fp_old:.4f} est_af={est_old:.5f}"
            )
        else:
            print(f"  OLD DUCKDB: no row found for these variants")

        hap_key = tuple(case_rids)
        final_row = final_by_hap.get(hap_key)
        if final_row is not None:
            mp_idx = final_row.max_pop
            print(
                f"  NEW HT:     max_pop={pops_legend[mp_idx]} "
                f"popmax_AC={final_row.max_empirical_AC} "
                f"popmax_AF={final_row.max_empirical_AF:.5f} "
                f"min_var_freq={final_row.min_variant_frequency:.5f} "
                f"fp={final_row.fraction_phased:.4f} "
                f"est_af={final_row.estimated_gnomad_AF:.5f}"
            )
        elif ht_final is not None:
            print(
                "  NEW HT:     not present (est_af below the threshold this run was filtered to)"
            )

        hap_key = tuple(case_rids)
        hap_row = hap_ac_by_hap.get(hap_key)
        if hap_row is None:
            print(
                "  HAP_AC: case haplotype NOT enumerated as a contiguous sub-fragment of "
                "any parent block (algorithmically absent — not just filtered)"
            )
        else:
            per_pop_AC = list(hap_row.per_pop_AC)
            vrows = [pa_to_v[pa] for pa in expected]
            rows_out = []
            for p in range(n_pops):
                ans = [v.frequencies_by_pop.get(p) for v in vrows]
                if any(an is None for an in ans):
                    rows_out.append(None)
                    continue
                min_an = min(a.AN for a in ans)
                if min_an == 0:
                    rows_out.append(None)
                    continue
                emp_af = per_pop_AC[p] / min_an
                min_local = min(a.AF[1] for a in ans)
                if min_local == 0:
                    rows_out.append((emp_af, None, None))
                    continue
                fp = emp_af / min_local
                gnomad_min = min(v.freq[p].AF for v in vrows)
                est_af = gnomad_min * fp
                rows_out.append((emp_af, fp, est_af))

            defined = [(i, r) for i, r in enumerate(rows_out) if r is not None]
            print(
                "  raw per_pop_AC vector: "
                + ", ".join(f"{pops_legend[i]}={ac}" for i, ac in enumerate(per_pop_AC) if ac > 0)
            )
            if defined:
                max_idx = max(defined, key=lambda kv: kv[1][0])[0]
                emp_af, fp, est_af = rows_out[max_idx]
                est_s = f"{est_af:.5f}" if est_af is not None else "n/a"
                fp_s = f"{fp:.4f}" if fp is not None else "n/a"
                print(
                    f"  inferred max_pop={pops_legend[max_idx]} (AC={per_pop_AC[max_idx]}, "
                    f"emp_AF={emp_af:.5f}, fp={fp_s}, est_af={est_s})"
                )
                if est_af is not None:
                    verdict = "PASS" if est_af >= 0.005 else "FAIL"
                    print(f"  est_af vs 0.005 threshold: {verdict}")

                print("  per-pop breakdown:")
                for i, row in enumerate(rows_out):
                    if row is None or per_pop_AC[i] == 0:
                        continue
                    emp_af, fp, est_af = row
                    fp_s = f"{fp:.4f}" if fp is not None else "n/a"
                    est_s = f"{est_af:.5f}" if est_af is not None else "n/a"
                    print(
                        f"    {pops_legend[i]}: AC={per_pop_AC[i]} emp_AF={emp_af:.5f} "
                        f"fp={fp_s} est_af={est_s}"
                    )
            else:
                print("  no pop has min_AN > 0 for these variants")

        case_rid_set = set(case_rids)
        touching = [
            p
            for p in parents_touched
            if any(c.row_idx in case_rid_set for c in p.parent_block)
        ]
        pattern_counts: Counter = Counter()
        for p in touching:
            rids = tuple(c.row_idx for c in p.parent_block)
            pattern_counts[(rids, p.pop_int)] += 1
        by_tuple: dict = {}
        for (rids, popi), c in pattern_counts.items():
            by_tuple.setdefault(rids, Counter())[popi] += c

        if by_tuple:
            print("  PARENT BLOCK PATTERNS touching at least one case variant:")
            n = len(case_rids)
            for rids in sorted(by_tuple, key=lambda t: -sum(by_tuple[t].values())):
                per_pop = [by_tuple[rids].get(i, 0) for i in range(n_pops)]
                labels = ", ".join(
                    f"{pops_legend[i]}={c}" for i, c in enumerate(per_pop) if c > 0
                )
                tot = sum(per_pop)
                is_pure = case_rid_set == set(rids)
                is_contig_sub = False
                for i in range(len(rids) - n + 1):
                    if set(rids[i : i + n]) == case_rid_set:
                        is_contig_sub = True
                        break
                marker = (
                    " [PURE]"
                    if is_pure
                    else (" [contiguous-sub]" if is_contig_sub else " [non-contig]")
                )
                print(
                    f"    [{fmt_block(rids)}]{marker}  total={tot}  ({labels})"
                )
        print()


if __name__ == "__main__":
    main()
