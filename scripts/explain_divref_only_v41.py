"""Explain why DivRef 1.1 chr22 variants are absent from the gnomAD v4.1 joint sites table
under the intersection filter.

Reads the `compare_divref_gnomad` workflow's `chr22.joint_41.divref_not_in_gnomad.tsv`
(DivRef variants with no match in v4.1 joint after both-filters PASS), looks each one up
in the v4.1 joint sites Hail table, and bins by reason for absence:

  - absent_from_v41:      not in the HT at all (e.g. v4.1 callset didn't call the site)
  - below_af_threshold:   in HT, both filter sets empty (= PASS), but `max(pop_AF) < 0.005`
                          across the selected populations
  - exome_filter_only:    in HT with `max(pop_AF) >= 0.005`, exomes.filters non-empty,
                          genomes.filters empty
  - genome_filter_only:   in HT with `max(pop_AF) >= 0.005`, exomes.filters empty,
                          genomes.filters non-empty
  - both_filters_nonempty: in HT with `max(pop_AF) >= 0.005`, both filter sets non-empty

For the three filtered buckets, also reports the top filter codes observed.

By default points at the GCS public-data v4.1 joint sites table. Run via `pixi run python`
so the Hail+Spark+gcs-connector environment is set up.
"""

import argparse
from collections import Counter
from pathlib import Path

import hail as hl

from divref import defaults
from divref.hail import hail_init


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--variants-tsv",
        default="data/analysis/compare_divref_gnomad/chr22.joint_41.divref_not_in_gnomad.tsv",
        help="TSV with a `variants` column of `chr:pos:ref:alt` strings to look up. "
        "Default: the compare_divref_gnomad workflow's chr22 joint_41 divref-only list.",
    )
    parser.add_argument(
        "--gnomad-joint-ht",
        default=(
            "gs://gcp-public-data--gnomad/release/4.1/ht/joint/gnomad.joint.v4.1.sites.ht"
        ),
        help="URI of the gnomAD v4.1 joint sites Hail table. Default points at the GCS "
        "public-data mirror.",
    )
    parser.add_argument(
        "--contig",
        default="chr22",
        help="Contig to restrict the HT scan to.",
    )
    parser.add_argument(
        "--af-threshold",
        type=float,
        default=0.005,
        help="Per-population AF threshold used in the original comparison.",
    )
    parser.add_argument(
        "--populations",
        nargs="+",
        default=defaults.POPULATIONS,
        help="Population codes whose `joint.freq[pop].AF` to take the max over for the "
        "`below_af_threshold` bucket.",
    )
    parser.add_argument(
        "--gcs-credentials-path",
        type=Path,
        default=Path("~/.config/gcloud/application_default_credentials.json").expanduser(),
        help="Path to GCP application default credentials JSON.",
    )
    parser.add_argument(
        "--top-n-filter-codes",
        type=int,
        default=15,
        help="How many distinct filter-code combinations to list per bucket.",
    )
    parser.add_argument(
        "--spark-driver-memory-gb",
        type=int,
        default=16,
        help="Spark driver memory in GB.",
    )
    parser.add_argument(
        "--spark-executor-memory-gb",
        type=int,
        default=16,
        help="Spark executor memory in GB.",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()

    # Load the variant list (TSV header: `variants` column with chr:pos:ref:alt strings).
    with Path(args.variants_tsv).open() as fh:
        header = fh.readline().strip()
        if header != "variants":
            raise ValueError(
                f"Expected first line of {args.variants_tsv} to be `variants`, got {header!r}"
            )
        target_variants = [line.strip() for line in fh if line.strip()]
    print(f"Loaded {len(target_variants)} target variants from {args.variants_tsv}")

    hail_init(
        gcs_credentials_path=args.gcs_credentials_path,
        spark_driver_memory_gb=args.spark_driver_memory_gb,
        spark_executor_memory_gb=args.spark_executor_memory_gb,
    )

    # Build a Hail set of `chr:pos:ref:alt` strings for fast filtering driver-side.
    target_set = hl.literal(set(target_variants), dtype=hl.tset(hl.tstr))

    # Read v4.1 joint sites HT, restrict to contig, build the canonical variant key.
    ht = hl.read_table(args.gnomad_joint_ht)
    ht = hl.filter_intervals(ht, [hl.parse_locus_interval(args.contig, reference_genome="GRCh38")])
    ht = ht.annotate(
        _variant=hl.format(
            "%s:%s:%s:%s",
            ht.locus.contig,
            hl.str(ht.locus.position),
            ht.alleles[0],
            ht.alleles[1],
        )
    )
    ht = ht.filter(target_set.contains(ht._variant))

    # Build a map population -> joint.freq index by walking the freq_meta global.
    # Match ONLY entries with exactly `{group: "adj", gen_anc: <pop>}` — additional
    # keys like `subset`, `downsampling`, `sex` indicate restricted freq slices that
    # would silently shadow the primary adj+gen_anc entry if the loop is not strict.
    freq_meta = hl.eval(ht.joint_globals.freq_meta)
    pop_to_idx = {}
    for i, meta in enumerate(freq_meta):
        if dict(meta) == {"group": "adj", "gen_anc": meta.get("gen_anc")} and meta.get("gen_anc") in args.populations:
            pop_to_idx[meta["gen_anc"]] = i
    missing = [p for p in args.populations if p not in pop_to_idx]
    if missing:
        raise ValueError(
            f"Could not find joint.freq entries for populations {missing}. "
            f"Found keys: {sorted(pop_to_idx)}"
        )
    print(f"Resolved population indices: {pop_to_idx}")

    # Find the index of the adj-overall (population-pooled, adj-filtered) entry in
    # exomes.freq. That entry's `AN` tells us whether the exome assay saw the position
    # at all: AN == 0 means no exome coverage, AN > 0 means coverage but no alt.
    exome_freq_meta = hl.eval(ht.exomes_globals.freq_meta)
    exome_adj_idx = None
    for i, meta in enumerate(exome_freq_meta):
        # Strict exact-match: only the primary `{group: "adj"}` entry, no other keys.
        if dict(meta) == {"group": "adj"}:
            exome_adj_idx = i
            break
    if exome_adj_idx is None:
        raise ValueError(
            "Could not find the adj-overall entry in exomes_globals.freq_meta. "
            f"First few entries: {exome_freq_meta[:5]}"
        )
    print(f"exomes.freq adj-overall index: {exome_adj_idx} ({exome_freq_meta[exome_adj_idx]})")

    pop_af_struct = hl.struct(
        **{
            pop: hl.coalesce(ht.joint.freq[idx].AF, 0.0)
            for pop, idx in pop_to_idx.items()
        }
    )
    ht = ht.annotate(
        _pop_afs=pop_af_struct,
        _max_pop_af=hl.max(list(pop_af_struct.values())),
        _exome_filters=hl.sorted(hl.array(ht.exomes.filters)),
        _genome_filters=hl.sorted(hl.array(ht.genomes.filters)),
        _exome_AC=hl.coalesce(ht.exomes.freq[exome_adj_idx].AC, 0),
        _exome_AN=hl.coalesce(ht.exomes.freq[exome_adj_idx].AN, 0),
    )

    rows = ht.select(
        "_variant", "_pop_afs", "_max_pop_af",
        "_exome_filters", "_genome_filters",
        "_exome_AC", "_exome_AN",
    ).collect()
    present_variants = {r._variant for r in rows}
    absent = [v for v in target_variants if v not in present_variants]

    bucket_counts: Counter = Counter()
    bucket_filter_codes: dict[str, Counter] = {
        "below_af_threshold": Counter(),
        "exome_filter_only": Counter(),
        "genome_filter_only": Counter(),
        "both_filters_nonempty": Counter(),
    }
    bucket_pop_afs: dict[str, dict[str, list[float]]] = {
        b: {p: [] for p in args.populations}
        for b in (
            "below_af_threshold",
            "exome_filter_only",
            "genome_filter_only",
            "both_filters_nonempty",
        )
    }
    variant_bucket: dict[str, str] = {}
    for r in rows:
        ef = list(r._exome_filters or [])
        gf = list(r._genome_filters or [])
        passes_both = (len(ef) == 0) and (len(gf) == 0)
        af = r._max_pop_af or 0.0
        if passes_both and af >= args.af_threshold:
            # Both filter sets empty AND AF >= threshold: should have appeared in the
            # comparison's intersection of v4.1 PASS + AF threshold. Flag as anomaly.
            bucket = "unexpected_pass"
        elif passes_both or af < args.af_threshold:
            bucket = "below_af_threshold"
        elif len(ef) > 0 and len(gf) == 0:
            bucket = "exome_filter_only"
        elif len(ef) == 0 and len(gf) > 0:
            bucket = "genome_filter_only"
        else:
            bucket = "both_filters_nonempty"

        bucket_counts[bucket] += 1
        variant_bucket[r._variant] = bucket
        if bucket in bucket_filter_codes:
            bucket_filter_codes[bucket][(tuple(ef), tuple(gf))] += 1
        if bucket in bucket_pop_afs:
            for pop in args.populations:
                bucket_pop_afs[bucket][pop].append(getattr(r._pop_afs, pop) or 0.0)

    bucket_counts["absent_from_v41"] = len(absent)
    for v in absent:
        variant_bucket[v] = "absent_from_v41"

    total = sum(bucket_counts.values())
    print()
    print(f"=== Bucket summary (n={total} of {len(target_variants)} input variants) ===")
    for name in (
        "absent_from_v41",
        "below_af_threshold",
        "exome_filter_only",
        "genome_filter_only",
        "both_filters_nonempty",
        "unexpected_pass",
    ):
        count = bucket_counts.get(name, 0)
        pct = (100.0 * count / len(target_variants)) if target_variants else 0.0
        print(f"  {name:<23} {count:>8}  ({pct:5.1f}%)")

    for bucket in ("exome_filter_only", "genome_filter_only", "both_filters_nonempty",
                   "below_af_threshold"):
        codes = bucket_filter_codes[bucket]
        if not codes:
            continue
        print()
        print(f"=== top filter-code combos in `{bucket}` (top {args.top_n_filter_codes}) ===")
        for (ef, gf), n in codes.most_common(args.top_n_filter_codes):
            ef_s = "[]" if not ef else "[" + ",".join(ef) + "]"
            gf_s = "[]" if not gf else "[" + ",".join(gf) + "]"
            print(f"  exomes={ef_s:<30} genomes={gf_s:<30} count={n}")

    def _quantiles(values: list[float]) -> str:
        if not values:
            return "n=0"
        s = sorted(values)
        n = len(s)
        return (
            f"n={n} min={s[0]:.4g} p10={s[n // 10]:.4g} p50={s[n // 2]:.4g} "
            f"p90={s[min(n - 1, n * 9 // 10)]:.4g} max={s[-1]:.4g}"
        )

    for bucket in ("exome_filter_only", "genome_filter_only", "both_filters_nonempty",
                   "below_af_threshold"):
        per_pop = bucket_pop_afs[bucket]
        if not any(per_pop.values()):
            continue
        print()
        print(f"=== per-pop AF distribution in `{bucket}` ===")
        for pop in args.populations:
            print(f"  {pop}: {_quantiles(per_pop[pop])}")

    # Exome coverage check for the AC0-only-failing variants in the
    # `exome_filter_only` bucket. AC0 means "no high-quality alt genotype in the
    # exomes"; on a whole-genome input that is expected for any variant outside the
    # exome capture target. `exomes.freq[adj].AN` is the number of called alleles in
    # the exomes track at this site -- AN == 0 means no exome data at all (position
    # outside capture or otherwise uncalled), AN > 0 means the assay saw the position
    # and genuinely observed no alts.
    ac0_an_values: list[int] = []
    for r in rows:
        ef = list(r._exome_filters or [])
        gf = list(r._genome_filters or [])
        if "AC0" not in ef or len(gf) != 0:
            continue
        ac0_an_values.append(int(r._exome_AN or 0))

    if ac0_an_values:
        ac0_total = len(ac0_an_values)
        n_zero = sum(1 for v in ac0_an_values if v == 0)
        n_under_100 = sum(1 for v in ac0_an_values if 0 < v < 100)
        n_100_to_10k = sum(1 for v in ac0_an_values if 100 <= v < 10_000)
        n_10k_to_100k = sum(1 for v in ac0_an_values if 10_000 <= v < 100_000)
        n_100k_plus = sum(1 for v in ac0_an_values if v >= 100_000)

        print()
        print(f"=== exome AN distribution for AC0-only exome failures (n={ac0_total}) ===")

        def _line(label: str, count: int) -> str:
            pct = (100.0 * count / ac0_total) if ac0_total else 0.0
            return f"  {label:<32} {count:>8}  ({pct:5.1f}%)"

        print(_line("AN == 0 (no exome coverage)", n_zero))
        print(_line("0 < AN < 100", n_under_100))
        print(_line("100 <= AN < 10,000", n_100_to_10k))
        print(_line("10,000 <= AN < 100,000", n_10k_to_100k))
        print(_line("AN >= 100,000 (in capture)", n_100k_plus))
        print(f"  AN quantiles: {_quantiles([float(v) for v in ac0_an_values])}")

    # Per-variant TSV next to the input file: original variant, bucket, per-pop AFs,
    # filter codes. Lets downstream consumers (the blog / spot-check scripts) read the
    # explanation row-by-row without re-running Hail.
    input_path = Path(args.variants_tsv)
    out_path = input_path.with_suffix("")
    out_path = out_path.with_name(out_path.name + ".explained.tsv")
    rows_by_variant = {r._variant: r for r in rows}
    pop_cols = list(args.populations)
    with out_path.open("w") as fh:
        header = (
            "variant\tbucket\tmax_pop_af\t"
            + "\t".join(f"AF_{p}" for p in pop_cols)
            + "\texomes_filters\tgenomes_filters\texomes_AC\texomes_AN"
        )
        fh.write(header + "\n")
        n_pad_cells = 1 + len(pop_cols) + 2 + 2  # max_pop_af + per-pop + filters + AC/AN
        for v in target_variants:
            bucket = variant_bucket.get(v, "absent_from_v41")
            r = rows_by_variant.get(v)
            if r is None:
                cells = ["NA"] * n_pad_cells
                fh.write(f"{v}\t{bucket}\t" + "\t".join(cells) + "\n")
                continue
            af_cells = "\t".join(
                f"{getattr(r._pop_afs, p):.5f}" if getattr(r._pop_afs, p) is not None else "NA"
                for p in pop_cols
            )
            ef = ",".join(r._exome_filters) if r._exome_filters else ""
            gf = ",".join(r._genome_filters) if r._genome_filters else ""
            fh.write(
                f"{v}\t{bucket}\t{r._max_pop_af:.5f}\t{af_cells}\t{ef}\t{gf}"
                f"\t{int(r._exome_AC or 0)}\t{int(r._exome_AN or 0)}\n"
            )
    print()
    print(f"wrote per-variant explanation TSV: {out_path}")


if __name__ == "__main__":
    main()
