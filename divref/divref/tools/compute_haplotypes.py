"""Tool to compute haplotypes from VCF files with gnomAD population frequency annotations."""

import logging
import os
from pathlib import Path

import hail as hl
from fgpyo.io import assert_directory_exists
from fgpyo.io import assert_path_is_readable
from fgpyo.io import assert_path_is_writable

from divref import defaults

logger = logging.getLogger(__name__)


def _form_parent_blocks(
    cols_ht: hl.Table,
    window_size: int,
) -> hl.Table:
    """
    Form per-(sample, strand) parent blocks via adjacency at `window_size`.

    For each sample and each chromosome strand, sorts that strand's alt-carrier variants by
    genomic position, walks the sorted list, cuts into blocks at any gap whose number of
    intervening reference bases is ≥ `window_size` (the gap accounts for ref allele length,
    matching `divref.haplotype.variant_distance`), and discards blocks of length < 2.

    Args:
        cols_ht: Hail table with one row per sample. Required fields:
            - `col_idx` (int): sample identifier.
            - `pop_int` (int): population index.
            - `left_carriers` (array of struct(locus, row_idx, ref_len)): variants on the left
              strand where the sample carries the alt allele. Order is not assumed.
            - `right_carriers` (same shape): variants on the right strand.
        window_size: Adjacency-gap threshold in bp. A gap of exactly `window_size` reference
            bases triggers a split.

    Returns:
        Hail table with one row per (sample, strand, parent block) and fields:
            - `col_idx` (int): inherited from input.
            - `pop_int` (int): inherited from input.
            - `strand` (int): 0 for left, 1 for right.
            - `parent_block` (array of struct(locus, row_idx, ref_len)): the block, sorted by
              `locus.position`, length ≥ 2.
        The output is unkeyed.
    """
    cols_ht = cols_ht.annotate(
        left_carriers=hl.sorted(cols_ht.left_carriers, key=lambda v: v.locus.position),
        right_carriers=hl.sorted(cols_ht.right_carriers, key=lambda v: v.locus.position),
    )

    def _blocks_for(carriers: hl.Expression) -> hl.Expression:
        """Build adjacency-cut blocks of length ≥ 2 from a position-sorted carrier array."""
        n = hl.len(carriers)
        breakpoints = hl.range(1, n).filter(
            lambda i: carriers[i].locus.position
            - carriers[i - 1].locus.position
            - carriers[i - 1].ref_len
            >= window_size
        )

        def get_range(i: hl.Expression) -> hl.Expression:
            start = hl.if_else(i == 0, 0, breakpoints[i - 1])
            end = hl.if_else(i == hl.len(breakpoints), n, breakpoints[i])
            return hl.range(start, end)

        block_ranges = (
            hl.range(0, hl.len(breakpoints) + 1).map(get_range).filter(lambda r: hl.len(r) >= 2)
        )
        return block_ranges.map(lambda r: r.map(lambda i: carriers[i]))

    cols_ht = cols_ht.annotate(
        left_blocks=_blocks_for(cols_ht.left_carriers),
        right_blocks=_blocks_for(cols_ht.right_carriers),
    ).key_by()

    left = cols_ht.select("col_idx", "pop_int", parent_block=cols_ht.left_blocks).annotate(strand=0)
    right = cols_ht.select("col_idx", "pop_int", parent_block=cols_ht.right_blocks).annotate(
        strand=1
    )
    combined = left.union(right)
    combined = combined.explode("parent_block")
    return combined.select("col_idx", "pop_int", "strand", "parent_block")


def _enumerate_subfragments(parents_ht: hl.Table) -> hl.Table:
    """
    Enumerate every adjacency-contiguous sub-fragment of length ≥ 2 from each parent block.

    For a parent block of length N, emits N(N-1)/2 sub-fragments — one per `(i, j)` index
    pair with `0 ≤ i < j < N`, namely `parent_block[i : j + 1]`. The full parent block is
    included as the case `i = 0, j = N - 1`.

    Args:
        parents_ht: Hail table with one row per (sample, strand, parent block) — the output
            of `_form_parent_blocks`. Required fields:
            - `col_idx` (int)
            - `pop_int` (int)
            - `strand` (int)
            - `parent_block` (array of struct(locus, row_idx, ref_len)), length ≥ 2.

    Returns:
        Hail table with one row per (sample, strand, parent block, sub-fragment) and fields:
            - `col_idx`, `pop_int`, `strand` inherited from input.
            - `sub_fragment` (array of struct(locus, row_idx, ref_len)): the sub-fragment,
              length ≥ 2, preserving the parent's variant order.
        Output is unkeyed.
    """
    parents_ht = parents_ht.key_by()
    n = hl.len(parents_ht.parent_block)
    sub_fragments = hl.range(0, n).flatmap(
        lambda i: hl.range(i + 1, n).map(lambda j: parents_ht.parent_block[i : j + 1])
    )
    parents_ht = parents_ht.annotate(sub_fragments=sub_fragments)
    parents_ht = parents_ht.explode("sub_fragments")
    return parents_ht.select("col_idx", "pop_int", "strand", sub_fragment=parents_ht.sub_fragments)


def _aggregate_containment_ac(subfragments_ht: hl.Table, n_pops: int) -> hl.Table:
    """
    Group sub-fragment rows by haplotype and produce a per-population containment AC vector.

    Each (sample, strand, parent block) contributes at most one row per unique sub-fragment
    haplotype string (sub-fragment positions within a parent are distinct), so the row count
    per (haplotype, `pop_int`) equals the number of parent blocks in that population
    containing the sub-fragment as adjacency-contiguous.

    Args:
        subfragments_ht: Hail table with one row per (sample, strand, parent block,
            sub-fragment). Required fields:
            - `pop_int` (int)
            - `sub_fragment` (array of struct(..., row_idx, ...)): the sub-fragment carriers.
        n_pops: Number of populations. The output `per_pop_AC` array has this length, indexed
            by `pop_int`.

    Returns:
        Hail table keyed by `haplotype` with one row per unique sub-fragment haplotype:
            - `haplotype` (array<int64>): row indices of the sub-fragment carriers.
            - `per_pop_AC` (array<int64>): length `n_pops`, AC at index `p` is the count of
              distinct parent blocks in pop `p` whose sub-fragment yielded this haplotype.
    """
    sf = subfragments_ht.key_by()
    sf = sf.select(
        haplotype=sf.sub_fragment.map(lambda v: v.row_idx),
        pop_int=sf.pop_int,
    )
    counts = sf.group_by("haplotype").aggregate(per_pop_counts=hl.agg.counter(sf.pop_int))
    counts = counts.annotate(
        per_pop_AC=hl.range(0, n_pops).map(lambda p: counts.per_pop_counts.get(p, hl.int64(0)))
    )
    return counts.select("per_pop_AC")


def _attach_component_info(hap_table: hl.Table, variants_ht: hl.Table) -> hl.Table:
    """
    Attach per-variant component information to each haplotype row.

    Looks up `(locus, alleles, freq, frequencies_by_pop)` for every `row_idx` in each row's
    `haplotype` array and produces three parallel-length arrays.

    Implementation: collects the variants table into a driver-side dict keyed by `row_idx`,
    broadcasts it as a Hail literal, and indexes per-haplotype. Driver memory is proportional
    to the number of variants in `variants_ht` (i.e., variants that pass `variant_freq_threshold`
    upstream); for typical chr1 inputs this is in the hundreds of thousands of small structs.

    Args:
        hap_table: Hail table with at least these fields:
            - `haplotype` (array<int64>): row indices into `variants_ht.row_idx`.
        variants_ht: Hail table with one row per variant. Must contain fields `row_idx`,
            `locus`, `alleles`, `freq`, and `frequencies_by_pop`.

    Returns:
        `hap_table` annotated with three new arrays parallel-length to `haplotype`:
            - `variants` (array<struct(locus, alleles)>).
            - `gnomad_freqs` (array<freq array>): the per-variant gnomAD frequency array
              (one inner array element per population).
            - `frequencies_by_pop` (array<dict<int, struct(AF, AC, AN)>>): the per-variant
              call-stats grouping. Used downstream for per-population min-AN computation.
    """
    components_dict = hl.literal(
        dict(
            zip(
                variants_ht.row_idx.collect(),
                variants_ht.aggregate(
                    hl.agg.collect(
                        hl.struct(
                            locus=variants_ht.locus,
                            alleles=variants_ht.alleles,
                            freq=variants_ht.freq,
                            frequencies_by_pop=variants_ht.frequencies_by_pop,
                        )
                    )
                ),
                strict=True,
            )
        )
    )
    hap_table = hap_table.annotate(
        _components=hap_table.haplotype.map(lambda idx: components_dict[hl.int32(idx)])
    )
    hap_table = hap_table.annotate(
        variants=hap_table._components.map(lambda c: hl.struct(locus=c.locus, alleles=c.alleles)),
        gnomad_freqs=hap_table._components.map(lambda c: c.freq),
        frequencies_by_pop=hap_table._components.map(lambda c: c.frequencies_by_pop),
    )
    return hap_table.drop("_components")


def _compute_metrics(hap_table: hl.Table, n_pops: int) -> hl.Table:
    """
    Annotate per-haplotype frequency / phasing summary fields.

    For each population `p`:
      - `min_AN_p = min over segment variants of frequencies_by_pop[p].AN`
      - `empirical_AF_p = per_pop_AC[p] / min_AN_p` (or missing when `min_AN_p == 0`)

    Then `max_pop = argmax(empirical_AF)`, and the remaining summary fields are derived from
    that population's component data. Sorted `all_pop_freqs` puts the largest AF first;
    populations with missing AF sort to the end.

    Args:
        hap_table: Hail table with `haplotype`, `per_pop_AC`, `variants`, `gnomad_freqs`, and
            `frequencies_by_pop` already attached (see `_attach_component_info`).
        n_pops: Number of populations. Drives the per-pop iteration and the length of the
            `all_pop_freqs` array.

    Returns:
        `hap_table` with these added fields:
            - `max_pop` (int)
            - `max_empirical_AF` (float64)
            - `max_empirical_AC` (int)
            - `min_variant_frequency` (float64): min gnomAD AF over segment variants for `max_pop`.
            - `all_pop_freqs` (array<struct(pop, empirical_AC, empirical_AF)>): sorted by AF desc.
            - `fraction_phased` (float64): `max_empirical_AF / min_variant_frequency`.
            - `estimated_gnomad_AF` (float64): min over segment variants of
              `gnomad_freqs[i][max_pop].AF * fraction_phased`.
    """
    pops_range = hl.range(0, n_pops)
    hap_table = hap_table.annotate(
        _per_pop_AF=pops_range.map(
            lambda p: hl.bind(
                lambda min_an: hl.if_else(
                    hl.is_defined(min_an) & (min_an > 0),
                    hl.float64(hap_table.per_pop_AC[p]) / hl.float64(min_an),
                    hl.missing(hl.tfloat64),
                ),
                hl.min(
                    hap_table.frequencies_by_pop.map(
                        lambda fbp: fbp.get(p, hl.missing(fbp.dtype.value_type)).AN
                    )
                ),
            )
        )
    )
    hap_table = hap_table.annotate(max_pop=hl.argmax(hap_table._per_pop_AF))
    hap_table = hap_table.annotate(
        max_empirical_AF=hap_table._per_pop_AF[hap_table.max_pop],
        max_empirical_AC=hap_table.per_pop_AC[hap_table.max_pop],
        min_variant_frequency=hl.min(hap_table.gnomad_freqs.map(lambda x: x[hap_table.max_pop].AF)),
        all_pop_freqs=hl.sorted(
            pops_range.map(
                lambda p: hl.struct(
                    pop=p,
                    empirical_AC=hap_table.per_pop_AC[p],
                    empirical_AF=hap_table._per_pop_AF[p],
                )
            ),
            key=lambda s: s.empirical_AF,
            reverse=True,
        ),
    )
    hap_table = hap_table.annotate(
        fraction_phased=hap_table.max_empirical_AF / hap_table.min_variant_frequency,
    )
    hap_table = hap_table.annotate(
        estimated_gnomad_AF=hl.min(
            hap_table.gnomad_freqs.map(
                lambda x: x[hap_table.max_pop].AF * hap_table.fraction_phased
            )
        ),
    )
    return hap_table.drop("_per_pop_AF")


def _is_contiguous_subarray(short: tuple[int, ...], long_: tuple[int, ...]) -> bool:
    """Return True if `short` appears as a strictly shorter contiguous slice of `long_`."""
    n, m = len(short), len(long_)
    if n >= m:
        return False
    for i in range(m - n + 1):
        if long_[i : i + n] == short:
            return True
    return False


def _find_subsumed_haplotypes(
    by_ac: dict[tuple[int, ...], list[tuple[int, ...]]],
) -> set[tuple[int, ...]]:
    """Within each AC group, return haplotypes properly contained in a longer group member."""
    drop: set[tuple[int, ...]] = set()
    for haps in by_ac.values():
        # Sort by length descending so each candidate is checked against earlier (longer) ones.
        sorted_haps = sorted(haps, key=len, reverse=True)
        for i, h_short in enumerate(sorted_haps):
            for h_long in sorted_haps[:i]:
                if _is_contiguous_subarray(h_short, h_long):
                    drop.add(h_short)
                    break
    return drop


def _apply_containment_dedup(hap_table: hl.Table) -> hl.Table:
    """
    Drop sub-fragment rows subsumed by a larger fragment with identical per-pop AC.

    For each pair of rows `X`, `Y`, drops `X` if:
      - `X.haplotype` is a strictly shorter adjacency-contiguous sub-array of `Y.haplotype`, AND
      - `X.per_pop_AC == Y.per_pop_AC` (element-wise across all populations).

    Implementation: collects `(haplotype, per_pop_AC)` to the driver, groups by per-pop AC
    tuple, scans each group for proper sub-array containment, and filters the Hail table by
    the resulting set of stringified haplotype keys. Driver memory is proportional to the
    number of unique haplotypes; pairwise scan within each AC group is O(G²) padded-string
    contains, which is bounded because most groups are small.

    Args:
        hap_table: Hail table with `haplotype` (array<int64>) and `per_pop_AC` (array<int64>).

    Returns:
        `hap_table` with subsumed rows removed.
    """
    rows = hap_table.aggregate(
        hl.agg.collect(hl.struct(haplotype=hap_table.haplotype, per_pop_AC=hap_table.per_pop_AC))
    )
    by_ac: dict[tuple[int, ...], list[tuple[int, ...]]] = {}
    for r in rows:
        by_ac.setdefault(tuple(r.per_pop_AC), []).append(tuple(r.haplotype))

    drop = _find_subsumed_haplotypes(by_ac)
    if not drop:
        return hap_table

    drop_strs = {",".join(str(x) for x in h) for h in drop}
    drop_lit = hl.literal(drop_strs, dtype=hl.tset(hl.tstr))
    hap_table = hap_table.annotate(_hap_str=hl.delimit(hap_table.haplotype.map(hl.str), ","))
    hap_table = hap_table.filter(~drop_lit.contains(hap_table._hap_str))
    return hap_table.drop("_hap_str")


def compute_haplotypes(
    *,
    vcfs_path: Path,
    gnomad_va_file: Path,
    gnomad_sa_file: Path,
    window_size: int,
    variant_freq_threshold: float,
    haplotype_freq_threshold: float,
    output_base: Path,
    temp_dir: Path = Path("/tmp"),
    spark_driver_memory_gb: int = 1,
    spark_executor_memory_gb: int = 1,
) -> None:
    """
    Compute population haplotypes from VCF files with gnomAD frequency annotations.

    Reads VCF files, annotates variants with gnomAD population allele frequencies, and
    extracts phased haplotypes via per-sample adjacency block formation: for each
    (sample, strand) pair, alt-carrier variants are walked in genomic order and grouped into
    parent blocks at gaps ≥ `window_size`. Every adjacency-contiguous sub-fragment of length
    ≥ 2 is enumerated from each parent, with `empirical_AC` counted as the number of parent
    blocks containing it. Sub-fragments are dropped when subsumed by a longer surviving
    fragment with identical per-population AC.

    Args:
        vcfs_path: Path or glob pattern to input VCF files.
        gnomad_va_file: Path to the gnomAD variant annotations Hail table
            (from extract_gnomad_afs).
        gnomad_sa_file: Path to the gnomAD sample metadata Hail table
            (from extract_sample_metadata).
        window_size: Adjacency-gap threshold in bp for both per-sample parent block formation
            and per-parent sub-fragment emission.
        variant_freq_threshold: Minimum gnomAD population allele frequency to retain a variant.
        haplotype_freq_threshold: Minimum estimated gnomAD allele frequency for the haplotype to
            be retained.
        output_base: Base output path; writes `{output_base}.variants.ht` (per-variant
            checkpoint) and the final `{output_base}.ht`.
        temp_dir: Local directory for Hail temporary files.
        spark_driver_memory_gb: Memory in GB to allocate to the Spark driver.
        spark_executor_memory_gb: Memory in GB to allocate to the Spark executor.
    """
    assert_path_is_readable(vcfs_path)
    assert_directory_exists(gnomad_va_file)
    assert_directory_exists(gnomad_sa_file)

    if spark_driver_memory_gb < 1:
        raise ValueError(
            f"Spark driver memory must be at least 1GB. Saw {spark_driver_memory_gb}GB."
        )
    if spark_executor_memory_gb < 1:
        raise ValueError(
            f"Spark executor memory must be at least 1GB. Saw {spark_executor_memory_gb}GB."
        )

    os.environ["PYSPARK_SUBMIT_ARGS"] = (
        f"--driver-memory {spark_driver_memory_gb}g "
        f"--executor-memory {spark_executor_memory_gb}g "
        "pyspark-shell"
    )
    hl.init(tmp_dir=str(temp_dir))

    gnomad_sa = hl.read_table(str(gnomad_sa_file))
    gnomad_va = hl.read_table(str(gnomad_va_file))
    gnomad_va = gnomad_va.filter(
        hl.max(gnomad_va.pop_freqs.map(lambda x: x.AF)) >= variant_freq_threshold
    )

    mt = hl.import_vcf(
        str(vcfs_path),
        reference_genome=defaults.REFERENCE_GENOME,
        min_partitions=64,
        force_bgz=True,
    )
    mt = mt.select_rows().select_cols()
    mt = mt.annotate_rows(freq=gnomad_va[mt.row_key].pop_freqs)
    mt = mt.filter_rows(hl.is_defined(mt.freq))

    pop_legend: list[str] = gnomad_va.globals.pops.collect()[0]
    n_pops: int = len(pop_legend)
    pop_ints = {pop: i for i, pop in enumerate(pop_legend)}
    mt = mt.annotate_cols(pop_int=hl.literal(pop_ints).get(gnomad_sa[mt.col_key].pop))
    mt = mt.filter_cols(hl.is_defined(mt.pop_int))
    mt = mt.add_row_index().add_col_index()
    mt = mt.filter_entries(mt.freq[mt.pop_int].AF >= variant_freq_threshold)

    mt = mt.annotate_rows(
        frequencies_by_pop=hl.agg.group_by(mt.pop_int, hl.agg.call_stats(mt.GT, 2)),
    )
    mt = mt.annotate_cols(
        left_carriers=hl.agg.filter(
            mt.GT[0] != 0,
            hl.agg.collect(
                hl.struct(
                    locus=mt.locus,
                    row_idx=mt.row_idx,
                    ref_len=hl.len(mt.alleles[0]),
                )
            ),
        ),
        right_carriers=hl.agg.filter(
            mt.GT[1] != 0,
            hl.agg.collect(
                hl.struct(
                    locus=mt.locus,
                    row_idx=mt.row_idx,
                    ref_len=hl.len(mt.alleles[0]),
                )
            ),
        ),
    )

    variants_ht = mt.rows().select("freq", "row_idx", "frequencies_by_pop")
    variants_ht = variants_ht.checkpoint(f"{str(output_base)}.variants.ht", overwrite=True)

    if variants_ht.head(1).count() == 0:
        raise ValueError(f"No variants found with minimum population AF {variant_freq_threshold}.")

    cols_ht = mt.cols().select("col_idx", "pop_int", "left_carriers", "right_carriers")
    cols_ht = cols_ht.checkpoint(f"{str(output_base)}.cols.ht", overwrite=True)

    parents_ht = _form_parent_blocks(cols_ht, window_size)
    parents_ht = parents_ht.checkpoint(f"{str(output_base)}.parents.ht", overwrite=True)
    subfragments_ht = _enumerate_subfragments(parents_ht)
    hap_table = _aggregate_containment_ac(subfragments_ht, n_pops)
    hap_table = hap_table.checkpoint(f"{str(output_base)}.hap_ac.ht", overwrite=True)
    hap_table = _attach_component_info(hap_table, variants_ht.key_by())
    hap_table = _compute_metrics(hap_table, n_pops)

    hap_table = hap_table.filter(
        (hap_table.min_variant_frequency > 0)
        & (hap_table.estimated_gnomad_AF >= haplotype_freq_threshold)
    )
    hap_table = _apply_containment_dedup(hap_table)
    hap_table = hap_table.drop("frequencies_by_pop")

    count_after_filter: int = hap_table.count()
    logger.info(
        f"{count_after_filter} haplotypes remaining after filtering and "
        f"containment dedup at window size {window_size}"
    )

    logger.info("Writing final %s.ht ...", output_base)
    hap_table = hap_table.annotate_globals(pops=hl.literal(pop_legend))
    hap_table.key_by("haplotype").naive_coalesce(64).write(f"{str(output_base)}.ht", overwrite=True)
