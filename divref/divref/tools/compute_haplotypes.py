"""Tool to compute haplotypes from VCF files with gnomAD population frequency annotations."""

import logging
import os
from pathlib import Path

import hail as hl
from fgpyo.io import assert_directory_exists
from fgpyo.io import assert_path_is_readable

from divref import defaults

logger = logging.getLogger(__name__)


def _haploid_adjusted_call(
    locus: hl.Expression,
    gt: hl.Expression,
    sex_karyotype: hl.Expression,
) -> hl.Expression:
    """
    Collapse chrX non-PAR male genotypes to haploid; pass every other call through unchanged.

    Males (`sex_karyotype == "XY"`) at chrX non-PAR loci are encoded as pseudo-homozygous diploid
    (0|0 / 1|1) in the SHAPEIT5 BCFs, but gnomAD reports chrX non-PAR allele numbers with males
    counted as haploid. Replacing such a call with `hl.call(gt[0])` makes the downstream
    `call_stats` AN/AC match the gnomAD convention. It is lossless given the verified encoding
    (`GT[0] == GT[1]` there) but would undercount a heterozygous male non-PAR call if one were
    ever emitted. Autosomes, PAR, and chrY are unaffected.

    Args:
        locus: The variant locus expression.
        gt: The diploid genotype call expression.
        sex_karyotype: The sample's sex-karyotype string expression.

    Returns:
        A call expression: haploid `call(gt[0])` for chrX non-PAR males, else `gt` unchanged.
    """
    is_male_nonpar = locus.in_x_nonpar() & (sex_karyotype == "XY")
    return hl.if_else(is_male_nonpar, hl.call(gt[0]), gt)


def _compute_locus_groups(
    variants_ht: hl.Table,
    window_size: int,
) -> dict[int, int]:
    """
    Assign each variant a global `locus_group` id by cutting at gaps ≥ `window_size`.

    Two variants share a `locus_group` iff a chain of variants in `variants_ht` connects them
    with every consecutive gap (number of intervening reference bases) `< window_size`. Variants
    on different contigs always get different groups. Because connectivity uses any variant
    (not just a single sample's carriers), two carriers in different locus groups for any sample
    are guaranteed to be ≥ `window_size` apart — i.e., partitioning carrier processing by
    `locus_group` cannot put variants into the same parent block that the per-sample
    adjacency-cut would have separated.

    Args:
        variants_ht: variants table with `row_idx` (int64), `locus`, and `ref_len` fields.
        window_size: adjacency-gap threshold in bp.

    Returns:
        `group_of`, mapping `row_idx → locus_group_id`.
    """
    rows = variants_ht.key_by().select("row_idx", "locus", "ref_len").collect()
    logger.info("compute_locus_groups: collected %d variants to the driver", len(rows))
    rows.sort(key=lambda v: (v.locus.contig, v.locus.position))
    group_of: dict[int, int] = {}
    current_group = -1
    prev_end = -1
    prev_contig = ""
    for v in rows:
        if v.locus.contig != prev_contig or v.locus.position - prev_end >= window_size:
            current_group += 1
        group_of[v.row_idx] = current_group
        prev_end = max(prev_end, v.locus.position + v.ref_len)
        prev_contig = v.locus.contig
    return group_of


def _multi_member_locus_groups(rows_with_groups: hl.Table) -> hl.Table:
    """
    Return the `locus_group` ids that have at least two member variants.

    A singleton locus group (one variant) can never form a parent block of ≥2 carriers, so its
    variant never reaches a haplotype. Filtering to these multi-member groups before the
    per-sample entries explosion is exact (identical downstream result) and shrinks that step —
    the dominant early-stage cost on large chromosomes, where most variants are isolated
    singletons. Returns the small group-id table (one row per surviving group) so callers can
    filter by `locus_group` membership directly, without a row-level re-key/join.

    Args:
        rows_with_groups: table with a `locus_group` field.

    Returns:
        Table keyed by `locus_group`, one row per group with ≥2 members.
    """
    sizes = rows_with_groups.group_by(rows_with_groups.locus_group).aggregate(n=hl.agg.count())
    return sizes.filter(sizes.n > 1)


def _form_parent_blocks(
    blocks_ht: hl.Table,
    window_size: int,
) -> hl.Table:
    """
    Form per-(sample, strand) parent blocks via adjacency at `window_size`.

    For each input row (one per `(col_idx, locus_group, strand)`), sorts that row's alt-carrier
    variants by genomic position, walks the sorted list, cuts into blocks at any gap whose
    number of intervening reference bases is ≥ `window_size` (the gap accounts for ref allele
    length, matching `divref.haplotype.variant_distance`), and discards blocks of length < 2.
    Per-sample cuts can still split a single locus group when intermediate variants in the
    group are not carried by this sample.

    Args:
        blocks_ht: Hail table with one row per `(sample, locus_group, strand)`. Required fields:
            - `col_idx` (int): sample identifier.
            - `pop_int` (int): population index.
            - `strand` (int): 0 for left, 1 for right.
            - `carriers` (array of struct(locus, row_idx, ref_len)): variants in this locus
              group where the sample carries the alt allele on this strand. Order is not
              assumed.
        window_size: Adjacency-gap threshold in bp. A gap of exactly `window_size` reference
            bases triggers a split.

    Returns:
        Hail table with one row per (sample, strand, parent block) and fields:
            - `col_idx`, `pop_int`, `strand` inherited from input.
            - `parent_block` (array of struct(locus, row_idx, ref_len)): the block, sorted by
              `locus.position`, length ≥ 2.
        The output is unkeyed.
    """
    blocks_ht = blocks_ht.key_by()
    blocks_ht = blocks_ht.annotate(
        carriers=hl.sorted(
            blocks_ht.carriers,
            key=lambda v: (v.locus.position, v.ref_len, v.row_idx),
        ),
    )
    carriers = blocks_ht.carriers
    n = hl.len(carriers)
    # Cut at index i when carriers[i] starts ≥ window_size beyond the max end of the active
    # block (carriers[0..i-1]), not just the immediately preceding variant. A shorter
    # overlapping carrier can otherwise pull the boundary back and force a false split.
    breakpoints = hl.range(1, n).filter(
        lambda i: carriers[i].locus.position
        - hl.max(hl.range(0, i).map(lambda k: carriers[k].locus.position + carriers[k].ref_len))
        >= window_size
    )

    def get_range(i: hl.Expression) -> hl.Expression:
        start = hl.if_else(i == 0, 0, breakpoints[i - 1])
        end = hl.if_else(i == hl.len(breakpoints), n, breakpoints[i])
        return hl.range(start, end)

    block_ranges = (
        hl.range(0, hl.len(breakpoints) + 1).map(get_range).filter(lambda r: hl.len(r) >= 2)
    )
    blocks_ht = blocks_ht.annotate(
        parent_blocks=block_ranges.map(lambda r: r.map(lambda i: blocks_ht.carriers[i]))
    )
    blocks_ht = blocks_ht.explode("parent_blocks")
    return blocks_ht.select("col_idx", "pop_int", "strand", parent_block=blocks_ht.parent_blocks)


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

    Implementation: a distributed semi-join first restricts `variants_ht` to only the variants
    referenced by some haplotype (most variants are isolated singletons that never enter a parent
    block, so they are dead weight in the broadcast); those are collected into a driver-side dict
    keyed by `row_idx`, broadcast as a Hail literal, and indexed per-haplotype. Driver memory is
    therefore proportional to the number of haplotype-referenced variants, not to all variants in
    `variants_ht` or those passing `variant_freq_threshold` upstream (measured ~6.5x fewer on
    chr2). Hail-side alternatives (explode+group_by, per-element table indexing inside `.map()`)
    were attempted but trigger Hail IR compiler bugs in 0.2.137 — revisit when upgrading.

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
    # Restrict the broadcast to variants that actually appear in a haplotype. Most variants are
    # isolated singletons (median locus-group size 1) that never enter a parent block, so they
    # are dead weight in the broadcast (measured ~6.5x fewer on chr2). The output is identical:
    # every row_idx the per-haplotype lookup below uses is, by construction, referenced here and
    # so retained. Done as a distributed semi-join -- no driver collect, no extra literal.
    referenced = hap_table.key_by()
    referenced = (
        referenced.select(row_idx=referenced.haplotype)
        .explode("row_idx")
        .key_by("row_idx")
        .distinct()
    )
    variants_ht = variants_ht.key_by("row_idx").semi_join(referenced)

    # INSTRUMENTATION: bracket the driver-side collect so we can attribute memory peaks.
    # `len(pairs)` reads the already-materialized Python list (no extra Hail evaluation).
    logger.info("attach_component_info: collecting variant components to driver ...")
    pairs = variants_ht.aggregate(
        hl.agg.collect(
            hl.tuple([
                variants_ht.row_idx,
                hl.struct(
                    locus=variants_ht.locus,
                    alleles=variants_ht.alleles,
                    freq=variants_ht.freq,
                    frequencies_by_pop=variants_ht.frequencies_by_pop,
                ),
            ])
        )
    )
    logger.info("attach_component_info: collected %d variant components for broadcast", len(pairs))
    # Keyed by int32 to match the `hl.int32(idx)` lookup below; the explicit dtype also keeps the
    # literal well-typed when `pairs` is empty (no haplotype-referenced variants), where
    # `hl.literal({})` would otherwise fail to infer the key/value types. The int64->int32 downcast
    # of `row_idx` at the lookup is safe: `row_idx` is a per-contig variant index, well below 2^31.
    components_dict = hl.literal(
        dict(pairs),
        dtype=hl.tdict(
            hl.tint32,
            hl.tstruct(
                locus=variants_ht.locus.dtype,
                alleles=variants_ht.alleles.dtype,
                freq=variants_ht.freq.dtype,
                frequencies_by_pop=variants_ht.frequencies_by_pop.dtype,
            ),
        ),
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
            - `min_variant_frequency` (float64): min over component variants of the *local*
              HGDP+1KG `call_stats` alt AF in `max_pop` — the rarest haplotype-component
              variant's frequency in the population where the haplotype is most common.
            - `all_pop_freqs` (array<struct(pop, empirical_AC, empirical_AF, fraction_phased,
              estimated_gnomad_AF)>): per-population view of all three frequency-derived
              metrics, sorted by `empirical_AF` descending. `pop` indexes into the table's
              `globals.pops` legend. Each entry's `fraction_phased` and `estimated_gnomad_AF`
              are computed using *that pop's own* `empirical_AF` and local
              `call_stats.AF[1]` — not `max_pop`'s. So
              `all_pop_freqs[p].fraction_phased = empirical_AF_p / min_local_AF_in_p` and
              `all_pop_freqs[p].estimated_gnomad_AF = min_i(gnomad_freqs[i][p].AF) *
              all_pop_freqs[p].fraction_phased`. The scalar fields below correspond to the
              entry whose `pop == max_pop` (note `all_pop_freqs` is sorted by `empirical_AF`,
              not pop-indexed, so `all_pop_freqs[max_pop]` would be a positional access — not
              the desired entry).
            - `fraction_phased` (float64): `max_empirical_AF / min_variant_frequency`.
              The proportion of chromosomes carrying the rarest component variant in
              `max_pop` (in HGDP+1KG) that also carry the full haplotype.
            - `estimated_gnomad_AF` (float64): min over segment variants of
              `gnomad_freqs[i][max_pop].AF * fraction_phased`. The haplotype's projected
              frequency in `max_pop` of the broader gnomAD population, extrapolated from
              the HGDP+1KG LD pattern. Equivalent to the `estimated_gnomad_AF` value in the
              `all_pop_freqs` entry whose `pop == max_pop`.
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
        ),
        # Per-population min over component variants of the *local* HGDP+1KG `call_stats` alt
        # AF — the rarest haplotype-component variant's frequency in that pop, measured in
        # the BCF we just processed (not the broader gnomAD sites table). Used as the
        # denominator of `fraction_phased`, which then scales gnomAD per-variant AFs into
        # per-haplotype frequency estimates. Using the local AF (rather than the gnomAD
        # sites AF) is essential: with gnomAD AFs on both sides of the formula the gnomAD
        # multiplier cancels and `estimated_gnomad_AF` collapses to `max_empirical_AF`.
        _min_local_variant_AF_by_pop=pops_range.map(
            lambda p: hl.min(
                hap_table.frequencies_by_pop.map(
                    lambda fbp: fbp.get(p, hl.missing(fbp.dtype.value_type)).AF[1]
                )
            )
        ),
    )
    # Per-pop fraction_phased and estimated_gnomad_AF, computed up front so we can bundle
    # all three frequency-derived metrics into `all_pop_freqs` and then pick out the scalar
    # values at `max_pop` for backward compatibility.
    hap_table = hap_table.annotate(
        _per_pop_fraction_phased=pops_range.map(
            lambda p: hl.bind(
                lambda emp_af, min_var_af: hl.if_else(
                    hl.is_defined(emp_af) & hl.is_defined(min_var_af) & (min_var_af > 0),
                    emp_af / min_var_af,
                    hl.missing(hl.tfloat64),
                ),
                hap_table._per_pop_AF[p],
                hap_table._min_local_variant_AF_by_pop[p],
            )
        ),
    )
    hap_table = hap_table.annotate(
        _per_pop_estimated_gnomad_AF=pops_range.map(
            lambda p: hl.bind(
                lambda fp: hl.if_else(
                    hl.is_defined(fp),
                    hl.min(hap_table.gnomad_freqs.map(lambda x: x[p].AF * fp)),
                    hl.missing(hl.tfloat64),
                ),
                hap_table._per_pop_fraction_phased[p],
            )
        ),
    )
    hap_table = hap_table.annotate(max_pop=hl.argmax(hap_table._per_pop_AF))
    hap_table = hap_table.annotate(
        max_empirical_AF=hap_table._per_pop_AF[hap_table.max_pop],
        max_empirical_AC=hap_table.per_pop_AC[hap_table.max_pop],
        min_variant_frequency=hap_table._min_local_variant_AF_by_pop[hap_table.max_pop],
        fraction_phased=hap_table._per_pop_fraction_phased[hap_table.max_pop],
        estimated_gnomad_AF=hap_table._per_pop_estimated_gnomad_AF[hap_table.max_pop],
        all_pop_freqs=hl.sorted(
            pops_range.map(
                lambda p: hl.struct(
                    pop=p,
                    empirical_AC=hap_table.per_pop_AC[p],
                    empirical_AF=hap_table._per_pop_AF[p],
                    fraction_phased=hap_table._per_pop_fraction_phased[p],
                    estimated_gnomad_AF=hap_table._per_pop_estimated_gnomad_AF[p],
                )
            ),
            key=lambda s: s.empirical_AF,
            reverse=True,
        ),
    )
    return hap_table.drop(
        "_per_pop_AF",
        "_min_local_variant_AF_by_pop",
        "_per_pop_fraction_phased",
        "_per_pop_estimated_gnomad_AF",
    )


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
    """
    Within each AC group, return haplotypes properly contained in a longer group member.

    Sorts each group by haplotype length descending and short-circuits the inner scan when
    the running candidate is no longer strictly longer than the current target — proper
    containment requires `len(h_long) > len(h_short)`, so once equal-length candidates are
    reached, no remaining candidate in the sorted list can subsume the target.
    """
    drop: set[tuple[int, ...]] = set()
    for haps in by_ac.values():
        sorted_haps = sorted(haps, key=len, reverse=True)
        for i, h_short in enumerate(sorted_haps):
            n_short = len(h_short)
            for h_long in sorted_haps[:i]:
                if len(h_long) <= n_short:
                    break
                if _is_contiguous_subarray(h_short, h_long):
                    drop.add(h_short)
                    break
    return drop


def _apply_containment_dedup(hap_table: hl.Table) -> tuple[hl.Table, int]:
    """
    Drop sub-fragment rows subsumed by a larger fragment with identical per-pop AC.

    For each pair of rows `X`, `Y`, drops `X` if:
      - `X.haplotype` is a strictly shorter adjacency-contiguous sub-array of `Y.haplotype`, AND
      - `X.per_pop_AC == Y.per_pop_AC` (element-wise across all populations).

    Implementation: restricts to rows whose `per_pop_AC` vector is shared by at least one
    other row (a drop requires an AC-equal pair, so rows with a unique AC are irrelevant),
    collects only those `(haplotype, per_pop_AC)` to the driver, groups by per-pop AC tuple,
    scans each group for proper sub-array containment, and filters the full Hail table by the
    resulting set of stringified haplotype keys. Driver memory is proportional to the number
    of haplotypes in multi-member AC groups, not the total; the per-AC group counts and the
    membership filter run distributed in Hail. Within each AC group, candidates are sorted by
    length descending and the inner scan exits as soon as the running candidate is no longer
    strictly longer than the target (proper containment requires `len(h_long) > len(h_short)`).
    This keeps the worst case at O(G²) but prunes most pairs in practice when groups span a
    wide length range.

    Args:
        hap_table: Hail table with `haplotype` (array<int64>) and `per_pop_AC` (array<int64>).

    Returns:
        A `(hap_table, n_dropped)` tuple: the table with subsumed rows removed, and the number of
        haplotypes dropped.
    """
    # Only haplotypes that SHARE a per_pop_AC vector with another row can participate in
    # containment (a drop requires `X.per_pop_AC == Y.per_pop_AC`). Rows whose per_pop_AC is
    # unique can neither be subsumed nor subsume anything, so they are irrelevant to the dedup
    # and need not be pulled to the driver. Restricting the collect to multi-member per_pop_AC
    # groups keeps this — the memory-dominant step on large chromosomes, where the collect
    # materializes every haplotype on the driver — proportional to the ambiguous subset only.
    # The per-AC group counts and the semi-join filter run distributed in Hail.
    ac_counts = hap_table.group_by(_ac=hap_table.per_pop_AC).aggregate(n=hl.agg.count())
    ambiguous = hap_table.filter(ac_counts[hap_table.per_pop_AC].n > 1)

    # INSTRUMENTATION: bracket the driver-side collect; `len(rows)`/`len(by_ac)` read
    # already-materialized Python objects (no extra Hail evaluation).
    logger.info("containment_dedup: collecting multi-AC haplotypes to driver ...")
    rows = ambiguous.aggregate(
        hl.agg.collect(hl.struct(haplotype=ambiguous.haplotype, per_pop_AC=ambiguous.per_pop_AC))
    )
    by_ac: dict[tuple[int, ...], list[tuple[int, ...]]] = {}
    for r in rows:
        by_ac.setdefault(tuple(r.per_pop_AC), []).append(tuple(r.haplotype))
    logger.info(
        "containment_dedup: collected %d haplotypes across %d multi-member AC groups",
        len(rows),
        len(by_ac),
    )

    drop = _find_subsumed_haplotypes(by_ac)
    logger.info("containment_dedup: dropping %d subsumed haplotypes", len(drop))
    if not drop:
        return hap_table, 0

    drop_strs = {",".join(str(x) for x in h) for h in drop}
    drop_lit = hl.literal(drop_strs, dtype=hl.tset(hl.tstr))
    hap_table = hap_table.annotate(_hap_str=hl.delimit(hap_table.haplotype.map(hl.str), ","))
    hap_table = hap_table.filter(~drop_lit.contains(hap_table._hap_str))
    return hap_table.drop("_hap_str"), len(drop)


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
    min_partitions: int = 64,
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

    Samples without an inferred population assignment are excluded from the haplotype
    computation (their `pop_int` is undefined, so `compute_haplotypes` can't credit
    their carriers to any population). In gnomAD HGDP+1KG v3.1.2 this drops every aneuploid
    sample (`X`, `XXY`, `XYY`, `ambiguous`) plus samples in populations the workflow doesn't
    configure. The aneuploid samples in particular have `pop = None` because gnomAD's
    PCA-based pop inference declines to assign a population to non-XX/XY karyotypes; this
    drop is therefore the de-facto guarantee that only XX and XY genotypes reach the algorithm.

    On chrX non-PAR loci, males (sex_karyotype == "XY") are treated as haploid. The
    SHAPEIT5 phased BCFs encode all chrX non-PAR genotypes as diploid pseudo (males
    appear as 0|0 / 1|1), but gnomAD reports chrX non-PAR allele numbers with males
    counted as haploid; without the correction, empirical AC/AN would not match the
    gnomAD HT track. Concretely, frequencies_by_pop uses hl.call(GT[0]) for these
    samples, and only the left strand is collected when building haplotype carrier sets.
    Autosomes, PAR1, PAR2, and chrY are unaffected.

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
        output_base: Base output path. Writes intermediate checkpoints
            `{output_base}.variants.ht`, `.blocks.ht`, `.parents.ht`, and `.hap_ac.ht`
            (the tool does not delete them; the Snakemake rule removes them post-run) and
            the final `{output_base}.ht`.
        temp_dir: Local directory for Hail temporary files.
        spark_driver_memory_gb: Memory in GB to allocate to the Spark driver.
        spark_executor_memory_gb: Memory in GB to allocate to the Spark executor.
        min_partitions: Minimum partitions for `import_vcf`. Higher values give finer map-side
            granularity, reducing per-task memory in the downstream entries->blocks shuffle.
            Default 64 (the prior hard-coded value).
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
        min_partitions=min_partitions,
        force_bgz=True,
    )
    mt = mt.select_rows().select_cols()
    mt = mt.annotate_rows(freq=gnomad_va[mt.row_key].pop_freqs)
    mt = mt.filter_rows(hl.is_defined(mt.freq))

    pop_legend: list[str] = gnomad_va.globals.pops.collect()[0]
    n_pops: int = len(pop_legend)
    pop_ints = {pop: i for i, pop in enumerate(pop_legend)}
    sa_row = gnomad_sa[mt.col_key]
    mt = mt.annotate_cols(
        pop_int=hl.literal(pop_ints).get(sa_row.pop),
        sex_karyotype=sa_row.sex_karyotype,
    )
    sample_counts = mt.aggregate_cols(
        hl.struct(total=hl.agg.count(), assigned=hl.agg.count_where(hl.is_defined(mt.pop_int)))
    )
    logger.info(
        "Dropped %d of %d samples whose population is not in the DivRef legend "
        "(a gnomAD ancestry outside %s, or an unassigned population, indicating an "
        "aneuploid karyotype)",
        sample_counts.total - sample_counts.assigned,
        sample_counts.total,
        pop_legend,
    )
    mt = mt.filter_cols(hl.is_defined(mt.pop_int))
    mt = mt.add_row_index().add_col_index()
    mt = mt.filter_entries(mt.freq[mt.pop_int].AF >= variant_freq_threshold)

    # Count chrX non-PAR males as haploid so empirical AC/AN match gnomAD's convention; see
    # `_haploid_adjusted_call`. Autosomes, PAR, and chrY keep the original diploid call.
    adjusted_gt = _haploid_adjusted_call(mt.locus, mt.GT, mt.sex_karyotype)
    mt = mt.annotate_rows(
        # call_stats n_alleles=2: gnomAD HGDP+1KG sites are biallelic (one alt per row).
        frequencies_by_pop=hl.agg.group_by(mt.pop_int, hl.agg.call_stats(adjusted_gt, 2)),
        ref_len=hl.len(mt.alleles[0]),
    )

    variants_ht = mt.rows().select("freq", "row_idx", "frequencies_by_pop", "ref_len")
    variants_ht = variants_ht.checkpoint(f"{str(output_base)}.variants.ht", overwrite=True)

    if variants_ht.head(1).count() == 0:
        raise ValueError(f"No variants found with minimum population AF {variant_freq_threshold}.")

    group_of = _compute_locus_groups(variants_ht, window_size)
    group_lit = hl.literal(group_of, dtype=hl.tdict(hl.tint64, hl.tint32))
    mt = mt.annotate_rows(locus_group=group_lit[mt.row_idx])

    # Drop variants alone in their locus group before the per-sample entries explosion below
    # (the dominant early-stage cost): a singleton group can never form a parent block of ≥2
    # carriers, so its variant never reaches a haplotype. Exact, and shrinks the explosion since
    # most variants are isolated singletons. Filter by `locus_group` membership against the small
    # multi-member-group table -- no row-level re-key shuffle.
    multi_member_groups = _multi_member_locus_groups(mt.rows())
    mt = mt.filter_rows(hl.is_defined(multi_member_groups[mt.locus_group]))

    # Skip the right strand for non-PAR males so each male carrier is counted once (via the
    # left strand). Without this, the pseudo-`1|1` encoding would double-count male
    # haplotypes in chrX non-PAR. Re-derived here rather than reused from the earlier
    # `is_male_nonpar` expression because `mt` was reassigned by `annotate_rows`/
    # `add_col_index` above, so the older expression refers to a stale matrix-table snapshot.
    is_male_nonpar_row = mt.locus.in_x_nonpar() & (mt.sex_karyotype == "XY")
    entries = mt.select_entries(
        is_left=mt.GT[0] != 0,
        is_right=(mt.GT[1] != 0) & ~is_male_nonpar_row,
    ).entries()
    entries = entries.filter(entries.is_left | entries.is_right)
    entries = entries.key_by()
    entries = entries.select(
        "col_idx",
        "pop_int",
        "locus_group",
        "is_left",
        "is_right",
        carrier=hl.struct(locus=entries.locus, row_idx=entries.row_idx, ref_len=entries.ref_len),
    )

    grouped = entries.group_by("col_idx", "locus_group").aggregate(
        pop_int=hl.agg.take(entries.pop_int, 1)[0],
        carriers_left=hl.agg.filter(entries.is_left, hl.agg.collect(entries.carrier)),
        carriers_right=hl.agg.filter(entries.is_right, hl.agg.collect(entries.carrier)),
    )
    grouped = grouped.key_by()
    left = grouped.select("col_idx", "pop_int", carriers=grouped.carriers_left).annotate(strand=0)
    right = grouped.select("col_idx", "pop_int", carriers=grouped.carriers_right).annotate(strand=1)
    blocks_ht = left.union(right)
    blocks_ht = blocks_ht.filter(hl.len(blocks_ht.carriers) >= 2)
    blocks_ht = blocks_ht.checkpoint(f"{str(output_base)}.blocks.ht", overwrite=True)

    parents_ht = _form_parent_blocks(blocks_ht, window_size)
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
    hap_table, n_dropped = _apply_containment_dedup(hap_table)
    hap_table = hap_table.drop("frequencies_by_pop")

    count_after_filter: int = hap_table.count()
    logger.info(
        f"{count_after_filter} haplotypes remaining after filtering and "
        f"containment dedup at window size {window_size} "
        f"({count_after_filter + n_dropped} total before dedup, {n_dropped} subsumed)"
    )

    logger.info("Writing final %s.ht ...", output_base)
    hap_table = hap_table.annotate_globals(pops=hl.literal(pop_legend))
    # Coalesce to a fixed output partition count (independent of the input-side `min_partitions`)
    # to bound the number of part-files written for the final table.
    hap_table.key_by("haplotype").naive_coalesce(64).write(f"{str(output_base)}.ht", overwrite=True)
