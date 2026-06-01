# Rebuilding DivRef: a customizable resource of human variation and haplotype sequences

Finding matches in the human genome to short sequences is a common problem in bioinformatics.
For example, CRISPR experiments which seek to edit the genome are guided by short sequences matching the desired editing site.
Knowing what possible off-target edits could occur is important: they could have unintended consequences.
The human reference genome assembly doesn't cover all the possible matches that are common in populations across the world. An alternative allele at a variant site, or possibly multiple alleles joined in a haplotype, could be sufficient to create a strong match with real biological effects.

The [DivRef](https://zenodo.org/records/14802613) resource bundle developed by E9 Genomics addresses this. It is a set of FASTA files containing common haplotypes from the HGDP+1KG gnomAD 3.1.2 phased subset of genomes, plus additional common single variants from gnomAD, both with a context of +/-25bp. Common in this case is at least 0.5% allele frequency in one of the five populations covered by the HGDP+1KG dataset: AFR, AMR, EAS, NFE, and SAS. There is also a DuckDB index database with additional information such as the variants in each FASTA record and the per-population allele frequencies. You can use these FASTA files directly with CRISPR off-target site finding tools.

> [!NOTE]
> If you use CALITAS for this purpose, the bundle contains a remapping script that can be used in conjunction with the index database to give you output files with coordinates in the GRCh38 reference genome assembly, making the output seamlessly integrate into workflows that already handle CALITAS output on the reference genome.

The [DivRef generation workflow](https://github.com/e9genomics/human-diversity-reference) is a set of standalone Python scripts plus a Makefile, with some inputs hard-coded and others unrecorded.
I wanted a bundle where I could trust the provenance of every haplotype and variant, and tune parameters for individual Fulcrum clients.
I [re-implemented](https://github.com/fg-labs/divref-wf) it as a Snakemake workflow wrapping a Python toolkit, with a configuration schema, GCS or AWS Open Data ingestion, and per-chromosome parallelism so a full rebuild takes less than a day on a sufficiently large VM.

*Note: end-to-end benchmark numbers (wall time, peak memory, instance specs) are pending a full whole-genome rebuild on a suitably provisioned EC2 instance and will be added before publishing.*

| Parameter | Default | Why you'd change it |
|---|---|---|
| `populations` | AFR, AMR, EAS, NFE, SAS | Restrict or extend ancestry coverage |
| `min_variant_allele_freq` | 0.005 | Tighter / looser variant inclusion |
| `min_haplotype_allele_freq` | 0.005 | Tighter / looser haplotype inclusion |
| `gnomad_variant_annotation_source` | gnomAD 4.1 joint exomes + genomes | Pick a different gnomAD release for the single-variants |
| `sequence_window_size` | 25 | Flanking reference context length around each variant/haplotype |
| `chromosomes` | chr1..chr22, chrX, chrY | Iterate on a subset during development |

The actual config schema separates the HGDP+1KG-derived and gnomAD-derived parameters (see [`workflows/config/config_schema.yml`](../workflows/config/config_schema.yml) for the exact names).

## DivRef 1.1 doesn't contain gnomAD 4.1 variants

The published bundle README states that DivRef's HGDP+1KG-derived haplotypes are "merged with the common variants in the gnomAD 4.1 release".
To verify this, I compared the DivRef 1.1 gnomAD variant sites on chr22 from the DuckDB index against the variants from three different gnomAD sources, all filtered to a minimum 0.5% allele frequency in at least one of the 5 populations.

| gnomAD release | gnomAD variants | In DivRef 1.1 | DivRef 1.1-only | gnomAD-only |
|---|---:|---:|---:|---:|
| 3.1.2 HGDP+1KG subset (4K) | 506,983 | 506,983 | 0 | 0 |
| 3.1.2 genomes (76K) | 511,006 | 506,983 | 0 | 4,023 |
| 4.1 joint (76K genomes, 730K exomes) | 497,031 | 490,931 | 16,052 | 6,100 |

The 3.1.2 HGDP+1KG row matches exactly: every DivRef single-variant entry traces back to the v3.1.2 HGDP+1KG subset.
The 3.1.2 genomes row adds ~4K chr22 variants which are common in the broader 76K-genome cohort.

For v4.1 joint, I keep variants with genome filter `PASS` and exome filter `PASS` or only `AC0`.
The `AC0` flag is gnomAD's "no sample had a high-quality genotype" filter.
On the v4.1 exome track, `AC0` almost always means the position is outside the exome capture footprint rather than that exome data exists and is low quality (looking up the `AC0`-failing variants from a strict-intersection extract shows median exome `AN` = 34 of a possible ~1.46M for a fully called site).
Treating "exome `PASS` or `AC0`-only" as passing on the exome side retains good genome-supported variants that would otherwise be discarded for purely coverage reasons.

16,052 chr22 variants in DivRef 1.1 are not in v4.1 joint above the 0.5% threshold:

| Bucket | Count | % of 16,052 |
|---|---:|---:|
| Fails v4.1 genome filter, passes exome | 9,072 | 56.5% |
| Fails v4.1 exome filter (not AC0-only), passes genome | 4,876 | 30.4% |
| `max(pop_AF) < 0.005` in v4.1 | 1,629 | 10.1% |
| Fails both v4.1 filter sets | 467 | 2.9% |
| Absent from v4.1 callset entirely | 8 | < 0.1% |

Most of the missing variants (89.9%) pass only one or neither of v4.1's two filter sets, or are simply not called at all.
The remaining 10.1% fall below the AF filter in the larger dataset.
A smaller set of ~6K variants are more common in v4.1 joint and now pass the 0.5% threshold.
Extrapolating to the whole genome, I estimate ~420K v4.1 common variants absent from DivRef 1.1 for these 5 populations, and ~1.1M DivRef variants that v4.1's filter sets no longer consider high enough quality or high enough AF to include.

## Improving haplotype computation

The original DivRef algorithm bins variants into staggered 100 bp windows, computes one haplotype-aggregation pass at offset 0 and another at offset 50 bp, and unions the two results. The haplotypes are later split if the gap between any two variants is > 25 bp. Haplotypes with only one variant remaining are dropped, and the rest are de-duplicated, dropping all but one representative of any identity groups. The intent is to avoid losing haplotypes that straddle a fixed bin boundary and later prune down to the desired sequence context size, but the process has two ways it can go wrong:

1. Per-bin grouping can split a true haplotype if it crosses both staggered boundaries, so the haplotype is reconstructed only partially in each pass.
2. When the same haplotype is emitted by both window passes with different AC values due to the inclusion of different samples, the downstream deduplication picks one arbitrarily rather than combining the evidence. This can potentially lead to either under- or over-counting.

I replaced it with a per-sample adjacency algorithm:

1. For each (sample, strand) pair, walk that sample's alt-carrier variants in genomic order. Variants where the sample is ref on this strand are skipped; a window with no alts at all is identical to the reference sequence and is already covered by the genome FASTA, so it doesn't earn a separate haplotype entry.
2. Cut into parent blocks at any ref-aware gap ≥ `window_size` (the gap accounts for ref-allele length, so a 3 bp deletion absorbs 3 bp of the apparent distance).
3. Enumerate every contiguous sub-fragment of length ≥ 2 from each parent block.
4. Per-population AC = count of distinct parent blocks containing each sub-fragment (no double-counting; one block contributes one observation regardless of which sub-fragments it spawns).
5. Containment dedup: drop any sub-fragment whose per-pop AC vector exactly matches a strictly longer enclosing fragment, because it carries no independent evidence.

The result accounts correctly for indel-driven gaps and produces deterministic AC counts.

### Estimating gnomAD allele frequency for a haplotype

Both algorithms filter their output by the same per-haplotype quantity, `estimated_gnomad_AF`.
The name is slightly misleading.
It's not the haplotype's empirical AF in the HGDP+1KG callset, it's a heuristic estimate of how often the haplotype occurs in the broader gnomAD population (remember that v3.1.2 has 76K genomes).
We compute it by scaling each component variant's gnomAD frequency by a measure of how tightly the haplotype's carriers stick together within the HGDP+1KG dataset:

For a given population `p` and haplotype `(v1, v2, ..., vn)`:

1. Count the number of distinct parent blocks in `p` that contain `(v1, v2, ..., vn)` as a contiguous sub-fragment:
  - `hgdp_1kg_haplotype_AC[p]`
2. For each variant `v`, look up:
  - `hgdp_1kg_AF[v,p]` = local HGDP+1KG AF for variant `v` in pop `p`, calculated on the fly from samples in `p`
  - `hgdp_1kg_AN[v,p]` = total HGDP+1KG alleles for `v` in `p` (`hgdp_1kg_AN[v,p]` collapses to 0 for variants where `gnomad_AF[v,p] < min_variant_allele_freq`, because we drop all population `p` entries for those variants)
  - `gnomad_AF[v,p]` = per-variant gnomAD-sites AF for `v` in pop `p`, looked up from gnomAD Hail table
3. Empirical HGDP+1KG AF for the haplotype in this population:
  - `hgdp_1kg_haplotype_AF[p] = hgdp_1kg_haplotype_AC[p] / min(hgdp_1kg_AN[v,p])`
4. Fraction phased. Intuitively, this value is close to 1 when most carriers of the rarest component variant also carry the full haplotype, and close to 0 when they do not.
  - `fraction_phased[p] = hgdp_1kg_haplotype_AF[p] / min(hgdp_1kg_AF[v,p])`
5. Estimated gnomAD AF:
  - `estimated_gnomad_AF[p] = min(gnomad_AF[v,p]) * fraction_phased[p]`

If `min(hgdp_1kg_AN[v,p]) == 0`, then `hgdp_1kg_haplotype_AF[p]` is missing.
If `min(hgdp_1kg_AF[v,p]) == 0` or missing, `fraction_phased[p]` is missing; either makes `estimated_gnomad_AF[p]` missing.

`max_pop` is the population with the highest empirical haplotype AF, so haplotypes common in some ancestries and rare in others aren't penalized by being averaged.

### Sex chromosome haplotypes

The original DivRef workflow doesn't compute haplotypes on chrX or chrY.
The new workflow extends haplotype computation to chrX, with two complications that don't arise on the autosomes.
chrY still only contains single gnomAD variants and their sequence contexts.

#### Three input BCFs

The gnomAD HGDP+1KG v3.1.2 release ships chrX as three separate SHAPEIT5 phased BCFs covering PAR1, non-PAR, and PAR2.
PAR1 and PAR2 use the SHAPEIT5 common-variant track; non-PAR uses the rare-variant track.
The workflow concatenates the three in genomic order into a single chrX VCF before feeding it into the haplotype computation tool.

#### Haploid males in non-PAR

The SHAPEIT5 BCFs encode every chrX non-PAR genotype as pseudo-diploid.
Males appear as `0|0` or `1|1` rather than as a single haploid call.
gnomAD's per-variant chrX non-PAR `AN` counts males as haploid, so without a correction `hgdp_1kg_AN[v, p]` on chrX non-PAR would be inflated by the number of XY samples in pop `p` relative to gnomAD's `AN`, and `hgdp_1kg_AF[v, p]` would be correspondingly deflated.
We apply the correction in two places on chrX non-PAR loci.
First, we only look at the first allele for male samples, so each male contributes a single haploid allele to the per-population aggregation.
This is lossless because both alleles are always the same.
Second, the per-sample adjacency walk only collects carriers from the left strand for non-PAR males; the right strand would double-count the same haploid call.
Autosomes, PAR1, PAR2, and chrY are unaffected.

The ploidy correction relies on every sample having a determinable sex karyotype.
gnomAD HGDP+1KG v3.1.2's PCA-based pop inference declines to assign a population to any sex-aneuploid sample (`X`, `XXY`, `XYY`, `ambiguous`), so those samples are dropped by the population filter and are never counted for any haplotypes on autosomes or chrX.

### Comparing against the original algorithm

Comparing the two algorithms on chr22 (≥0.5% for at least one per-population AF, 25 bp window):

![chr22 haplotype overlap between the original and new algorithms: 29,548 shared, 1,333 original-only, 1,680 new-only](../data/analysis/compute_haplotypes/algo_comparison.venn.png)

1,327 of the 1,333 original-only haplotypes are proper contiguous sub-fragments of some new haplotype with identical per-population AC vectors.
For 1,242 of those, both algorithms find the longer fragment, but only the new recognises the short one as redundant.
For the remaining 85 sub-fragments, the original algorithm emitted the short fragment but never accumulated enough AC for the longer haplotype to pass the AF filter, while the new algorithm's containment counting correctly credited the longer fragment with carriers from every parent block in which it appears.

Of the 1,680 new-only haplotypes, 1,590 contain only variants the original already emitted in some other haplotype, 54 mix shared and novel variants, and 36 are entirely novel.
Most (21) of those 36 are 2 variant haplotypes, one contains 8.

### A closer look at the six original-only haplotypes the new algorithm doesn't emit

There were six haplotypes that weren't subsumed by anything new.
Re-running haplotype computation with the threshold lowered to 0.002 lets us see what the new algorithm computed for each.
We'll label the variants in each haplotype `(v1, v2, ...)`.

| # | variants | gap pattern | description |
|---:|---|---|---|
| 1 | `19951135:G>A`, `19951136:A>G` | 1 bp | Two adjacent SNPs with a third common variant nearby (`19951156:A>G`); every carrier of `(v1, v2)` also carries the third variant. |
| 2 | `40457473:AGAAAGAAAGAAAGAAG>A`, `40457477:AGAAAGAAAGAAG>A` | 4 bp | Two short deletions inside a VNTR region with at least 12 other competing biallelic alleles at 40457457–40457477. |
| 3 | `32402562:G>A`, `32402565:TCAG>T` | 3 bp | Adjacent SNP plus a 3 bp deletion with no relevant intermediate carriers. |
| 4 | `22627350:T>C`, `22627366:G>A` | 16 bp | One common intermediate variant `22627351:G>T` is carried by 5,890 chromosomes, vs 6 carrying the `(v1, v2)` haplotype; the intermediate dominates `v1`'s haplotype landscape. |
| 5 | `24600987:A>G`, `24601005:T>C`, `24601022:A>C` | 18 bp and 17 bp consecutive; 35 bp `v1`→`v3` | Three-variant haplotype where most `v1` carriers (400 chromosomes) don't carry `v2`; their `v1`→`v3` walk gap is 35 bp ≥ 25, so the walk is cut. |
| 6 | `24626868:C>T`, `24626892:A>G` | 24 bp | Non-contiguous subset of a 3-variant new haplotype `(v1, v2, v3)`; `(v1, v3)` is the only old haplotype on chr22 whose variants don't form a contiguous sub-array of any new haplotype. |

| # | old `max_pop` | old `AC` | old `fp` | old `est AF` | new `max_pop` | new `AC` | new `fp` | new `est AF` |
|---:|---|---:|---:|---:|---|---:|---:|---:|
| 1 | afr | 7 | **1.0000** | 0.00505 | afr | 7 | 0.0864 | 0.00272 |
| 2 | amr | 4 | 0.4000 | 0.00572 | **eas** | **6** | 0.0909 | 0.00318 |
| 3 | nfe | 11 | 0.8461 | 0.00531 | **sas** | **14** | 0.5833 | 0.00462 |
| 4 | afr | 6 | **1.0000** | 0.00515 | afr | 6 | 0.8571 | 0.00443 |
| 5 | sas | 9 | **1.0000** | 0.00540 | sas | 9 | 0.9000 | 0.00486 |
| 6 | eas | 10 | 0.1099 | 0.00606 | **sas** | **13** | 0.4643 | 0.00446 |

Shorthand used in the table (all values are taken at the row's `max_pop`):
- `AC` = `hgdp_1kg_haplotype_AC[max_pop]`
- `fp` = `fraction_phased[max_pop]`
- `est AF` = `estimated_gnomad_AF[max_pop]`

**Mechanism 1: `fp = 1.0000` exactly in the original (cases 1, 4, 5).**
These three cases have the same `max_pop` and the same `AC` as new, but old's `fp` is exactly 1 and new's is much lower.
For `fp = hgdp_1kg_haplotype_AF[max_pop] / min(hgdp_1kg_AF[v, max_pop])` to land at exactly 1.0, the numerator must equal the denominator.
That happens repeatedly because the original pipeline's `create_duckdb_index` runs `split_haplotypes(25)` on longer haplotypes and the sub-haplotypes inherit `min(hgdp_1kg_AF[v, max_pop])` from the parent.
The parent's minimum can come from a variant that isn't in the resulting sub-haplotype at all.
If that extra variant happens to be the rarest in the parent's `max_pop` and has the same AC as the haplotype, `min(hgdp_1kg_AF[v, max_pop])` ends up equal to `hgdp_1kg_haplotype_AF[max_pop]` and `fp` collapses to 1.0.
The new algorithm computes `min(hgdp_1kg_AF[v, max_pop])` on the sub-haplotype's actual component variants, so its `fp` reflects only the relevant variants.

**Mechanism 2: `max_pop` shifts under containment counting (cases 2, 3, 6).**
The new single-pass containment counting redistributes `hgdp_1kg_haplotype_AC[p]` across populations differently from the original's per-bin tuple aggregation.
Containment AC includes chromosomes from longer parent blocks where the haplotype appears as a contiguous sub-fragment, not just chromosomes whose exact bin tuple is the haplotype.
In cases 2, 3, and 6 the redistribution moves `max_pop` to a population with a larger sample size (and therefore larger `hgdp_1kg_AN[v, max_pop]`), which lowers `hgdp_1kg_haplotype_AF[max_pop]` for the same `AC` and pulls `est AF` below 0.005.

### AC counts

The new algorithm's `popmax_empirical_AC` tends to be larger than the original's for shared haplotypes.
Among the 29,548 haplotypes found by both algorithms:

| | Count | % of shared |
|---|---:|---:|
| Same popmax AC   | 25,718 | 87.0% |
| New > old        |  3,785 | 12.8% |
| Old > new        |     45 |  0.2% |

Both columns report `popmax_empirical_AC` — the AC at whichever population has the highest empirical haplotype AF.
Two things move that value between the algorithms.
First, the new algorithm's containment counting credits a haplotype with every parent block where it appears as a contiguous sub-fragment, so within any single population the new AC is at least the per-bin tuple AC the original would have computed.
Second, `popmax` may shift to a different population: when it shifts to a population with more carriers, `popmax_AC` rises (the `new > old` rows); when it shifts to a population that happens to have a smaller AN but a higher empirical AF, the same haplotype's `popmax_AC` can fall (the 45 `old > new` rows).

## Summary

The new [divref-wf](https://github.com/fg-labs/divref-wf) produces a haplotype catalog whose frequencies are deterministic and reproducibly comparable across gnomAD releases.
You can rebuild it quickly and painlessly with whatever AF threshold or population mix your application requires.
