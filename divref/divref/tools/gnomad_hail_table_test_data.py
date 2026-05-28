"""Tool to fetch a small subset of a gnomAD sites hail table and HGDP sample information."""

import logging
import random
from pathlib import Path

import hail as hl
from fgpyo.io import assert_path_is_writable

from divref import defaults
from divref.alias import HailPath
from divref.hail import hail_init

logger = logging.getLogger(__name__)


# Per-(pop, sex_karyotype) subsample targets for the test fixture. Counts are rounded
# to 10 and approximate the per-pop XX/XY distribution in the full HGDP+1KG sample
# meta. Aneuploid samples are not included — they have `pop = None` in the upstream
# gnomAD data and would be dropped by `compute_haplotypes`'s population filter anyway.
_TEST_SUBSAMPLE_TYPED_PER_POP: dict[tuple[str, str], int] = {
    ("afr", "XX"): 50,
    ("afr", "XY"): 70,
    ("amr", "XX"): 40,
    ("amr", "XY"): 30,
    ("eas", "XX"): 40,
    ("eas", "XY"): 60,
    ("sas", "XX"): 40,
    ("sas", "XY"): 60,
    ("nfe", "XX"): 40,
    ("nfe", "XY"): 50,
}


def _select_test_subsample(sa: hl.Table, seed: int) -> list[str]:
    """
    Select a stratified subsample of the sample-metadata Hail table for test fixtures.

    Args:
        sa: Sample-metadata Hail table with at minimum the `gnomad_population_inference.pop` and
            `gnomad_sex_imputation.sex_karyotype` fields. Key is `s` (sample id).
        seed: Random seed for reproducibility. The same seed always selects the same sample IDs.

    Returns:
        The list of sample IDs to keep.

    Raises:
        ValueError: if any required `(pop, sex_karyotype)` bucket has fewer samples than the
            configured target — meaning the upstream meta has shrunk and the test plan needs
            updating.
    """
    rng = random.Random(seed)

    # Collect (sample_id, pop, sex_karyotype) for every sample. The table is small (~4151 rows),
    # so driver-side sampling is fine and avoids non-determinism from Hail's randomness.
    rows = sa.select(
        pop=sa.gnomad_population_inference.pop,
        sex_karyotype=sa.gnomad_sex_imputation.sex_karyotype,
    ).collect()

    typed_buckets: dict[tuple[str, str], list[str]] = {}
    for r in rows:
        if r.pop is not None and r.sex_karyotype in ("XX", "XY"):
            typed_buckets.setdefault((r.pop, r.sex_karyotype), []).append(r.s)

    # Sort each bucket so `rng.sample` is deterministic across runs (Hail row order can drift).
    for ids in typed_buckets.values():
        ids.sort()

    kept_ids: list[str] = []
    for (pop, sex_kt), n_target in _TEST_SUBSAMPLE_TYPED_PER_POP.items():
        ids = typed_buckets.get((pop, sex_kt), [])
        if len(ids) < n_target:
            raise ValueError(
                f"Typed bucket ({pop}, {sex_kt}) has {len(ids)} samples; target is {n_target}. "
                f"Upstream sample meta may have changed — update the test subsample plan."
            )
        kept_ids.extend(rng.sample(ids, n_target))

    return kept_ids


def gnomad_hail_table_test_data(
    *,
    in_gnomad_hgdp_variant_annotation_table: HailPath = defaults.GNOMAD_HGDP_1KG_VARIANT_ANNOTATION_HAIL_TABLE,  # noqa: E501
    in_gnomad_hgdp_sample_metadata: HailPath = defaults.GNOMAD_HGDP_1KG_SAMPLE_METADATA_HAIL_TABLE,  # noqa: E501
    loci: list[str],
    out_variant_annotation_dir: Path,
    out_sample_metadata: Path,
    out_samples_txt: Path,
    subsample_seed: int = 0,
    gcs_credentials_path: Path = Path("~/.config/gcloud/application_default_credentials.json"),
    spark_driver_memory_gb: int = 1,
    spark_executor_memory_gb: int = 1,
) -> None:
    """
    Extract subsets of gnomAD HGDP/1KG variant annotations and sample metadata for testing.

    The sample metadata subset is always reduced to the hardcoded test-fixture subsample plan:
    480 stratified-random samples across the five DivRef populations (see
    `_TEST_SUBSAMPLE_TYPED_PER_POP`). The variant annotation subsets cover the full input
    cohort.

    Per-locus variant annotation Hail tables are written as
    `{out_variant_annotation_dir}/{safe_locus_name}.ht`, where `safe_locus_name` is the
    locus string with `:` and `-` replaced by `_` (e.g. `chr1:100001-200000` becomes
    `chr1_100001_200000`).

    Args:
        in_gnomad_hgdp_variant_annotation_table: Path to the gnomAD HGDP/1KG variant annotation
            Hail table.
        in_gnomad_hgdp_sample_metadata: Path to the gnomAD HGDP/1KG sample metadata Hail table.
        loci: Locus intervals for variant filtering. One Hail table is written per locus.
            Sharing a single tool invocation across multiple loci avoids the per-locus Hail /
            Spark startup cost.
        out_variant_annotation_dir: Directory in which the per-locus variant annotation Hail
            tables are written. See the locus-to-filename convention above.
        out_sample_metadata: Output path for the sample metadata Hail table, stripped to key,
            `gnomad_population_inference`, and `gnomad_sex_imputation`, and reduced to the
            test-fixture subsample.
        out_samples_txt: Path to write the kept sample IDs as plain text, one per line. Used
            by downstream `bcftools view -S` to subset VCFs to the same sample set as the
            meta.
        subsample_seed: Random seed used by the subsample plan for reproducible selection.
        gcs_credentials_path: Path to GCS default credentials JSON file.
        spark_driver_memory_gb: Memory in GB to allocate to the Spark driver.
        spark_executor_memory_gb: Memory in GB to allocate to the Spark executor.
    """
    if not loci:
        raise ValueError("loci must contain at least one interval.")

    out_variant_annotation_dir.mkdir(parents=True, exist_ok=True)
    out_paths: list[Path] = []
    for locus in loci:
        safe_name = locus.replace(":", "_").replace("-", "_")
        out_path = out_variant_annotation_dir / f"{safe_name}.ht"
        assert_path_is_writable(out_path)
        out_paths.append(out_path)
    assert_path_is_writable(out_sample_metadata)
    assert_path_is_writable(out_samples_txt)

    hail_init(
        gcs_credentials_path.expanduser(),
        spark_driver_memory_gb=spark_driver_memory_gb,
        spark_executor_memory_gb=spark_executor_memory_gb,
    )

    va = hl.read_table(in_gnomad_hgdp_variant_annotation_table)
    for locus, out_path in zip(loci, out_paths, strict=True):
        roi = [hl.parse_locus_interval(locus, reference_genome=defaults.REFERENCE_GENOME)]
        logger.info(f"Filtering to {locus}.")
        va_subset = hl.filter_intervals(va, roi)
        logger.info(f"Writing {va_subset.count()} variants to {out_path}.")
        va_subset.write(str(out_path), overwrite=True)

    sa = hl.read_table(in_gnomad_hgdp_sample_metadata)
    sa = sa.select("gnomad_population_inference", "gnomad_sex_imputation").select_globals()

    kept_ids = _select_test_subsample(sa, seed=subsample_seed)
    logger.info(f"Test-subsample plan selected {len(kept_ids)} samples.")
    kept_ids_set = hl.literal(set(kept_ids), dtype=hl.tset(hl.tstr))
    sa = sa.filter(kept_ids_set.contains(sa.s))

    # `sorted(kept_ids)` for deterministic file content; the order doesn't matter for
    # `bcftools view -S` but stable order makes git diffs sensible.
    out_samples_txt.write_text("\n".join(sorted(kept_ids)) + "\n")
    logger.info(f"Wrote {len(kept_ids)} sample IDs to {out_samples_txt}.")

    logger.info(f"Writing {sa.count()} samples to {out_sample_metadata}.")
    sa.naive_coalesce(1).write(str(out_sample_metadata), overwrite=True)
