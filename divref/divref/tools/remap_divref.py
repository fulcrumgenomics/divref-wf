"""Tool to remap DivRef haplotype coordinates to reference genome coordinates."""

import csv
import json
import logging
import os
from pathlib import Path
from typing import Any
from typing import Optional

import duckdb
import pandas as pd
from pydantic import BaseModel
from pydantic import ConfigDict
from pydantic import Field
from tqdm import tqdm

logger = logging.getLogger(__name__)

_GNOMAD_AF_COLUMN_PREFIX = "gnomAD_AF_"
_ESTIMATED_GNOMAD_AF_COLUMN_PREFIX = "estimated_gnomad_AF_"


class Variant(BaseModel):
    """A genomic variant with chromosome, position, reference, and alternate alleles."""

    chromosome: str
    position: int
    reference: str
    alternate: str

    def render(self) -> str:
        """
        Return the variant in colon-delimited format.

        Returns:
            String in the form chromosome:position:reference:alternate.
        """
        return f"{self.chromosome}:{self.position}:{self.reference}:{self.alternate}"


class ReferenceMapping(BaseModel):
    """A mapped interval on the reference genome corresponding to a DivRef haplotype region."""

    chromosome: str
    start: int
    end: int
    variants_involved: list[Variant]
    first_variant_index: Optional[int]
    last_variant_index: Optional[int]
    population_frequencies: dict[str, list[float]]

    def variants_involved_str(self) -> str:
        """
        Return a comma-delimited string of all variants involved in this mapping.

        Returns:
            Comma-separated variant strings in chromosome:position:reference:alternate format.
        """
        return ",".join([v.render() for v in self.variants_involved])


class Haplotype(BaseModel):
    """A DivRef haplotype sequence with metadata and population frequency information."""

    # Field names use aliases to match DuckDB column names (which use mixedCase). Extra columns
    # (e.g. the per-pop `empirical_AC_{POP}` family) are ignored: this model only declares the
    # fields `remap_divref` actually reads.
    model_config = ConfigDict(populate_by_name=True, extra="ignore")

    sequence_id: str
    sequence: str
    sequence_length: int
    n_variants: int
    popmax_fraction_phased: float
    popmax_empirical_af: float = Field(alias="popmax_empirical_AF")
    popmax_empirical_ac: int = Field(alias="popmax_empirical_AC")
    popmax_estimated_gnomad_af: float = Field(alias="popmax_estimated_gnomad_AF")
    max_pop: str
    variants: str
    source: str
    # Per-pop comma-delimited gnomAD AF strings, keyed by pop label. The dict's iteration order
    # is the order of `joint_pops_legend` used when constructing the index.
    gnomad_afs: dict[str, str]
    # Per-pop scalar estimated gnomAD AF for the whole haplotype, keyed by pop label. None for
    # pops with no source data on this row. Iteration order matches `joint_pops_legend`.
    estimated_gnomad_af_per_pop: dict[str, Optional[float]]

    _variants: Optional[list[Variant]] = None

    @classmethod
    def from_row(cls, row: dict[str, Any], pops_legend: list[str]) -> "Haplotype":
        """
        Build a Haplotype from a DuckDB `sequences` row, separating per-pop columns.

        Args:
            row: Mapping from DuckDB column name to value for one `sequences` row.
            pops_legend: Ordered list of population labels expected as `gnomAD_AF_{pop}` and
                `estimated_gnomad_AF_{pop}` columns. The resulting `gnomad_afs` and
                `estimated_gnomad_af_per_pop` dicts preserve this ordering.

        Returns:
            Haplotype instance with `gnomad_afs` and `estimated_gnomad_af_per_pop` populated
            from the per-pop columns.
        """
        # A NULL value in the DuckDB column (e.g. `gnomAD_AF_mid` on an HGDP_haplotype row
        # whose HGDP source legend doesn't include the joint pop) reads through polars as
        # `None`. Substitute a comma-delimited string of `NA` tokens matching the haplotype's
        # variant count so `_parse_pop_freqs` produces a list-of-zeros of the right length.
        n_variants = int(row["n_variants"])
        missing_freqs = ",".join(["NA"] * n_variants)
        gnomad_afs: dict[str, str] = {
            pop: (row[f"{_GNOMAD_AF_COLUMN_PREFIX}{pop}"] or missing_freqs)
            for pop in pops_legend
        }
        estimated_gnomad_af_per_pop: dict[str, Optional[float]] = {
            pop: row[f"{_ESTIMATED_GNOMAD_AF_COLUMN_PREFIX}{pop}"] for pop in pops_legend
        }
        base: dict[str, Any] = {
            k: v
            for k, v in row.items()
            if not k.startswith(_GNOMAD_AF_COLUMN_PREFIX)
            and not k.startswith(_ESTIMATED_GNOMAD_AF_COLUMN_PREFIX)
        }
        return cls(
            **base,
            gnomad_afs=gnomad_afs,
            estimated_gnomad_af_per_pop=estimated_gnomad_af_per_pop,
        )

    def parsed_variants(self) -> list[Variant]:
        """
        Parse the comma-delimited variants string into Variant objects.

        Returns:
            List of Variant objects parsed from the variants field.
        """
        if self._variants is not None:
            return self._variants
        vs = []
        for v_str in self.variants.split(","):
            chrom, pos, ref, alt = v_str.strip().split(":")
            vs.append(Variant(chromosome=chrom, position=int(pos), reference=ref, alternate=alt))
        self._variants = vs
        return vs

    def contig(self) -> str:
        """
        Return the chromosome of the first variant in this haplotype.

        Returns:
            Chromosome name (e.g. 'chr1').
        """
        return self.parsed_variants()[0].chromosome

    def reference_mapping(self, start: int, end: int, context_size: int) -> ReferenceMapping:
        """
        Map a [start, end) interval in haplotype sequence space to reference genome coordinates.

        Accounts for insertions and deletions when translating coordinates. For positions
        within a variant interval, snaps to the variant boundary (start for the left edge,
        end for the right edge). For positions in reference-only sequence, translates
        relative to the nearest preceding variant.

        Args:
            start: Start position (0-indexed, inclusive) in haplotype sequence space.
            end: End position (0-indexed, exclusive) in haplotype sequence space.
            context_size: Number of flanking reference bases prepended to the haplotype sequence.

        Returns:
            ReferenceMapping with the corresponding reference genome interval and variant metadata.
        """
        vs = self.parsed_variants()

        # Build [start, end) intervals in 0-indexed haplotype sequence space for each variant.
        # index_translation converts locus positions to string indices: locus - translation = index.
        variant_intervals: list[tuple[int, int]] = []
        index_translation = vs[0].position - context_size
        for v in vs:
            v_start = v.position - index_translation
            v_end = v_start + len(v.alternate)
            index_translation += len(v.reference) - len(v.alternate)
            variant_intervals.append((v_start, v_end))

        first_variant_index: Optional[int] = None
        last_variant_index: Optional[int] = None
        for i, (v_start, v_end) in enumerate(variant_intervals):
            if _intervals_overlap(start, end, v_start, v_end):
                if first_variant_index is None:
                    first_variant_index = i
                last_variant_index = i

        reference_coord_start = _translate_coordinate_to_ref(start, -1, vs, variant_intervals)
        reference_coord_end = _translate_coordinate_to_ref(end, 1, vs, variant_intervals)

        all_pop_freqs = {pop: _parse_pop_freqs(encoded) for pop, encoded in self.gnomad_afs.items()}

        if first_variant_index is not None and last_variant_index is not None:
            variants_involved = vs[first_variant_index : last_variant_index + 1]
        else:
            variants_involved = []

        return ReferenceMapping(
            chromosome=self.contig(),
            start=reference_coord_start,
            end=reference_coord_end,
            variants_involved=variants_involved,
            first_variant_index=first_variant_index,
            last_variant_index=last_variant_index,
            population_frequencies=all_pop_freqs,
        )


def _intervals_overlap(start1: int, end1: int, start2: int, end2: int) -> bool:
    return start1 < end2 and start2 < end1


# Missing-AF tokens seen in the TSVs exported by Hail/`create_duckdb_index`: "null" and "NA" come
# from Hail's missing representation; the empty string appears when a pop has no entry at all.
_MISSING_AF_TOKENS = frozenset({"null", "NA", ""})


def _parse_pop_freqs(encoded: str) -> list[float]:
    return [0.0 if v in _MISSING_AF_TOKENS else float(v) for v in encoded.split(",")]


def _translate_coordinate_to_ref(
    coord: int,
    sign: int,
    vs: list[Variant],
    variant_intervals: list[tuple[int, int]],
) -> int:
    """
    Translate a haplotype sequence coordinate back to a reference genome position.

    If the coordinate falls before the first variant, it is translated relative to
    that variant's reference position. If it falls within a variant interval, it snaps
    to the variant's reference start (sign < 0) or end (sign > 0). Otherwise it is
    translated relative to the end of the last preceding variant on the reference.

    Args:
        coord: 0-indexed coordinate in haplotype sequence space.
        sign: Negative to snap to variant start, positive to snap to variant end.
        vs: List of Variant objects in order.
        variant_intervals: Corresponding [start, end) intervals in haplotype sequence space.

    Returns:
        Reference genome position (1-based locus coordinate).
    """
    first_variant_start = variant_intervals[0][0]
    if (coord < first_variant_start and sign == -1) or (
        coord - 1 < first_variant_start and sign == 1
    ):
        return vs[0].position - (first_variant_start - coord)

    last_smaller_variant = 0
    for i, (v_start, v_end) in enumerate(variant_intervals):
        if v_start <= coord < v_end:
            if sign < 0:
                return vs[i].position
            else:
                return vs[i].position + len(vs[i].reference)
        if v_start > coord:
            break
        last_smaller_variant = i

    v = vs[last_smaller_variant]
    v_end_ref = v.position + len(v.reference)
    return v_end_ref + (coord - variant_intervals[last_smaller_variant][1])


def _get_index_connection(index_path: Optional[Path]) -> duckdb.DuckDBPyConnection:
    if index_path is None:
        for root, _dirs, files in os.walk(Path.cwd()):
            for file in files:
                if file.endswith(".duckdb"):
                    index_path = Path(root) / file
                    break
    if index_path is None:
        raise RuntimeError(
            "Unable to find a DuckDB index file. Pass --index-path or run from the "
            "same directory as the index file."
        )
    return duckdb.connect(str(index_path))


def remap_divref(  # noqa: C901
    *,
    input_path: Path,
    output_path: Path,
    index_path: Optional[Path] = None,
    separator: str = "\t",
    batch_size: int = 25000,
) -> None:
    """
    Remap DivRef haplotype coordinates to reference genome coordinates for CALITAS output.

    Reads a CALITAS output TSV, looks up each haplotype sequence in the DivRef DuckDB
    index, translates the haplotype-space coordinates back to reference genome positions,
    and writes an augmented TSV with reference coordinates and variant metadata appended.

    Args:
        input_path: Path to the CALITAS output file.
        output_path: Path to write the remapped output file.
        index_path: Path to the DivRef DuckDB index file. If not provided, the tool
            searches the directory containing this script for a .duckdb file.
        separator: Field delimiter used in both input and output files.
        batch_size: Number of rows to process per database query batch.
    """
    conn = _get_index_connection(index_path)

    df: pd.DataFrame = pd.read_csv(input_path, sep=separator)
    chrom_field: str = "chromosome"
    start_field: str = "coordinate_start"
    end_field: str = "coordinate_end"
    strand_field: str = "strand"
    padded_target_field: str = "padded_target"
    unpadded_target_field: str = "unpadded_target_sequence"

    required_fields: list[str] = [
        chrom_field,
        start_field,
        end_field,
        strand_field,
        padded_target_field,
        unpadded_target_field,
    ]
    if not all(x in df.columns for x in required_fields):
        raise ValueError(f"Required fields not found in input file: {', '.join(required_fields)}")

    if df[chrom_field].dtype != object:
        df[chrom_field] = df[chrom_field].astype(str)

    version_row = conn.execute("SELECT * FROM VERSION").fetchone()
    if version_row is None:
        raise RuntimeError("Index is missing VERSION table — ensure this is a valid DivRef index.")
    version: str = version_row[0]

    window_size_row = conn.execute("SELECT * FROM window_size").fetchone()
    if window_size_row is None:
        raise RuntimeError(
            "Index is missing window_size table — ensure this is a valid DivRef index."
        )
    window_size: int = window_size_row[0]

    joint_pops_legend_row = conn.execute("SELECT * FROM joint_pops_legend").fetchone()
    if joint_pops_legend_row is None:
        raise RuntimeError(
            "Index is missing joint_pops_legend table — ensure this is a valid DivRef index."
        )
    joint_pops_legend: list[str] = json.loads(joint_pops_legend_row[0])

    contigs: list[str] = []
    starts: list[int] = []
    ends: list[int] = []
    variants_involved: list[str] = []
    all_variants: list[str] = []
    n_variants_involved: list[int] = []
    popmax_empirical_af: list[float] = []
    popmax_empirical_ac: list[int] = []
    max_pop: list[str] = []
    source: list[str] = []
    # One column per pop in the joint legend. `gnomad_af_per_pop` holds the comma-delimited
    # per-variant AF strings (matches DuckDB's `gnomAD_AF_{POP}` columns verbatim);
    # `estimated_gnomad_af_per_pop` holds the per-pop scalar haplotype-level estimated AF.
    gnomad_af_per_pop: dict[str, list[str]] = {pop: [] for pop in joint_pops_legend}
    estimated_gnomad_af_per_pop: dict[str, list[Optional[float]]] = {
        pop: [] for pop in joint_pops_legend
    }

    for batch_start in tqdm(range(0, len(df), batch_size)):
        batch_end = min(batch_start + batch_size, len(df))
        batch_df = df.iloc[batch_start:batch_end]
        batch_hap_ids = batch_df[chrom_field].tolist()

        results = conn.execute(
            """
            SELECT * FROM sequences
            WHERE sequences.sequence_id IN (SELECT unnest($1::STRING[]))
            """,
            [batch_hap_ids],
        ).fetchall()

        columns = [desc[0] for desc in conn.description]
        id_to_hap: dict[str, Haplotype] = {}
        for row in results:
            hap = Haplotype.from_row(dict(zip(columns, row, strict=True)), joint_pops_legend)
            id_to_hap[hap.sequence_id] = hap

        for _, df_row in batch_df.iterrows():
            start: int = df_row[start_field]
            end: int = df_row[end_field]
            hap_id: str = df_row[chrom_field]
            strand: str = df_row[strand_field]
            padded_target: str = df_row[padded_target_field]
            target: str = df_row[unpadded_target_field]

            padded_len_adj = len(padded_target.replace("-", "")) - len(target)
            if strand == "+":
                end += padded_len_adj
            else:
                start -= padded_len_adj

            found_hap = id_to_hap.get(hap_id)
            if found_hap is None:
                raise RuntimeError(
                    f"Unable to find haplotype for {hap_id} — ensure you are aligning against "
                    f"the same DivRef version as this index (DivRef-v{version})"
                )
            rm = found_hap.reference_mapping(start, end, window_size)

            contigs.append(rm.chromosome)
            starts.append(rm.start)
            ends.append(rm.end)
            all_variants.append(found_hap.variants)
            variants_involved.append(rm.variants_involved_str())
            n_variants_involved.append(len(rm.variants_involved))
            popmax_empirical_af.append(found_hap.popmax_empirical_af)
            popmax_empirical_ac.append(found_hap.popmax_empirical_ac)
            max_pop.append(found_hap.max_pop)
            source.append(found_hap.source)
            for pop in joint_pops_legend:
                gnomad_af_per_pop[pop].append(found_hap.gnomad_afs[pop])
                estimated_gnomad_af_per_pop[pop].append(found_hap.estimated_gnomad_af_per_pop[pop])

    df["divref_sequence_id"] = df[chrom_field]
    df["divref_start"] = df[start_field]
    df["divref_end"] = df[end_field]
    df[chrom_field] = contigs
    df[start_field] = starts
    df[end_field] = ends
    df["genome_build"] = f"DivRef-v{version}"
    df["all_variants"] = all_variants
    df["variants_involved"] = variants_involved
    df["n_variants_involved"] = n_variants_involved
    df["popmax_empirical_AF"] = popmax_empirical_af
    df["popmax_empirical_AC"] = popmax_empirical_ac
    df["max_pop"] = max_pop
    df["variant_source"] = source
    for pop in joint_pops_legend:
        df[f"{_GNOMAD_AF_COLUMN_PREFIX}{pop}"] = gnomad_af_per_pop[pop]
        df[f"{_ESTIMATED_GNOMAD_AF_COLUMN_PREFIX}{pop}"] = estimated_gnomad_af_per_pop[pop]

    df.to_csv(output_path, sep=separator, index=False, quoting=csv.QUOTE_MINIMAL)
    logger.info("Wrote remapped output to %s", output_path)
