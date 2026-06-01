"""Build the masked chr1+chrX GRCh38 test reference committed under divref/tests/data.

Hail's `ReferenceGenome.add_sequence` requires the FASTA to contain an entry for every
GRCh38 contig at its full declared length (it rejects a chr1+chrX-only FASTA with
"Contigs missing in FASTA ... present in reference genome 'GRCh38'"). So this writes all
contigs from the full reference's `.fai`, filled with `N` except the chr1 and chrX test
windows, which carry real sequence. bgzipped, the all-`N` contigs compress to ~12 MB, small
enough to commit (the `.fa.gz` is tracked via git-lfs).

The real-sequence windows cover the test loci (chr1:100001-200000, chrX:50000000-50025000)
plus a 5 kb margin, which is ample for the flanking reference context that
`compute-haplotypes` / `append-contig-to-duckdb-index` extract around each variant.

Requires the full GRCh38 reference (NOT committed) at
`data/work/inputs/Homo_sapiens_assembly38.fasta` (+ `.fai`).

Regenerate with:

    pixi run bash -c "python scripts/build_test_reference.py | bgzip \
        > divref/tests/data/test_reference.chr1_chrX.fa.gz"
    pixi run samtools faidx divref/tests/data/test_reference.chr1_chrX.fa.gz
    # Hail and append-contig-to-duckdb-index look for the index at the with_suffix('.fai') path:
    cp divref/tests/data/test_reference.chr1_chrX.fa.gz.fai \
       divref/tests/data/test_reference.chr1_chrX.fa.fai
    rm divref/tests/data/test_reference.chr1_chrX.fa.gz.fai
"""

import subprocess
import sys

FULL_FASTA = "data/work/inputs/Homo_sapiens_assembly38.fasta"
FULL_FAI = "data/work/inputs/Homo_sapiens_assembly38.fai"
# Real-sequence windows (1-based, inclusive): test loci + 5 kb margin.
WINDOWS: dict[str, tuple[int, int]] = {
    "chr1": (95001, 205000),
    "chrX": (49995001, 50030000),
}
LINE = 60


def fetch(contig: str, start: int, end: int) -> str:
    """Return the reference bases for `contig:start-end` (1-based inclusive) via samtools."""
    result = subprocess.run(
        ["samtools", "faidx", FULL_FASTA, f"{contig}:{start}-{end}"],
        capture_output=True,
        text=True,
        check=True,
    )
    return "".join(result.stdout.splitlines()[1:])


def main() -> None:
    """Write the masked all-contig FASTA to stdout (pipe into bgzip)."""
    real = {contig: (s, e, fetch(contig, s, e)) for contig, (s, e) in WINDOWS.items()}
    for contig, (s, e, seq) in real.items():
        if len(seq) != e - s + 1:
            raise ValueError(f"{contig}: fetched {len(seq)} bp, expected {e - s + 1}")

    write = sys.stdout.write
    with open(FULL_FAI) as fai:
        for line in fai:
            name, length_str = line.split("\t")[:2]
            length = int(length_str)
            write(f">{name}\n")
            if name in real:
                s, e, seq_real = real[name]
                seq = "N" * (s - 1) + seq_real + "N" * (length - e)
                write("\n".join(seq[i : i + LINE] for i in range(0, length, LINE)))
                write("\n")
            else:
                full_lines, remainder = divmod(length, LINE)
                write(("N" * LINE + "\n") * full_lines)
                if remainder:
                    write("N" * remainder + "\n")


if __name__ == "__main__":
    main()
