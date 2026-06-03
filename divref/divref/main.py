import logging
import sys
from collections.abc import Callable

import defopt

from divref.tools.append_contig_to_duckdb_index import append_contig_to_duckdb_index
from divref.tools.compute_haplotypes import compute_haplotypes
from divref.tools.create_divref_fasta import create_divref_fasta
from divref.tools.extract_gnomad_afs import extract_gnomad_afs
from divref.tools.extract_gnomad_single_afs import extract_gnomad_single_afs
from divref.tools.extract_sample_metadata import extract_sample_metadata
from divref.tools.finalize_duckdb_index import finalize_duckdb_index
from divref.tools.gnomad_hail_table_test_data import gnomad_hail_table_test_data
from divref.tools.init_duckdb_index import init_duckdb_index
from divref.tools.remap_divref import remap_divref

_tools: list[Callable[..., None]] = [
    append_contig_to_duckdb_index,
    compute_haplotypes,
    create_divref_fasta,
    extract_gnomad_afs,
    extract_gnomad_single_afs,
    extract_sample_metadata,
    finalize_duckdb_index,
    gnomad_hail_table_test_data,
    init_duckdb_index,
    remap_divref,
]


def setup_logging(level: str = "INFO") -> None:
    """
    Set up basic logging to print to the console.

    Args:
        level: Logging level string (e.g. "INFO", "DEBUG", "WARNING").
    """
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(name)s:%(funcName)s:%(lineno)s [%(levelname)s]: %(message)s",
    )


def run() -> None:
    """Set up logging, then hand over to defopt for running command line tools."""
    setup_logging()
    logger = logging.getLogger("divref")
    logger.info("Executing: " + " ".join(sys.argv))
    defopt.run(
        funcs=_tools,
        argv=sys.argv[1:],
    )
    logger.info("Finished executing successfully.")
