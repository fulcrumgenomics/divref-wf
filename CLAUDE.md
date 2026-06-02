# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**divref-wf** creates a DivRef-style resource bundle (FASTA sequences and DuckDB index) for common human variation, based on the [DivRef project](https://github.com/e9genomics/human-diversity-reference). It wraps original Python scripts with improved typing, parameterization, and unit tests, and adds Snakemake workflow orchestration.

## Environment & Commands

This project uses **two package managers**:

- **`uv`** — manages the `divref` Python package under `./divref/`
- **`pixi`** — manages the full workspace (Python toolkit + Snakemake + Hail)

### Python Toolkit (`./divref/`)

```bash
uv run --directory divref poe fix-and-check-all   # Fix then check everything (required before commit)
uv run --directory divref poe check-all            # Check format, lint, tests, and types
uv run --directory divref poe fix-all              # Auto-fix formatting and linting
uv run --directory divref pytest                   # Run tests
uv run --directory divref pytest tests/test_main.py::test_name  # Run a single test
uv run --directory divref mypy divref/             # Type-check only
```

### Workspace (Pixi)

```bash
pixi run fix-and-check-all   # Fix and check toolkit + Snakemake linting
pixi run lint --check        # Validate Snakemake files with snakefmt
pixi run download            # Run the Snakemake download workflow
pixi run setup-gcs           # Download GCS connector JAR (required once for Hail on GCS)
```

## Architecture

### Repository Layout

```
divref/                    # Python package (uv-managed)
  divref/
    main.py                # CLI entry point; registers tools in _tools list
    alias.py               # HailPath type alias (str; accepts local, gs://, hdfs://)
    defaults.py            # Package-wide constants: POPULATIONS, REFERENCE_GENOME, gnomAD HT URIs
    hail.py                # Hail initialization with GCS connector setup
    haplotype.py           # Shared Hail utilities for haplotype sequence/windowing
    haplotype_compat.py    # Pure variant-overlap classification for the haplotype_filter flag
    duckdb_index.py        # Shared helpers for the init/append/finalize DuckDB-index tools
    tools/                 # One module per CLI subcommand
  tests/                   # pytest tests
  pyproject.toml           # Package deps, ruff/mypy/pytest config
workflows/                 # Snakemake workflows
  generate_divref.smk      # Main workflow (extract → haplotypes → reference download)
  create_test_data.smk     # Generates gnomAD subset for unit tests
  config/config.yml        # Workflow configuration (chromosomes, populations, paths)
pixi.toml                  # Workspace config (snakemake + hail environments)
```

### CLI Pattern

Tools are plain functions registered in `main.py`; **defopt** auto-generates the CLI from their docstrings. To add a new tool:

1. Create `divref/tools/<name>.py` with a keyword-only function and Google docstring
2. Import and add it to `_tools` in `main.py`

```bash
divref <tool-name> --arg value   # Invokes the registered tool
```

### Tool Pipeline (execution order)

The tools implement a data pipeline:
1. `extract_gnomad_afs` / `extract_gnomad_single_afs` → per-population allele frequency Hail table
2. `extract_sample_metadata` → simplified sample→population mapping table
3. `compute_haplotypes` → groups phased variants into haplotype windows using Hail
4. `init_duckdb_index` → create the DuckDB and write the population-legend + version metadata; then `append_contig_to_duckdb_index` (once per chromosome, each in a fresh JVM) → append that contig's merged haplotype + gnomAD-sites rows with reference-context sequences; then `finalize_duckdb_index` → build the `sequence_id` index. The append step also computes the `haplotype_filter` column (VCF-style: `PASS`, else the `;`-joined sorted flags from `divref/divref/haplotype_compat.py`) flagging haplotypes whose component variants overlap (a phasing artifact, not dropped). `compatibility_flag` classifies every variant pair (not just adjacent ones) into a reason from the `REASONS` taxonomy, then adds the `end_extends_past_rightmost_variant` flag when an earlier long-reference variant reaches past the rightmost one. A co-located SNP+indel is the one genuinely compatible overlap and yields no flag. Correspondingly `haplo_coordinates` sets `end` from the maximum reference end over all variants so a deletion's span is never truncated. `haplotype_filter` flows through `remap_divref` into the final CALITAS TSV.
5. `create_divref_fasta` → per-chromosome FASTA files streamed from the DuckDB index (final deliverable)
6. `remap_divref` → maps haplotype coordinates back to reference genome (post-CALITAS step)

Steps 4 and 5 were split out of a single `create_fasta_and_index` tool (PR #39). The index builder (originally `create_duckdb_index`, with chunked writes via `--polars-chunk-size` from PR #42) was later split into per-chromosome `init_duckdb_index` / `append_contig_to_duckdb_index` / `finalize_duckdb_index` so each contig runs in a fresh JVM, bounding file-descriptor use (a single long-lived process exhausted the per-process limit on a whole-genome run); their shared helpers live in `divref/divref/duckdb_index.py`. `create_divref_fasta` streams FASTA output (PR #43).

`extract_gnomad_single_afs` is an alternative to `extract_gnomad_afs` supporting both gnomAD v4.1 (JOINT) and v3.1.2 (HGDP+1KG) table schemas; it is used when the workflow's `gnomad_variant_annotation_source` config selects a gnomAD source different from the haplotype source (the haplotypes themselves always come from gnomAD 3.1.2 HGDP+1KG phased genotypes).

`gnomad_hail_table_test_data` is a separate tool (registered in `main.py`, but not part of the pipeline) used by `workflows/create_test_data.smk` to generate test-data subsets of gnomAD Hail tables.

### Key Shared Modules

**`haplotype.py`** (Hail expressions)
- `get_haplo_sequence(context_size, variants)` — builds haplotype sequence strings with flanking reference context. Composes variants left to right onto the reference with a running cursor (not allele concatenation), so any overlap resolves to a defined sequence; the window spans through the maximum reference end across all variants (`_max_reference_end`), matching `haplo_coordinates`.
- `haplo_coordinates(window_size, variants)` — start/end window; `end` uses `_max_reference_end`, not the last-by-position variant.
- `variant_distance(v1, v2)` — reference bases between two variants (accounts for indel length)

**`haplotype_compat.py`** (pure functions over `(contig, pos, ref, alt)` tuples, no Hail)
- `compatibility_flag(variants_str)` — the VCF-style `haplotype_filter` value for a haplotype's `variants` string: `PASS` or `;`-joined sorted flags.
- `classify_haplotype` / `classify_pair` — map overlapping variant pairs to a reason in the `REASONS` taxonomy (e.g. `snp_in_deletion`, `overlapping_deletions`, `insertion_in_deletion`); a co-located SNP+indel composes and returns `None`.
- `end_extends_past_rightmost_variant` — True when an earlier variant's reference reaches past the rightmost one (drives the `end_extends_past_rightmost_variant` flag and always co-occurs with an overlap).

**`hail.py`**: `hail_init(gcs_credentials_path)` — sets `GOOGLE_APPLICATION_CREDENTIALS`, verifies GCS connector JAR (installed via `pixi run setup-gcs`), then calls `hl.init()` with Spark GCS config.

**`defaults.py`**: `POPULATIONS`, `REFERENCE_GENOME`, and gnomAD HGDP+1KG Hail-table URI defaults shared across tools.

### Data Models (`remap_divref.py`)

Pydantic `frozen=True` models: `Variant`, `ReferenceMapping`, `Haplotype` — used for type-safe coordinate remapping. `Haplotype` uses field aliases to match mixedCase column names in the DuckDB index produced by `append_contig_to_duckdb_index`.

### Snakemake Workflows

- `workflows/generate_divref.smk` — main workflow. Reads `workflows/config/config.yml` (validated against `config_schema.yml`). Per-chromosome rules feed into a single `create_divref_index` rule (one DuckDB, built by looping `init_duckdb_index` → per-contig `append_contig_to_duckdb_index` → `finalize_duckdb_index`) and a `create_divref_fasta` rule (per-chromosome FASTAs).
- `workflows/create_test_data.smk` — generates the gnomAD Hail-table subsets committed under `divref/tests/data/` (used by `pytest`). Run when test inputs need refreshing, not on every test run.
- `workflows/compare_divref_gnomad.smk` — analysis-only workflow comparing DivRef 1.1 against multiple gnomAD releases (requires the `analysis` pixi environment).

Workflow knobs worth knowing about (in `config.yml`):
- `gnomad_variant_annotation_source` — selects which gnomAD release the single-variant track is drawn from (drives `extract_gnomad_single_afs`); haplotype track always comes from HGDP+1KG.
- `polars_chunk_size` — chunk size for the streaming DuckDB index writer.
- `sequence_window_size` — flanking reference context size around each haplotype/variant in the FASTA output.

### Analysis Scripts (`scripts/`)

Ad-hoc analysis and plotting that support the [blog post](docs/blog.md) and its [reproduction notes](docs/blog_reproduction.md); see `scripts/README.md`. These are not production code (mostly Claude-generated) and are not part of the deliverable pipeline. Run them via `pixi run python` / `pixi run Rscript`. They read existing DuckDB indexes and emit machine-readable summary TSVs (downstream plots and the blog read counts from those TSVs, never hardcoded values).

- `compare_haplotypes.py` — compares the old (DivRef 1.1) vs new `compute_haplotypes` haplotype catalogs from two DuckDB indexes: partitions into shared / old-only / new-only, checks shared-haplotype sequence equality (partitioned by overlap, not the PASS flag), and compares AC. Pure Python + DuckDB, no Spark. Emits `algo_comparison.summary.tsv`, consumed by `compare_haplotypes_venn.R`.
- `compare_divref_gnomad.R` — plotting for the `compare_divref_gnomad.smk` analysis workflow.
- `build_test_reference.py`, `explain_divref_only_v41.py`, `inspect_threshold_edge_haplotypes.py` — one-off investigations.

## Git Workflow

### Commit Granularity

Commit after completing one of:
- A single function/method implementation
- One refactoring step (rename, extract, move)
- A bug fix with its regression test
- A documentation update

**Size guidelines:**
- Per commit: 100–300 lines preferred, 400 max
- Per PR: No hard limit, but consider splitting if >800 lines or >5 unrelated files

**Good commit scope examples:**
- `Add FastaIndex.validate() method`
- `Rename species_map → species_to_ref_fasta_map`
- `Fix off-by-one in BED coordinate parsing`

### Commit Messages

Use [Conventional Commits](https://www.conventionalcommits.org/) for commit messages and PR
titles+bodies. Common types: `feat`, `fix`, `chore`, `docs`, `refactor`, `test`.

```
<type>: <imperative description> (<72 chars total)

Detailed body explaining:
- What changed
- Why (link issues with "Closes #123" or "Related to #456")
- Any non-obvious implementation choices
```

### Commit Rules
- Run `uv run poe fix-and-check-all` before each commit; all checks must pass
- No merge commits
- Do not rebase without explicit user approval
- **Never mix formatting and functional changes.** If unavoidable, isolate formatting into separate commits at start or end of branch.

### Branch Hygiene
- Use `.gitignore` liberally
- Never commit: IDE files, personal test files, local debug data, commented-out code

## Coding Conventions

### Organization
- Extract logic into small–medium functions with clear inputs/outputs
- Scope variables tightly; limit visibility to where needed
- Use block comments for visual separation when function extraction isn't practical

### Naming
- Meaningful names, even if long: `species_to_ref_fasta_map` not `species_map`
- Short names only for tight scope (loop indices, single-line lambdas)
- Signal behavior in function names: `to_y()`, `is_valid()` → returns value; `update_x()` → side effect

## Testing

### Principles
- Generate test data programmatically; avoid committing test data files
- Test behavior, not implementation—tests should survive refactoring
- Cover: expected behavior, error conditions, boundary cases
- Scale rigor to code longevity: thorough for shared code, lighter for one-off scripts

### Coverage Expectations
- New public functions: at least one happy-path test + one error case
- Bug fixes: add a regression test that would have caught the bug
- Performance-critical code: include benchmark or explain in PR why not needed

## Documentation Maintenance

When modifying code, update as needed:
- [ ] Docstrings (if signature or behavior changed)
- [ ] README.md (if usage patterns changed)
- [ ] Migration notes (if breaking change)

## Python-Specific

### Pragmatism
- Balance functional, OOP, and imperative—use what's clearest
- When in doubt, prefer pure functions and immutable data
- Know your utility libraries; contribute upstream rather than writing one-offs

### Style
- Heavier use of classes and type annotations than typical Python
- Prefer `@dataclass(frozen=True)` and Pydantic models with `frozen=True`

### Functions
- Functions should have **either** returns **or** side effects, not both
- Exceptions: logging, caching (where side effect is performance-only)

### Documentation
- Google-style docstrings with `Args:`, `Returns:`, `Yields:`, and `Raises:` blocks
- Docstrings are required on all public functions/classes
- Code comments should explain non-obvious choices and complex logic

### Typing
- **Required:** Type annotations on all function parameters and returns
- **Parameters:** Accept the most general type practical (e.g., `Iterable` over `List`)
- **Returns:** Return the most specific type without exposing implementation details
- Annotate locals when: they become return values, or called function lacks hints
- Use type aliases or `NewType` for complex structures
- Avoid `Any`—prefer type alias or `TypeVar`
- Avoid `cast()` and `type: ignore`—prefer alternatives, but when unavoidable (e.g., incorrect upstream stubs), document the reason inline.

