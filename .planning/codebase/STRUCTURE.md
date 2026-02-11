# Codebase Structure

**Analysis Date:** 2026-02-11

## Directory Layout

```
viral-ngs/
├── src/viral_ngs/                   # Main source code
│   ├── __init__.py                  # Package init; exports core module
│   ├── core/                        # Core tool abstractions and utilities
│   ├── assemble/                    # Genome assembly tool wrappers
│   ├── classify/                    # Metagenomic classification tool wrappers
│   ├── phylo/                       # Phylogenetic analysis tool wrappers
│   ├── read_utils.py                # Read processing utilities and CLI
│   ├── assembly.py                  # Assembly workflow and CLI
│   ├── metagenomics.py              # Classification workflow and CLI
│   ├── interhost.py                 # Interhost phylogenetics CLI
│   ├── intrahost.py                 # Intrahost variant calling CLI
│   ├── taxon_filter.py              # Taxonomic filtering CLI
│   ├── reports.py                   # Report generation CLI
│   ├── illumina.py                  # Illumina-specific utilities CLI
│   ├── kmer_utils.py                # K-mer utilities CLI
│   ├── ncbi.py                      # NCBI utilities CLI
│   ├── file_utils.py                # File handling utilities CLI
│   └── broad_utils.py               # Broad-specific utilities CLI
├── tests/
│   ├── __init__.py                  # Test utilities and base classes
│   ├── conftest.py                  # pytest configuration and fixtures
│   ├── unit/                        # Unit tests
│   │   ├── core/                    # Core module tests
│   │   ├── assemble/                # Assembly module tests
│   │   ├── classify/                # Classification module tests
│   │   └── phylo/                   # Phylogenetics module tests
│   └── input/                       # Test input data
│       ├── TestDepleteHuman/        # Human read depletion test data
│       ├── TestMetagenomicsSimple/  # Metagenomics test data
│       ├── s3/                      # Mock S3 test data
│       └── ...                      # Other test datasets
├── docker/                          # Docker build configuration
├── docs/                            # Documentation
├── pyproject.toml                   # Package metadata and build config
├── AGENTS.md                        # Comprehensive development guidelines
├── CLAUDE.md                        # Claude Code quick reference
├── README.md                        # Project overview
└── .github/                         # GitHub workflows and CI/CD
```

## Directory Purposes

**`src/viral_ngs/core/`:**
- Purpose: Tool abstractions, base classes, and cross-cutting utilities
- Contains: Tool base class, tool wrappers for external executables, file/misc utilities
- Key files:
  - `__init__.py`: Tool base class, InstallMethod, submodule imports
  - `cmd.py`: CLI framework (argument parsing, logging, subcommand registration)
  - `file.py`: File I/O, compression, FASTA/BAM handling, temp files
  - `misc.py`: Threading, string processing, sequence utilities
  - `samtools.py`, `bwa.py`, `picard.py`, etc.: Individual tool wrappers (170+ lines each)

**`src/viral_ngs/assemble/`:**
- Purpose: Tool wrappers specific to genome assembly workflows
- Contains: SPAdes, Gap2Seq, MUMmer, MAFFT, Muscle, Skani, WGSim, VCF utilities
- Naming: Each tool has its own module: `spades.py`, `mummer.py`, etc.

**`src/viral_ngs/classify/`:**
- Purpose: Tool wrappers specific to metagenomic classification
- Contains: BLAST, BMTagger, Kraken2, Krona, Taxonomy, KMA, KB, LAST, KMC
- Naming: One module per classifier

**`src/viral_ngs/phylo/`:**
- Purpose: Tool wrappers and utilities for phylogenetic analysis
- Contains: snpEff, VPhaser2, GenBank utilities, feature table handling, VCF tools
- Large files: `feature_table_types.py` (1297 lines), `feature_table.py` (335 lines)

**`tests/unit/`:**
- Purpose: Unit tests organized to mirror source structure
- Organization: `core/`, `assemble/`, `classify/`, `phylo/` subdirectories
- Naming: `test_<module>.py` (e.g., `test_tools_bwa.py`, `test_assembly.py`)
- Base class: `TestCaseWithTmp` from `tests/__init__.py` provides tmpdir fixture

**`tests/input/`:**
- Purpose: Test data for unit and integration tests
- Organization: One directory per test class (e.g., `TestDepleteHuman/`, `TestMetagenomicsSimple/`)
- Purpose of subdirs:
  - `s3/`: Mock S3 bucket structure with test databases
  - `db/`: Taxonomic and sequence databases for testing
  - `expected/`, `aligned-expected/`: Expected output for validation tests

**`docker/`:**
- Purpose: Container definitions for Docker-centric development
- Structure:
  - `Dockerfile.assemble`: Build image with assembly tools
  - `Dockerfile.core`: Build image with core bioinformatics tools
  - `requirements/`: Conda environment files
    - `core.txt`: Core dependencies (samtools, picard, gatk, etc.)
    - `tests.txt`: Testing dependencies (pytest, pytest-xdist)

## Key File Locations

**Entry Points:**
- `src/viral_ngs/__init__.py`: Package exports `viral_ngs.core` module; defines `__version__`
- `src/viral_ngs/read_utils.py`: CLI entry point `read_utils` (defined in pyproject.toml)
- `src/viral_ngs/assembly.py`: CLI entry point `assembly`
- `src/viral_ngs/metagenomics.py`: CLI entry point `metagenomics`
- `src/viral_ngs/interhost.py`: CLI entry point `interhost`
- `src/viral_ngs/intrahost.py`: CLI entry point `intrahost`

**Configuration:**
- `pyproject.toml`: Package metadata, dependencies, entry points, pytest config
- `docker/requirements/core.txt`: Conda dependencies for all tools
- `pyproject.toml` [project.scripts]: Maps CLI command names to module:main functions

**Core Logic:**
- `src/viral_ngs/core/__init__.py`: Tool base class (Tool, InstallMethod, PrexistingUnixCommand)
- `src/viral_ngs/core/file.py`: File operations, FASTA/BAM/TSV handling, compression (1257 lines)
- `src/viral_ngs/core/misc.py`: Utilities, threading, string/seq helpers (1361 lines)
- `src/viral_ngs/core/cmd.py`: CLI framework, argument parsing, logging setup (321 lines)

**Testing:**
- `tests/__init__.py`: TestCaseWithTmp base class, assert helpers, test fixtures
- `tests/conftest.py`: pytest configuration, tmpdir scope fixtures
- `tests/unit/core/test_tools_*.py`: Tool wrapper tests (e.g., `test_tools_bwa.py`)
- `tests/unit/assemble/test_assembly.py`: Assembly workflow tests

## Naming Conventions

**Files:**
- Tool wrappers: `<tool_name>.py` (lowercase, single word or hyphenated) in appropriate subpackage
  - Example: `src/viral_ngs/core/bwa.py`, `src/viral_ngs/assemble/spades.py`
- Test files: `test_<module>.py` or `test_<tool>.py`
  - Example: `tests/unit/core/test_tools_bwa.py`, `tests/unit/assemble/test_assembly.py`
- CLI modules: `<feature>.py` (readable name) in `src/viral_ngs/`
  - Example: `read_utils.py`, `metagenomics.py`, `interhost.py`

**Directories:**
- Package: lowercase, matches module name
  - Example: `src/viral_ngs/core/`, `src/viral_ngs/assemble/`
- Test data: PascalCase, matches test class name
  - Example: `tests/input/TestDepleteHuman/`, `tests/input/TestMetagenomicsSimple/`

**Python Naming (within files):**
- Tool classes: PascalCase with "Tool" suffix
  - Example: `SamtoolsTool`, `Bwa`, `PicardTool`, `SpadesTool`
- Functions: snake_case
  - Example: `assemble_spades()`, `trim_rmdup_subsamp_reads()`
- Module constants: UPPERCASE
  - Example: `TOOL_NAME = 'samtools'`
- Private/internal functions: Leading underscore
  - Example: `_get_tool_version()`, `_attempt_install()`

## Where to Add New Code

**New Tool Wrapper (e.g., new bioinformatics tool):**
- Create file: `src/viral_ngs/<category>/<toolname>.py` where `<category>` is `core`, `assemble`, `classify`, or `phylo`
- Pattern:
  ```python
  # At top:
  from . import Tool, PrexistingUnixCommand
  import shutil

  TOOL_NAME = 'mytool'

  class MyToolTool(Tool):
      def __init__(self, install_methods=None):
          if install_methods is None:
              install_methods = [PrexistingUnixCommand(shutil.which(TOOL_NAME))]
          super().__init__(install_methods=install_methods)

      def _get_tool_version(self):
          # Implement version detection
          pass

      def my_method(self, arg1, arg2):
          # Implement tool-specific methods
          pass
  ```
- Test file: `tests/unit/<category>/test_tools_<toolname>.py`
- Pattern: Extend `TestCaseWithTmp`, use `self.input()` helper for test data

**New CLI Subcommand:**
- Add to existing module (e.g., `read_utils.py`) or create new module if workflow-specific
- Pattern in `src/viral_ngs/<module>.py`:
  ```python
  __commands__ = []

  def parser_mycommand(parser=argparse.ArgumentParser()):
      parser.add_argument('input_file', help='...')
      parser.add_argument('--option', help='...')
      util_cmd.common_args(parser, (('loglevel', None),))
      util_cmd.attach_main(parser, main_mycommand)
      return parser

  def main_mycommand(input_file, option=None):
      """Short description of command."""
      # Implementation
      return 0

  __commands__.append(('mycommand', parser_mycommand))
  ```
- Add entry point to `pyproject.toml` [project.scripts] if new top-level command

**New Domain Module (e.g., new workflow):**
- Create `src/viral_ngs/<workflow>.py` or expand existing domain subpackage
- If large: Consider creating `src/viral_ngs/<workflow>/` subpackage with submodules
- Must include `__commands__` registry for CLI integration
- Follow existing patterns in `metagenomics.py` or `assembly.py`

**Utilities:**
- Shared helper functions: Add to `src/viral_ngs/core/misc.py` or `src/viral_ngs/core/file.py`
- Domain-specific helpers: Keep in relevant domain module file
- Example: `trim_rmdup_subsamp_reads()` in `read_utils.py` (read-processing specific)

**Test Fixtures:**
- Parametrized/reusable test helpers: Add to `tests/__init__.py` as functions or base classes
- Module-specific fixtures: Add to test module file as pytest fixtures
- Test data: Place in `tests/input/<TestClassName>/` directory, matching test class name

## Special Directories

**`tests/input/s3/`:**
- Purpose: Mock S3 bucket structure for testing database downloads
- Generated: No, committed to repo
- Committed: Yes, contains test-sized database files
- Structure: `s3/sabeti-public-dbs/` mimics AWS S3 bucket layout

**`src/viral_ngs/core/__pycache__/`, `*.egg-info/`:**
- Purpose: Python cache and package installation metadata
- Generated: Yes, automatically during build/test
- Committed: No (in .gitignore)

**`docker/requirements/`:**
- Purpose: Conda environment specifications for Docker builds
- Generated: No, manually maintained
- Committed: Yes
- Files:
  - `core.txt`: Main bioinformatics tools and dependencies
  - `tests.txt`: Testing framework and plugins

---

*Structure analysis: 2026-02-11*
