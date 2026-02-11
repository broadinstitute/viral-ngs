# Architecture

**Analysis Date:** 2026-02-11

## Pattern Overview

**Overall:** Modular tool-wrapper architecture with layered design

**Key Characteristics:**
- Tool-oriented design: Each bioinformatics tool is wrapped in its own class (e.g., `SamtoolsTool`, `Bwa`, `SpadesTool`)
- Domain-driven organization: Code organized by biological workflow stages (assembly, classification, phylogenetics)
- Dependency injection for tools: Install methods are passed to Tool classes, enabling testability
- Command-line interface (CLI) pattern: Each top-level module exposes multiple subcommands via argparse

## Layers

**Core Tool Abstraction Layer:**
- Purpose: Base classes and installation mechanisms for external bioinformatics tools
- Location: `src/viral_ngs/core/`
- Contains: Tool base class, InstallMethod classes, generic tool wrappers (samtools, bwa, picard, gatk, etc.)
- Depends on: External executables (via PATH), subprocess, pysam
- Used by: All domain modules and high-level CLI modules

**Utility/Infrastructure Layer:**
- Purpose: Cross-cutting functionality for file operations, logging, threading, validation
- Location: `src/viral_ngs/core/file.py`, `src/viral_ngs/core/misc.py`, `src/viral_ngs/core/cmd.py`
- Contains:
  - `file.py` (1257 lines): FASTA/BAM/TSV handling, compression support, temp file management, sequence I/O
  - `misc.py` (1361 lines): Threading utilities, string processing, sequence utilities, YAML/JSON helpers
  - `cmd.py` (321 lines): Argument parsing framework, logging setup, CLI infrastructure
- Depends on: BioPython, pysam, subprocess
- Used by: All tool wrappers and domain modules

**Domain-Specific Tool Wrapper Layer:**
- Purpose: Group related tool wrappers by biological workflow
- Locations:
  - `src/viral_ngs/assemble/`: Genome assembly tools (SPAdes, Gap2Seq, MUMmer, MAFFT, Muscle, Skani, WGSim)
  - `src/viral_ngs/classify/`: Metagenomic classification tools (BLAST, BMTagger, Kraken2, Krona, Taxonomy, KMA, KB, LAST, KMC)
  - `src/viral_ngs/phylo/`: Phylogenetic analysis tools (snpEff, VPhaser2, GenBank utilities, feature table handling)
- Depends on: Core Tool classes, utility layer
- Used by: High-level CLI modules and user-facing functions

**Business Logic Layer:**
- Purpose: High-level workflow orchestration and read processing
- Locations:
  - `src/viral_ngs/read_utils.py`: Read validation, format conversion, duplicate removal, trimming, filtering
  - `src/viral_ngs/assembly.py`: De novo assembly workflows (SPAdes integration, contig filtering, refinement)
  - `src/viral_ngs/metagenomics.py`: Metagenomic classification workflows and report generation
  - `src/viral_ngs/taxon_filter.py`: Taxonomic filtering of reads
  - `src/viral_ngs/reports.py`: Report generation and visualization
  - `src/viral_ngs/illumina.py`: Illumina sequencing-specific utilities
- Depends on: Tool wrappers, utility layer, external libraries (BioPython, numpy, pandas, matplotlib)
- Used by: CLI entry points

**CLI Entry Point Layer:**
- Purpose: Command-line interface and argument parsing
- Location: Top-level modules with `__commands__` registry:
  - `src/viral_ngs/read_utils.py`: read_utils CLI
  - `src/viral_ngs/assembly.py`: assembly CLI
  - `src/viral_ngs/metagenomics.py`: metagenomics CLI
  - `src/viral_ngs/interhost.py`: interhost phylogenetics CLI
  - `src/viral_ngs/intrahost.py`: intrahost variant calling CLI
  - `src/viral_ngs/taxon_filter.py`: taxon_filter CLI
  - `src/viral_ngs/reports.py`: reports CLI
  - `src/viral_ngs/file_utils.py`: file_utils CLI
  - `src/viral_ngs/illumina.py`: illumina CLI
  - `src/viral_ngs/kmer_utils.py`: kmer_utils CLI
  - `src/viral_ngs/ncbi.py`: ncbi CLI
  - `src/viral_ngs/broad_utils.py`: broad_utils CLI
- Pattern: Each module contains `__commands__` list of (command_name, parser_function) tuples
- Depends on: Business logic layer, core utilities
- Entry points defined in: `pyproject.toml` [project.scripts]

## Data Flow

**Typical Read Processing Workflow:**

1. **Input**: Sequencing reads (FASTQ or BAM)
2. **Read QC & Processing** (`read_utils.py`):
   - Quality trimming (Trimmomatic)
   - Duplicate removal (Picard MarkDuplicates)
   - Format conversion (samtools, picard)
   - Subsampling for assembly
3. **Classification** (`metagenomics.py`):
   - K-mer based methods (Kraken2, KMC)
   - Sequence search (BLAST, BMTagger)
   - Taxonomy assignment
   - Report generation
4. **Assembly** (`assembly.py`):
   - De novo assembly (SPAdes)
   - Quality filtering
   - Reference-guided refinement (Gap2Seq)
5. **Phylogenetic Analysis** (`interhost.py`, `intrahost.py`):
   - Multiple sequence alignment (MAFFT, Muscle, MUMmer)
   - Variant calling (VPhaser2)
   - Variant annotation (snpEff)
   - VCF processing
6. **Output**: Assemblies, variant calls, taxonomic reports

**State Management:**
- Stateless functional approach: Each function takes input files and produces output files
- Temporary file management: Uses `viral_ngs.core.file.mkstempfname()` and context managers
- Logging: Module-level loggers with configurable levels (DEBUG, INFO, WARNING, ERROR)
- Error handling: Custom exceptions (QCError, DenovoAssemblyError, InvalidBamHeaderError) for workflow failures

## Key Abstractions

**Tool Base Class:**
- Purpose: Encapsulates external tool execution with install machinery
- Location: `src/viral_ngs/core/__init__.py`
- Pattern: Subclasses override `_get_tool_version()` and implement tool-specific methods
- Examples: `SamtoolsTool`, `Bwa`, `PicardTool`, `GatkTool`
- Usage: `tool = SamtoolsTool(); tool.install(); tool.view(args, inFile, outFile)`

**InstallMethod Interface:**
- Purpose: Pluggable installation strategies
- Location: `src/viral_ngs/core/__init__.py`
- Pattern: Base class `InstallMethod` with method `is_installed()` and `executable_path()`
- Implementations: `PrexistingUnixCommand` (search PATH for executable)
- Enables: Testing without tool installation, multiple install strategies

**Temporary File Context Managers:**
- Purpose: Safely manage temporary files with cleanup
- Location: `src/viral_ngs/core/file.py`
- Pattern: `@contextlib.contextmanager` decorated functions
- Examples: `tmp_dir()`, `tempfnames()`, `bam2fq_tmp()`, `execute_tmp()`

**CLI Subcommand Registration:**
- Purpose: Flexible CLI framework without hard-coded command list
- Location: Module-level `__commands__` list in each CLI module
- Pattern: Each command is `(name, parser_function)` tuple
- Parser functions: Take optional parser, configure arguments, call `util_cmd.attach_main()`, return parser
- Main functions: Wrapped with `@util_cmd.main_command` to map Namespace args to function params

## Entry Points

**Assembly CLI:**
- Location: `src/viral_ngs/assembly.py`
- Triggers: Command `assembly <subcommand>`
- Responsibilities:
  - De novo assembly orchestration (assemble_spades)
  - Reference-guided refinement (gap2seq_assemble, order_and_orient)
  - Assembly quality control and filtering
  - Consensus sequence generation

**Metagenomics CLI:**
- Location: `src/viral_ngs/metagenomics.py`
- Triggers: Command `metagenomics <subcommand>`
- Responsibilities:
  - Taxonomic classification of reads
  - K-mer and sequence-based classification
  - Taxonomic filtering
  - Multi-sample reports

**Read Utils CLI:**
- Location: `src/viral_ngs/read_utils.py`
- Triggers: Command `read_utils <subcommand>`
- Responsibilities:
  - Read format conversion (BAM/FASTQ)
  - Quality trimming and filtering
  - Duplicate removal
  - Subsampling

**Interhost Phylogenetics CLI:**
- Location: `src/viral_ngs/interhost.py`
- Triggers: Command `interhost <subcommand>`
- Responsibilities:
  - Multiple sequence alignment
  - Whole genome phylogenetics
  - Consensus calling

**Intrahost Phylogenetics CLI:**
- Location: `src/viral_ngs/intrahost.py`
- Triggers: Command `intrahost <subcommand>`
- Responsibilities:
  - Intra-sample variant detection
  - Variant annotation
  - SNP calling at consensus level

## Error Handling

**Strategy:** Fail-fast with descriptive exceptions and logging

**Patterns:**
- Custom exception hierarchy: `QCError` (quality failures), `DenovoAssemblyError` (assembly failures), `InvalidBamHeaderError` (input validation)
- Subprocess error handling: `subprocess.check_call()` and `subprocess.check_output()` raise `CalledProcessError` on non-zero exit
- Validation errors: Explicit checks for tool availability, input file existence, read group presence
- Logging: All modules use `log = logging.getLogger(__name__)` with structured debug/info/error messages

## Cross-Cutting Concerns

**Logging:**
- Configured via `util_cmd.setup_logger(log_level)` in each CLI module
- Format: `timestamp - module:line:function - level - message`
- Default level: INFO, overridable via `--loglevel` flag

**Validation:**
- Input validation: File existence, BAM header checks, FASTA format validation
- Data quality checks: Read count thresholds, assembly quality metrics, coverage requirements
- Implicit in tool-specific methods (e.g., `samtools.getReadGroups()` raises on invalid BAM)

**Authentication:**
- Environment variable based: Docker-centric workflow (see INTEGRATIONS.md)
- S3 access via boto3 (implicit in `viral_ngs.core.file` URL handling)
- No explicit authentication layer; relies on external tool credentials

**Resource Management:**
- Thread pooling: `concurrent.futures.ProcessPoolExecutor` for multi-read group alignment
- Memory limits: Configurable JVM memory for Java tools (Picard, GATK)
- Temporary storage: Base temp dir configurable via `--tmp_dir` flag

---

*Architecture analysis: 2026-02-11*
