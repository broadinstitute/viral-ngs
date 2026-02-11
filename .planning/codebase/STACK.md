# Technology Stack

**Analysis Date:** 2026-02-11

## Languages

**Primary:**
- Python 3.10+ (3.10, 3.11, 3.12 supported) - Core application and CLI tools

**Secondary:**
- Perl - Limited usage in FASTA processing scripts (via perl-bio-easel and sequip)
- Bash - Shell scripts in `docker/scripts/`

## Runtime

**Environment:**
- Docker (primary) - All images built on `mambaorg/micromamba:2.4.0-ubuntu24.04`
- micromamba 2.4.0+ - Conda-compatible package manager in containers
- Linux (x86_64 and ARM64 architectures supported)
- Python 3.12 in base Docker images

**Package Manager:**
- pip (Python) - For Python package installation (`requirements` via setuptools)
- micromamba (conda-compatible) - For bioinformatics tools and compiled dependencies
- conda/mamba - Symlinks provided for compatibility with scripts expecting conda
- Lockfile: Not explicit; conda channel resolution happens at build time

## Frameworks

**Core:**
- setuptools 61.0+ - Python package building and distribution
- setuptools-scm - Dynamic version management from git tags

**Testing:**
- pytest 7.0+ - Test runner
- pytest-cov 4.0+ - Coverage reporting
- pytest-xdist 3.0+ - Parallel test execution

**Build/Dev:**
- GitHub Actions - CI/CD pipeline (see `.github/workflows/docker.yml`)
- Docker - Multi-stage builds with layered images (baseimage → core → assemble/classify/phylo)
- Sphinx - Documentation generation (via Read the Docs)

## Key Dependencies

**Critical Python (from pyproject.toml):**
- arrow 0.12.1+ - Date/time handling
- biopython 1.72+ - Sequence analysis and GenBank access via Entrez API
- jinja2 2.11.3+ - Template rendering
- pandas 2.0.0+ - Data manipulation and analysis
- psutil 6.1.0+ - System and process utilities
- pysam 0.23.0+ - SAM/BAM file handling
- pyyaml 6.0+ - YAML configuration parsing
- scipy 1.10.0+ - Scientific computing
- lz4 2.2.1+ - Compression
- zstandard 0.23.0+ - Zstandard compression
- matplotlib 2.2.4+ - Data visualization
- seaborn 0.12.0+ - Statistical visualization
- pybedtools 0.7.10+ - BED file manipulation

**Infrastructure/Bioinformatics Tools (via conda, from docker/requirements/):**
- bbmap 39.10 - Read mapping and assembly evaluation
- bcftools 1.10+ - VCF/BCF file manipulation
- bedtools 2.29.2+ - Genomic interval operations
- bwa 0.7.17+ - Read alignment
- samtools 1.21+ - SAM/BAM manipulation
- picard 2.25.6 - SAM/BAM QC and manipulation
- gatk 3.8 - Variant calling
- fastqc 0.11.7+ - Read quality assessment
- trimmomatic 0.38+ - Read trimming
- spades 4.2.0+ - De novo genome assembly
- mafft 7.525+ - Multiple sequence alignment
- muscle 3.8.1551 - Sequence alignment
- minimap2 2.17+ - Sequence mapping
- kraken2 2.1.3+ - Metagenomic classification
- blast 2.15.0+ - Sequence similarity search
- snpeff 4.3.1t+ - Variant effect prediction
- lofreq 2.1.5+ - Variant calling
- mummer4 4.0.0+ - Sequence comparison
- cd-hit 4.6.8+ - Sequence clustering

**Cloud/Utilities (from docker/requirements/baseimage.txt):**
- miniwdl 1.11.1+ - WDL workflow execution (uses udocker for local containerization)
- udocker 1.3.16+ - Container runtime for non-privileged execution
- awscli 1.32.0+ - AWS command-line tools (S3, EC2, etc.)
- google-cloud-storage 2.14.0+ - Google Cloud Storage client
- jq 1.6+ - JSON query and manipulation
- parallel 20190922+ - GNU Parallel for concurrent execution
- pigz 2.4+ - Parallel gzip compression
- zstd 1.5.7+ - Zstandard compression utility

## Configuration

**Environment:**
- Configuration via command-line arguments (argparse)
- Optional environment variables:
  - `VIRAL_NGS_SOURCE_DIR` - Override project path (used in Docker)
  - `VIRAL_NGS_TMP_DIRKEEP` - Keep temporary files for debugging
  - `NOVOALIGN_PATH` - Path to novoalign binary (x86-only tool)
  - `NOVOALIGN_LICENSE_PATH` - Path to novoalign license file
  - `PYTEST_XDIST_WORKER_COUNT` - Override parallel test worker count
  - `PATH`, `PYTHONPATH` - Standard shell variables
  - `MINIWDL__SCHEDULER__CONTAINER_BACKEND` - Set to "udocker" in baseimage
- NCBI Entrez access: email address and optional `api_key` passed as function arguments

**Build:**
- `pyproject.toml` - Python package metadata and dependencies
- `docker/requirements/*.txt` - Conda package lists by functionality (baseimage, core, assemble, classify, phylo)
- `Dockerfile.baseimage`, `Dockerfile.core`, `Dockerfile.assemble`, `Dockerfile.classify`, `Dockerfile.phylo` - Multi-stage Docker builds
- `.readthedocs.yml` - Sphinx documentation build configuration
- `.codecov.yml` - Code coverage reporting
- `.coveragerc` - Coverage measurement configuration

## Platform Requirements

**Development:**
- Python 3.10+ (3.12 recommended)
- micromamba or conda for dependency installation
- Docker for testing (required per CLAUDE.md)
- Git for version management

**Production:**
- Docker images published to:
  - Quay.io: `quay.io/broadinstitute/viral-ngs` (primary registry)
  - GitHub Container Registry: `ghcr.io/broadinstitute/viral-ngs`
- Multi-architecture support: linux/amd64 (x86-64) and linux/arm64 (ARM64/Apple Silicon)
- x86-only tools (novoalign, mvicuna) excluded from ARM64 builds

## Image Variants

- `viral-ngs:main-baseimage` - Base with micromamba, Python 3.12, utilities
- `viral-ngs:main-core` - Core + bioinformatics tools (samtools, bwa, picard, etc.)
- `viral-ngs:main-assemble` - Core + assembly tools (SPAdes, MUMmer, MAFFT)
- `viral-ngs:main-classify` - Core + classification tools (Kraken2, BLAST, BMTagger)
- `viral-ngs:main-phylo` - Core + phylogenetics tools (SnpEff, LoFreq, VPhaser2)
- `viral-ngs:main-mega` - All tools combined

---

*Stack analysis: 2026-02-11*
