# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

viral-core is a set of scripts and tools for analyzing viral NGS (Next-Generation Sequencing) data. It provides utilities for demultiplexing Illumina data, working with sequence reads, read alignment, quality control, and various bioinformatics operations.

More detailed command line documentation: http://viral-core.readthedocs.org/
Developer documentation: [DEVELOPMENT_NOTES.md](DEVELOPMENT_NOTES.md)
Higher level pipelines: https://github.com/broadinstitute/viral-pipelines

## Development Environment

### Docker-Centric Development

The development paradigm is **intentionally docker-centric**. Developers need:
- Docker
- Git
- Text/code editor

Development workflow:

1. Check out the repository: `git clone https://github.com/broadinstitute/viral-core.git`
2. Pull and run the base image, mounting local checkout:
   ```bash
   docker run -it --rm -v `pwd`/viral-core:/opt/viral-ngs/source quay.io/broadinstitute/viral-core
   ```
3. If modifying conda dependencies in `requirements-conda.txt`, install them inside container:
   ```bash
   conda install <packages>
   ```
   Optionally snapshot the new image: `docker commit <image_hash> local/viral-core-dev`
4. Test code interactively inside container:
   ```bash
   cd /opt/viral-ngs/source
   pytest -rsxX -n auto test/unit
   ```
5. Push changes to GitHub for automated CI testing & builds

### Running Tests

```bash
# Run all unit tests in parallel
docker run --rm \
  -v $(pwd):/opt/viral-ngs/source \
  quay.io/broadinstitute/viral-core:$(git branch --show-current) \
  pytest -rsxX -n auto /opt/viral-ngs/source/test/unit
```

**Important test notes:**
- Uses `py.test` (not `nose` or `unittest`) for generative testing, fixtures, and parallelized execution
- Test input files in `test/input/<TestClassName>/`
- Access via `util.file.get_test_input_path(self)` for class-specific inputs
- Access parent directory via `util.file.get_test_path()`
- Custom pytest option: `--fixture-durations N` shows N slowest fixtures
- Custom pytest marker: `@pytest.mark.slow` for slow tests

### Documentation

Documentation is auto-built by Read the Docs for master branch. To test locally:

```bash
cd docs
make html
```

- Format: reStructuredText in `docs/`
- `sphinx-argparse` auto-generates argparse documentation
- CI tests doc builds but doesn't deploy (Read the Docs handles deployment separately)

## Code Architecture

### Module Structure

The codebase is modularized and layered:

**Top-level modules** (`illumina.py`, `read_utils.py`, `reports.py`, etc.):
- High-level command-line tools and workflows
- Define `__commands__` list
- Parser functions: `parser_<command_name>()`
- Main functions: `main_<command_name>(args)`

**`tools/` directory**: Wrappers for external bioinformatics tools
- Each tool has its own module: `tools/bwa.py`, `tools/samtools.py`, `tools/picard.py`, etc.
- Base `Tool` class in `tools/__init__.py` provides install machinery
- Tools track installation methods, executable paths, and versions
- When Python/binary dependencies are installed by conda, they use conda environment variables to specify paths

**`util/` directory**: Core utility functions
- `util/file.py`: File handling (gzip, bz2, lz4, zstd compression), FASTA/BAM/tab-text operations
- `util/misc.py`: General utilities (memoization, unique, timing, counting)
- `util/cmd.py`: Command-line argument parsing helpers
- `util/stats.py`: Statistical utilities
- `util/illumina_indices.py`: Illumina barcode/index reference data
- `util/version.py`: Version tracking

**Test structure** (`test/`):
- `test/unit/`: Unit test files (e.g., `test_tools_samtools.py`, `test_illumina.py`)
- `test/input/`: Static test input files organized by test class
- `test/__init__.py`: Test utilities
- `conftest.py`: Pytest configuration and fixtures

### Dependency Installation

When adding a new tool or dependency:

1. Check if conda package exists:
   ```bash
   conda search <package_name>              # default channel
   conda search -c bioconda <package_name>  # bioconda channel
   ```
2. If available, add to `requirements-conda.txt`
3. If unavailable, first create a recipe for bioconda (https://github.com/bioconda/bioconda-recipes)

Conda packages are installed in the active conda environment. The code accesses the active environment path via environment variables.

## CI/CD (GitHub Actions)

Configuration: `.github/workflows/build.yml`

**Build matrix** (triggered on every branch commit and pull request):

1. **Docker container build** (`build_docker` job):
   - Builds docker image with layer caching (pulls previous build parent)
   - Typical build: 10-20 seconds (only source code layer changes)
   - Dependency changes: +10 minutes (rebuilds heavy conda install layer)
   - Master branch: pushes to `quay.io/broadinstitute/viral-core:latest` + versioned tag
   - Non-master/PRs: pushes versioned tag with 10-week expiration
   - Requires Quay.io credentials in GitHub Actions secrets

2. **Test in docker** (`test_in_docker` job):
   - Pulls built docker image
   - Runs `pytest -n $(nproc) test/unit` inside container
   - Generates coverage reports
   - Sends coverage to coveralls.io

3. **Documentation test** (`test_docs` job):
   - Tests documentation build
   - Does NOT deploy (Read the Docs has separate auto-build process)
   - Completes in <1 minute

## Dependencies

### Conda Dependencies (`requirements-conda.txt`)

**Bioinformatics tools:**
- Aligners: bwa, minimap2, novoalign
- SAM/BAM: samtools, picard, gatk, fgbio
- Assembly: mvicuna
- QC/filtering: fastqc, prinseq, trimmomatic
- Utilities: bcftools, bedtools, cd-hit, bbmap, parallel, splitcode
- Data: ncbi-datasets-cli
- Compression: lbzip2, lz4-c, pigz, zstd, unzip

**Python packages:**
- Core: biopython, pysam, pybedtools, psutil
- Data: pandas, matplotlib, arrow
- Cloud: firecloud
- Compression: lz4, zstandard
- Templating: jinja2, csvkit, jq

### Test Dependencies (`requirements-conda-tests.txt`)

pytest, pytest-cov, pytest-xdist, pytest-mock, coveralls, lxml, mock, coverage

## Docker Image

The Dockerfile builds on a viral-baseimage (see `FROM` line in Dockerfile for current version):
- Python 3.12
- Conda environment: `viral-ngs-env`
- Install path: `/opt/viral-ngs`
- Source code: `/opt/viral-ngs/source`
- PYTHONPATH: `/opt/viral-ngs/source`

Note: Container includes R/readline symlink workaround needed for Picard's CollectInsertSizeMetrics
