# AGENTS.md

This document provides guidance for AI assistants (Claude Code, GitHub Copilot, etc.) working on this repository.

## Overview

viral-ngs is a consolidated monorepo for viral NGS (Next-Generation Sequencing) analysis tools. It provides:

- **Core utilities**: Read manipulation, Illumina demultiplexing, file handling, QC
- **Assembly**: Genome assembly, scaffolding, gap filling
- **Classification**: Metagenomic classification, taxonomy filtering, k-mer analysis
- **Phylogenetics**: Variant calling, consensus generation, annotation

**Related resources:**
- Command-line documentation: https://viral-ngs.readthedocs.org/
- Higher-level pipelines: https://github.com/broadinstitute/viral-pipelines

---

## Development Environment

### Docker-Centric Development

Development is **intentionally docker-centric**. Developers need:
- Docker
- Git
- Text/code editor

### Development Workflow

1. Clone the repository:
   ```bash
   git clone https://github.com/broadinstitute/viral-ngs.git
   ```

2. Run the container with local checkout mounted:
   ```bash
   docker run -it --rm \
     -v $(pwd):/opt/viral-ngs/source \
     quay.io/broadinstitute/viral-ngs:main-core
   ```

3. If modifying conda dependencies, install them inside the container:
   ```bash
   conda install <packages>
   ```

4. Test code interactively:
   ```bash
   cd /opt/viral-ngs/source
   pytest -rsxX -n auto tests/unit
   ```

5. Push changes to GitHub for automated CI testing

### Running Tests

```bash
# Run all unit tests in the core image
docker run --rm \
  -v $(pwd):/opt/viral-ngs/source \
  quay.io/broadinstitute/viral-ngs:main-core \
  pytest -rsxX -n auto /opt/viral-ngs/source/tests/unit

# Run specific module tests
docker run --rm \
  -v $(pwd):/opt/viral-ngs/source \
  quay.io/broadinstitute/viral-ngs:main-classify \
  pytest -rsxX -n auto /opt/viral-ngs/source/tests/unit/classify
```

**Test conventions:**
- Uses pytest (not nose or unittest)
- Test files in `tests/unit/<module>/`
- Test input files in `tests/input/<TestClassName>/`
- Access via `viral_ngs.core.file.get_test_input_path(self)`
- Custom marker: `@pytest.mark.slow` for slow tests

---

## Code Architecture

### Directory Structure

```
viral-ngs/
├── pyproject.toml              # Package configuration
├── src/viral_ngs/
│   ├── __init__.py             # Version detection
│   ├── py.typed                # PEP 561 marker
│   │
│   ├── # Command modules (CLI entry points)
│   ├── illumina.py             # Illumina demux commands
│   ├── read_utils.py           # Read manipulation commands
│   ├── assembly.py             # Assembly commands
│   ├── metagenomics.py         # Classification commands
│   ├── interhost.py            # Phylo commands
│   │
│   ├── core/                   # Core library (shared utilities + tool wrappers)
│   │   ├── __init__.py         # Tool/InstallMethod classes
│   │   ├── samtools.py         # Tool wrapper
│   │   ├── picard.py           # Tool wrapper
│   │   ├── file.py             # File utilities
│   │   ├── misc.py             # General utilities
│   │   └── ...
│   │
│   ├── assemble/               # Assembly tool wrappers
│   │   ├── __init__.py
│   │   ├── spades.py
│   │   └── ...
│   │
│   ├── classify/               # Classification tool wrappers
│   │   ├── __init__.py
│   │   ├── kraken2.py
│   │   └── ...
│   │
│   └── phylo/                  # Phylogenetics tool wrappers
│       ├── __init__.py
│       ├── mafft.py
│       └── ...
│
├── docker/
│   ├── Dockerfile.baseimage    # Base with conda/python
│   ├── Dockerfile.core         # Core tools
│   ├── Dockerfile.assemble     # + assembly tools
│   ├── Dockerfile.classify     # + classification tools
│   ├── Dockerfile.phylo        # + phylo tools
│   ├── Dockerfile.mega         # All tools combined
│   ├── install-conda-deps.sh
│   └── requirements/
│       ├── baseimage.txt
│       ├── core.txt
│       ├── core-x86.txt        # x86-only packages
│       ├── assemble.txt
│       ├── classify.txt
│       ├── phylo.txt
│       └── phylo-x86.txt       # x86-only phylo packages
│
├── tests/
│   ├── conftest.py
│   ├── unit/
│   │   ├── core/
│   │   ├── assemble/
│   │   ├── classify/
│   │   └── phylo/
│   └── input/
│
├── scripts/                    # Utility scripts
├── .github/workflows/
│   └── docker.yml              # CI/CD workflow
└── docs/
```

### Command Module Pattern

Command modules define CLI entry points:

```python
__commands__ = []

def parser_<command_name>(parser=argparse.ArgumentParser()):
    # Define arguments
    return parser

def main_<command_name>(args):
    # Implementation
    pass

__commands__.append(('command_name', parser_command_name))
```

### Tool Wrapper Pattern

Tool wrappers in `core/`, `assemble/`, `classify/`, `phylo/`:

```python
import viral_ngs.core as core

class SamtoolsTool(core.Tool):
    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [core.PrexistingUnixCommand('samtools')]
        super().__init__(install_methods=install_methods)

    def execute(self, command, *args):
        # Run samtools with arguments
        pass
```

---

## Import Patterns

### Standard imports

```python
# Within command modules (illumina.py, assembly.py, etc.)
import viral_ngs.core as core
import viral_ngs.core.file as util_file
import viral_ngs.core.misc as util_misc

# Using tools
samtools = core.samtools.SamtoolsTool()
bwa = core.bwa.BwaTool()

# Using utilities
util_file.mkstempfname()
util_misc.available_cpu_count()
```

### Within core/ modules (use relative imports)

```python
from . import samtools, picard
from .file import mkstempfname
from .misc import available_cpu_count
```

### Within subpackages (assemble/, classify/, phylo/)

```python
import viral_ngs.core as core
import viral_ngs.core.file as util_file

# For other tools in same subpackage
from . import mummer, mafft
```

### Key rules

1. **Prefer full imports**: `import viral_ngs.core.samtools` over `from viral_ngs.core import samtools`
2. **Use relative imports within packages**: `from . import X` inside core/, assemble/, etc.
3. **No backward compat stubs**: `viral_ngs.tools` and `viral_ngs.util` don't exist

---

## Dependencies

### Conda-First Approach

ALL runtime dependencies are installed via conda for speed and binary compatibility.
The `pyproject.toml` has empty dependencies - conda handles everything.

### Adding Dependencies

1. Check conda availability:
   ```bash
   conda search <package>              # default channel
   conda search -c bioconda <package>  # bioconda channel
   ```

2. Add to appropriate requirements file:
   - `docker/requirements/core.txt` - core dependencies
   - `docker/requirements/assemble.txt` - assembly-specific
   - `docker/requirements/classify.txt` - classification-specific
   - `docker/requirements/phylo.txt` - phylo-specific

3. For x86-only packages (no ARM64 build), add to `core-x86.txt`

### Dependency Resolution

When building derivative images, ALL dependencies must be installed in a single resolver call:

```bash
/tmp/install-conda-deps.sh /tmp/requirements/core.txt /tmp/requirements/assemble.txt
```

This prevents version regressions. **Never install incrementally.**

---

## Docker Images

### Image Hierarchy

```
baseimage (conda/python)
    └── core (core tools)
        ├── assemble (+ assembly tools)
        ├── classify (+ classification tools)
        ├── phylo (+ phylo tools)
        └── mega (all tools)
```

### Tag Format

```
quay.io/broadinstitute/viral-ngs:2.6.0-core
quay.io/broadinstitute/viral-ngs:2.6.0-classify
quay.io/broadinstitute/viral-ngs:2.6.0              # mega (no suffix)
quay.io/broadinstitute/viral-ngs:main-core          # main branch
quay.io/broadinstitute/viral-ngs:latest             # alias for main mega
```

### Building Locally

```bash
# Build baseimage
docker build -t viral-ngs:baseimage -f docker/Dockerfile.baseimage .

# Build core (needs baseimage)
docker build --build-arg BASEIMAGE=viral-ngs:baseimage \
  -t viral-ngs:core -f docker/Dockerfile.core .

# Build derivatives (need core)
docker build --build-arg BASEIMAGE=viral-ngs:core \
  -t viral-ngs:classify -f docker/Dockerfile.classify .
```

---

## CI/CD

### GitHub Actions Workflow

The `.github/workflows/docker.yml` workflow handles building and testing:

**Build Jobs:**
1. **paths-filter**: Detect which code paths changed (using `dorny/paths-filter`)
2. **get-version**: Extract version from git describe
3. **build-baseimage**: Build base image with conda/python
4. **build-core**: Build core image (depends on baseimage)
5. **build-derivatives**: Build assemble, classify, phylo in parallel (depend on core)
6. **build-mega**: Build all-in-one image (depends on core)

**Test Jobs:**
- **test-core**: Runs after build-core, tests `tests/unit/core/`
- **test-assemble**: Runs after build-derivatives, tests `tests/unit/assemble/`
- **test-classify**: Runs after build-derivatives, tests `tests/unit/classify/`
- **test-phylo**: Runs after build-derivatives, tests `tests/unit/phylo/`

**Smart Test Scoping:**
Tests only run when relevant code changes:
- Core tests: `src/viral_ngs/*.py`, `core/**`, `util/**`, `tests/unit/core/**`
- Assemble tests: `assemble/**`, `assembly.py`, or core changes
- Classify tests: `classify/**`, `metagenomics.py`, `taxon_filter.py`, or core changes
- Phylo tests: `phylo/**`, `interhost.py`, `intrahost.py`, `ncbi.py`, or core changes
- Docker changes trigger all tests

**Coverage:**
Each test job uploads coverage to Codecov with flavor-specific flags.

### Multi-Architecture Support

- Images built for `linux/amd64` and `linux/arm64`
- x86-only packages (novoalign, mvicuna, table2asn) handled gracefully on ARM
- Cache stored on Quay.io registry (GHA 10GB limit too small)

### Feature Branch Images

Feature branch images get `quay.expires-after=10w` label for automatic cleanup.

---

## Coding Guidelines

### Test-Driven Development

1. Write tests first
2. Verify tests fail
3. Implement feature
4. Verify tests pass
5. Refactor if needed

### Avoid Over-Engineering

- Only make changes directly requested
- Don't add features beyond what was asked
- Don't add comments/docstrings to unchanged code
- Don't create abstractions for one-time operations

### Security

- Never introduce command injection, XSS, SQL injection vulnerabilities
- Validate at system boundaries (user input, external APIs)
- Trust internal code and framework guarantees

### Code Style

- Use explicit imports (prefer `import x.y.z` over `from x.y import z`)
- Follow existing patterns in the codebase
- Run tests before committing

---

## Common Operations

### Verify Python imports

```bash
docker run --rm viral-ngs:core python -c "
import viral_ngs.core
import viral_ngs.core.samtools
import viral_ngs.core.picard
import viral_ngs.core.file
import viral_ngs.core.misc
print('Core imports OK')
"
```

### Run syntax check on all files

```bash
find src tests -name "*.py" -exec python -m py_compile {} \;
```

### Check ARM64 package availability

```bash
conda search -c bioconda <package> --subdir linux-aarch64
```

---

## Troubleshooting

### Circular Import Errors

- Use relative imports within packages (`from . import X`)
- Check import order in `__init__.py` files
- Don't import parent package from child modules

### Missing Tool Errors

- Verify tool is in the appropriate requirements file
- Check if tool is x86-only (add to `core-x86.txt`)
- Tool wrappers should fail gracefully at runtime, not import time

### Docker Build Failures

- Check that BASEIMAGE arg points to existing image
- Verify all requirements files exist
- Check for conda resolution conflicts (install all deps together)
