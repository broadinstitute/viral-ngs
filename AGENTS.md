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
   micromamba install <packages>
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

**Important: Testing source code changes requires re-installing the package.**
The `-v` mount makes your local files visible on disk, but `viral_ngs` is already installed as a package inside the container image. Python imports resolve to the *installed* copy, not your mounted source files. If you've modified files under `src/viral_ngs/`, you must re-install before running tests:

```bash
# Run tests with local source changes applied
docker run --rm \
  -v $(pwd):/opt/viral-ngs/source \
  quay.io/broadinstitute/viral-ngs:main-core \
  bash -c "pip install -e /opt/viral-ngs/source --quiet && pytest -rsxX -n auto /opt/viral-ngs/source/tests/unit"
```

Changes to test files (`tests/`) and test inputs (`tests/input/`) are picked up automatically via the volume mount — the re-install is only needed when modifying the `src/viral_ngs/` package code.

Running pytest directly on the host will generally not work — most dependencies (bioinformatics tools, conda packages) are only available inside the Docker containers. Always test inside Docker.

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
│       ├── core-x86.txt        # x86-only core packages
│       ├── assemble.txt
│       ├── assemble-x86.txt    # x86-only assembly packages
│       ├── classify.txt
│       ├── classify-x86.txt    # x86-only classify packages
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
   micromamba search <package>              # default channel
   micromamba search -c bioconda <package>  # bioconda channel
   ```

2. Add to appropriate requirements file:
   - `docker/requirements/core.txt` - core dependencies
   - `docker/requirements/assemble.txt` - assembly-specific
   - `docker/requirements/classify.txt` - classification-specific
   - `docker/requirements/phylo.txt` - phylo-specific

3. For x86-only packages (no ARM64 build), add to the appropriate `-x86.txt` file:
   - `core-x86.txt` - novoalign, mvicuna
   - `classify-x86.txt` - bmtagger, kallisto, kb-python
   - `phylo-x86.txt` - table2asn

### Dependency Resolution

When building derivative images, ALL dependencies (including x86-only) must be installed in a **single resolver call** using the `--x86-only:` prefix:

```bash
# Single resolver call - x86-only files skipped on ARM64
/tmp/install-conda-deps.sh \
  /tmp/requirements/baseimage.txt \
  /tmp/requirements/core.txt \
  /tmp/requirements/classify.txt \
  --x86-only:/tmp/requirements/classify-x86.txt
```

This prevents version regressions. **Never install incrementally.**

The `install-conda-deps.sh` script:
- On x86: Includes all files in one micromamba call
- On ARM64: Skips files tagged with `--x86-only:` but includes others

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

**Build Architecture:**
Each image flavor is built using 3 parallel jobs for native multi-arch support:
1. `build-{flavor}-amd64` - runs on `ubuntu-latest`
2. `build-{flavor}-arm64` - runs on `ubuntu-24.04-arm` (native ARM runner)
3. `create-manifest-{flavor}` - combines arch-specific images into multi-arch manifest

This approach is 3-5x faster than QEMU emulation for ARM builds.

**Build Job Flow:**
```
paths-filter + get-version (parallel)
         ↓
build-baseimage-amd64  ←→  build-baseimage-arm64  (parallel)
         ↓                          ↓
    create-manifest-baseimage
         ↓
build-core-amd64  ←→  build-core-arm64  (parallel)
         ↓                    ↓
    create-manifest-core
         ↓
build-{assemble,classify,phylo,mega}-amd64  ←→  build-{...}-arm64  (parallel)
         ↓
    create-manifest-{flavor}
         ↓
    test-{flavor} + test-{flavor}-arm64 (ARM64 tests only on PRs with docker changes)
         ↓
    deploy-to-quay (push/tag events only)
```

**Test Jobs:**
- **test-core**: Runs on x86, tests `tests/unit/core/`
- **test-assemble**: Runs on x86, tests `tests/unit/assemble/`
- **test-classify**: Runs on x86, tests `tests/unit/classify/`
- **test-phylo**: Runs on x86, tests `tests/unit/phylo/`
- **test-{flavor}-arm64**: Runs on native ARM, only on PRs when docker files change

**Smart Test Scoping:**
Tests only run when relevant code changes:
- Core tests: `src/viral_ngs/*.py`, `core/**`, `util/**`, `tests/unit/core/**`
- Assemble tests: `assemble/**`, `assembly.py`, or core changes
- Classify tests: `classify/**`, `metagenomics.py`, `taxon_filter.py`, or core changes
- Phylo tests: `phylo/**`, `interhost.py`, `intrahost.py`, `ncbi.py`, or core changes
- Docker changes trigger all tests (including ARM64 tests on PRs)

**Coverage:**
Each x86 test job uploads coverage to Codecov with flavor-specific flags.

### Multi-Architecture Support

- Images built natively for `linux/amd64` and `linux/arm64` using parallel runners
- Multi-arch manifests created with OCI annotations using `docker buildx imagetools create`
- x86-only packages (novoalign, mvicuna, bmtagger, kallisto, kb-python, table2asn) handled via `--x86-only:` prefix in `install-conda-deps.sh`
- Python tool wrappers still importable on ARM64; only runtime execution fails for missing binaries
- Tests using x86-only tools have `@unittest.skipIf(IS_ARM, ...)` decorators
- Architecture-specific caches prevent cross-arch cache pollution

### ARM Test Skipping

Tests that use x86-only bioconda packages must be decorated to skip on ARM:

```python
from tests import IS_ARM

SKIP_X86_ONLY_REASON = "tool requires x86-only bioconda package (not available on ARM)"

@unittest.skipIf(IS_ARM, SKIP_X86_ONLY_REASON)
class TestSomeTool(TestCaseWithTmp):
    ...

# Or at method level:
@unittest.skipIf(IS_ARM, SKIP_X86_ONLY_REASON)
def test_specific_tool(self):
    ...
```

### Documentation Build

The `docs.yml` workflow builds Sphinx documentation. Key points:
- Uses `mock` to stub heavy dependencies (`Bio`, `pysam`, `scipy`, etc.) in `docs/conf.py`
- When adding new imports to source code, add corresponding mocks to `MOCK_MODULES` in `docs/conf.py`
- Runs `sphinx-build -W` (warnings as errors)

### Registry Strategy

- **GHCR (ghcr.io)**: Primary build registry, images pushed during CI for all events including PRs
- **Quay.io**: Production registry, images copied from GHCR after tests pass (push/tag events only)
- Feature branch images should be cleaned up periodically from Quay.io

---

## Coding Guidelines

### Agent Attribution

**Commit messages**: By default, do NOT include agent/model credits (e.g., "Co-Authored-By: Claude") in commit messages. This reduces noise in the git history.

**Code review comments**: DO include notes about agent/model involvement when writing code review comments (e.g., PR reviews, inline comments). This provides useful context about how the review was conducted.

**Explicit requests**: Include agent attribution in commits or elsewhere when explicitly requested by a human reviewer or contributor.

**Avoid amending pushed commits**: Do not use `git commit --amend` after a commit has been pushed to the remote. Amending pushed commits causes problems for collaboration. Instead, create a new commit with the fix. Amending is fine for local commits that haven't been pushed yet.

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
micromamba search -c bioconda <package> --subdir linux-aarch64
```

---

## Container Vulnerability Management

### Scanning

Container images are scanned for vulnerabilities using [Trivy](https://aquasecurity.github.io/trivy/):

- **On every PR/push**: `docker.yml` scans each image flavor after build (SARIF → GitHub Security tab, JSON → artifact)
- **Weekly schedule**: `container-scan.yml` scans the latest published images
- Scans filter to **CRITICAL/HIGH** severity, **ignore-unfixed**, and apply a Rego policy (`.trivy-ignore-policy.rego`)
- Per-CVE exceptions go in `.trivyignore` with mandatory justification comments

### Rego Policy (`.trivy-ignore-policy.rego`)

The Rego policy filters CVEs that are architecturally inapplicable to ephemeral batch containers:

- **AV:P** (Physical access required) — containers are cloud-hosted
- **AV:A** (Adjacent network required) — no attacker on same network segment
- **AV:L + UI:R** (Local + user interaction) — no interactive sessions
- **AV:L + PR:H** (Local + high privileges) — containers run non-root
- **AV:L + S:U** (Local + scope unchanged) — attacker already has code execution and impact stays within the ephemeral container

Changes to this policy should be reviewed carefully. The comments in the file explain the rationale and risk for each rule.

### Common Vulnerability Sources

**Python transitive deps**: Pin minimum versions in `docker/requirements/*.txt`. Prefer conda packages over pip. Check conda-forge availability before assuming a version exists — conda-forge often lags PyPI by days/weeks.

**Java fat JARs** (picard, gatk, snpeff, fgbio): Bioinformatics Java tools are distributed as uber JARs with all dependencies bundled inside. Trivy detects vulnerable libraries (log4j, commons-compress, etc.) baked into these JARs. Version bumps can cause ARM64 conda solver conflicts because Java tools pull in openjdk → harfbuzz → icu version chains that clash with other packages (r-base, boost-cpp, pyicu). Always check:
1. Whether the tool is actually flagged by Trivy (don't bump versions unnecessarily)
2. Whether the CVE applies (e.g., log4j 1.x is NOT vulnerable to Log4Shell)
3. Whether the desired version resolves on ARM64 before pushing

**Go binaries**: Some conda packages bundle compiled Go binaries (e.g., mafft's `dash_client`, google-cloud-sdk's `gcloud-crc32c`). If the binary is unused, delete it in the Dockerfile. Delete from **both** the installed location and `/opt/conda/pkgs/*/` (conda package cache) — Trivy scans the full filesystem.

**Vendored copies**: Packages like google-cloud-sdk and setuptools bundle their own copies of Python libraries that may be older than what's in the conda environment. Trivy flags these vendored copies separately. Options: delete the vendored directory (if not needed at runtime), or accept the risk in `.trivyignore` with justification.

### ARM64 Solver Conflicts

The conda solver on ARM64 (linux-aarch64) is more constrained than amd64 because fewer package builds exist. Common conflict patterns:

- **icu version conflicts**: Many packages (openjdk, r-base, boost-cpp, pyicu) pin specific icu version ranges. Bumping one package can make the entire environment unsolvable.
- **libdeflate/htslib conflicts**: lofreq 2.1.5 pins old htslib/libdeflate versions that conflict with newer pillow/libtiff.
- **openjdk version escalation**: snpeff 5.2+ requires openjdk>=11, 5.3+ requires openjdk>=21. Higher openjdk versions pull in harfbuzz→icu chains that conflict with everything.

When a solver conflict occurs: revert the change, check what version the solver was picking before, and pin to that exact version if it already addresses the CVE.

### Mitigation Decision Process

When triaging a CVE:

1. **Check the CVSS vector** — does the Rego policy already filter it?
2. **Identify the source package** — use Trivy JSON output (`PkgName`, `PkgPath`, `InstalledVersion`)
3. **Check if a fix version exists on conda-forge/bioconda** — not just on PyPI
4. **Test on ARM64** — solver conflicts are the most common failure mode
5. **If the fix version conflicts**: consider whether the CVE is exploitable in your deployment model. Document the risk assessment in `.trivyignore` or `vulnerability-mitigation-status.md`.
6. **If the vulnerable code is unused**: delete the binary/file inline in the Dockerfile (same RUN layer as install to avoid bloating images)

### Key Files

| File | Purpose |
|------|---------|
| `.trivy-ignore-policy.rego` | Rego policy for class-level CVE filtering |
| `.trivyignore` | Per-CVE exceptions with justifications |
| `.github/workflows/docker.yml` | Build-time scanning (SARIF + JSON) |
| `.github/workflows/container-scan.yml` | Weekly scheduled scanning |
| `vulnerability-mitigation-status.md` | Local-only tracking doc (not committed) |

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
