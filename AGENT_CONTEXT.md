# Agent Context: viral-ngs Monorepo Migration

> This document provides background context for the Claude Code agent that will implement the monorepo migration. Read this alongside `MONOREPO_IMPLEMENTATION_PLAN.md`.

---

## CURRENT STATUS

| Phase | Status | Notes |
|-------|--------|-------|
| Phase 0: Prepare repo | ✅ COMPLETE | Legacy archived, branches cleaned, secrets configured |
| Phase 1: Foundation | ✅ COMPLETE | pyproject.toml, baseimage, CI workflow created |
| Phase 2: Migrate viral-core | 🔲 NOT STARTED | Next step |
| Phase 3: Migrate derivatives | 🔲 NOT STARTED | assemble, phylo, classify |
| Phase 4: Finalize | 🔲 NOT STARTED | mega image, docs, badges |

**Next action:** Push to origin, then start Phase 2 (migrate viral-core with git history preservation).

---

## Overview

You are helping migrate 5 separate git repositories into a single monorepo at `github.com/broadinstitute/viral-ngs`. This is a substantial refactoring project that will modernize the codebase's CI/CD, Docker builds, and Python packaging.

## Source Repositories

All repositories are located at `/Users/dpark/dev/`:

| Repository | Path | Description |
|------------|------|-------------|
| viral-baseimage | NOT LOCAL (use GitHub API) | Base Docker image with conda/python |
| viral-core | `/Users/dpark/dev/viral-core` | Core utilities, illumina, read_utils, tools/ |
| viral-classify | `/Users/dpark/dev/viral-classify` | Metagenomics, taxonomy filtering |
| viral-assemble | `/Users/dpark/dev/viral-assemble` | Genome assembly tools |
| viral-phylo | `/Users/dpark/dev/viral-phylo` | Phylogenetic analysis |

**Target repository:** `/Users/dpark/dev/viral-ngs` (this directory - monorepo in progress)

## Codebase Architecture

### Current Structure (per repo)

Each repository follows this pattern:
```
repo/
├── *.py                    # Top-level command modules (illumina.py, assembly.py, etc.)
├── util/                   # Utility functions (file.py, cmd.py, misc.py)
├── tools/                  # Bioinformatics tool wrappers (samtools.py, picard.py)
├── docker/                 # Docker build scripts
├── test/
│   ├── unit/               # pytest tests
│   └── input/              # Test data files
├── docs/                   # Sphinx documentation
├── requirements-conda.txt  # Conda dependencies
├── Dockerfile
├── CLAUDE.md               # AI assistant guidance
└── DEVELOPMENT_NOTES.md    # Developer documentation
```

### Target Structure (monorepo)

```
viral-ngs/
├── pyproject.toml
├── src/viral_ngs/
│   ├── __init__.py
│   ├── illumina.py, read_utils.py, ...  (from viral-core)
│   ├── util/                             (from viral-core)
│   ├── tools/                            (from viral-core)
│   ├── metagenomics.py, taxon_filter.py  (from viral-classify)
│   ├── classify/                         (from viral-classify)
│   ├── assembly.py                       (from viral-assemble)
│   ├── assemble/                         (from viral-assemble)
│   ├── interhost.py, intrahost.py, ncbi.py (from viral-phylo)
│   └── phylo/                            (from viral-phylo)
├── docker/
│   ├── Dockerfile.baseimage, Dockerfile.core, etc.
│   └── requirements/
├── tests/
│   ├── conftest.py
│   ├── unit/
│   └── input/
├── .github/workflows/
├── docs/
├── README.md
├── CLAUDE.md               (points to AGENTS.md)
└── AGENTS.md               (consolidated guidance)
```

## Key Technical Details

### Python Command Pattern

Each command module follows this pattern:
```python
# At module top
__commands__ = []

def parser_<command_name>(parser=argparse.ArgumentParser()):
    # Define arguments
    return parser

def main_<command_name>(args):
    # Implementation
    pass

__commands__.append(('command_name', parser_command_name))
```

This is used by the CLI framework to auto-discover commands.

### Tool Wrappers

The `tools/` directory contains wrappers for external bioinformatics tools:
```python
# tools/samtools.py
class SamtoolsTool(tools.Tool):
    def __init__(self, install_methods=None):
        # ...

    def execute(self, command, *args):
        # Run samtools with arguments
```

Tools use conda environment variables to find executables.

### Conda Dependency Resolution Pattern

**IMPORTANT:** When building derivative images (classify, assemble, phylo), dependencies from BOTH the derivative AND core must be installed in a single resolver call:

```bash
# Current pattern in viral-assemble/Dockerfile line 10:
$VIRAL_NGS_PATH/docker/install-conda-dependencies.sh \
  $VIRAL_ASSEMBLE_PATH/requirements-conda.txt \
  $VIRAL_NGS_PATH/requirements-conda.txt
```

This prevents the resolver from making version regressions. **Preserve this pattern.**

### viral-classify Multi-Environment Situation

Currently uses 4 conda environments due to dependency conflicts:
- **Main env**: Most tools
- **env2**: kaiju → **RETIRE** (remove tool wrapper)
- **env3**: diamond → **RETIRE** (remove tool wrapper)
- **env4**: kraken2, krona → **TRY TO CONSOLIDATE** back into main

The goal is to simplify to a single environment if possible.

## Import Changes Required

All imports must change from flat to package style. Here are the patterns:

```python
# util imports
import util.file          →  from viral_ngs.util import file
from util.file import X   →  from viral_ngs.util.file import X
import util.cmd           →  from viral_ngs.util import cmd

# tools imports
import tools.samtools     →  from viral_ngs.tools import samtools
from tools import picard  →  from viral_ngs.tools import picard
import tools              →  from viral_ngs import tools

# Cross-module imports (in classify referencing core)
import read_utils         →  from viral_ngs import read_utils
from illumina import X    →  from viral_ngs.illumina import X
```

Use find/replace carefully. Test imports after each module migration.

## GitHub Actions Modernization

### Old Pattern (Travis-era, custom scripts)

```yaml
# Current viral-core/.github/workflows/build.yml
- name: Build Docker
  run: |
    ./github_actions_ci/deploy-docker.sh
```

### New Pattern (modern actions)

```yaml
# Use off-the-shelf actions
- uses: docker/setup-buildx-action@v3
- uses: docker/login-action@v3
- uses: docker/metadata-action@v5
- uses: docker/build-push-action@v6
- uses: dorny/paths-filter@v3  # For smart change detection
- uses: codecov/codecov-action@v4  # For coverage
```

## Docker Registry Strategy

### Single Repo with Flavor Tags

```
quay.io/broadinstitute/viral-ngs:2.6.0-core
quay.io/broadinstitute/viral-ngs:2.6.0-classify
quay.io/broadinstitute/viral-ngs:2.6.0              # mega (no suffix)
quay.io/broadinstitute/viral-ngs:main               # main branch mega
quay.io/broadinstitute/viral-ngs:latest             # alias for main mega
```

### Caching Strategy

Use registry cache on Quay.io (not GHA - 10GB limit too small for bioinformatics images):

```yaml
cache-from: type=registry,ref=quay.io/broadinstitute/viral-ngs:cache-core-amd64
cache-to: type=registry,ref=quay.io/broadinstitute/viral-ngs:cache-core-amd64,mode=max
```

## Micromamba Migration Notes

From viral-baseimage PR #29 and read-qc-tools repo, key lessons:

1. **API differences**: `conda config --add X` → `micromamba config append X`
2. **No default Python**: Must explicitly install Python in base image
3. **Symlink trick**: Create `conda` and `mamba` symlinks pointing to `/usr/bin/micromamba`
4. **Activation**: Set `ARG MAMBA_DOCKERFILE_ACTIVATE=1` for RUN commands
5. **PATH for non-interactive**: Must explicitly add `/opt/conda/bin` to PATH with `ENV` statement - the mambaorg entrypoint only sets PATH for interactive shells
6. **Base image**: Use `mambaorg/micromamba:2.4.0-ubuntu24.04` as base (simpler than building from scratch)
7. **Symlink location**: In mambaorg image, micromamba is at `/usr/bin/micromamba`, NOT `/usr/local/bin/`

## Test Strategy

### Current Pattern

Tests use pytest with these conventions:
- Test files in `test/unit/test_*.py`
- Test input files in `test/input/<TestClassName>/`
- Access via `util.file.get_test_input_path(self)`
- Run with: `pytest -n auto test/unit`

### New Pattern

- Reorganize tests by module in `tests/unit/`
- Use `dorny/paths-filter` for smart test selection
- Core changes → run all tests
- classify changes → run classify + core tests (cascade)

## Git History Preservation

Use `git filter-repo` to rewrite paths while preserving history:

```bash
pip install git-filter-repo

git clone <repo> <repo>-rewrite
cd <repo>-rewrite
git filter-repo --path-rename 'old/:new/' --force

# Then merge with --allow-unrelated-histories
```

After migration, `git log --follow <file>` will show full history.

## Files to Read for Context

If you need more context, these files are particularly useful:

| File | Purpose |
|------|---------|
| `/Users/dpark/dev/viral-core/CLAUDE.md` | Current AI guidance for viral-core |
| `/Users/dpark/dev/viral-core/DEVELOPMENT_NOTES.md` | Developer documentation |
| `/Users/dpark/dev/viral-core/Dockerfile` | Current Docker pattern |
| `/Users/dpark/dev/viral-core/.github/workflows/build.yml` | Current CI workflow |
| `/Users/dpark/dev/viral-classify/Dockerfile` | Multi-env pattern example |
| `/Users/dpark/dev/viral-core/requirements-conda.txt` | Core dependencies |
| `/Users/dpark/dev/viral-core/docker/install-conda-dependencies.sh` | Dep installer |

## Decision Log

These decisions have been made and should not be revisited:

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Python packaging | `viral-ngs` pip package | Modern, cleaner imports |
| Import style | `from viral_ngs.util import file` | Standard package pattern |
| Docker registry | Quay.io + ghcr.io | Redundancy, both registries |
| Multi-arch | QEMU via docker/setup-qemu-action | Simpler than native runners for now |
| Cache backend | Registry on Quay.io | GHA 10GB limit too small |
| Default branch | `main` | Modern convention |
| Code coverage | CodeCov | Replacing Coveralls |
| Retire tools | kaiju, diamond | Reduce complexity |
| Consolidate | kraken2, krona to main env | Try, keep multi-env if fails |
| Feature branch cleanup | `quay.expires-after=10w` label | No API token needed |
| License | MIT | Replacing BSD-3-Clause |
| Base image | `mambaorg/micromamba:2.4.0-ubuntu24.04` | Simpler than custom build |

## Common Pitfalls to Avoid

1. **Don't break imports**: Test after each module migration
2. **Don't lose git history**: Always use `git filter-repo`, never copy files manually
3. **Don't ignore conda resolution**: Always install deps together, not incrementally
4. **Don't forget ARM64**: Check bioconda package availability for each tool
5. **Don't hardcode paths**: Use `$VIRAL_NGS_PATH` environment variable
6. **Don't skip verification**: Run tests after each phase

## Communication

When working on this migration:

1. **Ask before major decisions**: If you encounter an architectural choice not covered here, ask
2. **Report blockers immediately**: Especially ARM64 package issues or conda conflicts
3. **Test incrementally**: Don't try to do everything at once
4. **Document surprises**: If something doesn't work as expected, note it for future reference

## Next Steps

1. Read `MONOREPO_IMPLEMENTATION_PLAN.md` for the detailed task list
2. **Phase 2**: Migrate viral-core with git history preservation using `git filter-repo`
3. Work through phases sequentially
4. Verify each phase before moving to the next

## What Was Done in Phase 0 & 1

### Phase 0 (Complete)
- Created `archive/legacy-monolith` tag and `archive/legacy` branch
- Deleted 148 stale feature branches (ct-*, is-*, sy-*, etc.)
- Deleted all legacy files (783 files, 83K+ lines)
- Replaced BSD-3-Clause license with MIT
- Renamed `master` to `main`, set as default branch
- Closed all 53 legacy issues with migration notice
- Configured secrets: `QUAY_USERNAME`, `QUAY_TOKEN`, `CODECOV_TOKEN`
- Updated repo description and topics

### Phase 1 (Complete)
- Created `pyproject.toml` with setuptools-scm
- Created `src/viral_ngs/__init__.py` and `py.typed`
- Created `docker/Dockerfile.baseimage` using mambaorg/micromamba
- Created `docker/install-conda-deps.sh` for unified dependency installation
- Created `.github/workflows/docker.yml` for multi-arch builds
- Verified baseimage builds locally with Python 3.12 and conda/mamba symlinks

## Reference Repositories

When implementing Docker patterns, refer to:
- `broadinstitute/read-qc-tools` - Simple micromamba-based Dockerfile pattern
- `broadinstitute/viral-baseimage` branch `ct-reduce-google-sdk-install-size` - PR #29 for micromamba migration ideas
- `broadinstitute/viral-core` - Current Dockerfile and install scripts for reference
