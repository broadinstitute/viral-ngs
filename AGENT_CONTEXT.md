# Agent Context: viral-ngs Monorepo Migration

> This document provides background context for the Claude Code agent that will implement the monorepo migration. Read this alongside `MONOREPO_IMPLEMENTATION_PLAN.md`.

---

## CURRENT STATUS

| Phase | Status | Notes |
|-------|--------|-------|
| Phase 0: Prepare repo | ✅ COMPLETE | Legacy archived, branches cleaned, secrets configured |
| Phase 1: Foundation | ✅ COMPLETE | pyproject.toml, baseimage, CI workflow created |
| Phase 2: Migrate viral-core | ✅ COMPLETE | Consolidated to `core/`, all deps via conda, x86-only handling |
| Phase 3a: Migrate viral-assemble | ✅ COMPLETE | 3,448 commits, tests in tests/unit/assemble/ |
| Phase 3b-c: Migrate phylo, classify | 🔲 NOT STARTED | phylo and classify remaining |
| Phase 4: Finalize | 🔲 NOT STARTED | mega image, docs, badges |

**Next action:** Start Phase 3b (migrate viral-phylo with git history preservation).

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
├── pyproject.toml              # With project.scripts for CLI entry points
├── src/viral_ngs/
│   ├── __init__.py             # Just imports core
│   ├── py.typed
│   ├── illumina.py             # Core command module (moved from core/)
│   ├── read_utils.py           # Core command module (moved from core/)
│   ├── reports.py              # Core command module (moved from core/)
│   ├── file_utils.py           # Core command module (moved from core/)
│   ├── broad_utils.py          # Core command module (moved from core/)
│   ├── assembly.py             # Assembly command module (Phase 3a)
│   ├── interhost.py, intrahost.py, ncbi.py  # Phylo modules (Phase 3b)
│   ├── metagenomics.py, taxon_filter.py     # Classify modules (Phase 3c)
│   ├── core/                   # Core library modules (utilities + tool wrappers)
│   │   ├── __init__.py         # Tool/InstallMethod classes + submodule imports
│   │   ├── samtools.py, picard.py, bwa.py, ...  # Tool wrappers
│   │   ├── file.py, misc.py, cmd.py, ...        # Utilities
│   │   └── errors.py, priorities.py, ...
│   ├── assemble/               # Assembly tool wrappers (Phase 3a)
│   │   ├── __init__.py
│   │   ├── spades.py, mummer.py, mafft.py, ...
│   │   └── skani.py, muscle.py, gap2seq.py, ...
│   ├── phylo/                  # Phylo tool wrappers (Phase 3b)
│   └── classify/               # Classify tool wrappers (Phase 3c)
├── docker/
│   ├── Dockerfile.baseimage, Dockerfile.core, etc.
│   └── requirements/
│       ├── core.txt            # All deps (Python + bioinformatics)
│       ├── core-x86.txt        # x86-only packages
│       └── ...
├── tests/
│   ├── conftest.py
│   ├── unit/
│   └── input/
├── scripts/                    # Migration/maintenance scripts
├── .github/workflows/
├── docs/
├── README.md
├── CLAUDE.md                   # (points to AGENTS.md)
└── AGENTS.md                   # (consolidated guidance)
```

**Note:** The original `tools/` and `util/` directories were consolidated into `core/` for simplicity. There are no backward compatibility stubs - all code uses `viral_ngs.core.*` imports directly. Non-core command modules (assembly.py, interhost.py, etc.) go directly under `src/viral_ngs/`, while their tool wrappers go in subpackages (`assemble/`, `phylo/`, `classify/`).

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

All imports use the consolidated `viral_ngs.core.*` pattern.

**Style preference:** Use full absolute imports (`import x.y.z`) over `from` imports (`from x.y import z`). This makes imports explicit and easier to trace.

```python
# PREFERRED: Full absolute imports
import viral_ngs.core.samtools
import viral_ngs.core.picard
import viral_ngs.core.file

# Then use as:
viral_ngs.core.samtools.SamtoolsTool()
viral_ngs.core.file.mkstempfname()

# ACCEPTABLE: from imports when accessing many items from same module
from viral_ngs.core.misc import available_cpu_count

# Within core/ modules, relative imports are required to avoid circular imports
from . import samtools, picard
from .file import mkstempfname

# For command modules at viral_ngs/ level (illumina.py, read_utils.py, etc.)
from .core import samtools, picard
from .core import file as util_file, misc as util_misc, cmd as util_cmd

# Old patterns (no longer used)
import util.file          →  import viral_ngs.core.file
from tools import samtools →  import viral_ngs.core.samtools
import read_utils         →  import viral_ngs.read_utils  # (read_utils is at viral_ngs level)
```

**Key points:**
- Prefer `import x.y.z` over `from x.y import z` for explicitness
- NO backward compatibility stubs - `viral_ngs.tools` and `viral_ngs.util` don't exist
- Command modules (illumina.py, read_utils.py, etc.) are at `viral_ngs/` level, NOT in `core/`
- Within `core/` modules, use relative imports (`from . import X`) to avoid circular imports
- Command modules import from core using `from .core import X`
- Test imports after each module migration

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
| Import style | `import x.y.z` (not `from x.y import z`) | Explicit, easier to trace |
| Module structure | Single `core/` directory | Simpler than separate tools/util |
| Backward compat | None (no stubs) | Clean codebase, fresh start |
| Python deps | All via conda | Faster install, binary compat |
| x86-only handling | Separate requirements + --x86-only flag | Graceful ARM support |
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

## Lessons Learned from Phase 2

### Circular Import Issues

When consolidating modules, be careful of circular imports:

1. **Import order matters**: In `core/__init__.py`, import `version` before `cmd` because `cmd` uses `version` at module level
2. **Use relative imports in core/**: Within `core/` modules, always use `from . import X` not `from viral_ngs.core import X`
3. **Avoid importing parent package**: Don't use `from viral_ngs import util` inside core modules - this creates cycles

### Conda Package Names

- Python `lz4` package is just `lz4` in conda (not `python-lz4`)
- `lz4-c` is the C library (separate from Python package)
- Always verify package names with `conda search <package>` before adding

### git-filter-repo Usage

Install in a standalone venv (not project dependency):
```bash
python -m venv ~/venvs/git-filter-repo
~/venvs/git-filter-repo/bin/pip install git-filter-repo
source ~/venvs/git-filter-repo/bin/activate
git filter-repo --path-rename ...
```

### x86-Only Packages

For tools without ARM64 support (novoalign, mvicuna):
1. Put in separate requirements file (`core-x86.txt`)
2. Use `install-conda-deps.sh --x86-only` flag
3. Script auto-detects architecture and skips on ARM
4. Tool wrapper handles missing tool at runtime, not import time

These packages are in bioconda but only have x86_64 builds. Use exact versions (e.g., `novoalign=3.09.04`, `mvicuna=1.0`) since newer versions may not exist.

### udocker Quirks

- udocker refuses to run as root by default. Use `--allow-root` flag for Docker build verification
- In Dockerfile, use `udocker --allow-root version` to verify installation
- The `MINIWDL__SCHEDULER__CONTAINER_BACKEND=udocker` environment variable tells miniwdl to use udocker

### qsv is x86-Only

The `qsv` package (fast CSV tool) is only available for linux-64, osx-64, win-64 - no ARM64 build. For multi-arch Docker builds, either skip it or add to an x86-only requirements file.

### PrexistingUnixCommand is the Only InstallMethod

Since all tools are installed via conda, `PrexistingUnixCommand` is the only `InstallMethod` subclass needed. It just checks if the tool exists in PATH.

### Avoid Merging GitHub Actions Workflows

When using `git filter-repo` to import repositories with full history, **exclude `.github/workflows/`** from the merge. Legacy workflows from source repos will trigger on the new repo and fail (wrong paths, missing secrets, different structure).

**Solution:** After merging, immediately delete any imported workflow files:
```bash
rm -rf .github/workflows/*.yml  # Delete imported workflows
git checkout HEAD -- .github/workflows/docker.yml  # Restore our workflows
```

Or use `--path-rename` to move them to a non-executable location during filter-repo.

### Be Selective with git filter-repo Imports

When using `git filter-repo` to import repositories, be **very selective** about what gets included. Don't just import everything and delete later - use explicit `--path` options to only include what you need:

```bash
# BAD: Import everything, delete unwanted files after
git filter-repo --path-rename 'util/:src/viral_ngs/core/'

# GOOD: Only import specific paths you actually want
git filter-repo \
    --path util/ \
    --path tools/ \
    --path '*.py' \
    --path-rename 'util/:src/viral_ngs/core/' \
    --path-rename 'tools/:src/viral_ngs/core/' \
    --force
```

**Files commonly NOT wanted from source repos:**
- `.coveralls.yml`, `.travis.yml` - obsolete CI configs
- `.flake8`, `.pylintrc`, `.style.yapf` - linter configs (consolidate to pyproject.toml if needed)
- `.landscape.yaml` - defunct Landscape.io service
- `github_actions_ci/` - legacy CI scripts (replaced by modern actions)
- `.github/workflows/*.yml` - legacy workflows that will fail in new context

**Cleanup applied in Phase 2:**
- Deleted 5 obsolete dotfiles
- Deleted 13 legacy CI scripts in `github_actions_ci/`

## Lessons Learned from Phase 3a

### Module Placement in Monorepo

Non-core modules have a specific structure:
- **Command modules** (e.g., `assembly.py`) go at `src/viral_ngs/` level, NOT inside `core/`
- **Tool wrappers** go in subpackages: `src/viral_ngs/assemble/`, `src/viral_ngs/phylo/`, etc.
- This keeps `core/` for truly core shared functionality

Example structure for assemble:
```
src/viral_ngs/
├── assembly.py              # Command module (at viral_ngs level)
└── assemble/                # Tool wrappers subpackage
    ├── __init__.py
    ├── spades.py
    ├── mummer.py
    └── ...
```

### Bulk Import Updates with sed

For updating many imports across files, use sed:
```bash
# Update imports in all Python files
sed -i '' 's/from util import/import viral_ngs.core./g' src/viral_ngs/**/*.py
sed -i '' 's/import util\./import viral_ngs.core./g' src/viral_ngs/**/*.py
```

**Key points:**
- Use `-i ''` on macOS (or `-i` on Linux) for in-place editing
- Always review changes with `git diff` after bulk updates
- Some complex patterns may need manual fixes

### Test Organization for Efficient CI

Tests are organized by module subdirectory:
```
tests/unit/
├── core/           # Core module tests
├── assemble/       # Assembly module tests
├── phylo/          # Phylo module tests
└── classify/       # Classify module tests
```

**Cascading test dependencies:**
- Core changes → run ALL tests
- assemble changes → run assemble + core tests
- phylo changes → run phylo + core tests

### Docker Build Tiers

The GitHub Actions workflow uses a three-tier build structure:
1. **build-baseimage** - Base image with conda/micromamba
2. **build-derivatives** - Images that depend on baseimage (core)
3. **build-on-core** - Images that depend on core (assemble, classify, phylo)

Each tier `needs:` the previous one to ensure proper dependency order.

### git filter-repo Path Rewrites for Non-Core Modules

When importing non-core modules, use path rewrites that reflect the correct structure:
```bash
git filter-repo \
    --path-rename 'assembly.py:src/viral_ngs/assembly.py' \
    --path-rename 'assemble/:src/viral_ngs/assemble/' \
    --path-rename 'test/unit/:tests/unit/assemble/' \
    --tag-rename '':'assemble-' \
    --force
```

**Important:** Command modules go to `src/viral_ngs/`, NOT `src/viral_ngs/core/`.

### Creating Package __init__.py for Tool Wrapper Subpackages

Each tool wrapper subpackage needs an `__init__.py` that exports its modules:
```python
"""
viral_ngs.assemble - Tool wrappers for genome assembly.
"""
from . import gap2seq
from . import mafft
from . import mummer
from . import muscle
from . import skani
from . import spades
from . import vcf
from . import wgsim
```

## Communication

When working on this migration:

1. **Ask before major decisions**: If you encounter an architectural choice not covered here, ask
2. **Report blockers immediately**: Especially ARM64 package issues or conda conflicts
3. **Test incrementally**: Don't try to do everything at once
4. **Document surprises**: If something doesn't work as expected, note it for future reference

## Next Steps

1. Read `MONOREPO_IMPLEMENTATION_PLAN.md` for the detailed task list
2. **Phase 3b**: Migrate viral-phylo with git history preservation using `git filter-repo`
   - Command modules: `interhost.py`, `intrahost.py`, `ncbi.py` → `src/viral_ngs/`
   - Tool wrappers: `phylo/` → `src/viral_ngs/phylo/`
   - Tests: `test/unit/` → `tests/unit/phylo/`
3. **Phase 3c**: Migrate viral-classify with git history preservation
   - Command modules: `metagenomics.py`, `taxon_filter.py`, `kmer_utils.py` → `src/viral_ngs/`
   - Tool wrappers: `classify/` → `src/viral_ngs/classify/`
   - Special: Retire kaiju/diamond, try to consolidate kraken2/krona to main env
4. Work through phases sequentially
5. Verify each phase before moving to the next

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

### Phase 2 (Complete)
- Imported viral-core with full git history (3,716 commits, all tags with `core-` prefix)
- Consolidated `tools/` and `util/` directories into single `core/` directory
- Updated all imports to `viral_ngs.core.*` pattern (no backward compat stubs)
- Moved ALL Python dependencies to conda (pyproject.toml has empty deps)
- Created `docker/requirements/core-x86.txt` for x86-only packages (novoalign, mvicuna)
- Updated `install-conda-deps.sh` with `--x86-only` flag for architecture-specific packages
- Created `docker/Dockerfile.core` with verification checks
- Updated all test imports to use `viral_ngs.core.*`
- Deleted legacy viral-core CI workflow (build.yml) that was accidentally merged with git history
- Docker build verified with all module imports working

### Baseimage Enhancements (Post-Phase 2)
- Created `docker/requirements/baseimage.txt` with general utilities
- Added to baseimage: miniwdl, udocker, awscli, google-cloud-storage, csvkit, jq, parallel, pigz, unzip, zstd
- Set `MINIWDL__SCHEDULER__CONTAINER_BACKEND=udocker` environment variable
- Moved general utilities from core.txt to baseimage.txt
- Added seaborn to core.txt for data visualization
- Updated pyproject.toml with Python dependencies for pip-based installs
- Note: qsv is x86-only (no ARM64 build), excluded for multi-arch support

### Phase 3a: viral-assemble (Complete)
- Imported viral-assemble with full git history (3,448 commits, tags with `assemble-` prefix)
- Structure: `assembly.py` at `src/viral_ngs/assembly.py`, tool wrappers in `src/viral_ngs/assemble/`
- Updated all imports to `viral_ngs.core.*` and `viral_ngs.assemble.*` pattern
- Created `docker/requirements/assemble.txt` with assembly tools (SPAdes, MUMmer, MAFFT, etc.)
- Created `docker/Dockerfile.assemble` building on core image
- Tests organized in `tests/unit/assemble/`, core tests moved to `tests/unit/core/`
- Added `build-on-core` job to GitHub Actions for images depending on core
- Imported `scripts/fasta-trim-terminal-ambigs.pl` Perl script

## Reference Repositories

When implementing Docker patterns, refer to:
- `broadinstitute/read-qc-tools` - Simple micromamba-based Dockerfile pattern
- `broadinstitute/viral-baseimage` branch `ct-reduce-google-sdk-install-size` - PR #29 for micromamba migration ideas
- `broadinstitute/viral-core` - Current Dockerfile and install scripts for reference
