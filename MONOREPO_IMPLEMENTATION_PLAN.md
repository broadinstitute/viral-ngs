# Monorepo Implementation Plan: Consolidating viral-* Repositories

> This document is the implementation plan for consolidating viral-core, viral-classify, viral-assemble, viral-phylo, and viral-baseimage into a single monorepo at `github.com/broadinstitute/viral-ngs`.

## Executive Summary

| Aspect | Decision |
|--------|----------|
| **Target repo** | Reuse dormant `github.com/broadinstitute/viral-ngs` |
| **Default branch** | `main` |
| **Python packaging** | Single pip package (`viral-ngs`) with `viral_ngs.*` imports |
| **Docker registry** | Single repo with flavor tags: `quay.io/broadinstitute/viral-ngs:2.6.0-core` |
| **Multi-arch** | Native parallel builds (ubuntu-24.04 + ubuntu-24.04-arm) |
| **Docker cache** | Registry cache on Quay.io (not GHA - 10GB limit too small) |
| **Code coverage** | CodeCov (replacing Coveralls) |
| **CI scoping** | `dorny/paths-filter` with cascading dependencies |
| **Git history** | Preserve via `git filter-repo` + merge |

---

## Phase 0: Prepare Existing viral-ngs Repository ✅ COMPLETE

**Goal:** Clean up the dormant repo for reuse.

### Tasks

1. **Audit current state** ✅
   - [x] Review any open issues/PRs (closed 53 issues with migration notice)
   - [x] Check branch structure (deleted 148 stale feature branches)
   - [x] Review any existing GitHub Actions workflows (none existed, Travis CI removed)
   - [x] Check for secrets/settings that need updating

2. **Archive old content** ✅
   ```bash
   git tag archive/legacy-monolith
   git branch archive/legacy
   git push origin archive/legacy-monolith archive/legacy
   ```

3. **Update repository settings** ✅
   - [x] Set default branch to `main`
   - [x] Add secrets: `QUAY_USERNAME`, `QUAY_TOKEN`, `CODECOV_TOKEN`
   - [x] Update repository description/topics
   - [x] Close all legacy issues with migration notice
   - [ ] Branch protection rules (deferred to Phase 4)
   - [x] GitHub Packages/ghcr.io (auto-enables on first push, uses GITHUB_TOKEN)
   - [x] Replaced BSD-3-Clause license with MIT

---

## Phase 1: Set Up Monorepo Foundation ✅ COMPLETE

**Goal:** Establish fresh structure with modern tooling.

### Completed Items
- [x] Created `pyproject.toml` with setuptools-scm versioning
- [x] Created `src/viral_ngs/__init__.py` with version detection
- [x] Created `src/viral_ngs/py.typed` marker for type hints
- [x] Created `docker/Dockerfile.baseimage` using `mambaorg/micromamba:2.4.0-ubuntu24.04`
- [x] Created `docker/install-conda-deps.sh` for unified dependency installation
- [x] Created `.github/workflows/docker.yml` for multi-arch builds to Quay.io + ghcr.io
- [x] Verified baseimage builds locally (Python 3.12, conda/mamba symlinks work)

### Directory Structure

```
viral-ngs/
├── pyproject.toml
├── src/viral_ngs/
│   ├── __init__.py
│   ├── py.typed
│   ├── assembly.py             # Assembly commands (Phase 3a)
│   ├── interhost.py, intrahost.py, ncbi.py  # Phylo commands (Phase 3b)
│   ├── metagenomics.py, taxon_filter.py     # Classify commands (Phase 3c)
│   ├── core/                   # Core modules consolidated here
│   │   ├── __init__.py         # Tool/InstallMethod classes + imports
│   │   ├── samtools.py, picard.py, bwa.py, ...  # Tool wrappers
│   │   ├── file.py, misc.py, cmd.py, ...        # Utilities
│   │   ├── read_utils.py, illumina.py, ...      # Command modules
│   │   └── errors.py, priorities.py, ...
│   ├── assemble/               # Assembly tool wrappers (Phase 3a)
│   │   ├── __init__.py
│   │   └── spades.py, mummer.py, mafft.py, ...
│   ├── phylo/                  # Phylo tool wrappers (Phase 3b)
│   └── classify/               # Classify tool wrappers (Phase 3c)
├── docker/
│   ├── Dockerfile.baseimage
│   ├── Dockerfile.core
│   ├── Dockerfile.classify
│   ├── Dockerfile.assemble
│   ├── Dockerfile.phylo
│   ├── Dockerfile.mega
│   ├── install-conda-deps.sh
│   └── requirements/
│       ├── core.txt            # All deps (Python + bioinformatics)
│       ├── core-x86.txt        # x86-only packages (novoalign, mvicuna)
│       ├── classify.txt
│       ├── assemble.txt
│       └── phylo.txt
├── tests/
│   ├── conftest.py
│   ├── unit/
│   │   ├── core/              # Core module tests
│   │   ├── assemble/          # Assembly module tests (Phase 3a)
│   │   ├── phylo/             # Phylo module tests (Phase 3b)
│   │   └── classify/          # Classify module tests (Phase 3c)
│   └── input/
├── scripts/                    # Migration/maintenance scripts
├── .github/workflows/
│   ├── docker.yml
│   └── test.yml
├── docs/
├── README.md
├── CLAUDE.md
├── AGENTS.md
└── SKILLS.md (optional)
```

### Tasks

1. **Create pyproject.toml**
   ```toml
   [build-system]
   requires = ["setuptools>=61.0", "setuptools-scm"]
   build-backend = "setuptools.build_meta"

   [project]
   name = "viral-ngs"
   dynamic = ["version"]
   description = "Tools for viral NGS data analysis"
   readme = "README.md"
   requires-python = ">=3.10"
   dependencies = [
       "biopython",
       "pysam",
       "pandas",
       # ... other runtime deps
   ]

   [project.optional-dependencies]
   dev = ["pytest", "pytest-xdist", "pytest-cov"]

   [project.scripts]
   illumina = "viral_ngs.illumina:main"
   read_utils = "viral_ngs.read_utils:main"
   # ... other CLI entry points

   [tool.setuptools_scm]

   [tool.setuptools.packages.find]
   where = ["src"]
   ```

2. **Create Dockerfile.baseimage**
   ```dockerfile
   FROM mambaorg/micromamba:2.4.0-ubuntu24.04

   ENV MAMBA_DOCKERFILE_ACTIVATE=1
   ENV INSTALL_PATH=/opt/viral-ngs
   ENV VIRAL_NGS_PATH=$INSTALL_PATH/source
   ENV PYTHONPATH=$VIRAL_NGS_PATH

   USER root

   # Symlinks for conda/mamba compatibility
   RUN ln -s /usr/local/bin/micromamba /usr/local/bin/conda && \
       ln -s /usr/local/bin/micromamba /usr/local/bin/mamba

   # Explicit Python install (not included by default in micromamba image)
   RUN micromamba install -y -n base python=3.12 pip && \
       micromamba clean --all --yes

   WORKDIR $VIRAL_NGS_PATH
   CMD ["/bin/bash"]
   ```

3. **Create .github/workflows/docker.yml** (see full workflow in Appendix A)

4. **Feature branch image expiration**

   Instead of a cleanup workflow requiring API tokens, use Quay.io's built-in expiration labels.
   In the Docker workflow, for non-main/non-release builds, append to Dockerfile:
   ```bash
   # For feature branch builds only
   echo "LABEL quay.expires-after=10w" >> docker/Dockerfile.${{ matrix.flavor }}
   ```

   This auto-expires feature branch images after 10 weeks while keeping release and main images permanent.

5. **Verify baseimage builds on both architectures**

---

## Phase 2: Migrate viral-core ✅ COMPLETE

**Goal:** Get core functionality working in new repo.

### What Was Done

1. **Git History Preservation** ✅
   - Installed git-filter-repo in standalone venv (`~/venvs/git-filter-repo`)
   - Cloned viral-core and rewrote paths to `src/viral_ngs/`
   - Merged with `--allow-unrelated-histories` (3,716 commits preserved)
   - All tags renamed with `core-` prefix

2. **Consolidated Structure** ✅
   - **Key change:** Merged `tools/` and `util/` into single `core/` directory
   - All imports now use `viral_ngs.core.*` pattern (e.g., `import viral_ngs.core.samtools`)
   - No backward compatibility stubs - clean imports only
   - Tool base classes (`Tool`, `InstallMethod`, `PrexistingUnixCommand`) in `core/__init__.py`

3. **Conda-Only Dependencies** ✅
   - ALL Python runtime dependencies installed via conda for speed
   - `pyproject.toml` has empty `dependencies = []`
   - All deps listed in `docker/requirements/core.txt`
   - Added scipy (required by priorities.py)

4. **x86-Only Tool Handling** ✅
   - Created `docker/requirements/core-x86.txt` for novoalign/mvicuna
   - Updated `install-conda-deps.sh` with `--x86-only` flag
   - On ARM, x86-only packages are skipped gracefully

5. **Docker Image** ✅
   - Created `docker/Dockerfile.core`
   - Verified build passes with all import checks

### Final Import Pattern

```python
# All modules are in viral_ngs.core.*
import viral_ngs.core.samtools
import viral_ngs.core.picard
import viral_ngs.core.file
import viral_ngs.core.misc

# Within core/ modules, use relative imports
from . import samtools, picard
from .file import mkstempfname
```

### Files Created/Modified

| File | Action |
|------|--------|
| `src/viral_ngs/core/` | New - all modules consolidated here |
| `src/viral_ngs/core/__init__.py` | Tool/InstallMethod classes + submodule imports |
| `docker/Dockerfile.core` | New |
| `docker/requirements/core.txt` | All Python + bioinformatics deps |
| `docker/requirements/core-x86.txt` | New - x86-only packages |
| `docker/install-conda-deps.sh` | Added `--x86-only` flag |
| `pyproject.toml` | Empty dependencies (all via conda) |
| `tests/` | Updated imports to viral_ngs.core.* |

### Removed (No Backward Compat Stubs)

- `src/viral_ngs/tools/` - removed, use `viral_ngs.core.*`
- `src/viral_ngs/util/` - removed, use `viral_ngs.core.*`

---

## Phase 3: Migrate Derivative Modules

### Phase 3a: viral-assemble ✅ COMPLETE

**What Was Done:**

1. **Git History Preservation**
   - Cloned viral-assemble and rewrote paths with git filter-repo
   - Merged with `--allow-unrelated-histories` (3,448 commits preserved)
   - All tags renamed with `assemble-` prefix

2. **Module Structure**
   - `assembly.py` → `src/viral_ngs/assembly.py` (at viral_ngs level, NOT in core/)
   - `assemble/` → `src/viral_ngs/assemble/` (tool wrappers subpackage)
   - Created `src/viral_ngs/assemble/__init__.py` with module exports

3. **Import Updates**
   - All imports updated to `viral_ngs.core.*` and `viral_ngs.assemble.*` pattern
   - Avoided `from x import y` where possible, using full imports instead

4. **Docker Infrastructure**
   - Created `docker/requirements/assemble.txt` (SPAdes, MUMmer, MAFFT, skani, etc.)
   - Created `docker/Dockerfile.assemble` building on core image
   - Added `build-on-core` job to `.github/workflows/docker.yml`

5. **Test Organization**
   - Tests imported to `tests/unit/assemble/`
   - Core tests moved from `tests/unit/test_*.py` to `tests/unit/core/`
   - Imported Perl script: `scripts/fasta-trim-terminal-ambigs.pl`

Tasks:
- [x] Import with history (3,448 commits, tags prefixed with `assemble-`)
- [x] Update imports (`viral_ngs.core.*`, `viral_ngs.assemble.*`)
- [x] Create `docker/requirements/assemble.txt`
- [x] Create `docker/Dockerfile.assemble`
- [x] Migrate tests (to `tests/unit/assemble/`)
- [x] Add `build-on-core` job to GitHub Actions

### Phase 3b: viral-phylo

```bash
git clone https://github.com/broadinstitute/viral-phylo.git viral-phylo-rewrite
cd viral-phylo-rewrite
git filter-repo --path-rename 'interhost.py:src/viral_ngs/interhost.py' \
                --path-rename 'intrahost.py:src/viral_ngs/intrahost.py' \
                --path-rename 'ncbi.py:src/viral_ngs/ncbi.py' \
                --path-rename 'phylo/:src/viral_ngs/phylo/' \
                --tag-rename '':'phylo-' \
                --force
```

Tasks:
1. [ ] Import with history
2. [ ] Update imports
3. [ ] Create `docker/requirements/phylo.txt`
4. [ ] Create `docker/Dockerfile.phylo`
5. [ ] Migrate tests
6. [ ] Verify build and tests

### Phase 3c: viral-classify

**Special considerations:**
- Retire kaiju (env2) and diamond (env3)
- Attempt to consolidate kraken2/krona (env4) into main environment

```bash
git clone https://github.com/broadinstitute/viral-classify.git viral-classify-rewrite
cd viral-classify-rewrite

# Exclude retired tool wrappers if desired
git filter-repo --path-rename 'metagenomics.py:src/viral_ngs/metagenomics.py' \
                --path-rename 'taxon_filter.py:src/viral_ngs/taxon_filter.py' \
                --path-rename 'kmer_utils.py:src/viral_ngs/kmer_utils.py' \
                --path-rename 'classify/:src/viral_ngs/classify/' \
                --tag-rename '':'classify-' \
                --force
```

Tasks:
1. [ ] Import with history (excluding kaiju/diamond wrappers)
2. [ ] Update imports
3. [ ] Create `docker/requirements/classify.txt`
4. [ ] **Test kraken2/krona consolidation** into main env
5. [ ] If consolidation fails, keep minimal multi-env pattern
6. [ ] Create `docker/Dockerfile.classify`
7. [ ] Migrate tests
8. [ ] Verify build and tests

---

## Phase 4: Finalize and Transition

### Documentation Consolidation

1. **Create AGENTS.md** - merge content from:
   - viral-core/CLAUDE.md
   - viral-core/DEVELOPMENT_NOTES.md
   - New monorepo-specific guidance

2. **Create minimal CLAUDE.md**
   ```markdown
   # CLAUDE.md
   See @AGENTS.md for project context and guidelines.
   ```

3. **Consolidate README.md** from all 5 repos

4. **Update README.md badges**
   - GitHub Actions build status badge
   - CodeCov coverage badge
   - Quay.io image badges (for each flavor)
   - License badge (MIT)

5. **Optional: Create SKILLS.md** for custom Claude Code commands

### Add Mega Image

Create `docker/Dockerfile.mega`:
```dockerfile
FROM quay.io/broadinstitute/viral-ngs:${VERSION}-core

COPY docker/requirements/ /tmp/requirements/
COPY docker/install-conda-deps.sh /tmp/

# Install ALL deps together
RUN /tmp/install-conda-deps.sh \
    /tmp/requirements/core.txt \
    /tmp/requirements/classify.txt \
    /tmp/requirements/assemble.txt \
    /tmp/requirements/phylo.txt

RUN python -c "from viral_ngs import assembly, metagenomics, interhost; print('mega OK')"
```

### Update viral-pipelines

1. [ ] Update WDL docker image references to new tag format
2. [ ] Test workflows with new images
3. [ ] Consider transition period with both old/new available

### Configure Branch Protection

1. [ ] Enable branch protection rules on `main`
2. [ ] Require pull request reviews before merging
3. [ ] Require status checks to pass (CI builds, tests)

### Archive Old Repositories

1. [ ] Add deprecation notices to old repos
2. [ ] Point to new monorepo
3. [ ] Archive (make read-only) after transition period

---

## Docker Image Tags

### Tag Format

```
quay.io/broadinstitute/viral-ngs:2.6.0-baseimage
quay.io/broadinstitute/viral-ngs:2.6.0-core
quay.io/broadinstitute/viral-ngs:2.6.0-classify
quay.io/broadinstitute/viral-ngs:2.6.0-assemble
quay.io/broadinstitute/viral-ngs:2.6.0-phylo
quay.io/broadinstitute/viral-ngs:2.6.0              # mega (no suffix)

quay.io/broadinstitute/viral-ngs:main               # main branch mega
quay.io/broadinstitute/viral-ngs:main-core          # main branch core
quay.io/broadinstitute/viral-ngs:latest             # alias for main mega

quay.io/broadinstitute/viral-ngs:feature-branch     # branch mega
quay.io/broadinstitute/viral-ngs:feature-branch-core
```

### Cache Tags

```
quay.io/broadinstitute/viral-ngs:cache-core-amd64
quay.io/broadinstitute/viral-ngs:cache-core-arm64
quay.io/broadinstitute/viral-ngs:cache-classify-amd64
# etc.
```

### Tag Lifecycle

- **Push on:** Tagged releases, branch pushes
- **Skip on:** PRs, merge queue builds
- **Auto-expire:** Feature branch images expire after 10 weeks via `quay.expires-after` label
- **Permanent:** Release tags and main branch images (no expiration label)

---

## Verification Checklist

### Per-Phase

| Phase | Verification |
|-------|-------------|
| 0 | Old content archived, settings updated, secrets configured |
| 1 | baseimage builds on amd64+arm64, pushes to quay.io+ghcr.io |
| 2 | core builds, all core tests pass, `pip install viral-ngs` works |
| 3a | ✅ assemble builds, assembly tests pass |
| 3b | phylo builds, phylo tests pass |
| 3c | classify builds (ideally single env), classify tests pass |
| 4 | mega image works, viral-pipelines WDLs work |

### End-to-End

- [ ] All 6 docker images build on linux/amd64 and linux/arm64
- [ ] `pytest tests/` passes in each image
- [ ] `pip install viral-ngs` works, imports succeed
- [ ] Run sample viral-pipelines workflow with new images
- [ ] CI build times are comparable or faster
- [ ] Images are not significantly larger

---

## Appendix A: Full Docker Workflow

```yaml
# .github/workflows/docker.yml
name: Build Docker Images

on:
  push:
    branches: [main]
    tags: ['v*']
  pull_request:
    branches: [main]

env:
  QUAY_REPO: quay.io/broadinstitute/viral-ngs
  GHCR_REPO: ghcr.io/broadinstitute/viral-ngs

jobs:
  changes:
    runs-on: ubuntu-latest
    outputs:
      core: ${{ steps.filter.outputs.core }}
      classify: ${{ steps.filter.outputs.classify }}
      assemble: ${{ steps.filter.outputs.assemble }}
      phylo: ${{ steps.filter.outputs.phylo }}
    steps:
      - uses: actions/checkout@v4
      - uses: dorny/paths-filter@v3
        id: filter
        with:
          filters: |
            core:
              - 'docker/requirements/core.txt'
              - 'docker/requirements/core-x86.txt'
              - 'src/viral_ngs/core/**'
              - 'docker/Dockerfile.core'
              - 'docker/Dockerfile.baseimage'
            classify:
              - 'docker/requirements/classify.txt'
              - 'src/viral_ngs/classify/**'
              - 'docker/Dockerfile.classify'
              - 'docker/requirements/core.txt'
              - 'src/viral_ngs/core/**'
            assemble:
              - 'docker/requirements/assemble.txt'
              - 'src/viral_ngs/assemble/**'
              - 'docker/Dockerfile.assemble'
              - 'docker/requirements/core.txt'
              - 'src/viral_ngs/core/**'
            phylo:
              - 'docker/requirements/phylo.txt'
              - 'src/viral_ngs/phylo/**'
              - 'docker/Dockerfile.phylo'
              - 'docker/requirements/core.txt'
              - 'src/viral_ngs/core/**'

  build:
    needs: changes
    if: github.event_name == 'push'  # Skip on PRs
    runs-on: ${{ matrix.runner }}
    strategy:
      fail-fast: false
      matrix:
        arch:
          - platform: linux/amd64
            runner: ubuntu-24.04
            suffix: amd64
          - platform: linux/arm64
            runner: ubuntu-24.04-arm
            suffix: arm64
        flavor: [baseimage, core, classify, assemble, phylo, mega]
    steps:
      - uses: actions/checkout@v4
      - uses: docker/setup-buildx-action@v3

      - uses: docker/login-action@v3
        with:
          registry: quay.io
          username: ${{ secrets.QUAY_USERNAME }}
          password: ${{ secrets.QUAY_TOKEN }}

      - uses: docker/login-action@v3
        with:
          registry: ghcr.io
          username: ${{ github.actor }}
          password: ${{ secrets.GITHUB_TOKEN }}

      - name: Build and push by digest
        uses: docker/build-push-action@v6
        id: build
        with:
          context: .
          file: docker/Dockerfile.${{ matrix.flavor }}
          platforms: ${{ matrix.arch.platform }}
          outputs: type=image,name=${{ env.QUAY_REPO }},push-by-digest=true,name-canonical=true,push=true
          cache-from: type=registry,ref=${{ env.QUAY_REPO }}:cache-${{ matrix.flavor }}-${{ matrix.arch.suffix }}
          cache-to: type=registry,ref=${{ env.QUAY_REPO }}:cache-${{ matrix.flavor }}-${{ matrix.arch.suffix }},mode=max

      - name: Export digest
        run: |
          mkdir -p /tmp/digests/${{ matrix.flavor }}
          digest="${{ steps.build.outputs.digest }}"
          touch "/tmp/digests/${{ matrix.flavor }}/${digest#sha256:}"

      - uses: actions/upload-artifact@v4
        with:
          name: digests-${{ matrix.flavor }}-${{ matrix.arch.suffix }}
          path: /tmp/digests/${{ matrix.flavor }}/*

  merge:
    needs: build
    runs-on: ubuntu-latest
    strategy:
      matrix:
        flavor: [baseimage, core, classify, assemble, phylo, mega]
    steps:
      - uses: actions/download-artifact@v4
        with:
          path: /tmp/digests
          pattern: digests-${{ matrix.flavor }}-*
          merge-multiple: true

      - uses: docker/setup-buildx-action@v3

      - uses: docker/login-action@v3
        with:
          registry: quay.io
          username: ${{ secrets.QUAY_USERNAME }}
          password: ${{ secrets.QUAY_TOKEN }}

      - uses: docker/metadata-action@v5
        id: meta
        with:
          images: ${{ env.QUAY_REPO }}
          flavor: |
            # mega image has no suffix
            suffix=${{ matrix.flavor != 'mega' && format('-{0}', matrix.flavor) || '' }}
          tags: |
            type=semver,pattern={{version}}
            type=ref,event=branch
            type=raw,value=${{ matrix.flavor != 'mega' && matrix.flavor || 'latest' }},enable=${{ github.ref == 'refs/heads/main' }}

      - name: Create manifest and push
        working-directory: /tmp/digests
        run: |
          docker buildx imagetools create \
            $(jq -cr '.tags | map("-t " + .) | join(" ")' <<< "$DOCKER_METADATA_OUTPUT_JSON") \
            $(printf '${{ env.QUAY_REPO }}@sha256:%s ' *)

  test:
    needs: merge
    runs-on: ubuntu-latest
    strategy:
      matrix:
        flavor: [core, classify, assemble, phylo]
    steps:
      - uses: actions/checkout@v4
      - name: Run tests in container
        run: |
          docker run --rm \
            -v ${{ github.workspace }}:/workspace \
            ${{ env.QUAY_REPO }}:${{ github.ref_name }}-${{ matrix.flavor }} \
            pytest -n auto /workspace/tests/unit/

  coverage:
    needs: test
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: codecov/codecov-action@v4
        with:
          token: ${{ secrets.CODECOV_TOKEN }}
```

---

## Appendix B: ARM64 Package Availability

Check before building:
```bash
conda search -c bioconda <package> --subdir linux-aarch64
```

| Package | ARM64 Status | Notes |
|---------|-------------|-------|
| novoalign | ❌ No | Commercial, x86-only, in `core-x86.txt` |
| mvicuna | ❌ No | x86-only, in `core-x86.txt` |
| gatk=3.8 | ✅ Yes | Legacy version, works on ARM |
| samtools | ✅ Yes | |
| bwa | ✅ Yes | |
| minimap2 | ✅ Yes | |
| picard | ✅ Yes | |
| kraken2 | ✅ Yes | |
| spades | ✅ Yes | |

### x86-Only Package Handling

For packages that lack ARM64 support:
1. Add to `docker/requirements/core-x86.txt` (or flavor-specific x86 file)
2. Install with `install-conda-deps.sh --x86-only <file>`
3. Script auto-detects architecture and skips on ARM
4. Tool wrapper should handle missing tool gracefully (fail on install, not import)
