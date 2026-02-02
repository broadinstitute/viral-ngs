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

## Phase 0: Prepare Existing viral-ngs Repository

**Goal:** Clean up the dormant repo for reuse.

### Tasks

1. **Audit current state**
   - [ ] Review any open issues/PRs (close or migrate if relevant)
   - [ ] Check branch structure (identify branches to archive/delete)
   - [ ] Review any existing GitHub Actions workflows
   - [ ] Check for secrets/settings that need updating

2. **Archive old content**
   ```bash
   git tag archive/legacy-monolith
   git branch archive/legacy
   git push origin archive/legacy-monolith archive/legacy
   ```

3. **Update repository settings**
   - [ ] Set default branch to `main`
   - [ ] Update branch protection rules
   - [ ] Add secrets: `QUAY_USERNAME`, `QUAY_TOKEN`, `QUAY_API_TOKEN`
   - [ ] Enable GitHub Packages (ghcr.io)
   - [ ] Update repository description/topics

---

## Phase 1: Set Up Monorepo Foundation

**Goal:** Establish fresh structure with modern tooling.

### Directory Structure

```
viral-ngs/
├── pyproject.toml
├── src/viral_ngs/
│   ├── __init__.py
│   └── py.typed
├── docker/
│   ├── Dockerfile.baseimage
│   ├── Dockerfile.core
│   ├── Dockerfile.classify
│   ├── Dockerfile.assemble
│   ├── Dockerfile.phylo
│   ├── Dockerfile.mega
│   ├── install-conda-deps.sh
│   └── requirements/
│       ├── core.txt
│       ├── classify.txt
│       ├── assemble.txt
│       └── phylo.txt
├── tests/
│   ├── conftest.py
│   ├── unit/
│   └── input/
├── .github/workflows/
│   ├── docker.yml
│   ├── test.yml
│   └── cleanup.yml
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

4. **Create .github/workflows/cleanup.yml**
   ```yaml
   name: Cleanup Docker Tags

   on:
     delete:
       branches: ['**']

   jobs:
     cleanup:
       if: github.event.ref_type == 'branch'
       runs-on: ubuntu-latest
       steps:
         - name: Delete Quay.io tags
           env:
             BRANCH: ${{ github.event.ref }}
             QUAY_TOKEN: ${{ secrets.QUAY_API_TOKEN }}
           run: |
             for flavor in "" "-baseimage" "-core" "-classify" "-assemble" "-phylo"; do
               TAG="${BRANCH}${flavor}"
               curl -X DELETE \
                 -H "Authorization: Bearer $QUAY_TOKEN" \
                 "https://quay.io/api/v1/repository/broadinstitute/viral-ngs/tag/$TAG" \
                 || true
             done
   ```

5. **Verify baseimage builds on both architectures**

---

## Phase 2: Migrate viral-core

**Goal:** Get core functionality working in new repo.

### Git History Preservation

```bash
# Install git-filter-repo
pip install git-filter-repo

# Clone viral-core and rewrite paths
git clone https://github.com/broadinstitute/viral-core.git viral-core-rewrite
cd viral-core-rewrite

# Rewrite paths: *.py -> src/viral_ngs/*.py, util/ -> src/viral_ngs/util/, etc.
git filter-repo --path-rename 'util/:src/viral_ngs/util/' \
                --path-rename 'tools/:src/viral_ngs/tools/' \
                --tag-rename '':'core-' \
                --force

# Merge into monorepo
cd /path/to/viral-ngs
git remote add viral-core ../viral-core-rewrite
git fetch viral-core
git merge viral-core/master --allow-unrelated-histories -m "Import viral-core with history"
git remote remove viral-core
```

### Import Changes Required

All imports must change from flat to package style:
```python
# Before
import util.file
from tools import samtools

# After
from viral_ngs.util import file
from viral_ngs.tools import samtools
```

Use this regex for find/replace:
```
# For "import util.X"
s/import util\./from viral_ngs.util import /g

# For "from util.X import Y"
s/from util\./from viral_ngs.util./g

# For "import tools.X"
s/import tools\./from viral_ngs.tools import /g
```

### Tasks

1. [ ] Import viral-core with history preservation
2. [ ] Update all imports to `viral_ngs.*` style
3. [ ] Create `docker/requirements/core.txt` from `requirements-conda.txt`
4. [ ] Create `docker/Dockerfile.core`
5. [ ] Migrate tests to `tests/unit/`
6. [ ] Verify: build image, run all core tests
7. [ ] Set up CodeCov integration

---

## Phase 3: Migrate Derivative Modules

### Phase 3a: viral-assemble

```bash
git clone https://github.com/broadinstitute/viral-assemble.git viral-assemble-rewrite
cd viral-assemble-rewrite
git filter-repo --path-rename 'assembly.py:src/viral_ngs/assembly.py' \
                --path-rename 'assemble/:src/viral_ngs/assemble/' \
                --tag-rename '':'assemble-' \
                --force
```

Tasks:
1. [ ] Import with history
2. [ ] Update imports
3. [ ] Create `docker/requirements/assemble.txt`
4. [ ] Create `docker/Dockerfile.assemble`
5. [ ] Migrate tests
6. [ ] Verify build and tests

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

4. **Optional: Create SKILLS.md** for custom Claude Code commands

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
- **Delete:** When branch is deleted on GitHub

---

## Verification Checklist

### Per-Phase

| Phase | Verification |
|-------|-------------|
| 0 | Old content archived, settings updated, secrets configured |
| 1 | baseimage builds on amd64+arm64, pushes to quay.io+ghcr.io |
| 2 | core builds, all core tests pass, `pip install viral-ngs` works |
| 3a | assemble builds, assembly tests pass |
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
              - 'src/viral_ngs/util/**'
              - 'src/viral_ngs/tools/**'
              - 'src/viral_ngs/*.py'
              - 'docker/Dockerfile.core'
              - 'docker/Dockerfile.baseimage'
            classify:
              - 'docker/requirements/classify.txt'
              - 'src/viral_ngs/classify/**'
              - 'src/viral_ngs/metagenomics.py'
              - 'src/viral_ngs/taxon_filter.py'
              - 'docker/Dockerfile.classify'
              - 'docker/requirements/core.txt'
              - 'src/viral_ngs/util/**'
              - 'src/viral_ngs/tools/**'
            assemble:
              - 'docker/requirements/assemble.txt'
              - 'src/viral_ngs/assemble/**'
              - 'src/viral_ngs/assembly.py'
              - 'docker/Dockerfile.assemble'
              - 'docker/requirements/core.txt'
              - 'src/viral_ngs/util/**'
              - 'src/viral_ngs/tools/**'
            phylo:
              - 'docker/requirements/phylo.txt'
              - 'src/viral_ngs/phylo/**'
              - 'src/viral_ngs/interhost.py'
              - 'src/viral_ngs/intrahost.py'
              - 'docker/Dockerfile.phylo'
              - 'docker/requirements/core.txt'
              - 'src/viral_ngs/util/**'
              - 'src/viral_ngs/tools/**'

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
| novoalign | ⚠️ No | Commercial, x86-only |
| gatk=3.8 | ⚠️ Check | Legacy version |
| samtools | ✅ Yes | |
| bwa | ✅ Yes | |
| minimap2 | ✅ Yes | |
| picard | ✅ Yes | |
| kraken2 | ✅ Yes | |
| spades | ✅ Yes | |

If a package lacks ARM64 support:
1. Document it as x86-only
2. Consider building from source in Dockerfile
3. Or skip that tool on ARM builds
