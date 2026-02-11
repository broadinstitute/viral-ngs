---
phase: 01-pre-implementation-validation
plan: 01
subsystem: infra
tags: [conda, genomad, docker, bioconda, tensorflow, multi-arch]

# Dependency graph
requires: []
provides:
  - geNomad conda dependency validated in classify.txt (x86_64 definitive, ARM64 theoretical)
  - Dependency placement decision documented (cross-platform in classify.txt)
  - Docker build proof that geNomad 1.11.2 + TF 2.19.1 + existing classify deps resolve together
affects: [02-tool-wrapper, 04-docker-integration]

# Tech tracking
tech-stack:
  added: [genomad>=1.11.0 (resolves to 1.11.2), tensorflow 2.19.1, keras 3.12.0, xgboost 3.1.3, mmseqs2 18.8cc5c]
  patterns: [cross-platform conda dependency placement in classify.txt]

key-files:
  created: []
  modified:
    - docker/requirements/classify.txt
    - .planning/PROJECT.md
    - .planning/STATE.md

key-decisions:
  - "geNomad placed in classify.txt (cross-platform), not classify-x86.txt; x86_64 validated via Docker build, ARM64 packages available but untested"
  - "ARM64 test skip decorators not needed -- all key dependencies (genomad, tensorflow, mmseqs2, xgboost) have linux-aarch64 packages"

patterns-established:
  - "Single conda dependency line for tools with complex transitive deps (no explicit tensorflow/xgboost/mmseqs2 lines)"

# Metrics
duration: 6min
completed: 2026-02-11
---

# Phase 1 Plan 1: Pre-Implementation Validation Summary

**geNomad 1.11.2 conda dependency validated on x86_64 via Docker build with TF 2.19.1, Keras 3.12.0, xgboost 3.1.3; ARM64 packages confirmed available via pixi search**

## Performance

- **Duration:** 6 min
- **Started:** 2026-02-11T22:29:11Z
- **Completed:** 2026-02-11T22:35:43Z
- **Tasks:** 2
- **Files modified:** 3

## Accomplishments

- Added `genomad>=1.11.0` to `docker/requirements/classify.txt` and validated x86_64 Docker build succeeds
- Verified geNomad 1.11.2 is installed and executable in the container (`genomad --version` outputs `geNomad, version 1.11.2`)
- Confirmed ARM64 package availability for all key dependencies (genomad, tensorflow, mmseqs2, xgboost) via pixi search
- Documented dependency placement decision in PROJECT.md and resolved ARM64 availability blocker in STATE.md

## Task Commits

Each task was committed atomically:

1. **Task 1: Add geNomad dependency and validate x86_64 Docker build** - `cbd0e7f3` (feat)
2. **Task 2: Document dependency placement decision** - `f5e48593` (docs)

## Files Created/Modified

- `docker/requirements/classify.txt` - Added `genomad>=1.11.0` in alphabetical order between blast and kma
- `.planning/PROJECT.md` - Added geNomad dependency placement decision and ARM64 decision to Key Decisions table
- `.planning/STATE.md` - Updated to reflect Phase 1 completion, resolved ARM64 blocker, advanced to Phase 2

## Decisions Made

- **geNomad in classify.txt (cross-platform):** All ARM64 dependencies confirmed available. No need for x86-only placement.
- **ARM64 test skip decorators deferred:** Not needed since ARM64 packages exist. Deferred final decision to Phase 2/4 when actual ARM64 testing can be evaluated.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Used pre-built remote Docker images for baseimage and core**
- **Found during:** Task 1 (Docker build chain)
- **Issue:** Local baseimage build failed due to pre-existing `miniwdl` / `pkg_resources` incompatibility in Dockerfile.baseimage (unrelated to geNomad)
- **Fix:** Pulled `quay.io/broadinstitute/viral-ngs:main-baseimage` and `quay.io/broadinstitute/viral-ngs:main-core` from registry instead of building locally
- **Files modified:** None (runtime Docker pull only)
- **Verification:** Classify image built successfully on top of remote core image
- **Committed in:** N/A (no file changes)

**2. [Rule 3 - Blocking] Added VERSION build-arg to classify Docker build**
- **Found during:** Task 1 (Docker build chain)
- **Issue:** `pip install .` failed inside Docker due to missing `SETUPTOOLS_SCM_PRETEND_VERSION` env var (no git history in build context)
- **Fix:** Added `--build-arg VERSION=0.0.0-dev` to the docker build command
- **Files modified:** None (build command change only)
- **Verification:** Build completed successfully with version arg set
- **Committed in:** N/A (no file changes)

**3. [Rule 3 - Blocking] Used pixi instead of conda for ARM64 package search**
- **Found during:** Task 1 (ARM64 availability check)
- **Issue:** `conda` is not installed on the host machine
- **Fix:** Used `pixi search` (available on host) to check ARM64 package availability, achieving the same result
- **Files modified:** None
- **Verification:** All packages found with ARM64 builds
- **Committed in:** N/A (no file changes)

---

**Total deviations:** 3 auto-fixed (all Rule 3 - blocking issues)
**Impact on plan:** All deviations were environmental workarounds. No scope change. The core validation (Docker build + ARM64 package check) was completed as planned with equivalent tools.

## Issues Encountered

- Pre-existing `miniwdl`/`pkg_resources` issue in Dockerfile.baseimage prevents local baseimage builds. Not blocking for this plan (used remote images), but may need fixing for CI/CD in Phase 4.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- classify.txt has geNomad dependency, Docker build validated on x86_64
- Ready to implement tool wrapper in Phase 2 (`src/viral_ngs/classify/genomad.py`)
- No blockers for Phase 2

## Self-Check: PASSED

- All 4 files exist (classify.txt, PROJECT.md, STATE.md, SUMMARY.md)
- Both task commits found (cbd0e7f3, f5e48593)
- All expected content verified (genomad>=1.11.0, x86_64 validated decision, Phase 1 completion)

---
*Phase: 01-pre-implementation-validation*
*Completed: 2026-02-11*
