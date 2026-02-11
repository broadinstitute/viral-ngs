# viral-ngs: geNomad Integration

## What This Is

Adding geNomad support to the viral-ngs bioinformatics toolkit. geNomad is a tool for identifying viruses and plasmids in genomic and metagenomic data. This integration will follow the existing classify tool patterns (Kraken2, BLAST, diamond) to provide end-to-end viral discovery capabilities in metagenomic assemblies.

## Core Value

Researchers can identify and classify viral sequences in metagenomic assemblies using geNomad's integrated annotation and taxonomic assignment, with results preserved in viral-ngs standard formats.

## Requirements

### Validated

<!-- Existing viral-ngs capabilities -->

- ✓ Metagenomic classification with Kraken2 — existing
- ✓ Sequence similarity search with BLAST — existing
- ✓ Tool wrapper pattern (Tool base class, InstallMethod) — existing
- ✓ CLI command registration via __commands__ list — existing
- ✓ Docker-centric development and testing — existing
- ✓ Conda-based dependency management — existing
- ✓ pytest-based unit testing framework — existing

### Active

<!-- geNomad integration scope -->

- [ ] geNomad tool wrapper class in `src/viral_ngs/classify/genomad.py`
- [ ] Single 'run' command for end-to-end geNomad execution
- [ ] Support for assembled contigs (FASTA input)
- [ ] Batch processing support for multiple FASTA files
- [ ] Preserve all geNomad outputs (virus/plasmid summaries, annotated sequences, taxonomy)
- [ ] User-provided database path (external database management)
- [ ] Sensible parameter defaults (minimal parameter exposure initially)
- [ ] Unit tests in `tests/unit/classify/test_genomad.py` with test data
- [x] Conda dependency in `docker/requirements/classify.txt` (genomad>=1.11.0, x86_64 validated)
- [ ] CLI integration in `src/viral_ngs/metagenomics.py` or dedicated module

### Out of Scope

- Custom database downloading/setup — users manage databases externally
- Advanced parameter tuning — start with defaults, can expand later
- Integration with other viral-ngs pipelines — focus on standalone functionality first
- Real-time processing or streaming — batch processing only
- GUI or web interface — CLI only

## Context

**Existing Architecture:**
- viral-ngs uses a modular tool-wrapper architecture with domain-driven organization
- Classification tools live in `src/viral_ngs/classify/` (kraken2, blast, diamond, etc.)
- Tool wrappers inherit from `core.Tool` base class with `InstallMethod` injection
- CLI commands registered via module-level `__commands__` list
- All dependencies installed via conda for reproducibility
- Docker-centric workflow with multi-arch support (x86_64, ARM64)

**geNomad Background:**
- End-to-end tool for virus/plasmid identification in genomic data
- Requires pre-downloaded database (user-managed in this integration)
- Produces multiple output files: summaries (TSV), annotated sequences (FASTA), taxonomy assignments
- Primary use case: viral discovery in metagenomic assemblies
- Documentation: https://portal.nersc.gov/genomad/quickstart.html

**Development Environment:**
- Tests run inside Docker containers: `quay.io/broadinstitute/viral-ngs:main-classify`
- Import pattern: `import viral_ngs.core as core` then `core.genomad.GenomadTool()`
- Test location: `tests/unit/classify/`
- Test framework: pytest with parallel execution (`pytest -rsxX -n auto`)

## Constraints

- **Tech Stack**: Python 3.10+, must use conda for all dependencies
- **Docker-Centric**: All testing must run in Docker containers, not on host
- **Pattern Consistency**: Follow existing tool wrapper patterns (see kraken2.py, blast.py)
- **Import Conventions**: Use `import viral_ngs.core as core` pattern, relative imports within packages
- **No pip Dependencies**: ALL runtime dependencies via conda (check `conda search -c bioconda genomad`)
- **Test Coverage**: Tests must pass in Docker with test data
- **Multi-arch**: Consider ARM64 compatibility (mark x86-only tests if needed)
- **Incremental**: Focus on minimal working implementation first

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| Add geNomad to classify.txt (cross-platform) | x86_64 validated via Docker build (geNomad 1.11.2, TF 2.19.1, xgboost 3.1.3); ARM64 packages available via pixi search (untested) | ✓ x86_64 validated |
| ARM64 test skip decorators for geNomad | Not needed - ARM64 packages exist (though untested due to x86_64-only base image) | Deferred to Phase 2/4 |
| Single 'run' command vs granular commands | User wants simplicity for end-to-end viral discovery workflow | — Pending |
| User-managed database path | External database management is simpler than auto-download, consistent with other tools | — Pending |
| Start with default parameters | Minimize complexity initially, can expose more controls later based on feedback | — Pending |
| Add to metagenomics.py CLI module | geNomad is classification tool, fits existing metagenomics command structure | — Pending |
| Use feature branch | Isolate development from main branch, enable PR-based review | ✓ Good |

---
*Last updated: 2026-02-11 after Phase 1 validation*
