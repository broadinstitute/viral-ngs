# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

See [AGENTS.md](AGENTS.md) for comprehensive project context and development guidelines.

## Quick Reference

- **Docker-centric development**: Run tests inside containers, not on host
- **Import pattern**: `from viral_ngs import core` then `core.samtools.SamtoolsTool()`
- **Test location**: `tests/unit/<module>/`
- **Dependencies**: ALL via conda, not pip (see `docker/requirements/*.txt`)

## Running Tests

```bash
docker run --rm \
  -v $(pwd):/opt/viral-ngs/source \
  quay.io/broadinstitute/viral-ngs:main-core \
  pytest -rsxX -n auto /opt/viral-ngs/source/tests/unit
```

## Key Files

| File | Purpose |
|------|---------|
| [AGENTS.md](AGENTS.md) | Full AI assistant guidance |
| [pyproject.toml](pyproject.toml) | Package configuration |
| [docker/](docker/) | Dockerfiles and requirements |
| [src/viral_ngs/](src/viral_ngs/) | Source code |
| [tests/](tests/) | Test files |
