# viral-ngs

[![Build Docker Images](https://github.com/broadinstitute/viral-ngs/actions/workflows/docker.yml/badge.svg)](https://github.com/broadinstitute/viral-ngs/actions/workflows/docker.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Consolidated tools for viral NGS (Next-Generation Sequencing) data analysis.

## Overview

viral-ngs provides a comprehensive suite of bioinformatics tools for viral genomics:

- **Core utilities**: Read manipulation, Illumina demultiplexing, file handling, quality control
- **Assembly**: Genome assembly, scaffolding, gap filling (SPAdes, MUMmer, MAFFT)
- **Classification**: Metagenomic classification, taxonomy filtering (Kraken2, BLAST, BMTagger)
- **Phylogenetics**: Variant calling, consensus generation, annotation (LoFreq, SnpEff, MUSCLE)

## Docker Images

Pre-built Docker images are available on [Quay.io](https://quay.io/repository/broadinstitute/viral-ngs) and [GitHub Container Registry](https://github.com/broadinstitute/viral-ngs/pkgs/container/viral-ngs):

| Image | Description |
|-------|-------------|
| `viral-ngs:latest` | All tools (mega image) |
| `viral-ngs:core` | Core utilities only |
| `viral-ngs:assemble` | Core + assembly tools |
| `viral-ngs:classify` | Core + classification tools |
| `viral-ngs:phylo` | Core + phylogenetics tools |

### Quick Start

```bash
# Pull the latest image
docker pull quay.io/broadinstitute/viral-ngs:latest

# Run interactively
docker run -it --rm quay.io/broadinstitute/viral-ngs:latest

# Mount local data directory
docker run -it --rm -v $(pwd)/data:/data quay.io/broadinstitute/viral-ngs:latest
```

### Version Tags

```bash
# Specific version
docker pull quay.io/broadinstitute/viral-ngs:2.6.0-core

# Main branch
docker pull quay.io/broadinstitute/viral-ngs:main-classify

# GitHub Container Registry (alternative)
docker pull ghcr.io/broadinstitute/viral-ngs:latest
```

## Installation

### Docker (Recommended)

The recommended way to use viral-ngs is via Docker images, which include all bioinformatics tool dependencies pre-configured.

### pip (Python package only)

For use as a Python library without bioinformatics tools:

```bash
pip install viral-ngs
```

Note: pip installation does not include external tools (samtools, bwa, etc.). Use Docker for the complete toolset.

## Documentation

- **Command-line documentation**: https://viral-ngs.readthedocs.org/
- **Developer guide**: [AGENTS.md](AGENTS.md)
- **Higher-level pipelines**: https://github.com/broadinstitute/viral-pipelines

## Development

### Running Tests

```bash
docker run --rm \
  -v $(pwd):/opt/viral-ngs/source \
  quay.io/broadinstitute/viral-ngs:main-core \
  pytest -rsxX -n auto /opt/viral-ngs/source/tests/unit
```

### Building Docker Images Locally

```bash
# Build baseimage
docker build -t viral-ngs:baseimage -f docker/Dockerfile.baseimage .

# Build core (requires baseimage)
docker build --build-arg BASEIMAGE=viral-ngs:baseimage \
  -t viral-ngs:core -f docker/Dockerfile.core .

# Build derivatives (require core)
docker build --build-arg BASEIMAGE=viral-ngs:core \
  -t viral-ngs:classify -f docker/Dockerfile.classify .
```

See [AGENTS.md](AGENTS.md) for comprehensive development documentation.

## Architecture

### Multi-Architecture Support

Images are built for both `linux/amd64` and `linux/arm64` (Apple Silicon, ARM servers).

Some tools (novoalign, mvicuna) are x86-only and will be skipped on ARM builds.

### Module Structure

```
src/viral_ngs/
├── core/           # Core utilities and tool wrappers
├── assemble/       # Assembly tool wrappers
├── classify/       # Classification tool wrappers
└── phylo/          # Phylogenetics tool wrappers
```

## Related Projects

- [viral-pipelines](https://github.com/broadinstitute/viral-pipelines) - WDL workflows using viral-ngs
- [Terra](https://terra.bio/) - Cloud platform for running viral-pipelines

## License

MIT License - see [LICENSE](LICENSE)

## Contributing

Contributions are welcome! Please see [AGENTS.md](AGENTS.md) for development guidelines.

## Citation

If you use viral-ngs in your research, please cite:

> Broad Institute viral-ngs: Tools for viral NGS data analysis.
> https://github.com/broadinstitute/viral-ngs
