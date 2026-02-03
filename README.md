# viral-ngs

[![Build and Test](https://github.com/broadinstitute/viral-ngs/actions/workflows/docker.yml/badge.svg)](https://github.com/broadinstitute/viral-ngs/actions/workflows/docker.yml)
[![codecov](https://codecov.io/gh/broadinstitute/viral-ngs/graph/badge.svg)](https://codecov.io/gh/broadinstitute/viral-ngs)
[![Documentation Status](https://readthedocs.org/projects/viral-ngs/badge/?version=latest)](https://viral-ngs.readthedocs.io/en/latest/?badge=latest)
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

### Conda + pip (Development)

For local development with bioinformatics tools:

1. **Clone the repository:**
   ```bash
   git clone https://github.com/broadinstitute/viral-ngs.git
   cd viral-ngs
   ```

2. **Create and activate a conda environment:**
   ```bash
   conda create -n viral-ngs python=3.12
   conda activate viral-ngs
   ```

   Or using micromamba (faster):
   ```bash
   micromamba create -n viral-ngs python=3.12
   micromamba activate viral-ngs
   ```

3. **Install bioinformatics tools via conda:**
   ```bash
   # Core tools only
   conda install -c conda-forge -c bioconda \
     --file docker/requirements/baseimage.txt \
     --file docker/requirements/core.txt

   # Or for all tools (mega)
   conda install -c conda-forge -c bioconda \
     --file docker/requirements/baseimage.txt \
     --file docker/requirements/core.txt \
     --file docker/requirements/assemble.txt \
     --file docker/requirements/classify.txt \
     --file docker/requirements/phylo.txt
   ```

4. **Install the viral-ngs Python package:**
   ```bash
   pip install -e .
   ```

   Note: Installing into an activated conda environment is safe - pip installs into the conda environment, not system Python.

5. **Verify installation:**
   ```bash
   read_utils --version
   python -c "from viral_ngs.core import samtools; print(samtools.SamtoolsTool().version())"
   ```

### pip (Python package only)

For use as a Python library without bioinformatics tools:

```bash
pip install viral-ngs
```

Note: pip installation does not include external tools (samtools, bwa, etc.). Use Docker or Conda for the complete toolset.

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
