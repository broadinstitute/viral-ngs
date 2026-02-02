"""
viral-ngs: Consolidated tools for viral NGS data analysis.

This package provides utilities for viral genome sequencing data analysis,
including read processing, assembly, classification, and phylogenetics.

All modules are in the `core` subpackage:
    import viral_ngs.core.samtools
    import viral_ngs.core.picard
    import viral_ngs.core.file
    import viral_ngs.core.misc
"""

try:
    from importlib.metadata import version, PackageNotFoundError
    try:
        __version__ = version("viral-ngs")
    except PackageNotFoundError:
        __version__ = "0.0.0.dev0"
except ImportError:
    __version__ = "0.0.0.dev0"

# Import core module - all functionality is in viral_ngs.core.*
from . import core
