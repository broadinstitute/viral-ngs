"""
viral-ngs: Consolidated tools for viral NGS data analysis.

This package provides utilities for viral genome sequencing data analysis,
including read processing, assembly, classification, and phylogenetics.
"""

try:
    from importlib.metadata import version, PackageNotFoundError
    try:
        __version__ = version("viral-ngs")
    except PackageNotFoundError:
        __version__ = "0.0.0.dev0"
except ImportError:
    __version__ = "0.0.0.dev0"
