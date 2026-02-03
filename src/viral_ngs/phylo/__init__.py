"""
viral_ngs.phylo - Phylogenetic analysis tools and utilities.

This package contains tool wrappers and utilities for:
- Multiple sequence alignment (MAFFT, MUMmer, MUSCLE)
- Variant calling and annotation (V-Phaser2, SnpEff)
- GenBank feature table manipulation
- VCF file utilities
"""
from . import feature_table
from . import feature_table_types
from . import genbank
from . import mafft
from . import mummer
from . import muscle
from . import snpeff
from . import vcf
from . import vphaser2
