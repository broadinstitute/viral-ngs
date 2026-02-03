"""
viral_ngs.assemble - Tool wrappers for genome assembly.

Contains wrappers for assembly-related bioinformatics tools:
- SPAdes, Gap2Seq, MAFFT, MUMmer, Muscle, etc.
"""

# Import tool wrappers for convenient access
from . import gap2seq
from . import mafft
from . import mummer
from . import muscle
from . import skani
from . import spades
from . import vcf
from . import wgsim
