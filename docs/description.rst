Description of the methods
==========================

Much more documentation to come...

TO DO: here we will put a high level description of the various tools that
exist here, perhaps with some pictures and such. We will describe why we
used certain tools and approaches / how other approaches fell short / what
kinds of problems certain steps are trying to solve.  Perhaps some links to
papers and such.  Kind of a mini-methods paper here.

Taxonomic read filtration
-------------------------


Human, contaminant, and duplicate read removal
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

BMTAGGER

BLAST

M-Vicuna


Taxonomic selection
~~~~~~~~~~~~~~~~~~~

LASTAL


Viral genome analysis
---------------------

Viral genome assembly
~~~~~~~~~~~~~~~~~~~~~

*de novo* genome assembly with Trinity_.  Reference-assisted
assembly improvements (scaffolding, orienting, etc) with
VFAT_ (which relies on MUSCLE_).

We then do two rounds of assembly improvement (Novoalign_ and GATK_).

.. _Trinity: http://trinityrnaseq.github.io/
.. _VFAT: http://www.broadinstitute.org/scientific-community/science/projects/viral-genomics/v-fat
.. _Novoalign: http://www.novocraft.com/products/novoalign/
.. _GATK: https://www.broadinstitute.org/gatk/


Intrahost variant identification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Intrahost variants (iSNVs) are identified from deep sequence coverage
using `V-Phaser2 <http://dx.doi.org/10.1186/1471-2164-14-674>`_.
For each sample, reads are first aligned to their own consensus
genome with Novoalign, followed by duplicate read removal with Picard and
local realignment with GATK IndelRealigner. V-Phaser2 is called on
each sample to produce a set of iSNV calls.

(then stuff about strand bias filter, then stuff about library
counts)

(then stuff about remapping all calls back to the reference assembly's
coordinate space and alleles using MUSCLE_, and merging calls across all
samples together, emitting in VCF format)

iSNVs are then annotated with snpEff_ and provided in both VCF and tabular
text formats.

.. _snpEff: http://snpeff.sourceforge.net/
.. _MUSCLE: http://www.drive5.com/muscle/




Taxonomic read identification
-----------------------------

Nothing here at the moment. That comes later, but we will later
integrate it when it's ready.

