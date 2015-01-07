[![Build Status](https://travis-ci.org/broadinstitute/viral-ngs.svg)](https://travis-ci.org/broadinstitute/viral-ngs)
[![Coverage Status](https://coveralls.io/repos/broadinstitute/viral-ngs/badge.png)](https://coveralls.io/r/broadinstitute/viral-ngs)

viral-ngs
=========

A set of scripts and tools for the analysis of viral NGS data.  More 
details to come.


Scripts:
 1. read_utils.py - various utilities for operating on reads
 2. taxon_filter.py - filter reads for or against certain species, genera,
   or other taxonomic classifications
 3. assembly.py - assemble reads, clean contigs, re-map and refine
   consensus assemblies
 4. interhost.py - align consensus sequences, call SNPs, annotate, produce
   trees
 5. intrahost.py - align reads back to consensus and call intrahost
   variants
 6. reports.py - various statistics
 7. broad_utils.py - tools specifically for handling Broad Institute
   walk-up sequencing outputs
