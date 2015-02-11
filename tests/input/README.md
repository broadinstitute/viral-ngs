Description of input files:

- ebola.fasta is a collection of Ebolavirus genomes from many Ebolavirus species
(EBOV, SUDV, etc)

- ebov-makona.fasta is our typical EBOV 2014 reference genome (the earliest Baize et al
sequence from Guinea)

- G5012.3.testreads.bam is a stripped down set of reads for one particular Sierra Leone
patient. These reads are unaligned. This only contains 9179 read pairs that map to
EBOV, none of which are PCR duplicates according to Picard after alignment. An additional
200 read pairs that fail to map to EBOV are included here as well.  This contains reads
from 12 read groups: six from one library prep and six from another independent prep.
There are no non-EBOV reads here.

- G5012.3.fasta is the consensus assembly that is created from the above reads.
