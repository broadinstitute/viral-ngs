#!/usr/bin/env python
''' This script contains a number of utilities for SNP calling, multi-alignment,
    phylogenetics, etc.
'''

__author__ = "PLACEHOLDER"
__commands__ = []

import Bio.AlignIO
from Bio import SeqIO
import argparse, logging, os
try :
    from itertools import zip_longest
except ImportError :
    from itertools import izip_longest as zip_longest
import tools.muscle
import util.cmd, util.file, util.vcf, util.misc
from collections import OrderedDict

log = logging.getLogger(__name__)


class CoordMapper :
    """ Maps coordinates between genome A and genome B.
        Coordinates are 1-based.
        Indels are handled as follows after corresponding sequences are aligned:
            Raises IndexError if base is past either end of the other sequence.
            If a base maps to a gap in the other species, returns the index
                of the closest upstream non-gap base.
            If a base is followed by a gap then instead of returning an integer,
                returns a two-element list representing the interval in the
                other species that aligns to this base and the subsequent gap.
    """
    def __init__(self, fastaA, fastaB, alignerTool = tools.muscle.MuscleTool) :
        """ The two genomes are described by fasta files with the same number of 
            chromosomes, and corresponding chromosomes must be in same order.
        """
        self.AtoB = OrderedDict() # {chrA : [chrB, mapper],...}
        self.BtoA = OrderedDict() # {chrB : [chrA, mapper],...}
        self._align(fastaA, fastaB, alignerTool())
    
    def mapAtoB(self, fromChrom, pos=None) :
        toChrom, mapper = self.AtoB[fromChrom]
        if pos != None:
            pos = mapper(pos)
        return (toChrom, pos)
    
    def mapBtoA(self, fromChrom, pos=None) :
        toChrom, mapper = self.BtoA[fromChrom]
        if pos != None:
            pos = mapper(pos)
        return (toChrom, pos)

    def _align(self, fastaA, fastaB, aligner) :
        alignInFileName = util.file.mkstempfname('.fasta')
        alignOutFileName = util.file.mkstempfname('.fasta')
        with util.file.open_or_gzopen(fastaA) as infA :
            with util.file.open_or_gzopen(fastaB) as infB :
                numSeqs = 0
                for recA, recB in zip_longest(SeqIO.parse(infA, 'fasta'),
                                               SeqIO.parse(infB, 'fasta')) :
                    assert (recA != None and recB != None), 'CoordMapper '\
                            'input files must have same number of sequences.'
                    numSeqs += 1
                    chrA = recA.id
                    chrB = recB.id
                    with open(alignInFileName, 'wt') as alignInFile :
                        SeqIO.write(recA, alignInFile, 'fasta')
                        SeqIO.write(recB, alignInFile, 'fasta')
                    aligner.execute(alignInFileName, alignOutFileName)
                    with open(alignOutFileName) as alignOutFile :
                        seqParser = SeqIO.parse(alignOutFile, 'fasta')
                        seqs = [seqRec.seq for seqRec in seqParser]
                        self.AtoB[chrA] = [chrB, self._Mapper(seqs[0], seqs[1])]
                        self.BtoA[chrB] = [chrA, self._Mapper(seqs[1], seqs[0])]
        assert numSeqs > 0, 'CoordMapper: no input sequences.'
        assert len(self.AtoB) == len(self.BtoA) == numSeqs, \
               'CoordMapper: duplicate sequence name.'
        os.unlink(alignInFileName)
        os.unlink(alignOutFileName)

    class _Mapper(object) :
        """Class that maps a 1-based position on one sequence to a 1-based 
              position or interval on another, raising IndexError if beyond end.
           Gap handling is described in CoordMapper main comment string.
        """
        # Current implementation is a space hog. If space utilitization needs to
        # be improved all that needs to be done is to replace the list that
        # currently holds the map information with a new class that implements
        # __len__, __getitem__, __setitem__, and append, and that efficiently
        # stores long runs of consecutive integers. The __getitem__ and
        # __setitem__ methods need to handle negative indices, but needn't
        # handle slices. Once the class is defined, change "Map = []" to
        # "Map = NewClass()".
        def __init__(self, seqFrom, seqTo) :
            """Input sequences are BioPython Seqs representing aligner output:
               sequences of equal length with gaps represented by dashes."""
            self.map = self._make_map(seqFrom, seqTo)
            #       map is a list whose ith element is result of mapping
            #       position i+1 of seqFrom to seqTo,
            #       or -1 if beyond either end of seqTo.

        def __call__(self, fromPos) :
            if fromPos < 1 or fromPos > len(self.map) :
                raise IndexError
            result = self.map[fromPos - 1]
            if result == -1 :
                raise IndexError
            return result

        @staticmethod
        def _make_map(seqFrom, seqTo) :
            Map = [] # Use uppercase to avoid conflict with built in map
            baseCountFrom = 0  # Number of real From bases upstream of cur pos
            baseCountTo = 0    # Number of real To bases upstream of cur pos
            numFromPastEnd = 0 # Number of bases in From past end of To
            numToPastEnd = 0   # Number of bases in To past end of From
            for bFrom, bTo in zip(seqFrom, seqTo) :
                realFromBase = bFrom != '-'
                realToBase   =   bTo != '-'
                assert(realFromBase or realToBase) # This simplifies the rest
                if realFromBase :
                    if realToBase :
                        result = baseCountTo + 1
                    elif baseCountTo == 0 :
                        result = -1 # Before start of To sequence
                    else :
                        result = Map[-1] # Aligned to gap; map to previous base
                    Map.append(result) # Sometimes changed to interval later
                    baseCountFrom += 1
                    numToPastEnd = 0
                else : # From base is a gap; extend map of previous real base.
                    if baseCountFrom != 0 :
                        if not type(Map[-1]) == list :
                            Map[-1] = [Map[-1], Map[-1]]
                        Map[-1][1] += 1
                    numToPastEnd += 1
                if realToBase :
                    baseCountTo += 1
                    numFromPastEnd = 0
                else : # May assume realFromBase
                    numFromPastEnd += 1
            for ii in range(numFromPastEnd) :
                Map[-ii - 1] = -1 # Mark positions past end of other sequence
            if type(Map[-1]) == list :
                Map[-1][1] -= numToPastEnd # Don't map to these; past end of seq
                if Map[-1][0] == Map[-1][1] :
                    Map[-1] = Map[-1][0]
            return Map


# modified version of rachel's call_snps_3.py follows
def call_snps_3(inFasta, outVcf, REF="KJ660346.2"):
    a=Bio.AlignIO.read(inFasta, "fasta")
    ref_idx = find_ref(a, REF)
    with open(outVcf, 'wt') as outf:
        outf.write(vcf_header(a))
        for row in make_vcf(a, ref_idx, REF):
            outf.write('\t'.join(map(str, row))+'\n')
def find_ref(a, ref):
    for i in range(len(a)):
        if a[i].id == ref:
            return i
    return -1
def vcf_header(a):
    header = "##fileformat=VCFv4.1\n"
    header += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
    header += "##contig=<ID=\"KM034562\",length=18957>\n"
    header += '#' + '\t'.join(['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'] + [x.id for x in a]) + '\n'
    return header
def make_vcf(a, ref_idx, chrom):
    bases=set(["A", "C", "G", "T"])
    for i in range(len(a[0])):
        alt = []
        for j in range(len(a)):
            if (a[j][i] != a[ref_idx][i]) and ((a[ref_idx][i] in bases) and (a[j][i] in bases)) and a[j][i] not in alt:
                alt.append(a[j][i])
        if len(alt) > 0:
            row = [chrom, i+1, '.', a[ref_idx][i], ','.join(alt), '.', '.', '.', 'GT']
            vars = []
            for k in range(len(a)):
                if a[k][i] == a[ref_idx][i]:
                    vars.append(0)
                elif a[k][i] not in bases:
                    vars.append(".")
                else:
                    for m in range(0, len(alt)):
                        if a[k][i] == alt[m]:
                            vars.append(m+1)
            yield row+vars
            

def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)
if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)
