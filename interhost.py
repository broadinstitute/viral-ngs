#!/usr/bin/env python
''' This script contains a number of utilities for SNP calling, multi-alignment,
    phylogenetics, etc.
'''

__author__ = "PLACEHOLDER"
__commands__ = []

import Bio.AlignIO
from Bio import SeqIO
import argparse, logging, os, array, bisect
try :
    from itertools import zip_longest
except ImportError :
    from itertools import izip_longest as zip_longest
import tools.muscle
import util.cmd, util.file, util.vcf, util.misc
from collections import OrderedDict

log = logging.getLogger(__name__)

# =========== CoordMapper =================

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
        Assumption: the aligner tool will never align two gaps, and will never
            put gaps in opposite species adjacent to each other without aligning
            a pair of real bases in between.
    """
    def __init__(self, fastaA, fastaB, alignerTool = tools.muscle.MuscleTool) :
        """ The two genomes are described by fasta files with the same number of 
            chromosomes, and corresponding chromosomes must be in same order.
        """
        self.AtoB = OrderedDict() # {chrA : [chrB, mapperAB], chrC : [chrD, mapperCD], ...}
        self.BtoA = OrderedDict() # {chrB : [chrA, mapperAB], chrD : [chrC, mapperCD], ...}
        
        self._align(fastaA, fastaB, alignerTool())
    
    def mapAtoB(self, fromChrom, pos=None, side=0) :
        ''' Map coordinates from genome A to genome B.
            If pos is None, map only the chromosome name
            If side is:
                < 0, return the left-most position on B
                ==0, return either the unique position on B or a (left,right) tuple
                > 0, return the right-most position on B
        '''
        toChrom, mapper = self.AtoB[fromChrom]
        if pos == None:
            return toChrom
        pos = mapper(pos, 0)
        if type(pos) != int and side != 0:
            pos = pos[0] if side<0 else pos[1]
        return (toChrom, pos)
    
    def mapBtoA(self, fromChrom, pos=None, side=0) :
        ''' Map coordinates from genome B to genome A.
            If pos is None, map only the chromosome name
            If side is:
                < 0, return the left-most position on A
                ==0, return either the unique position on A or a (left,right) tuple
                > 0, return the right-most position on A
        '''
        toChrom, mapper = self.BtoA[fromChrom]
        if pos == None:
            return toChrom
        pos = mapper(pos, 1)
        if type(pos) != int and side != 0:
            pos = pos[0] if side<0 else pos[1]
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
                        mapper = CoordMapper2Seqs(seqs[0], seqs[1])
                        self.AtoB[chrA] = [chrB, mapper]
                        self.BtoA[chrB] = [chrA, mapper]
        assert numSeqs > 0, 'CoordMapper: no input sequences.'
        assert len(self.AtoB) == len(self.BtoA) == numSeqs, \
               'CoordMapper: duplicate sequence name.'
        os.unlink(alignInFileName)
        os.unlink(alignOutFileName)

class CoordMapper2Seqs(object) :
    """ Map 1-based coordinates between two aligned sequences.
        Result is a coordinate or an interval, as described in CoordMapper main 
            comment string.
        Raise IndexError if beyond end.
        Input sequences must be already-aligned iterators through bases with
            gaps represented by dashes and all other characters assumed to be
            real bases. 
        Assumptions:
            - Sequences (including gaps) are same length.
            - Each sequence has at least one real base.
            - A gap is never aligned to a gap.
            - A gap in one sequence is never adjacent to a gap in the other;
                there must always be an intervening real base between two gaps.
    """
    """
    Implementation:
        mapArrays is a pair of arrays of equal length such that
        (mapArrays[0][n], mapArrays[1][n]) are the coordinates of a pair of
        aligned real bases on the two sequences. The only pairs that are 
        included are the first, the last, and the pair immediately following 
        any gap. Pairs are in increasing order. Coordinate mapping
        requires binary search in one of the arrays.
        Total space required, in bytes, is const + 8 * (number of indels).
        Time for a map in either direction is O(log(number of indels)).
    """
    
    def __init__(self, seq0, seq1) :
        self.mapArrays = [array.array('I'), array.array('I')]
        baseCount0 = 0  # Number of real bases in seq0 up to and including cur pos
        baseCount1 = 0  # Number of real bases in seq1 up to and including cur pos
        beforeStart = True # Haven't yet reached first pair of aligned real bases
        gapSinceLast = False # Have encounted a gap since last pair in mapArrays
        prevRealBase0 = prevRealBase1 = True
        for b0, b1 in zip_longest(seq0, seq1) :
            assert b0 != None and b1 != None, 'CoordMapper2Seqs: sequences '\
                'must be same length.'
            realBase0 = b0 != '-'
            realBase1 = b1 != '-'
            assert realBase0 or realBase1, 'CoordMapper2Seqs: gap aligned to gap.'
            assert (realBase0 or prevRealBase1) and (realBase1 or prevRealBase0),\
                 'CoordMapper2Seqs: gap in one sequence adjacent to gap in other.'
            prevRealBase0 = realBase0
            prevRealBase1 = realBase1
            baseCount0 += realBase0
            baseCount1 += realBase1
            if realBase0 and realBase1 :
                if beforeStart or gapSinceLast :
                    self.mapArrays[0].append(baseCount0)
                    self.mapArrays[1].append(baseCount1)
                    gapSinceLast = False
                    beforeStart = False
                finalPos0 = baseCount0 # Last pair of aligned real bases so far
                finalPos1 = baseCount1 # Last pair of aligned real bases so far
            else :
                gapSinceLast = True
        assert len(self.mapArrays[0]) != 0, 'CoordMapper2Seqs: no aligned bases.'
        if self.mapArrays[0][-1] != finalPos0 :
            self.mapArrays[0].append(finalPos0)
            self.mapArrays[1].append(finalPos1)

    def __call__(self, fromPos, fromWhich) :
        """ fromPos: 1-based coordinate
            fromWhich: if 0, map from 1st sequence to 2nd, o.w. 2nd to 1st."""
        if fromPos != int(fromPos) :
            raise TypeError('CoordMapper2Seqs: pos %s is not an integer' % fromPos)
        fromArray = self.mapArrays[fromWhich]
        toArray = self.mapArrays[1 - fromWhich]
        if fromPos < fromArray[0] or fromPos > fromArray[-1] :
            raise IndexError
        if fromPos == fromArray[-1] :
            result = toArray[-1]
        else :
            insertInd = bisect.bisect(fromArray, fromPos)
            prevFromPos = fromArray[insertInd - 1]
            nextFromPos = fromArray[insertInd]
            prevToPos = toArray[insertInd - 1]
            nextToPos = toArray[insertInd]
            assert(prevFromPos <= fromPos < nextFromPos)
            prevPlusOffset = prevToPos + (fromPos - prevFromPos)
            if fromPos == nextFromPos - 1 and prevPlusOffset < nextToPos - 1 :
                result = [prevPlusOffset, nextToPos - 1]
            else :
                result = min(prevPlusOffset, nextToPos - 1)
        return result

# ============================


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
