#!/usr/bin/env python
''' This script contains a number of utilities for SNP calling, multi-alignment,
    phylogenetics, etc.
'''

__author__ = "PLACEHOLDER"
__commands__ = []

# built-ins
import argparse
import logging
import os
import array
import bisect
import json
from itertools import permutations
from collections import OrderedDict, Sequence

from itertools import zip_longest # pylint: disable=E0611

from collections import MutableMapping as DictMixin

# third-party libraries
import Bio.AlignIO
from Bio import SeqIO

# module-specific
import tools.muscle
import tools.snpeff
import tools.mafft
import util.cmd
import util.file
import util.vcf

log = logging.getLogger(__name__)

# =========== CoordMapper =================

# CoordMapper extends DictMixin so that after the basic dict dict() interface methods are defined,
# we get higher level dictionary intervace methods for free


class CoordMapperError(Exception):
    def __init___(self, *args, **kwargs):
        super(CoordMapperError, self).__init__(self, *args, **kwargs)


class CoordMapper(DictMixin):
    """ Map (chrom, coordinate) between genome A and genome B.
        Coordinates are 1-based.
        Indels are handled as follows after corresponding sequences are aligned:
            Return (chrom, None) if base is past either end of other sequence.
            If a base maps to a gap in the other species, return the index
                of the closest upstream non-gap base.
            If a base is followed by a gap then instead of returning an integer,
                return a two-element list representing the interval in the
                other species that aligns to this base and the subsequent gap.
        Assumption: the aligner tool will never align two gaps, and will never
            put gaps in opposite species adjacent to each other without aligning
            a pair of real bases in between.
    """

    def __init__(self, alignerTool=tools.muscle.MuscleTool):
        """ The two genomes are described by fasta files with the same number of
            chromosomes, and corresponding chromosomes must be in same order.
        """

        # {
        #    chrA : {chrB: mapperAB, chrC: mapperAC},
        #    chrB : {chrA: mapperBA, chrC: mapperBC},
        #    chrC : {chrA: mapperCA, chrB: mapperCB}
        # }
        self.chrMaps = OrderedDict()
        self.chrMapsUngapped = OrderedDict()
        self.alignerTool = alignerTool()

    def __getitem__(self, key):
        return self.chrMaps[key]

    def __setitem__(self, key, value):
        raise TypeError("'%s' object does not support item assignment" % self.__class__.__name__)

    def __delitem__(self, key):
        raise TypeError("'%s' object does not support item deletion" % self.__class__.__name__)

    def __len__(self):
        return len(self.chrMaps)

    def __iter__(self):
        for i in self.chrMaps:
            yield i

    def __contains__(self, key):
        if key in self.chrMaps:
            return True
        else:
            return False

    def keys(self):
        return self.chrMaps.keys()

    def mapAtoB(self, fromChrom, fromPos=None, side=0):
        """ Map (chrom, coordinate) from genome A to genome B.
            If fromPos is None, map only the chromosome name
            If side is:
                < 0, return the left-most position on B
                ==0, return either the unique position on B or a [left,right] list
                > 0, return the right-most position on B
        """
        if len(self.chrMaps.keys()) != 4:
            raise LookupError(
                "CoordMapper.mapAtoB expects two input sequences and is provided only as a legacy function")

        return self.mapChr(fromChrom, list(self.chrMaps[fromChrom].keys())[0], fromPos, side)

    def mapBtoA(self, fromChrom, fromPos=None, side=0):
        """ Map (chrom, coordinate) from genome B to genome A.
            If fromPos is None, map only the chromosome name
            If side is:
                < 0, return the left-most position on A
                ==0, return either the unique position on A or a [left,right] list
                > 0, return the right-most position on A
        """
        if len(self.chrMaps.keys()) != 4:
            raise LookupError(
                "CoordMapper.mapBtoA expects two input sequences and is provided only as a legacy function")

        return self.mapChr(fromChrom, list(self.chrMaps[fromChrom].keys())[0], fromPos, side)

    def mapChr(self, fromChrom, toChrom, fromPos=None, side=0):
        """ Map (chrom, coordinate) from seq "fromChrom" to seq "toChrom".
            If fromPos is None, map only the chromosome name
            If side is:
                < 0, return the left-most position on A
                ==0, return either the unique position on A or a [left,right] list
                > 0, return the right-most position on A
        """

        if fromChrom not in self.chrMaps:
            raise KeyError("chr '%s' not found in CoordMapper relation map" % fromChrom)
        if toChrom not in self.chrMaps[fromChrom].keys():
            raise KeyError("chr '%s' not found in CoordMapper relation map" % toChrom)

        mapper = self.chrMaps[fromChrom][toChrom]

        if fromPos is None:
            return toChrom
        toPos = mapper(fromPos, 0)
        if isinstance(toPos, Sequence) and side != 0:
            toPos = toPos[0] if side < 0 else toPos[1]
        return (toChrom, toPos)

    def load_alignments(self, aligned_files, a_idx=None, b_idx=None):
        """ Loads aligned sequences into a CoordMapper instance.
            Any number of sequences >1 may be read in.
            Mappers may be accessed via CoordMapper.chrMaps where chrMaps may look like:
            ```
            {
                chrA : {chrB: mapperAB, chrC: mapperAC},
                chrB : {chrA: mapperBA, chrC: mapperBC},
                chrC : {chrA: mapperCA, chrB: mapperCB}
            }
            ```
        """
        for alignOutFileName in aligned_files:
            with open(alignOutFileName, 'rt') as alignOutFile:
                seqs = list(SeqIO.parse(alignOutFile, 'fasta'))

                # if len(list(seqs)) <2:
                #    raise Exception("Each aligned input file must contain >1 sequence.")

                # if mapping between specific sequences is specified
                if a_idx is not None and b_idx is not None:
                    assert a_idx >= 0 and b_idx >= 0
                    assert a_idx < len(seqs) and b_idx < len(seqs)

                    mapper = CoordMapper2Seqs(seqs[a_idx].seq, seqs[b_idx].seq)
                    self.chrMaps.setdefault(seqs[a_idx].id, OrderedDict())
                    mapDict = OrderedDict()
                    mapDict[seqs[b_idx].id] = mapper
                    self.chrMaps[seqs[a_idx].id] = mapDict

                    mapper = CoordMapper2Seqs(seqs[b_idx].seq, seqs[a_idx].seq)
                    self.chrMaps.setdefault(seqs[b_idx].id, OrderedDict())
                    mapDict = OrderedDict()
                    mapDict[seqs[a_idx].id] = mapper
                    self.chrMaps[seqs[b_idx].id] = mapDict
                # otherwise, map all possible pairwise permutations
                else:
                    for (seq1, seq2) in permutations(seqs, 2):
                        if (seq1.id == seq2.id):
                            raise KeyError("duplicate sequence names '%s', '%s'" % (seq1.id, seq2.id))

                        self.chrMaps.setdefault(seq1.id, OrderedDict())
                        self.chrMapsUngapped.setdefault(seq1.id, OrderedDict())
                        # if the sequence we are mapping onto is already in the map
                        # raise an error
                        # (could occur if same sequence is read in from multiple files)
                        if (seq2.id in self.chrMaps[seq1.id]):
                            raise KeyError(
                                "duplicate sequence name '%s' already in chrMap for %s" % (seq2.id, seq1.id))

                        mapper = CoordMapper2Seqs(seq1.seq, seq2.seq)
                        mapDict = self.chrMaps[seq1.id]
                        mapDict[seq2.id] = mapper
                        self.chrMaps[seq1.id] = mapDict

                        # ungapped strings
                        #longerSeqLen = max( len(seq1.seq.ungap("-")), len(seq2.seq.ungap("-")) )
                        #seq1UngappedPadded = str(seq1.seq.ungap("-")).ljust(longerSeqLen, "N")
                        #seq2UngappedPadded = str(seq2.seq.ungap("-")).ljust(longerSeqLen, "N")
                        #mapper = CoordMapper2Seqs(seq1UngappedPadded, seq2UngappedPadded)
                        #mapDict = self.chrMapsUngapped[seq1.id]
                        #mapDict[seq2.id] = mapper
                        #self.chrMapsUngapped[seq1.id] = mapDict

    def align_and_load_sequences(self, unaligned_fasta_files, aligner=None):
        aligner = self.alignerTool if aligner is None else aligner

        # transpose
        per_chr_fastas = transposeChromosomeFiles(unaligned_fasta_files)
        if not per_chr_fastas:
            raise Exception('no input sequences')
        # align
        alignOutFileNames = []
        for alignInFileName in per_chr_fastas:
            alignOutFileName = util.file.mkstempfname('.fasta')
            aligner.execute(alignInFileName, alignOutFileName)
            alignOutFileNames.append(alignOutFileName)
            os.unlink(alignInFileName)
        # read in
        self.load_alignments(alignOutFileNames)
        # clean up
        for f in alignOutFileNames:
            os.unlink(f)


class CoordMapper2Seqs(object):
    """ Map 1-based coordinates between two aligned sequences.
        Result is a coordinate or an interval, as described in CoordMapper main
            comment string.
        Return None if beyond end.
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
    #
    # Implementation:
    #     mapArrays is a pair of arrays of equal length such that
    #     (mapArrays[0][n], mapArrays[1][n]) are the coordinates of a pair of
    #     aligned real bases on the two sequences. The only pairs that are
    #     included are the first, the last, and the pair immediately following
    #     any gap. Pairs are in increasing order. Coordinate mapping
    #     requires binary search in one of the arrays.
    #     Total space required, in bytes, is const + 8 * (number of indels).
    #     Time for a map in either direction is O(log(number of indels)).
    #

    def __init__(self, seq0, seq1):
        self.mapArrays = [array.array('I'), array.array('I')]
        baseCount0 = 0  # Number of real bases in seq0 up to and including cur pos
        baseCount1 = 0  # Number of real bases in seq1 up to and including cur pos
        beforeStart = True  # Haven't yet reached first pair of aligned real bases
        gapSinceLast = False  # Have encounted a gap since last pair in mapArrays
        for b0, b1 in zip_longest(seq0, seq1):
            if b0 is None or b1 is None:
                raise Exception('CoordMapper2Seqs: sequences must be same length.')
            realBase0 = b0 != '-'
            realBase1 = b1 != '-'
            baseCount0 += realBase0
            baseCount1 += realBase1
            if realBase0 and realBase1:
                if beforeStart or gapSinceLast:
                    self.mapArrays[0].append(baseCount0)
                    self.mapArrays[1].append(baseCount1)
                    gapSinceLast = False
                    beforeStart = False
                finalPos0 = baseCount0  # Last pair of aligned real bases so far
                finalPos1 = baseCount1  # Last pair of aligned real bases so far
            else:
                gapSinceLast = True
        if len(self.mapArrays[0]) > 0:
            if self.mapArrays[0][-1] != finalPos0:
                self.mapArrays[0].append(finalPos0)
                self.mapArrays[1].append(finalPos1)

    def __call__(self, fromPos, fromWhich):
        """ fromPos: 1-based coordinate
            fromWhich: if 0, map from 1st sequence to 2nd, o.w. 2nd to 1st."""
        if len(self.mapArrays[0]) == 0:
            raise Exception('CoordMapper2Seqs: no aligned bases.')
        if fromPos != int(fromPos):
            raise TypeError('CoordMapper2Seqs: pos %s is not an integer' % fromPos)
        fromArray = self.mapArrays[fromWhich]
        toArray = self.mapArrays[1 - fromWhich]
        if fromPos < fromArray[0] or fromPos > fromArray[-1]:
            result = None
        elif fromPos == fromArray[-1]:
            result = toArray[-1]
        else:
            insertInd = bisect.bisect(fromArray, fromPos)
            prevFromPos = fromArray[insertInd - 1]
            nextFromPos = fromArray[insertInd]
            prevToPos = toArray[insertInd - 1]
            nextToPos = toArray[insertInd]
            assert (prevFromPos <= fromPos < nextFromPos)
            prevPlusOffset = prevToPos + (fromPos - prevFromPos)
            if fromPos == nextFromPos - 1 and prevPlusOffset < nextToPos - 1:
                result = [prevPlusOffset, nextToPos - 1]
            else:
                result = min(prevPlusOffset, nextToPos - 1)
        return result

# ========== snpEff annotation of VCF files ==================


def parser_snpEff(parser=argparse.ArgumentParser()):
    parser.add_argument("inVcf", help="Input VCF file")
    parser.add_argument("genomes", nargs='+', help="genome name (snpEff db name, or NCBI accessions)")
    parser.add_argument("outVcf", help="Output VCF file")
    parser.add_argument("--emailAddress",
                        help="""Your email address. To access the Genbank CoreNucleotide database,
        NCBI requires you to specify your email address with each request.
        In case of excessive usage of the E-utilities, NCBI will attempt to contact
        a user at the email address provided before blocking access.""")
    util.cmd.common_args(parser, (('tmp_dir', None), ('loglevel', None), ('version', None)))
    util.cmd.attach_main(parser, tools.snpeff.SnpEff().annotate_vcf, split_args=True)
    return parser


__commands__.append(('snpEff', parser_snpEff))

# =======================
# ***  align_mafft  ***
# =======================


def parser_general_mafft(parser=argparse.ArgumentParser()):
    parser.add_argument('inFastas', nargs='+', help='Input FASTA files.')

    group = parser.add_mutually_exclusive_group()
    group.add_argument('--localpair',
                       default=None,
                       action='store_true',
                       help='All pairwise alignments are computed with the Smith-Waterman algorithm.')
    group.add_argument('--globalpair',
                       default=None,
                       action='store_true',
                       help='All pairwise alignments are computed with the Needleman-Wunsch algorithm.')

    parser.add_argument('--preservecase',
                        default=None,
                        action='store_true',
                        help='Preserve base or aa case, as well as symbols.')
    parser.add_argument('--reorder',
                        default=None,
                        action='store_true',
                        help='Output is ordered aligned rather than in the order of the input (default: %(default)s).')
    parser.add_argument('--gapOpeningPenalty',
                        default=1.53,
                        type=float,
                        help='Gap opening penalty (default: %(default)s).')
    parser.add_argument('--ep', type=float, help='Offset (works like gap extension penalty).')
    parser.add_argument('--verbose', default=False, action='store_true', help='Full output (default: %(default)s).')
    parser.add_argument('--outputAsClustal',
                        default=None,
                        action='store_true',
                        help='Write output file in Clustal format rather than FASTA')
    parser.add_argument(
        '--maxiters',
        default=0,
        type=int,
        help="""Maximum number of refinement iterations (default: %(default)s).
                Note: if "--localpair" or "--globalpair" is specified this defaults to 1000.""")
    return parser


def parser_align_mafft(parser):
    parser = parser_general_mafft(parser)

    parser.add_argument('outFile', help='Output file containing alignment result (default format: FASTA)')

    util.cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, main_align_mafft)
    return parser


def main_align_mafft(args):
    ''' Run the mafft alignment on the input FASTA file.'''

    tools.mafft.MafftTool().execute(
        inFastas=args.inFastas,
        outFile=args.outFile,
        localpair=args.localpair,
        globalpair=args.globalpair,
        preservecase=args.preservecase,
        reorder=args.reorder,
        gapOpeningPenalty=args.gapOpeningPenalty,
        offset=args.ep,
        verbose=args.verbose,
        outputAsClustal=args.outputAsClustal,
        maxiters=args.maxiters,
        threads=args.threads)

    return 0


__commands__.append(('align_mafft', parser_align_mafft))

# =======================
# ***  multichr_mafft  ***
# =======================


def parser_multichr_mafft(parser):
    parser = parser_general_mafft(parser)

    parser.add_argument('outDirectory', help='Location for the output files (default is cwd: %(default)s)')
    parser.add_argument('--outFilePrefix',
                        default="aligned",
                        help='Prefix for the output file name (default: %(default)s)')
    parser.add_argument('--sampleRelationFile',
                        default=None,
                        help="""If the parameter sampleRelationFile is specified
        (as a file path), a JSON file will be written mapping
        sample name to sequence position in the output.""")
    parser.add_argument('--sampleNameListFile',
                        default=None,
                        help="""If the parameter sampleRelationFile is specified
        (as a file path), a file will be written mapping
        sample names in the order of their sequence
        positions in the output.""")

    util.cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, multichr_mafft)
    return parser


def multichr_mafft(args):
    ''' Run the mafft alignment on a series of chromosomes provided in sample-partitioned FASTA files. Output as FASTA.
        (i.e. file1.fasta would contain chr1, chr2, chr3; file2.fasta would also contain chr1, chr2, chr3)'''

    # get the absolute path to the output directory in case it has been specified as a relative path,
    # since MAFFT relies on its CWD for path resolution
    absoluteOutDirectory = os.path.abspath(args.outDirectory)

    # make the output directory if it does not exist
    if not os.path.isdir(absoluteOutDirectory):
        os.makedirs(absoluteOutDirectory)

    # prefix for output files
    prefix = "" if args.outFilePrefix is None else args.outFilePrefix

    # reorder the data into new FASTA files, where each FASTA file has only variants of its respective chromosome
    transposedFiles = transposeChromosomeFiles(args.inFastas, args.sampleRelationFile, args.sampleNameListFile)

    # since the FASTA files are
    for idx, filePath in enumerate(transposedFiles):

        # execute MAFFT alignment. The input file is passed within a list, since argparse ordinarily
        # passes input files in this way, and the MAFFT tool expects lists,
        # but in this case we are creating the input file ourselves
        tools.mafft.MafftTool().execute(
            inFastas=[os.path.abspath(filePath)],
            outFile=os.path.join(absoluteOutDirectory, "{}_{}.fasta".format(prefix, idx + 1)),
            localpair=args.localpair,
            globalpair=args.globalpair,
            preservecase=args.preservecase,
            reorder=args.reorder,
            gapOpeningPenalty=args.gapOpeningPenalty,
            offset=args.ep,
            verbose=args.verbose,
            outputAsClustal=args.outputAsClustal,
            maxiters=args.maxiters,
            threads=args.threads)

    return 0


__commands__.append(('multichr_mafft', parser_multichr_mafft))

# ============================

# modified version of rachel's call_snps_3.py follows


def call_snps_3(inFasta, outVcf, REF="KJ660346.2"):
    a = Bio.AlignIO.read(inFasta, "fasta")
    ref_idx = find_ref(a, REF)
    with open(outVcf, 'wt') as outf:
        outf.write(vcf_header(a))
        for row in make_vcf(a, ref_idx, REF):
            outf.write('\t'.join(map(str, row)) + '\n')


def find_ref(a, ref):
    for i in range(len(a)):
        if a[i].id == ref:
            return i
    return -1


def vcf_header(a):
    header  = "##fileformat=VCFv4.1\n"
    header += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
    header += "##contig=<ID=\"KM034562\",length=18957>\n"
    header += '#' + '\t'.join(['CHROM', 'POS', 'ID', 'REF', 'ALT',
                               'QUAL', 'FILTER', 'INFO', 'FORMAT'] + [x.id for x in a]) + '\n' # pylint: disable=E1101

    return header


def make_vcf(a, ref_idx, chrom):
    bases = set(["A", "C", "G", "T"])
    for i in range(len(a[0])):
        alt = []
        for j in range(len(a)):
            if (a[j][i] != a[ref_idx][i]) and ((a[ref_idx][i] in bases) and (a[j][i] in bases)) and a[j][i] not in alt:
                alt.append(a[j][i])
        if len(alt) > 0:
            row = [chrom, i + 1, '.', a[ref_idx][i], ','.join(alt), '.', '.', '.', 'GT']
            genos = []
            for k in range(len(a)):
                if a[k][i] == a[ref_idx][i]:
                    genos.append(0)
                elif a[k][i] not in bases:
                    genos.append(".")
                else:
                    for m in range(0, len(alt)):
                        if a[k][i] == alt[m]:
                            genos.append(m + 1)
            yield row + genos

class TranspositionError(Exception):
    def __init___(self, *args, **kwargs):
        super(TranspositionError, self).__init__(self, *args, **kwargs)

def transposeChromosomeFiles(inputFilenamesList, sampleRelationFile=None, sampleNameListFile=None):
    ''' Input:  a list of FASTA files representing a genome for each sample.
                Each file contains the same number of sequences (chromosomes, segments,
                etc) in the same order.
                If the parameter sampleRelationFile is specified (as a file path),
                a JSON file will be written mapping sample name to sequence position
                in the output.
        Output: a list of FASTA files representing all samples for each
                chromosome/segment for input to a multiple sequence aligner.
                The number of FASTA files corresponds to the number of chromosomes
                in the genome.  Each file contains the same number of samples
                in the same order.  Each output file is a tempfile.
    '''
    outputFilenames = []

    # open all files
    inputFilesList = [util.file.open_or_gzopen(x, 'r') for x in inputFilenamesList]
    # get BioPython iterators for each of the FASTA files specified in the input
    fastaFiles = [SeqIO.parse(x, 'fasta') for x in inputFilesList]

    # write out json file containing relation of
    # sample name to position in output
    if sampleRelationFile:
        with open(os.path.realpath(sampleRelationFile), "w") as outFile:
            # dict mapping sample->index, zero indexed
            sampleIdxMap = dict((os.path.basename(v).replace(".fasta", ""), k)
                                for k, v in enumerate(inputFilenamesList))
            json.dump(sampleIdxMap, outFile, sort_keys=True, indent=4, separators=(',', ': '))

    if sampleNameListFile:
        with open(os.path.realpath(sampleNameListFile), "w") as outFile:
            sampleNameList = [os.path.basename(v).replace(".fasta", "\n") for v in inputFilenamesList]
            outFile.writelines(sampleNameList)

    # for each interleaved record
    for chrRecordList in zip_longest(*fastaFiles):
        if any(rec is None for rec in chrRecordList):
            raise TranspositionError("input fasta files must all have the same number of sequences")

        outputFilename = util.file.mkstempfname('.fasta')
        outputFilenames.append(outputFilename)
        with open(outputFilename, "w") as outf:
            # write the corresonding records to a new FASTA file
            SeqIO.write(chrRecordList, outf, 'fasta')

    # close all input files
    for x in inputFilesList:
        x.close()

    return outputFilenames


def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)


if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)
