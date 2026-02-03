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
import json
from collections import OrderedDict

# third-party libraries
import Bio.AlignIO, Bio.SeqIO

# module-specific
from .core import cmd
from .core import file
from .core.misc import CoordMapperError, CoordMapper, CoordMapper2Seqs
from .phylo import muscle
from .phylo import snpeff
from .phylo import mafft
from .phylo import vcf

log = logging.getLogger(__name__)

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
    cmd.common_args(parser, (('tmp_dir', None), ('loglevel', None), ('version', None)))
    cmd.attach_main(parser, snpeff.SnpEff().annotate_vcf, split_args=True)
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

    cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    cmd.attach_main(parser, main_align_mafft)
    return parser


def main_align_mafft(args):
    ''' Run the mafft alignment on the input FASTA file.'''

    mafft.MafftTool().execute(
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

    cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    cmd.attach_main(parser, multichr_mafft)
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
    transposedFiles = file.transposeChromosomeFiles(args.inFastas, args.sampleRelationFile, args.sampleNameListFile)

    # since the FASTA files are
    for idx, filePath in enumerate(transposedFiles):

        # execute MAFFT alignment. The input file is passed within a list, since argparse ordinarily
        # passes input files in this way, and the MAFFT tool expects lists,
        # but in this case we are creating the input file ourselves
        mafft.MafftTool().execute(
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



def full_parser():
    return cmd.make_parser(__commands__, __doc__)


def main():
    cmd.main_argparse(__commands__, __doc__)


if __name__ == '__main__':
    main()
