#!/usr/bin/env python
'''This script contains a number of utilities for submitting our analyses
to NCBI's Genbank and SRA databases, as well as retreiving records from Genbank.
'''

__author__ = "tomkinsc@broadinstitute.org"
__commands__ = []

import argparse
import logging
import collections
import shutil
import os
import os.path

import Bio.SeqIO

import util.cmd
import util.file
import util.version
import util.genbank
import tools.tbl2asn
import interhost
import util.feature_table

log = logging.getLogger(__name__)


def fasta_chrlens(fasta):
    out = collections.OrderedDict()
    with open(fasta, 'rt') as inf:
        for seq in Bio.SeqIO.parse(inf, 'fasta'):
            out[seq.id] = len(seq)
    return out

def tbl_transfer_common(cmap, ref_tbl, out_tbl, alt_chrlens, oob_clip=False, ignore_ambig_feature_edge=False):
    """
        This function is the feature transfer machinery used by tbl_transfer()
        and tbl_transfer_prealigned(). cmap is an instance of CoordMapper.
    """

    ft = util.feature_table.FeatureTable(ref_tbl)
    remapped_ft = util.feature_table.FeatureTable()

    # sequence identifier
    refSeqID = [x for x in cmap.keys() if ft.refID in x][0]
    altid = list(set(cmap.keys()) - set([refSeqID]))[0]
    ft.refID = altid

    # feature with numeric coordinates (map them)
    def remap_function(start, end, feature):
        """
            start/end are SeqPositions from util.feature_table
        """
        strand = None
        if end.position >= start.position:
            strand = '+'
            start.position = cmap.mapChr(refSeqID, altid, start.position, -1)[1]
            end.position = cmap.mapChr(refSeqID, altid, end.position, 1)[1]
        else:
            # negative strand feature
            strand = '-'
            start.position = cmap.mapChr(refSeqID, altid, start.position, 1)[1]
            end.position = cmap.mapChr(refSeqID, altid, end.position, -1)[1]

        if ignore_ambig_feature_edge:
            start.location_operator = None
            end.location_operator = None
        
        if start.position and end.position:
            # feature completely within bounds
            return (start, end)
        elif start.position == None and end.position == None:
            # feature completely out of bounds
            return (None, None)
        else:
            # feature overhangs end of sequence
            if oob_clip:
                if feature.type == "CDS":
                    feature.add_note("sequencing did not capture complete CDS")
                if strand == '+':
                    # clip pos strand feature
                    if start.position == None:
                        # clip the beginning
                        if feature.type == 'CDS':
                            # for CDS features, clip in multiples of 3
                            r = (end.position if end.position is not None else alt_chrlens[altid])
                            start.position = '{}'.format((r % 3) + 1)
                            start.location_operator = '<'
                        else:
                            start.position = '1'
                            start.location_operator = '<'
                    if end.position == None:
                        # clip the end
                        end.position = '{}'.format(alt_chrlens[altid])
                        end.location_operator = '>'
                else:
                    # clip neg strand feature
                    if start.position == None:
                        # clip the beginning (right side)
                        r = alt_chrlens[altid]
                        if feature.type == 'CDS':
                            # for CDS features, clip in multiples of 3
                            l = (end.position if end.position is not None else 1) # new left
                            r = r - ((r-l+1) % 3) # create new right in multiples of 3 from left
                            if (r-l) < 3:
                                # less than a codon remains, drop it
                                return (None, None)
                        start.position = '{}'.format(r)
                        start.location_operator = '<'
                    if end.position == None:
                        # clip the end (left side)
                        end.position = '1'
                        end.location_operator = '<'
            else:
                return (None, None)
        return (start, end)

    ft.remap_locations(remap_function)

    with open(out_tbl, 'wt') as outf:
        exclude_patterns = [
            # regexp, matched anywhere in line
            r"protein_id"
        ]
        for line in ft.lines(exclude_patterns=exclude_patterns):
            outf.write(str(line) + '\n')
        # extra newline at the end
        outf.write('\n')


def tbl_transfer(ref_fasta, ref_tbl, alt_fasta, out_tbl, oob_clip=False, ignore_ambig_feature_edge=False):
    ''' This function takes an NCBI TBL file describing features on a genome
        (genes, etc) and transfers them to a new genome.
    '''
    cmap = interhost.CoordMapper()
    cmap.align_and_load_sequences([ref_fasta, alt_fasta])
    alt_chrlens = fasta_chrlens(alt_fasta)

    tbl_transfer_common(cmap, ref_tbl, out_tbl, alt_chrlens, oob_clip, ignore_ambig_feature_edge)


def parser_tbl_transfer(parser=argparse.ArgumentParser()):
    parser.add_argument("ref_fasta", help="Input sequence of reference genome")
    parser.add_argument("ref_tbl", help="Input reference annotations (NCBI TBL format)")
    parser.add_argument("alt_fasta", help="Input sequence of new genome")
    parser.add_argument("out_tbl", help="Output file with transferred annotations")
    parser.add_argument('--oob_clip',
                        default=False,
                        action='store_true',
                        help='''Out of bounds feature behavior.
        False: drop all features that are completely or partly out of bounds
        True:  drop all features completely out of bounds
               but truncate any features that are partly out of bounds''')
    parser.add_argument('--ignoreAmbigFeatureEdge',
                        dest="ignore_ambig_feature_edge",
                        default=False,
                        action='store_true',
                        help='''Ambiguous feature behavior.
        False: features specified as ambiguous ("<####" or ">####") are mapped, 
               where possible
        True:  features specified as ambiguous ("<####" or ">####") are interpreted
               as exact values''')
    util.cmd.common_args(parser, (('tmp_dir', None), ('loglevel', None), ('version', None)))
    util.cmd.attach_main(parser, tbl_transfer, split_args=True)
    return parser


__commands__.append(('tbl_transfer', parser_tbl_transfer))


def tbl_transfer_prealigned(inputFasta, refFasta, refAnnotTblFiles, outputDir, oob_clip=False, ignore_ambig_feature_edge=False):
    """
        This breaks out the ref and alt sequences into separate fasta files, and then
        creates unified files containing the reference sequence first and the alt second. Each of these unified files
        is then passed as a cmap to tbl_transfer_common.

        This function expects to receive one fasta file containing a multialignment of a single segment/chromosome along
        with the respective reference sequence for that segment/chromosome. It also expects a reference containing all
        reference segments/chromosomes, so that the reference sequence can be identified in the input file by name. It
        also expects a list of reference tbl files, where each file is named according to the ID present for its
        corresponding sequence in the refFasta. For each non-reference sequence present in the inputFasta, two files are
        written: a fasta containing the segment/chromosome for the same, along with its corresponding feature table as
        created by tbl_transfer_common.
    """

    ref_tbl = ""  # must be identified in list of tables
    ref_fasta_filename = ""
    matchingRefSeq = None

    if not os.path.exists(outputDir):
        os.makedirs(outputDir)

    # identify which of the sequences in the multialignment is the reference,
    # matching by ID to one of the sequences in the refFasta
    with util.file.compressed_open(inputFasta, 'rt') as inf:
        for seq in Bio.SeqIO.parse(inf, 'fasta'):
            with util.file.compressed_open(refFasta, 'rt') as reff:
                for refSeq in Bio.SeqIO.parse(reff, 'fasta'):
                    if seq.id == refSeq.id:
                        ref_fasta_filename = util.file.mkstempfname('.fasta')
                        matchingRefSeq = seq
                        break
            if matchingRefSeq:
                break

    if ref_fasta_filename == "":
        raise KeyError("No reference was found in the input file %s" % (inputFasta))

    # identify the correct feature table source based on its filename,
    # which should correspond to a unique component of the ref sequence ID (i.e. the genbank accession)
    for tblFilename in refAnnotTblFiles:
        # identify the correct feature table as the one that has an ID that is
        # part of the ref seq ID
        fileAccession = util.genbank.get_feature_table_id(tblFilename)
        if fileAccession == matchingRefSeq.id.split('|')[0]:
            ref_tbl = tblFilename
            break
    if ref_tbl == "":
        raise KeyError("No reference table was found for the reference %s" % (matchingRefSeq.id))

    # write out the desired sequences to separate fasta files
    with util.file.compressed_open(inputFasta, 'rt') as inf:
        for seq in Bio.SeqIO.parse(inf, 'fasta'):
            # if we are looking at the reference sequence in the multialignment,
            # continue to the next sequence
            if seq.id == matchingRefSeq.id:
                continue

            combined_fasta_filename = ""

            combined_fasta_filename = util.file.mkstempfname('.fasta')
            # write ref and alt sequences to a combined fasta file, sourced from the
            # alignment so gaps are present for the CoordMapper instance, cmap
            with open(combined_fasta_filename, 'wt') as outf:
                Bio.SeqIO.write([matchingRefSeq, seq], outf, "fasta")

            # create a filepath for the output table
            out_tbl = os.path.join(outputDir, seq.id + ".tbl")

            cmap = interhost.CoordMapper()
            cmap.load_alignments([combined_fasta_filename])
            # sequences in the fasta file here should NOT include gaps
            # since alt_chrlens is only used in the case where features would 
            # extend beyond the genome (for reporting >{seq.len})
            alt_chrlens = {}#fasta_chrlens(combined_fasta_filename)
            alt_chrlens[seq.id] = len(seq.seq.ungap("-"))
            alt_chrlens[matchingRefSeq.id] = len(matchingRefSeq.seq.ungap("-"))

            tbl_transfer_common(cmap, ref_tbl, out_tbl, alt_chrlens, oob_clip, ignore_ambig_feature_edge)

def parser_tbl_transfer_prealigned(parser=argparse.ArgumentParser()):
    parser.add_argument("inputFasta",
                        help="""FASTA file containing input sequences,
        including pre-made alignments and reference sequence""")
    parser.add_argument("refFasta", help="FASTA file containing the reference genome")
    parser.add_argument("refAnnotTblFiles",
                        nargs='+',
                        help="""Name of the reference feature tables,
        each of which should have a filename comrised of [refId].tbl
        so they can be matched against the reference sequences""")
    parser.add_argument("outputDir", help="The output directory")
    parser.add_argument('--oob_clip',
                        default=False,
                        action='store_true',
                        help='''Out of bounds feature behavior.
        False: drop all features that are completely or partly out of bounds
        True:  drop all features completely out of bounds
               but truncate any features that are partly out of bounds''')
    parser.add_argument('--ignoreAmbigFeatureEdge',
                        dest="ignore_ambig_feature_edge",
                        default=False,
                        action='store_true',
                        help='''Ambiguous feature behavior.
        False: features specified as ambiguous ("<####" or ">####") are mapped, 
               where possible
        True:  features specified as ambiguous ("<####" or ">####") are interpreted
               as exact values''')
    util.cmd.common_args(parser, (('tmp_dir', None), ('loglevel', None), ('version', None)))
    util.cmd.attach_main(parser, tbl_transfer_prealigned, split_args=True)
    return parser


__commands__.append(('tbl_transfer_prealigned', parser_tbl_transfer_prealigned))


def fetch_fastas(accession_IDs, destinationDir, emailAddress, forceOverwrite, combinedFilePrefix, fileExt,
                 removeSeparateFiles, chunkSize, api_key=None):
    '''
        This function downloads and saves the FASTA files
        from the Genbank CoreNucleotide database given a given list of accession IDs.
    '''
    util.genbank.fetch_fastas_from_genbank(accession_IDs,
                                           destinationDir,
                                           emailAddress,
                                           forceOverwrite,
                                           combinedFilePrefix,
                                           removeSeparateFiles,
                                           fileExt,
                                           "fasta",
                                           chunkSize=chunkSize,
                                           api_key=api_key)


def fetch_feature_tables(accession_IDs, destinationDir, emailAddress, forceOverwrite, combinedFilePrefix, fileExt,
                         removeSeparateFiles, chunkSize, api_key=None):
    '''
        This function downloads and saves
        feature tables from the Genbank CoreNucleotide database given a given list of accession IDs.
    '''
    util.genbank.fetch_feature_tables_from_genbank(accession_IDs,
                                                   destinationDir,
                                                   emailAddress,
                                                   forceOverwrite,
                                                   combinedFilePrefix,
                                                   removeSeparateFiles,
                                                   fileExt,
                                                   "ft",
                                                   chunkSize=chunkSize,
                                                   api_key=api_key)


def fetch_genbank_records(accession_IDs, destinationDir, emailAddress, forceOverwrite, combinedFilePrefix, fileExt,
                          removeSeparateFiles, chunkSize, api_key=None):
    '''
        This function downloads and saves
        full flat text records from Genbank CoreNucleotide database given a given list of accession IDs.
    '''
    util.genbank.fetch_full_records_from_genbank(accession_IDs,
                                                 destinationDir,
                                                 emailAddress,
                                                 forceOverwrite,
                                                 combinedFilePrefix,
                                                 removeSeparateFiles,
                                                 fileExt,
                                                 "gb",
                                                 chunkSize=chunkSize,
                                                 api_key=api_key)


def parser_fetch_reference_common(parser=argparse.ArgumentParser()):
    parser.add_argument("emailAddress",
                        help="""Your email address. To access Genbank databases,
        NCBI requires you to specify your email address with each request.
        In case of excessive usage of the E-utilities, NCBI will attempt to contact
        a user at the email address provided before blocking access. This email address should
        be registered with NCBI. To register an email address, simply send
        an email to eutilities@ncbi.nlm.nih.gov including your email address and
        the tool name (tool='https://github.com/broadinstitute/viral-ngs').""")
    parser.add_argument("--api_key",
                        dest="api_key",
                        help="""Your NCBI API key. If an API key is not provided, NCBI 
                        requests are limited to 3/second. If an API key is provided, 
                        requests may be submitted at a rate up to 10/second. 
                        For more information, see: https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/
                        """)
    parser.add_argument("destinationDir", help="Output directory with where .fasta and .tbl files will be saved")
    parser.add_argument("accession_IDs", nargs='+', help="List of Genbank nuccore accession IDs")
    parser.add_argument('--forceOverwrite',
                        default=False,
                        action='store_true',
                        help='''Overwrite existing files, if present.''')
    parser.add_argument('--combinedFilePrefix',
                        help='''The prefix of the file containing the combined concatenated
                 results returned by the list of accession IDs, in the order provided.''')
    parser.add_argument('--fileExt', default=None, help='''The extension to use for the downloaded files''')
    parser.add_argument('--removeSeparateFiles',
                        default=False,
                        action='store_true',
                        help='''If specified, remove the individual files and leave only the combined file.''')
    parser.add_argument('--chunkSize',
                        default=1,
                        type=int,
                        help='''Causes files to be downloaded from GenBank in chunks of N accessions.
        Each chunk will be its own combined file, separate from any combined
        file created via --combinedFilePrefix (default: %(default)s). If chunkSize is
        unspecified and >500 accessions are provided, chunkSize will be set to 500 to
        adhere to the NCBI guidelines on information retreival.''')
    return parser


def parser_fetch_fastas(parser):
    parser = parser_fetch_reference_common(parser)

    util.cmd.common_args(parser, (('tmp_dir', None), ('loglevel', None), ('version', None)))
    util.cmd.attach_main(parser, fetch_fastas, split_args=True)
    return parser


__commands__.append(('fetch_fastas', parser_fetch_fastas))


def parser_fetch_feature_tables(parser):
    parser = parser_fetch_reference_common(parser)

    util.cmd.common_args(parser, (('tmp_dir', None), ('loglevel', None), ('version', None)))
    util.cmd.attach_main(parser, fetch_feature_tables, split_args=True)
    return parser


__commands__.append(('fetch_feature_tables', parser_fetch_feature_tables))


def parser_fetch_genbank_records(parser):
    parser = parser_fetch_reference_common(parser)

    util.cmd.common_args(parser, (('tmp_dir', None), ('loglevel', None), ('version', None)))
    util.cmd.attach_main(parser, fetch_genbank_records, split_args=True)
    return parser


__commands__.append(('fetch_genbank_records', parser_fetch_genbank_records))


def fasta2fsa(infname, outdir, biosample=None):
    ''' copy a fasta file to a new directory and change its filename to end in .fsa
        for NCBI's sake. '''
    outfname = os.path.basename(infname)
    if outfname.endswith('.fasta'):
        outfname = outfname[:-6]
    elif outfname.endswith('.fa'):
        outfname = outfname[:-3]
    if not outfname.endswith('.fsa'):
        outfname = outfname + '.fsa'
    outfname = os.path.join(outdir, outfname)
    with open(infname, 'rU') as inf:
        with open(outfname, 'wt') as outf:
            for line in inf:
                if line.startswith('>') and biosample:
                    line = line.rstrip('\n\r')
                    line = '{} [biosample={}]\n'.format(line, biosample)
                outf.write(line)
    return outfname


def make_structured_comment_file(cmt_fname, name=None, seq_tech=None, coverage=None):
    with open(cmt_fname, 'wt') as outf:
        outf.write("StructuredCommentPrefix\t##Genome-Assembly-Data-START##\n")
        # note: the <tool name> v. <version name> format is required by NCBI, don't remove the " v. "
        outf.write("Assembly Method\tgithub.com/broadinstitute/viral-ngs v. {}\n".format(util.version.get_version()))
        if name:
            outf.write("Assembly Name\t{}\n".format(name))
        if coverage:
            coverage = str(coverage)
            if not coverage.endswith('x') or coverage.endswith('X'):
                coverage = coverage + 'x'
            outf.write("Genome Coverage\t{}\n".format(coverage))
        if seq_tech:
            outf.write("Sequencing Technology\t{}\n".format(seq_tech))
        outf.write("StructuredCommentSuffix\t##Genome-Assembly-Data-END##\n")


def prep_genbank_files(templateFile, fasta_files, annotDir,
                       master_source_table=None, comment=None, sequencing_tech=None,
                       coverage_table=None, biosample_map=None, organism=None, mol_type=None):
    ''' Prepare genbank submission files.  Requires .fasta and .tbl files as input,
        as well as numerous other metadata files for the submission.  Creates a
        directory full of files (.sqn in particular) that can be sent to GenBank.
    '''
    # get coverage map
    coverage = {}
    if coverage_table:
        for row in util.file.read_tabfile_dict(coverage_table):
            if row.get('sample') and row.get('aln2self_cov_median'):
                coverage[row['sample']] = row['aln2self_cov_median']

    # get biosample id map
    biosample = {}
    if biosample_map:
        for row in util.file.read_tabfile_dict(biosample_map):
            if row.get('sample') and row.get('BioSample'):
                biosample[row['sample']] = row['BioSample']

    # make output directory
    util.file.mkdir_p(annotDir)
    for fn in fasta_files:
        if not fn.endswith('.fasta'):
            raise Exception("fasta files must end in .fasta")
        sample_base = os.path.basename(fn)[:-6]

        # for each segment/chromosome in the fasta file,
        # create a separate new *.fsa file
        with open(fn, "r") as inf:
            asm_fasta = Bio.SeqIO.parse(inf, 'fasta')
            for idx, seq_obj in enumerate(asm_fasta):
                sample = sample_base + "-" + str(idx+1)

                # write the segment to a temp .fasta file
                # in the same dir so fasta2fsa functions as expected
                out_file_name = os.path.join(os.path.dirname(fn),sample+".fasta")
                with open(out_file_name, "w") as out_chr_fasta:
                    Bio.SeqIO.write(seq_obj, out_chr_fasta, "fasta")

                # make .fsa files
                fasta2fsa(out_file_name, annotDir, biosample=biosample.get(sample_base))
                # remove the .fasta file
                os.unlink(out_file_name)

                # make .src files
                if master_source_table:
                    out_src_fname = os.path.join(annotDir, sample + '.src')
                    with open(master_source_table, 'rt') as inf:
                        with open(out_src_fname, 'wt') as outf:
                            outf.write(inf.readline())
                            for line in inf:
                                row = line.rstrip('\n').split('\t')
                                if row[0] == sample_base or row[0] == sample:
                                    row[0] = sample
                                    outf.write('\t'.join(row) + '\n')

                # make .cmt files
                make_structured_comment_file(os.path.join(annotDir, sample + '.cmt'),
                                             name=sample,
                                             coverage=coverage.get(sample),
                                             seq_tech=sequencing_tech)

    # run tbl2asn (relies on filesnames matching by prefix)
    tbl2asn = tools.tbl2asn.Tbl2AsnTool()
    source_quals = []
    if organism:
        source_quals.append(('organism', organism))
    if mol_type:
        source_quals.append(('mol_type', mol_type))
    tbl2asn.execute(templateFile, annotDir, comment=comment,
        per_genome_comment=True, source_quals=source_quals)


def parser_prep_genbank_files(parser=argparse.ArgumentParser()):
    parser.add_argument('templateFile', help='Submission template file (.sbt) including author and contact info')
    parser.add_argument("fasta_files", nargs='+', help="Input fasta files")
    parser.add_argument("annotDir",
                        help="Output directory with genbank submission files (.tbl files must already be there)")
    parser.add_argument('--comment', default=None, help='comment field')
    parser.add_argument('--sequencing_tech', default=None, help='sequencing technology (e.g. Illumina HiSeq 2500)')
    parser.add_argument('--master_source_table', default=None, help='source modifier table')
    parser.add_argument('--organism', default=None, help='species name')
    parser.add_argument('--mol_type', default=None, help='molecule type')
    parser.add_argument("--biosample_map",
                        help="""A file with two columns and a header: sample and BioSample.
        This file may refer to samples that are not included in this submission.""")
    parser.add_argument('--coverage_table',
                        default=None,
                        help='''A genome coverage report file with a header row.  The table must
        have at least two columns named sample and aln2self_cov_median.  All other
        columns are ignored. Rows referring to samples not in this submission are
        ignored.''')
    util.cmd.common_args(parser, (('tmp_dir', None), ('loglevel', None), ('version', None)))
    util.cmd.attach_main(parser, prep_genbank_files, split_args=True)
    return parser


__commands__.append(('prep_genbank_files', parser_prep_genbank_files))


def prep_sra_table(lib_fname, biosampleFile, md5_fname, outFile):
    ''' This is a very lazy hack that creates a basic table that can be
        pasted into various columns of an SRA submission spreadsheet.  It probably
        doesn't work in all cases.
    '''
    metadata = {}

    with util.file.compressed_open(biosampleFile, 'rt') as inf:
        header = inf.readline()
        for line in inf:
            row = line.rstrip('\n\r').split('\t')
            metadata.setdefault(row[0], {})
            metadata[row[0]]['biosample_accession'] = row[1]

    with util.file.compressed_open(md5_fname, 'rt') as inf:
        for line in inf:
            row = line.rstrip('\n\r').split()
            s = os.path.basename(row[1])
            assert s.endswith('.cleaned.bam')
            s = s[:-len('.cleaned.bam')]
            metadata.setdefault(s, {})
            metadata[s]['filename'] = row[1]
            metadata[s]['MD5_checksum'] = row[0]

    with open(outFile, 'wt') as outf:
        header = ['biosample_accession', 'sample_name', 'library_ID', 'filename', 'MD5_checksum']
        outf.write('\t'.join(header) + '\n')
        with open(lib_fname, 'rU') as inf:
            for line in inf:
                lib = line.rstrip('\n\r')
                parts = lib.split('.')
                assert len(parts) > 1 and parts[-1].startswith('l')
                s = '.'.join(parts[:-1])
                metadata.setdefault(s, {})
                metadata[s]['library_ID'] = lib
                metadata[s]['sample_name'] = s
                outf.write('\t'.join(metadata[s].get(h, '') for h in header) + '\n')


def parser_prep_sra_table(parser=argparse.ArgumentParser()):
    parser.add_argument('lib_fname',
                        help='A file that lists all of the library IDs that will be submitted in this batch')
    parser.add_argument("biosampleFile",
                        help="""A file with two columns and a header: sample and BioSample.
        This file may refer to samples that are not included in this submission.""")
    parser.add_argument("md5_fname",
                        help="""A file with two columns and no header.  Two columns are MD5 checksum and filename.
        Should contain an entry for every bam file being submitted in this batch.
        This is typical output from "md5sum *.cleaned.bam".""")
    parser.add_argument("outFile",
                        help="Output table that contains most of the variable columns needed for SRA submission.")
    util.cmd.common_args(parser, (('loglevel', None), ('version', None)))
    util.cmd.attach_main(parser, prep_sra_table, split_args=True)
    return parser


__commands__.append(('prep_sra_table', parser_prep_sra_table))


def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)


if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)
