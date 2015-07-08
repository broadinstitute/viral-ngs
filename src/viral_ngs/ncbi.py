#!/usr/bin/env python
'''This script contains a number of utilities for submitting our analyses
to NCBI's Genbank and SRA databases.
'''

__author__ = "PLACEHOLDER"
__commands__ = []

import argparse, logging, collections, shutil, os, os.path
import Bio.SeqIO
import util.cmd, util.file, util.version
import tools.tbl2asn
import interhost

log = logging.getLogger(__name__)



def fasta_chrlens(fasta):
    out = collections.OrderedDict()
    with open(fasta, 'rt') as inf:
        for seq in Bio.SeqIO.parse(inf, 'fasta'):
            out[seq.id] = len(seq)
    return out

def tbl_transfer_common(cmap, ref_tbl, out_tbl, alt_chrlens, oob_clip=False):
    """
        This function is the feature transfer machinery used by tbl_transfer() 
        and tbl_transfer_prealigned(). cmap is an instance of CoordMapper.
    """

    with open(ref_tbl, 'rt') as inf:
        with open(out_tbl, 'wt') as outf:
            for line in inf:
                line = line.rstrip('\r\n')
                if not line:
                    pass
                elif line.startswith('>'):
                    # sequence identifier
                    if not line.startswith('>Feature '):
                        raise Exception("not sure how to handle a non-Feature record")
                    seqid = line[len('>Feature '):].strip()
                    if not (seqid.startswith('gb|') and seqid.endswith('|') and len(seqid)>4):
                        raise Exception("reference annotation does not refer to a GenBank accession")
                    seqid = seqid[3:-1]
                    altid = cmap.mapAtoB(seqid)
                    line = '>Feature ' + altid
                    feature_keep = True
                elif line[0] != '\t':
                    # feature with numeric coordinates (map them)
                    row = line.split('\t')
                    if not len(row)>=2:
                        raise Exception("this line has only one column?")
                    row[0] = int(row[0])
                    row[1] = int(row[1])
                    if row[1] >= row[0]:
                        row[0] = cmap.mapAtoB(seqid, row[0], -1)[1]
                        row[1] = cmap.mapAtoB(seqid, row[1],  1)[1]
                    else:
                        # negative strand feature
                        row[0] = cmap.mapAtoB(seqid, row[0],  1)[1]
                        row[1] = cmap.mapAtoB(seqid, row[1], -1)[1]
                    
                    if row[0] and row[1]:
                        feature_keep = True
                    elif row[0]==None and row[1]==None:
                        # feature completely out of bounds
                        feature_keep = False
                        continue
                    else:
                        # feature overhangs end of sequence
                        if oob_clip:
                            if row[0]==None:
                                row[0] = '<1'
                            if row[1]==None:
                                row[1] = '>{}'.format(alt_chrlens[altid])
                        else:
                            feature_keep = False
                            continue
                    line = '\t'.join(map(str,row))
                else:
                    # feature notes
                    if not feature_keep:
                        # skip any lines that follow a skipped feature
                        continue
                    elif 'protein_id' in line:
                        # skip any lines that refer to an explicit protein_id
                        continue
                outf.write(line+'\n')

def tbl_transfer(ref_fasta, ref_tbl, alt_fasta, out_tbl, oob_clip=False):
    ''' This function takes an NCBI TBL file describing features on a genome
        (genes, etc) and transfers them to a new genome.
    '''
    cmap = interhost.CoordMapper()
    cmap.align_and_load_sequences([ref_fasta, alt_fasta])
    alt_chrlens = fasta_chrlens(alt_fasta)
    
    tbl_transfer_common(cmap, ref_tbl, out_tbl, alt_chrlens, oob_clip)

def parser_tbl_transfer(parser=argparse.ArgumentParser()):
    parser.add_argument("ref_fasta", help="Input sequence of reference genome")
    parser.add_argument("ref_tbl", help="Input reference annotations (NCBI TBL format)")
    parser.add_argument("alt_fasta", help="Input sequence of new genome")
    parser.add_argument("out_tbl", help="Output file with transferred annotations")
    parser.add_argument('--oob_clip', default=False, action='store_true',
        help='''Out of bounds feature behavior.
        False: drop all features that are completely or partly out of bounds
        True:  drop all features completely out of bounds
               but truncate any features that are partly out of bounds''')
    util.cmd.common_args(parser, (('tmpDir',None), ('loglevel',None), ('version',None)))
    util.cmd.attach_main(parser, tbl_transfer, split_args=True)
    return parser
__commands__.append(('tbl_transfer', parser_tbl_transfer))

def tbl_transfer_prealigned(input_fasta, ref_seq_name, alt_seq_name, ref_tbl, out_tbl, oob_clip=False):
    """
        This breaks out the ref and alt sequences into separate fasta files, and then
        creates a unified files containing the reference sequence first and the alt second
    """
    ref_fasta_filename = ""
    alt_fasta_filename = ""
    combined_fasta_filename = ""

    # write out the desired sequences to separate fasta files
    with util.file.open_or_gzopen(input_fasta, 'r') as inf:
        for seq in Bio.SeqIO.parse(inf, 'fasta'):
            if seq.id == ref_seq_name:
                ref_fasta_filename = util.file.mkstempfname('.fasta')
                with open(ref_fasta_filename, 'wt') as outf:
                    for line in util.file.fastaMaker(seq):
                        outf.write(line)
            if seq.id == alt_seq_name:
                alt_fasta_filename = util.file.mkstempfname('.fasta')
                with open(alt_fasta_filename, 'wt') as outf:
                    for line in util.file.fastaMaker(seq):
                        outf.write(line)
    if ref_fasta_filename == "": 
        raise KeyError("The ref sequence '%s' was not found in the file %s" % (ref_seq_name, input_fasta) )
    if alt_fasta_filename == "":
        raise KeyError("The alt sequence '%s' was not found in the file %s" % (alt_seq_name, input_fasta) )

    with open(combined_fasta_filename, 'wt') as outf:
        with open(ref_fasta_filename, 'wt') as reff:
            for line in reff:
                outf.write(line)
        with open(alt_fasta_filename, 'wt') as altf:
            for line in altf:
                outf.write(line)

    cmap = interhost.CoordMapper()
    cmap.load_alignments([combined_fasta_filename])
    alt_chrlens = fasta_chrlens(alt_fasta_filename)
    
    tbl_transfer_common(cmap, ref_tbl, out_tbl, alt_chrlens, oob_clip)

def parser_tbl_transfer_prealigned(parser=argparse.ArgumentParser()):
    parser.add_argument("input_fasta", help="FASTA file containing input sequences, including pre-made alignments and reference sequence")
    parser.add_argument("ref_seq_name", help="Name of the reference sequence within the input sequences")
    parser.add_argument("alt_seq_name", help="Name of the alternate sequence within the input sequences")
    parser.add_argument("ref_tbl", help="Input reference annotations (NCBI TBL format)")
    parser.add_argument("out_tbl", help="Output file with transferred annotations")
    parser.add_argument('--oob_clip', default=False, action='store_true',
        help='''Out of bounds feature behavior.
        False: drop all features that are completely or partly out of bounds
        True:  drop all features completely out of bounds
               but truncate any features that are partly out of bounds''')
    util.cmd.common_args(parser, (('tmpDir',None), ('loglevel',None), ('version',None)))
    util.cmd.attach_main(parser, tbl_transfer_prealigned, split_args=True)
    return parser
__commands__.append(('tbl_transfer_prealigned', parser_tbl_transfer_prealigned))    

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
        outf.write("Assembly Method\tgithub.com/broadinstitute/viral-ngs v. {}\n".format(
            util.version.get_version()))
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
        coverage_table=None, biosample_map=None):
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
        sample = os.path.basename(fn)[:-6]
        # make .fsa files
        fasta2fsa(fn, annotDir, biosample=biosample.get(sample))
        # make .src files
        if master_source_table:
            shutil.copy(master_source_table, os.path.join(annotDir, sample+'.src'))
        # make .cmt files
        make_structured_comment_file(os.path.join(annotDir, sample+'.cmt'),
            name=sample, coverage=coverage.get(sample), seq_tech=sequencing_tech)
    
    # run tbl2asn
    tbl2asn = tools.tbl2asn.Tbl2AsnTool()
    tbl2asn.execute(templateFile, annotDir, comment=comment, per_genome_comment=True)

def parser_prep_genbank_files(parser=argparse.ArgumentParser()):
    parser.add_argument('templateFile',
        help='Submission template file (.sbt) including author and contact info')
    parser.add_argument("fasta_files", nargs='+',
        help="Input fasta files")
    parser.add_argument("annotDir",
        help="Output directory with genbank submission files (.tbl files must already be there)")
    parser.add_argument('--comment', default=None,
        help='comment field')
    parser.add_argument('--sequencing_tech', default=None,
        help='sequencing technology (e.g. Illumina HiSeq 2500)')
    parser.add_argument('--master_source_table', default=None,
        help='source modifier table')
    parser.add_argument("--biosample_map",
        help="""A file with two columns and a header: sample and BioSample.
        This file may refer to samples that are not included in this submission.""")
    parser.add_argument('--coverage_table', default=None,
        help='''A genome coverage report file with a header row.  The table must
        have at least two columns named sample and aln2self_cov_median.  All other
        columns are ignored. Rows referring to samples not in this submission are
        ignored.''')
    util.cmd.common_args(parser, (('tmpDir',None), ('loglevel',None), ('version',None)))
    util.cmd.attach_main(parser, prep_genbank_files, split_args=True)
    return parser
__commands__.append(('prep_genbank_files', parser_prep_genbank_files))


def prep_sra_table(lib_fname, biosampleFile, md5_fname, outFile):
    ''' This is a very lazy hack that creates a basic table that can be
        pasted into various columns of an SRA submission spreadsheet.  It probably
        doesn't work in all cases.
    '''
    metadata = {}

    with open(biosampleFile, 'rU') as inf:
        header = inf.readline()
        for line in inf:
            row = line.rstrip('\n\r').split('\t')
            metadata.setdefault(row[0], {})
            metadata[row[0]]['biosample_accession'] = row[1]

    with open(md5_fname, 'rU') as inf:
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
        outf.write('\t'.join(header)+'\n')
        with open(lib_fname, 'rU') as inf:
            for line in inf:
                lib = line.rstrip('\n\r')
                parts = lib.split('.')
                assert len(parts)>1 and parts[-1].startswith('l')
                s = '.'.join(parts[:-1])
                metadata.setdefault(s, {})
                metadata[s]['library_ID'] = lib
                metadata[s]['sample_name'] = s
                outf.write('\t'.join(metadata[s].get(h,'') for h in header)+'\n')

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
    util.cmd.common_args(parser, (('loglevel',None), ('version',None)))
    util.cmd.attach_main(parser, prep_sra_table, split_args=True)
    return parser
__commands__.append(('prep_sra_table', parser_prep_sra_table))


def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)
if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)
