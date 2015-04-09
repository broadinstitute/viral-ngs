#!/usr/bin/env python
'''This script contains a number of utilities for submitting our analyses
to NCBI's Genbank and SRA databases.
'''

__author__ = "PLACEHOLDER"
__commands__ = []

import argparse, logging, collections
import Bio.SeqIO
import util.cmd, util.file, util.vcf, util.misc
import interhost

log = logging.getLogger(__name__)


def tbl_reader(text_lines):
    ''' This function reads an iterator of unformatted text lines that represent
        features in NCBI's "tbl" format.  It returns them as an iterator of
        features:
            (seqid, start, stop, type, child_records)
        child_records is an ordered list of tuples:
            (key, value)
        This function accepts files that describe multiple sequences.
    '''
    (cur_seqid, start, stop, ftype, children) = (None, None, None, None, None)
    for line_num, line in enumerate(text_lines):
        try:
            line = line.rstrip('\n\r')
            if not line:
                pass
            elif line.startswith('>'):
                # new sequence
                if not line.startswith('>Feature '):
                    raise Exception("not sure how to handle a non-Feature record")
                else:
                    if cur_seqid and start and stop:
                        yield (cur_seqid, start, stop, ftype, children)
                    cur_seqid = line[len('>Feature '):].strip()
                    start, stop = (None, None)
            elif line[0] not in (' ', '\t'):
                # new feature
                if not cur_seqid:
                    raise Exception("feature started before the sequence id was specified")
                if cur_seqid and start and stop:
                    yield (cur_seqid, start, stop, ftype, children)
                start, stop, ftype = line.split('\t')
                start = int(start)
                stop = int(stop)
                children = []
            else:
                # child info
                k, v = line.strip().split('\t')
                children.append((k,v))
        except:
            log.exception("error parsing tbl file at line {}, contents: '{}'".format(
                line_num, line))
            raise
    # at EOF, emit last feature
    if cur_seqid and start and stop:
        yield (cur_seqid, start, stop, ftype, children)


def fasta_chrlens(fasta):
    out = collections.OrderedDict()
    with open(fasta, 'rt') as inf:
        for seq in Bio.SeqIO.parse(inf, 'fasta'):
            out[seq.id] = len(seq)
    return out

def tbl_transfer(ref_fasta, ref_tbl, alt_fasta, out_tbl, oob_clip=False):
    ''' This function takes an NCBI TBL file describing features on a genome
        (genes, etc) and transfers them to a new genome.
    '''
    cmap = interhost.CoordMapper(ref_fasta, alt_fasta)
    alt_chrlens = fasta_chrlens(alt_fasta)
    
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
                    # feature with numeric coordinates (translate them)
                    row = line.split('\t')
                    if not len(row)>=2:
                        raise Exception("this line has only one column?")
                    row[0] = int(row[0])
                    row[1] = int(row[1])
                    
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
                                row[0] = 1
                            if row[1]==None:
                                row[1] = alt_chrlens[altid]
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



def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)
if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)
