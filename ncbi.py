#!/usr/bin/env python
'''This script contains a number of utilities for submitting our analyses
to NCBI's Genbank and SRA databases.
'''

__author__ = "PLACEHOLDER"
__commands__ = []

import argparse, logging
import util.cmd, util.file, util.vcf, util.misc

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




def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)
if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)
