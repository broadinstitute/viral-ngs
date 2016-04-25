#!/usr/bin/env python
''' This script contains a number of utilities for metagenomic analyses.
'''

__author__ = "yesimon@broadinstitute.org"

import argparse
import gzip
import shutil
import util.file
import tools.kraken
import tools.krona
import tools.diamond
import tools.picard
__commands__ = []


def kraken(inBam, db, outReport=None, outReads=None,
           filterThreshold=None, numThreads=1):
    assert outReads or outReport, (
        'Either --outReads or --outReport must be specified.')

    tmp_fastq1 = util.file.mkstempfname('.1.fastq')
    tmp_fastq2 = util.file.mkstempfname('.2.fastq')
    picard = tools.picard.SamToFastqTool()
    picard_opts = {
        'CLIPPING_ATTRIBUTE': tools.picard.SamToFastqTool.illumina_clipping_attribute,
        'CLIPPING_ACTION': 'X'
    }
    picard.execute(inBam, tmp_fastq1, tmp_fastq2,
                   picardOptions=tools.picard.PicardTools.dict_to_picard_opts(picard_opts),
                   JVMmemory=picard.jvmMemDefault)

    kraken_tool = tools.kraken.Kraken()
    tmp_reads = util.file.mkstempfname('.kraken')
    opts = {
        '--paired': None,
        '--threads': numThreads,
    }
    # Could be optimized in 3.5 piping directly to kraken-filter.
    kraken_tool.classify(db, [tmp_fastq1, tmp_fastq2], tmp_reads, options=opts)

    if filterThreshold:
        opts = {
            '--threshold': filterThreshold,
        }

        tmp_filtered_reads = util.file.mkstempfname('.filtered-kraken')
        kraken_tool.execute('kraken-filter', db, tmp_filtered_reads, args=[tmp_reads],
                            options=opts)
    else:
        tmp_filtered_reads = tmp_reads

    if outReads:
        with open(tmp_filtered_reads, 'rb') as f_in:
            with gzip.open(outReads, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

    if outReport:
        kraken_tool.execute('kraken-report', db, outReport,
                            args=[tmp_filtered_reads])

def krona(inTsv, outHtml, queryColumn=None, taxidColumn=None,
          scoreColumn=None, noHits=None, noRank=None):

    krona_tool = tools.krona.Krona()
    if inTsv.endswith('.gz'):
        tmp_tsv = util.file.mkstempfname('.tsv')
        with gzip.open(inTsv, 'rb') as f_in:
            with open(tmp_tsv, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
                to_import = [tmp_tsv]
    else:
        to_import = [inTsv]

    krona_tool.import_taxonomy(
        to_import, outHtml, query_column=queryColumn, taxid_column=taxidColumn,
        score_column=scoreColumn, no_hits=noHits, no_rank=noRank)


def parser_kraken(parser=argparse.ArgumentParser()):
    parser.add_argument('inBam', help='Input unaligned reads, BAM format.')
    parser.add_argument('db', help='Kraken database directory.')
    parser.add_argument('--outReport', help='Kraken report output file.')
    parser.add_argument('--outReads', help='Kraken per read output file.')
    parser.add_argument('--filterThreshold',
                        default=0.05,
                        type=float,
                        help='Kraken filter threshold (default %(default)s)')
    parser.add_argument('--numThreads', default=1, help='Number of threads to run. (default %(default)s)')
    util.cmd.common_args(parser, (('loglevel', None), ('version', None),
                                  ('tmp_dir', None)))
    util.cmd.attach_main(parser, kraken, split_args=True)
    return parser


def parser_krona(parser=argparse.ArgumentParser()):
    parser.add_argument('inTsv', help='Input tab delimited file.')
    parser.add_argument('outHtml', help='Output html report.')
    parser.add_argument('--queryColumn', help='Column of query id. (default %(default)s)',
                        type=int, default=2)
    parser.add_argument('--taxidColumn', help='Column of taxonomy id. (default %(default)s)',
                        type=int, default=3)
    parser.add_argument('--scoreColumn', help='Column of score. (default %(default)s)',
                        type=int)
    parser.add_argument('--noHits', help='Include wedge for no hits.',
                        action='store_true')
    parser.add_argument('--noRank', help='Include no rank assignments.',
                        action='store_true')
    util.cmd.common_args(parser, (('loglevel', None), ('version', None)))
    util.cmd.attach_main(parser, krona, split_args=True)
    return parser

def diamond(inBam, db, outM8, numThreads=1):
    tmp_fastq = util.file.mkstempfname('.fastq')
    tmp_fastq2 = util.file.mkstempfname('.fastq')
    picard = tools.picard.SamToFastqTool()
    picard_opts = {
        'CLIPPING_ATTRIBUTE': tools.picard.SamToFastqTool.illumina_clipping_attribute,
        'CLIPPING_ACTION': 'X'
    }
    picard.execute(inBam, tmp_fastq, tmp_fastq2,
                   picardOptions=tools.picard.PicardTools.dict_to_picard_opts(picard_opts),
                   JVMmemory=picard.jvmMemDefault)


    diamond_tool = tools.diamond.Diamond()
    diamond_tool.install()
    tmp_alignment = util.file.mkstempfname('.daa')
    tmp_m8 = util.file.mkstempfname('.diamond.m8')
    diamond_tool.blastx(db, [tmp_fastq, tmp_fastq2], tmp_alignment,
                        options={'--threads': numThreads})
    diamond_tool.view(tmp_alignment, tmp_m8,
                      options={'--threads': numThreads})
    with open(tmp_m8, 'rb') as f_in:
        with gzip.open(outM8, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)


def parser_diamond(parser=argparse.ArgumentParser()):
    parser.add_argument('inBam', help='Input unaligned reads, BAM format.')
    parser.add_argument('db', help='Diamond database directory.')
    parser.add_argument('outM8', help='Blast m8 formatted output file.')
    parser.add_argument('--numThreads', default=1, help='Number of threads (default: %(default)s)')
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, diamond, split_args=True)
    return parser

__commands__.append(('kraken', parser_kraken))
__commands__.append(('diamond', parser_diamond))
__commands__.append(('krona', parser_krona))

def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)

if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)
