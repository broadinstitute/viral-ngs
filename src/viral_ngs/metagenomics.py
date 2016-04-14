#!/usr/bin/env python
import argparse
import gzip
import shutil
import util.file
import tools.kraken
import tools.krona
import tools.diamond
__commands__ = []


def kraken(inFastq1, inFastq2, db, outReport=None, outReads=None,
           filterThreshold=None, threads=1):
    assert outReads or outReport, (
        'Either --outReads or --outReport must be specified.')

    kraken_tool = tools.kraken.Kraken()
    tmp_reads = util.file.mkstempfname('.kraken')
    tmp_filtered_reads = util.file.mkstempfname('.filtered-kraken')
    opts = {
        '--paired': None,
    }
    # Could be optimized in 3.5 piping directly to kraken-filter.
    kraken_tool.classify(db, [inFastq1, inFastq2], tmp_reads, options=opts)
    opts = {
        '--threshold': filterThreshold,
    }

    kraken_tool.execute('kraken-filter', db, tmp_filtered_reads, args=[tmp_reads],
                        options=opts)

    if outReads:
        with open(tmp_filtered_reads, 'rb') as f_in:
            with gzip.open(outReads, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

    if outReport:
        kraken_tool.execute('kraken-report', db, outReport, args=[outReads],
                            options=opts)

def krona(inTsv, outHtml, queryColumn=None, taxidColumn=None,
          scoreColumn=None, noHits=None, noRank=None):
    krona_tool = tools.krona.Krona()
    krona_tool.import_taxonomy(
        [inTsv], outHtml, query_column=queryColumn, taxid_column=taxidColumn,
        score_column=scoreColumn, no_hits=noHits, no_rank=noRank)


def parser_kraken(parser=argparse.ArgumentParser()):
    parser.add_argument('inFastq1', help='First end of paired fastq.')
    parser.add_argument('inFastq2', help='First end of paired fastq.')
    parser.add_argument('db', help='Kraken database directory.')
    parser.add_argument('--outReport', help='Kraken report output file.')
    parser.add_argument('--outReads', help='Kraken per read output file.')
    parser.add_argument('--filterThreshold',
                        default=0.05,
                        type=float,
                        help='Kraken filter threshold')
    util.cmd.common_args(parser, (('loglevel', None), ('version', None),
                                  ('tmp_dir', None)))
    util.cmd.attach_main(parser, kraken, split_args=True)


def parser_krona(parser=argparse.ArgumentParser()):
    parser.add_argument('inTsv', help='Input tab delimited file.')
    parser.add_argument('outHtml', help='Output html report.')
    parser.add_argument('--queryColumn', help='Column of query id.',
                        type=int)
    parser.add_argument('--taxidColumn', help='Column of taxonomy id.',
                        type=int)
    parser.add_argument('--scoreColumn', help='Column of score.',
                        type=int)
    parser.add_argument('--noHits', help='Include wedge for no hits.',
                        action='store_true')
    parser.add_argument('--noRank', help='Include no rank assignments.',
                        action='store_true')
    util.cmd.common_args(parser, (('loglevel', None), ('version', None)))
    util.cmd.attach_main(parser, krona, split_args=True)
    return parser


__commands__.append(('kraken', parser_kraken))
__commands__.append(('krona', parser_krona))

if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)
