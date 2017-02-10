#!/usr/bin/env python
''' Functions to create reports from genomics pipeline data.
'''

__author__ = "dpark@broadinstitute.org"
__commands__ = []

import argparse
import logging
import glob
import os
import time
from collections import OrderedDict
import csv
import math
import shutil

import pysam
from pybedtools import BedTool
import Bio.SeqIO
import Bio.AlignIO
from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import util.cmd
import util.file
import util.misc
import tools.samtools
import tools.bwa
import assembly
import interhost
from util.stats import mean, median

log = logging.getLogger(__name__)


def get_assembly_stats(sample,
                       cov_thresholds=(1, 5, 20, 100),
                       assembly_dir='data/02_assembly', assembly_tmp='tmp/02_assembly',
                       align_dir='data/02_align_to_self', reads_dir='data/01_per_sample',
                       raw_reads_dir='data/00_raw'):
    ''' Fetch assembly-level statistics for a given sample '''
    out = {'sample': sample}
    samtools = tools.samtools.SamtoolsTool()
    header = ['sample',
              'reads_raw',
              'reads_cleaned',
              'reads_taxfilt',
              'assembled_trinity',
              'trinity_in_reads',
              'n_contigs',
              'contig_len',
              'unambig_bases',
              'pct_unambig',
              'aln2self_reads_tot',
              'aln2self_reads_aln',
              'aln2self_reads_rmdup',
              'aln2self_pct_nondup',
              'aln2self_cov_median',
              'aln2self_cov_mean',
              'aln2self_cov_mean_non0',] + ['aln2self_cov_%dX' % t for t in cov_thresholds]

    # per-sample unaligned read stats
    for adj in ('cleaned', 'taxfilt'):
        reads_bam = os.path.join(reads_dir, '.'.join((sample, adj, 'bam')))
        if os.path.isfile(reads_bam):
            out['reads_' + adj] = samtools.count(reads_bam)
    if os.path.isdir(raw_reads_dir):
        out['reads_raw'] = sum(samtools.count(bam)
            # correct issue where sample names containing other sample names as substrings leads
            # to extra files being included in the count
            #
            # add a dot before the wildcard, and assume the sample name is found before the dot.
            # this works for now since dots are the filename field separators
            # and leading/trailing dots are stripped from sample names in util.file.string_to_file_name()
            # TODO: replace this with better filtering?
            for bam in glob.glob(os.path.join(raw_reads_dir, sample + ".*.bam")))
        sample_raw_fname = os.path.join(raw_reads_dir, sample + ".bam")
        if os.path.isfile(sample_raw_fname):
            # if "00_raw/sample.bam" exists, these were not demuxed by snakemake
            if out['reads_raw']:
                # if sample.bam AND sample.library.flowcell.lane.bam exist, we have a problem!
                out['reads_raw'] = 'ambiguous filenames in raw reads directory!'
            else:
                # just count the sample.bam reads
                out['reads_raw'] = samtools.count(sample_raw_fname)

    # pre-assembly stats
    out['assembled_trinity'] = os.path.isfile(os.path.join(assembly_tmp, sample +
                                                           '.assembly1-trinity.fasta')) and 1 or 0
    sub_bam = os.path.join(assembly_tmp, sample + '.subsamp.bam')
    if os.path.isfile(sub_bam):
        out['trinity_in_reads'] = samtools.count(sub_bam)

    # assembly stats
    assembly_fname = os.path.join(assembly_dir, sample + '.fasta')
    if not os.path.isfile(assembly_fname):
        assembly_fname = os.path.join(assembly_tmp, sample + '.assembly2-scaffolded.fasta')
        if not os.path.isfile(assembly_fname):
            out['n_contigs'] = 0
    if os.path.isfile(assembly_fname):
        with open(assembly_fname, 'rt') as inf:
            counts = [(len(s), assembly.unambig_count(s.seq)) for s in Bio.SeqIO.parse(inf, 'fasta') if len(s) > 0]
        out['n_contigs'] = len(counts)
        out['contig_len'] = ','.join(str(x) for x, y in counts)
        out['unambig_bases'] = ','.join(str(y) for x, y in counts)
        out['pct_unambig'] = ','.join(str(float(y) / x) for x, y in counts)

    # read counts from align-to-self
    bam_fname = os.path.join(align_dir, sample + '.bam')
    if os.path.isfile(bam_fname):
        out['aln2self_reads_tot'] = samtools.count(bam_fname)
        out['aln2self_reads_aln'] = samtools.count(bam_fname, opts=['-F', '4'])
        out['aln2self_reads_rmdup'] = samtools.count(bam_fname, opts=['-F', '1028'])
        if out['aln2self_reads_aln']:
            out['aln2self_pct_nondup'] = float(out['aln2self_reads_rmdup']) / out['aln2self_reads_aln']

    # genome coverage stats
    bam_fname = os.path.join(align_dir, sample + '.mapped.bam')
    if os.path.isfile(bam_fname):
        with pysam.AlignmentFile(bam_fname, 'rb') as bam:
            coverages = list([pcol.nsegments for pcol in bam.pileup()])
        if coverages:
            out['aln2self_cov_median'] = median(coverages)
            out['aln2self_cov_mean'] = "%0.3f" % mean(coverages)
            out['aln2self_cov_mean_non0'] = "%0.3f" % mean([n for n in coverages if n > 0])
            for thresh in cov_thresholds:
                out['aln2self_cov_%dX' % thresh] = sum(1 for n in coverages if n >= thresh)

    return (header, out)


def assembly_stats(samples, outFile, cov_thresholds, assembly_dir, assembly_tmp, align_dir, reads_dir, raw_reads_dir):
    ''' Fetch assembly-level statistics for a given sample '''
    header_written = False
    with open(outFile, 'wt') as outf:
        for sample in samples:
            log.info("fetching stats on " + sample)
            header, out = get_assembly_stats(sample,
                                             cov_thresholds=cov_thresholds,
                                             assembly_dir=assembly_dir,
                                             assembly_tmp=assembly_tmp,
                                             align_dir=align_dir,
                                             reads_dir=reads_dir,
                                             raw_reads_dir=raw_reads_dir)
            if not header_written:
                outf.write('\t'.join(map(str, header)) + '\n')
                header_written = True
            outf.write('\t'.join([str(out.get(h, '')) for h in header]) + '\n')
            outf.flush()


def parser_assembly_stats(parser=argparse.ArgumentParser()):
    parser.add_argument('samples', nargs='+', help='Sample names.')
    parser.add_argument('outFile', help='Output report file.')
    parser.add_argument('--cov_thresholds',
                        nargs='+',
                        type=int,
                        default=(1, 5, 20, 100),
                        help='Genome coverage thresholds to report on. (default: %(default)s)')
    parser.add_argument('--assembly_dir',
                        default='data/02_assembly',
                        help='Directory with assembly outputs. (default: %(default)s)')
    parser.add_argument('--assembly_tmp',
                        default='tmp/02_assembly',
                        help='Directory with assembly temp files. (default: %(default)s)')
    parser.add_argument('--align_dir',
                        default='data/02_align_to_self',
                        help='Directory with reads aligned to own assembly. (default: %(default)s)')
    parser.add_argument('--reads_dir',
                        default='data/01_per_sample',
                        help='Directory with unaligned filtered read BAMs. (default: %(default)s)')
    parser.add_argument('--raw_reads_dir',
                        default='data/00_raw',
                        help='Directory with unaligned raw read BAMs. (default: %(default)s)')
    util.cmd.attach_main(parser, assembly_stats, split_args=True)
    return parser


__commands__.append(('assembly_stats', parser_assembly_stats))


def alignment_summary(inFastaFileOne, inFastaFileTwo, outfileName=None, printCounts=False):
    """ Write or print pairwise alignment summary information for sequences in two FASTA
        files, including SNPs, ambiguous bases, and indels.
    """
    gap = '-'
    ambiguous = 'N'
    aligner = tools.muscle.MuscleTool()

    per_chr_fastas = interhost.transposeChromosomeFiles([inFastaFileOne, inFastaFileTwo])

    results = OrderedDict()
    results["same_unambig"]  = 0
    results["snp_unambig"]   = 0
    results["indel_unambig"] = 0
    results["indel_ambig"]   = 0
    results["ambig_one"]     = 0
    results["ambig_two"]     = 0
    results["ambig_both"]    = 0
    results["unambig_both"]  = 0

    for chr_fasta in per_chr_fastas:
        same_unambig  = 0
        snp_unambig   = 0
        indel_unambig = 0
        indel_ambig   = 0
        ambig_one     = 0
        ambig_two     = 0
        ambig_both    = 0
        unambig_both  = 0

        alignOutFileName = util.file.mkstempfname('.fasta')
        aligner.execute(chr_fasta, alignOutFileName, fmt="clw")

        with open(alignOutFileName, "r") as f:
            alignment = Bio.AlignIO.read(f, "clustal")

            for col_idx in range(0, alignment.get_alignment_length()):
                col = alignment[:, col_idx]
                c1 = col[0]
                c2 = col[1]

                if (c1 in ambiguous
                   and c2 in ambiguous):
                    ambig_both +=1
                elif c1 in ambiguous:
                    ambig_one += 1
                elif c2 in ambiguous:
                    ambig_two += 1

                if (c1 in IUPACUnambiguousDNA().letters
                   and c2 in IUPACUnambiguousDNA().letters):
                    unambig_both += 1
                    if c1 == c2:
                        same_unambig += 1
                    else:
                        snp_unambig += 1

                if ((c1 == gap and
                    c2 in IUPACUnambiguousDNA().letters) or
                   (c2 == gap and
                    c1 in IUPACUnambiguousDNA().letters)):
                    indel_unambig += 1

                if ((c1 == gap and
                    c2 in ambiguous) or
                   (c2 == gap and
                    c1 in ambiguous)):
                    indel_ambig += 1

        if printCounts:
            print("Counts for this segment/chromosome:")
            print("same_unambig ", same_unambig)
            print("snp_unambig  ", snp_unambig)
            print("indel_unambig", indel_unambig)
            print("indel_ambig  ", indel_ambig)
            print("ambig_one    ", ambig_one)
            print("ambig_two    ", ambig_two)
            print("ambig_both   ", ambig_both)
            print("unambig_both ", unambig_both)

        results["same_unambig"]  += same_unambig
        results["snp_unambig"]   += snp_unambig
        results["indel_unambig"] += indel_unambig
        results["indel_ambig"]   += indel_ambig
        results["ambig_one"]     += ambig_one
        results["ambig_two"]     += ambig_two
        results["ambig_both"]    += ambig_both
        results["unambig_both"]  += unambig_both

    if printCounts:
        print("\nCounts for this sample:")
        print("same_unambig ", results["same_unambig"])
        print("snp_unambig  ", results["snp_unambig"])
        print("indel_unambig", results["indel_unambig"])
        print("indel_ambig  ", results["indel_ambig"])
        print("ambig_one    ", results["ambig_one"])
        print("ambig_two    ", results["ambig_two"])
        print("ambig_both   ", results["ambig_both"])
        print("unambig_both ", results["unambig_both"])

    if outfileName:
        with open(outfileName, "wt") as of:
            csvout = csv.writer(of, delimiter='\t')
            csvout.writerow(list(results.keys()))
            csvout.writerow(list(results.values()))

def parser_alignment_summary(parser=argparse.ArgumentParser()):
    parser.add_argument('inFastaFileOne', help='First fasta file for an alignment')
    parser.add_argument('inFastaFileTwo', help='First fasta file for an alignment')
    parser.add_argument('--outfileName', help='Output file for counts in TSV format')
    parser.add_argument('--printCounts', help='', action='store_true')
    util.cmd.attach_main(parser, alignment_summary, split_args=True)
    return parser
__commands__.append(('alignment_summary', parser_alignment_summary))


def consolidate_fastqc(inDirs, outFile):
    '''Consolidate multiple FASTQC reports into one.'''
    with util.file.open_or_gzopen(outFile, 'wt') as outf:
        header = ['Sample']
        out_n = 0
        for sdir in inDirs:
            out = {}
            with open(os.path.join(sdir, 'summary.txt'), 'rt') as inf:
                for line in inf:
                    v, k, fn = line.strip().split('\t')
                    out[k] = v
                    if out_n == 0:
                        header.append(k)
                    if not fn.endswith('.bam'):
                        raise TypeError("%s not a bam file" % fn)
                    out['Sample'] = fn[:-len('.bam')]
            if out_n == 0:
                outf.write('\t'.join(header) + '\n')
            outf.write('\t'.join([out.get(h, '') for h in header]) + '\n')
            out_n += 1


def parser_consolidate_fastqc(parser=argparse.ArgumentParser()):
    parser.add_argument('inDirs', help='Input FASTQC directories.', nargs='+')
    parser.add_argument('outFile', help='Output report file.')
    util.cmd.attach_main(parser, consolidate_fastqc, split_args=True)
    return parser


__commands__.append(('consolidate_fastqc', parser_consolidate_fastqc))


def get_bam_info(bamstats_dir):
    libs = {}
    for fn in glob.glob(os.path.join(bamstats_dir, "*.txt")):
        with util.file.open_or_gzopen(fn, 'rt') as inf:
            bam = {}
            for line in inf:
                k, v = line.rstrip('\n').split('\t')
                bam[k] = v
        libs.setdefault(bam['Sample'], {})
        libs[bam['Sample']][bam['BAM']] = bam['Total reads']
    return libs


def get_lib_info(runfile):
    libs = {}
    for lane in util.file.read_tabfile_dict(runfile):
        for well in util.file.read_tabfile_dict(lane['barcode_file']):
            libname = well['sample'] + '.l' + well['library_id_per_sample']
            libs.setdefault(libname, [])
            plate = well['Plate']
            if plate.lower().startswith('plate'):
                plate = plate[5:]
            well_id = well['Well'][0].upper() + "%02d" % int(well['Well'][1:])
            dat = [well['sample'], lane['flowcell'] + '.' + lane['lane'], well['barcode_1'] + '-' + well['barcode_2'],
                   plate.strip() + ':' + well_id, get_earliest_date(lane['bustard_dir']), well.get('Tube_ID', '')]
            libs[libname].append(dat)
    return libs


def get_earliest_date(inDir):
    fnames = [inDir] + [os.path.join(inDir, x) for x in os.listdir(inDir)]
    earliest = min(os.path.getmtime(fn) for fn in fnames)
    return time.strftime("%Y-%m-%d", time.localtime(earliest))


def consolidate_spike_count(inDir, outFile):
    '''Consolidate multiple spike count reports into one.'''
    with open(outFile, 'wt') as outf:
        for fn in os.listdir(inDir):
            fn = os.path.join(inDir, fn)
            s = os.path.basename(fn)
            if not s.endswith('.spike_count.txt'):
                raise Exception()
            s = s[:-len('.spike_count.txt')]
            with open(fn, 'rt') as inf:
                for line in inf:
                    if not line.startswith('Input bam') and not line.startswith('*'):
                        spike, count = [line.strip().split('\t')[i] for i in [0,2]]
                        outf.write('\t'.join([s, spike, count]) + '\n')


def parser_consolidate_spike_count(parser=argparse.ArgumentParser()):
    parser.add_argument('inDir', help='Input spike count directory.')
    parser.add_argument('outFile', help='Output report file.')
    util.cmd.attach_main(parser, consolidate_spike_count, split_args=True)
    return parser


__commands__.append(('consolidate_spike_count', parser_consolidate_spike_count))


# =========================


def parser_plot_coverage_common(parser=argparse.ArgumentParser()):    # parser needs add_help=False?
    parser.add_argument('in_bam', help='Input reads, BAM format.')
    parser.add_argument('out_plot_file', help='The generated chart file')
    parser.add_argument(
        '--plotFormat',
        dest="plot_format",
        default=None,
        type=str,
        choices=list(plt.gcf().canvas.get_supported_filetypes().keys()),
        metavar='',
        help="File format of the coverage plot. By default it is inferred from the file extension of out_plot_file, but it can be set explicitly via --plotFormat. Valid formats include: "
        + ", ".join(list(plt.gcf().canvas.get_supported_filetypes().keys()))
    )
    parser.add_argument(
        '--plotDataStyle',
        dest="plot_data_style",
        default="filled",
        type=str,
        choices=["filled", "line", "dots"],
        metavar='',
        help="The plot data display style. Valid options: " + ", ".join(["filled", "line", "dots"]) +
        " (default: %(default)s)"
    )
    parser.add_argument(
        '--plotStyle',
        dest="plot_style",
        default="ggplot",
        type=str,
        choices=plt.style.available,
        metavar='',
        help="The plot visual style. Valid options: " + ", ".join(plt.style.available) + " (default: %(default)s)"
    )
    parser.add_argument(
        '--plotWidth',
        dest="plot_width",
        default=880,
        type=int,
        help="Width of the plot in pixels (default: %(default)s)"
    )
    parser.add_argument(
        '--plotHeight',
        dest="plot_height",
        default=680,
        type=int,
        help="Width of the plot in pixels (default: %(default)s)"
    )
    parser.add_argument(
        '--plotDPI',
        dest="plot_dpi",
        default=plt.gcf().get_dpi(),
        type=int,
        help="dots per inch for rendered output, more useful for vector modes (default: %(default)s)"
    )
    parser.add_argument(
        '--plotTitle',
        dest="plot_title",
        default="Coverage Plot",
        type=str,
        help="The title displayed on the coverage plot (default: '%(default)s')"
    )
    parser.add_argument(
        '-q', dest="base_q_threshold",
        default=None, type=int,
        help="The minimum base quality threshold"
    )
    parser.add_argument(
        '-Q', dest="mapping_q_threshold",
        default=None,
        type=int, help="The minimum mapping quality threshold"
    )
    parser.add_argument(
        '-m',
        dest="max_coverage_depth",
        default=None,
        type=int,
        help="The max coverage depth (default: %(default)s)"
    )
    parser.add_argument('-l', dest="read_length_threshold", default=None, type=int, help="Read length threshold")
    parser.add_argument(
        '--outSummary',
        dest="out_summary",
        default=None,
        type=str,
        help="Coverage summary TSV file. Default is to write to temp."
    )
    return parser


def plot_coverage(
    in_bam,
    out_plot_file,
    plot_format,
    plot_data_style,
    plot_style,
    plot_width,
    plot_height,
    plot_dpi,
    plot_title,
    base_q_threshold,
    mapping_q_threshold,
    max_coverage_depth,
    read_length_threshold,
    plot_only_non_duplicates=False,
    out_summary=None
    ):
    ''' 
        Generate a coverage plot from an aligned bam file
    '''

    # TODO: remove this:
    #coverage_tsv_file = "/Users/tomkinsc/Downloads/plottest/test_multisegment.tsv"

    samtools = tools.samtools.SamtoolsTool()

    # check if in_bam is aligned, if not raise an error
    num_mapped_reads = samtools.count(in_bam, opts=["-F", "4"])
    if num_mapped_reads == 0:
        raise Exception(
            """The bam file specified appears to have zero mapped reads. 'plot_coverage' requires an aligned bam file. You can try 'align_and_plot_coverage' if you don't mind a simple bwa alignment. \n File: %s"""
            % in_bam
        )

    if out_summary is None:
        coverage_tsv_file = util.file.mkstempfname('.summary.tsv')
    else:
        coverage_tsv_file = out_summary

    bam_dupe_processed = util.file.mkstempfname('.dupe_processed.bam')
    if plot_only_non_duplicates:
        # TODO: this is probably not necessary since "samtools depth" does not count marked duplicates
        # write a new bam file; exclude reads with the 1024 flag set (PCR or optical duplicates)
        samtools.view(["-F", "1024"], in_bam, bam_dupe_processed)
    else:
        bam_dupe_processed = in_bam

    # call samtools sort
    bam_sorted = util.file.mkstempfname('.sorted.bam')
    samtools.sort(bam_dupe_processed, bam_sorted, args=["-O", "bam"])

    if plot_only_non_duplicates:
        os.unlink(bam_dupe_processed)

    # call samtools index
    samtools.index(bam_sorted)

    # call samtools depth
    opts = []
    opts += ['-aa']    # report coverate at "absolutely all" positions
    if base_q_threshold:
        if not plot_only_non_duplicates:
            # Note: "bedtools genomecov" will count depth including duplicates, but does
            # not expose options for filtering by quality. When duplicates
            # are excluded, "samtools depth" is used which does support quality filtering
            # We use either samtools or bedtools, because the former ignores marked duplicates
            # from its depth count while bedtools includes them. 
            log.warning("'-q' ignored since --plotOnlyNonDuplicates is absent")
        opts += ["-q", str(base_q_threshold)]
    if mapping_q_threshold:
        if not plot_only_non_duplicates:
            log.warning("'-Q' ignored since --plotOnlyNonDuplicates is absent")
        opts += ["-Q", str(mapping_q_threshold)]
    if max_coverage_depth:
        if not plot_only_non_duplicates:
            log.warning("'-m' ignored since --plotOnlyNonDuplicates is absent")
        opts += ["-m", str(max_coverage_depth)]
    if read_length_threshold:
        if not plot_only_non_duplicates:
            log.warning("'-l' ignored since --plotOnlyNonDuplicates is absent")
        opts += ["-l", str(read_length_threshold)]

    # add option here for bedtools to report coverage w/ duplicates 
    # (and then samtools for no-dups)
    #
    # Ex.
    #   samtools depth -aa mapped-to-ref.with-dups.tmp.bam
    #   bedtools genomecov -ibam mapped-to-ref.with-dups.tmp.bam -d
    if not plot_only_non_duplicates:
        bt = BedTool(bam_sorted)
        # "d=True" is the equivalent of passing "-d" to the bedtools CLI
        bt.genome_coverage(d=True).saveas(coverage_tsv_file)
    else:
        samtools.depth(bam_sorted, coverage_tsv_file, opts)
    os.unlink(bam_sorted)

    # ---- create plot based on coverage_tsv_file ----

    segment_depths = OrderedDict()
    domain_max = 0
    with open(coverage_tsv_file, "r") as tabfile:
        for row in csv.reader(tabfile, delimiter='\t'):
            segment_depths.setdefault(row[0], []).append(int(row[2]))
            domain_max += 1

    domain_max = 0
    with plt.style.context(plot_style):
        fig = plt.gcf()
        DPI = plot_dpi or fig.get_dpi()
        fig.set_size_inches(float(plot_width) / float(DPI), float(plot_height) / float(DPI))

        font_size = (2.5 * plot_height) / float(DPI)

        ax = plt.subplot()    # Defines ax variable by creating an empty plot

        # Set the tick labels font
        for label in (ax.get_xticklabels() + ax.get_yticklabels()):
            label.set_fontsize(font_size)

        for segment_num, (segment_name, position_depths) in enumerate(segment_depths.items()):
            prior_domain_max = domain_max
            domain_max += len(position_depths)

            colors = list(plt.rcParams['axes.prop_cycle'].by_key()['color'])    # get the colors for this style
            segment_color = colors[segment_num % len(colors)]    # pick a color, offset by the segment index

            if plot_data_style == "filled":
                plt.fill_between(
                    range(prior_domain_max, domain_max),
                    position_depths, [0] * len(position_depths),
                    linewidth=0,
                    antialiased=True,
                    color=segment_color
                )
            elif plot_data_style == "line":
                plt.plot(range(prior_domain_max, domain_max), position_depths, antialiased=True, color=segment_color)
            elif plot_data_style == "dots":
                plt.plot(
                    range(prior_domain_max, domain_max),
                    position_depths,
                    'ro',
                    antialiased=True,
                    color=segment_color
                )

        plt.title(plot_title, fontsize=font_size * 1.2)
        plt.xlabel("bp", fontsize=font_size * 1.1)
        plt.ylabel("read depth", fontsize=font_size * 1.1)

        # to squash a backend renderer error on OSX related to tight layout
        if plt.get_backend().lower() in ['agg', 'macosx']:
            fig.set_tight_layout(True)
        else:
            fig.tight_layout()

        plt.savefig(out_plot_file, format=plot_format, dpi=DPI)    #, bbox_inches='tight')
        log.info("Coverage plot saved to: " + out_plot_file)

    if not out_summary:
        os.unlink(coverage_tsv_file)


def parser_plot_coverage(parser=argparse.ArgumentParser()):
    parser = parser_plot_coverage_common(parser)
    parser.add_argument(
        '--plotOnlyNonDuplicates',
        dest="plot_only_non_duplicates",
        action="store_true",
        help="Plot only non-duplicates (samtools -F 1024), coverage counted by bedtools rather than samtools."
    )
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, plot_coverage, split_args=True)
    return parser


__commands__.append(('plot_coverage', parser_plot_coverage))


def align_and_plot_coverage(
    out_plot_file,
    plot_format,
    plot_data_style,
    plot_style,
    plot_width,
    plot_height,
    plot_dpi,
    plot_title,
    base_q_threshold,
    mapping_q_threshold,
    max_coverage_depth,
    read_length_threshold,
    out_summary,
    in_bam,
    ref_fasta,
    out_bam=None,
    sensitive=False,
    excludeDuplicates=False,
    JVMmemory=None,
    picardOptions=None,
    min_score_to_output=None,
    aligner="bwa",
    aligner_options='',
    novoalign_license_path=None
):
    ''' 
        Take reads, align to reference with BWA-MEM, and generate a coverage plot
    '''

    # TODO: use read_utils.py::align_and_fix in place of the duplicated alignment code here
    # The main difference is the presence/absence of GATK's local_realign

    if out_bam is None:
        bam_aligned = util.file.mkstempfname('.aligned.bam')
    else:
        bam_aligned = out_bam

    assert aligner in ["bwa", "novoalign"]
    if aligner_options is None:
        if aligner=="novoalign":
            aligner_options = '-r Random -l 40 -g 40 -x 20 -t 100 -k'
        elif aligner=='bwa':
            aligner_options = '-T 30' # quality threshold

    samtools = tools.samtools.SamtoolsTool()

    ref_indexed = util.file.mkstempfname('.reference.fasta')
    shutil.copyfile(ref_fasta, ref_indexed)

    aln_bam = util.file.mkstempfname('.bam')
    if aligner=="bwa":
        bwa = tools.bwa.Bwa()
        

        bwa.index(ref_indexed)

        bwa_opts = aligner_options.split()
        if sensitive:
            bwa_opts + "-k 12 -A 1 -B 1 -O 1 -E 1".split()

        # get the quality threshold from the opts
        # for downstream filtering
        bwa_map_threshold = min_score_to_output or 30
        if '-T' in bwa_opts:
            if bwa_opts.index("-T")+1 <= len(bwa_opts):
                bwa_map_threshold = int(bwa_opts[bwa_opts.index("-T")+1])

        bwa.align_mem_bam(in_bam, ref_indexed, aln_bam, options=bwa_opts, min_qual=bwa_map_threshold)
    elif aligner=="novoalign":
        
        tools.novoalign.NovoalignTool(license_path=novoalign_license_path).index_fasta(ref_indexed)

        tools.novoalign.NovoalignTool(license_path=novoalign_license_path).execute(
            in_bam, ref_indexed, aln_bam,
            options=aligner_options.split(),
            JVMmemory=JVMmemory
        )

    aln_bam_dupe_processed = util.file.mkstempfname('.filtered_dupe_processed.bam')
    if excludeDuplicates:
        opts = list(picardOptions)
        dupe_removal_out_metrics = util.file.mkstempfname('.metrics')
        tools.picard.MarkDuplicatesTool().execute(
            [aln_bam], aln_bam_dupe_processed,
            dupe_removal_out_metrics, picardOptions=opts,
            JVMmemory=JVMmemory
        )
    else:
        aln_bam_dupe_processed = aln_bam

    samtools.sort(aln_bam_dupe_processed, bam_aligned)
    os.unlink(aln_bam)
    
    if excludeDuplicates:
        os.unlink(aln_bam_dupe_processed)

    samtools.index(bam_aligned)

    # -- call plot function --
    plot_coverage(
        bam_aligned, out_plot_file, plot_format, plot_data_style, plot_style, plot_width, plot_height, plot_dpi, plot_title,
        base_q_threshold, mapping_q_threshold, max_coverage_depth, read_length_threshold, excludeDuplicates, out_summary
    )

    # remove the output bam, unless it is needed
    if out_bam is None:
        os.unlink(bam_aligned)

    # remove the files created by bwa index. 
    # The empty extension causes the original fasta file to be removed
    for ext in [".amb", ".ann", ".bwt", ".bwa", ".pac", ".sa", ""]:
        file_to_remove = ref_indexed + ext
        if os.path.isfile(file_to_remove):
            os.unlink(file_to_remove)


def parser_align_and_plot_coverage(parser=argparse.ArgumentParser()):
    parser = parser_plot_coverage_common(parser)
    parser.add_argument('ref_fasta', default=None, help='Reference genome, FASTA format.')
    parser.add_argument(
        '--outBam',
        dest="out_bam",
        default=None,
        help='Output aligned, indexed BAM file. Default is to write to temp.'
    )
    parser.add_argument(
        '--sensitive', action="store_true",
        help="Equivalent to giving bwa: '-k 12 -A 1 -B 1 -O 1 -E 1'. Only relevant if the bwa aligner is selected (the default). "
    )
    parser.add_argument(
        '--excludeDuplicates', action="store_true",
        help="MarkDuplicates with Picard and only plot non-duplicates"
    )
    parser.add_argument(
        '--JVMmemory',
        default=tools.picard.MarkDuplicatesTool.jvmMemDefault,
        help='JVM virtual memory size (default: %(default)s)'
    )
    parser.add_argument(
        '--picardOptions',
        default=[],
        nargs='*',
        help='Optional arguments to Picard\'s MarkDuplicates, OPTIONNAME=value ...'
    )
    parser.add_argument(
        '-T',
        dest="min_score_to_output",
        default=30,
        type=int,
        help="The min score to output during alignment (default: %(default)s)"
    )
    parser.add_argument('--aligner', choices=['novoalign', 'bwa'], default='bwa', help='aligner (default: %(default)s)')
    parser.add_argument('--aligner_options', default=None, help='aligner options (default for novoalign: "-r Random -l 40 -g 40 -x 20 -t 100 -k", bwa: "-T 30"')
    parser.add_argument(
        '--NOVOALIGN_LICENSE_PATH',
        default=None,
        dest="novoalign_license_path",
        help='A path to the novoalign.lic file. This overrides the NOVOALIGN_LICENSE_PATH environment variable. (default: %(default)s)'
    )

    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, align_and_plot_coverage, split_args=True)
    return parser


__commands__.append(('align_and_plot_coverage', parser_align_and_plot_coverage))


# =======================


def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)


if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)
