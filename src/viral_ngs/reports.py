#!/usr/bin/env python3
''' Functions to create reports from genomics pipeline data.
'''

__author__ = "dpark@broadinstitute.org"
__commands__ = []

import argparse
import logging
import glob
import os
import time
from collections import OrderedDict, defaultdict
import csv
import math
import shutil

import pysam
from pybedtools import BedTool
import Bio.SeqIO
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot

import util.cmd
import util.file
import util.misc
import tools.samtools
import tools.bwa
import tools.fastqc
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
            counts = [(len(s), util.misc.unambig_count(s.seq)) for s in Bio.SeqIO.parse(inf, 'fasta') if len(s) > 0]
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


def genome_coverage_stats_only(mapped_bam, chr_name=None, cov_thresholds=(1, 5, 20, 100)):
    out = {}
    with pysam.AlignmentFile(mapped_bam, 'rb') as bam:
        coverages = list([pcol.nsegments for pcol in bam.pileup(chr_name)])
    if coverages:
        out['aln2self_cov_median'] = median(coverages)
        out['aln2self_cov_mean'] = "%0.3f" % mean(coverages)
        out['aln2self_cov_mean_non0'] = "%0.3f" % mean([n for n in coverages if n > 0])
        for thresh in cov_thresholds:
            out['aln2self_cov_%dX' % thresh] = sum(1 for n in coverages if n >= thresh)
    return out


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
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, assembly_stats, split_args=True)
    return parser


__commands__.append(('assembly_stats', parser_assembly_stats))


def _get_samples_from_bam(bam):
    with pysam.AlignmentFile(bam) as af:
        return set(rg['SM'] for rg in af.header['RG'])
def _get_chrs_from_bam(bam):
    with pysam.AlignmentFile(bam) as af:
        return list(sq['SN'] for sq in af.header['SQ'])

def parser_coverage_only(parser=argparse.ArgumentParser()):
    parser.add_argument('mapped_bams', nargs='+', help='Aligned-to-self mapped bam files.')
    parser.add_argument('out_report', help='Output report file.')
    parser.add_argument('--cov_thresholds',
                        nargs='+',
                        type=int,
                        default=(1, 5, 20, 100),
                        help='Genome coverage thresholds to report on. (default: %(default)s)')
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, coverage_only, split_args=True)
    return parser
def coverage_only(mapped_bams, out_report, cov_thresholds=(1, 5, 20, 100)):
    header = ['sample','aln2self_cov_median', 'aln2self_cov_mean', 'aln2self_cov_mean_non0']
    header += ['aln2self_cov_%dX' % thresh for thresh in cov_thresholds]
    with open(out_report, 'wt') as outf:
        outf.write('\t'.join(header)+'\n')
        for bam in mapped_bams:
            # check for index and auto-create if needed
            with pysam.AlignmentFile(bam) as af:
                is_indexed = af.has_index()
            if not is_indexed:
                pysam.index(bam)
            # get unique sample name
            samples = _get_samples_from_bam(bam)
            if len(samples) != 1:
                raise Exception("input bam file {} has {} unique samples: {} (require one unique sample)".format(bam, len(samples), str(samples)))
            sample_name = samples.pop()
            # get and write coverage stats
            row = genome_coverage_stats_only(bam, cov_thresholds=cov_thresholds)
            row['sample'] = sample_name
            outf.write('\t'.join([str(row.get(h,'')) for h in header])+'\n')
            # for multi-seg genomes, also do per-chr stats
            chrs = _get_chrs_from_bam(bam)
            if len(chrs) > 1:
                for i in range(len(chrs)):
                    row = genome_coverage_stats_only(bam, chr_name=chrs[i], cov_thresholds=cov_thresholds)
                    row['sample'] = "{}-{}".format(sample_name, i+1)
                    outf.write('\t'.join([str(row.get(h,'')) for h in header])+'\n')
__commands__.append(('coverage_only', parser_coverage_only))



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
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
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


def consolidate_spike_count(in_dir, out_file):
    '''Consolidate multiple spike count reports into one.'''
    with open(out_file, 'wt') as outf:
        for fn in os.listdir(in_dir):
            fn = os.path.join(in_dir, fn)
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
    parser.add_argument('in_dir', metavar="inDir", help='Input spike count directory.')
    parser.add_argument('out_file', metavar="outFile", help='Output report file.')
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, consolidate_spike_count, split_args=True)
    return parser


__commands__.append(('consolidate_spike_count', parser_consolidate_spike_count))


def aggregate_spike_count(in_dir, out_file):
    '''aggregate multiple spike count reports into one.'''
    spike_in_sample_counts = defaultdict(dict) # For a given spikein ID, map to sample name and corresponding count
    samples_seen = []
    with open(out_file, 'wt') as outf:
        for fn in glob.glob(os.path.realpath(in_dir)+"/*.spike_count.txt"):# os.listdir():
            #fn = os.path.join(in_dir, fn)
            s = os.path.basename(fn)
            if not s.endswith('.spike_count.txt'):
                raise Exception()
            if s.find('.spike_count.txt'):
                s = s[:-len('.spike_count.txt')]
            if s not in samples_seen:
                samples_seen.append(s)
            with open(fn, 'rt') as inf:
                for line in inf:
                    if not line.startswith('Input bam') and not line.startswith('*'):
                        spike, count = [line.strip().split('\t')[i] for i in [0,2]]
                        spike_in_sample_counts[spike][s] = count
                        #outf.write('\t'.join([s, spike, count]) + '\n')
        outf.write("\t".join(["spike-in"]+samples_seen)+"\n")
        for spike in sorted(spike_in_sample_counts.keys()):
            row = []
            row.append(spike)
            for s in samples_seen:
                if s in spike_in_sample_counts[spike]:
                    row.append(spike_in_sample_counts[spike][s])
                else:
                    row.append("0")
            outf.write("\t".join(row)+"\n")


def parser_aggregate_spike_count(parser=argparse.ArgumentParser()):
    parser.add_argument('in_dir', metavar="inDir", help='Input spike count directory.')
    parser.add_argument('out_file', metavar="outFile", help='Output report file.')
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, aggregate_spike_count, split_args=True)
    return parser
__commands__.append(('aggregate_spike_count', parser_aggregate_spike_count))


def aggregate_alignment_counts(in_reports, out_file):
    '''aggregate multiple reports from read_utils.py bwamem_idxstats into one report.'''
    seq_in_sample_counts = defaultdict(dict) # For a given ref sequence ID, map to input file and corresponding count
    input_files_seen = []
    with open(out_file, 'wt') as outf:
        for in_report in in_reports:
            short_name = os.path.basename(in_report)
            for suffix in ['.txt','.tsv']:
                if short_name.endswith(suffix):
                    short_name = short_name[:-len(suffix)]   
            if short_name not in input_files_seen:
                input_files_seen.append(short_name)
            with open(in_report, 'rt') as inf:
                for line in inf:
                    if not line.startswith('Input bam') and not line.startswith('*'):
                        seq_mapped_to, count = [line.strip().split('\t')[i] for i in [0,2]]
                        seq_in_sample_counts[seq_mapped_to][short_name] = count
        outf.write("\t".join(["seq_mapped_to"]+sorted(input_files_seen))+"\n")
        for seq_mapped_to in sorted(seq_in_sample_counts.keys()):
            row = []
            row.append(seq_mapped_to)
            for s in sorted(input_files_seen):
                if s in seq_in_sample_counts[seq_mapped_to]:
                    row.append(seq_in_sample_counts[seq_mapped_to][s])
                else:
                    row.append("0")
            outf.write("\t".join(row)+"\n")

def parser_aggregate_alignment_counts(parser=argparse.ArgumentParser()):
    parser.add_argument('in_reports', nargs="+", metavar="in_reports", help='tsv reports with alignment counts from read_utils.py bwamem_idxstats')
    parser.add_argument('out_file', metavar="outFile", help='Output report file.')
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, aggregate_alignment_counts, split_args=True)
    return parser
__commands__.append(('aggregate_alignment_counts', parser_aggregate_alignment_counts))


# =========================


def parser_plot_coverage_common(parser=argparse.ArgumentParser()):    # parser needs add_help=False?
    parser.add_argument('in_bam', help='Input reads, BAM format.')
    parser.add_argument('out_plot_file', help='The generated chart file')
    parser.add_argument(
        '--plotFormat',
        dest="plot_format",
        default=None,
        type=str,
        choices=list(matplotlib.pyplot.gcf().canvas.get_supported_filetypes().keys()),
        metavar='',
        help="File format of the coverage plot. By default it is inferred from the file extension of out_plot_file, but it can be set explicitly via --plotFormat. Valid formats include: "
        + ", ".join(list(matplotlib.pyplot.gcf().canvas.get_supported_filetypes().keys()))
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
        choices=matplotlib.pyplot.style.available,
        metavar='',
        help="The plot visual style. Valid options: " + ", ".join(matplotlib.pyplot.style.available) + " (default: %(default)s)"
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
        default=matplotlib.pyplot.gcf().get_dpi(),
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
        '--plotXLimits',
        dest="plot_x_limits",
        nargs=2,
        default=None,
        type=int,
        help="Limits on the x-axis of the coverage plot; args are '<min> <max>'"
    )
    parser.add_argument(
        '--plotYLimits',
        dest="plot_y_limits",
        nargs=2,
        default=None,
        type=int,
        help="Limits on the y-axis of the coverage plot; args are '<min> <max>'"
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
    parser.add_argument('-l',
        dest="read_length_threshold",
        default=None,
        type=int,
        help="Read length threshold"
    )
    parser.add_argument(
        '--binLargePlots',
        dest="bin_large_plots",
        action="store_true",
        help="Plot summary read depth in one-pixel-width bins for large plots."
    )
    parser.add_argument(
        '--binningSummaryStatistic',
        dest="binning_summary_statistic",
        choices=["max", "min"],
        type=str,
        default="max",
        help="Statistic used to summarize each bin (max or min)."
    )
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
    plot_x_limits,
    plot_y_limits,
    base_q_threshold,
    mapping_q_threshold,
    max_coverage_depth,
    read_length_threshold,
    plot_only_non_duplicates=False,
    bin_large_plots=False,
    binning_summary_statistic="max",
    out_summary=None
    ):
    ''' 
        Generate a coverage plot from an aligned bam file
    '''
    samtools = tools.samtools.SamtoolsTool()

    # check if in_bam is aligned, if not raise an error
    num_mapped_reads = samtools.count(in_bam, opts=["-F", "4"])
    if num_mapped_reads == 0:
        raise Exception(
            """The bam file specified appears to have zero mapped reads. 'plot_coverage' requires an aligned bam file. You can try 'align_and_plot_coverage' if the plot input bam file contains reads and you don't mind a simple bwa alignment. \n File: %s"""
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
        samtools.view(["-F", "1024", '-@', '3'], in_bam, bam_dupe_processed)
    else:
        bam_dupe_processed = in_bam

    # only sort if not sorted
    bam_sorted = util.file.mkstempfname('.sorted.bam')
    should_remove_sorted = True
    if not util.file.bam_is_sorted(bam_dupe_processed):
        samtools.sort(bam_dupe_processed, bam_sorted, args=["-O", "bam"])
        if plot_only_non_duplicates:
            os.unlink(bam_dupe_processed)
    else:
        bam_sorted = bam_dupe_processed
        if not plot_only_non_duplicates:
            # in this case we are passing through the original in_bam directly
            should_remove_sorted = False

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

    # only remove the sorted bam if it is not the original input bam
    # which we use directly in some casess 
    if should_remove_sorted:
        os.unlink(bam_sorted)

    # ---- create plot based on coverage_tsv_file ----

    segment_depths = OrderedDict()
    domain_max = 0
    with open(coverage_tsv_file, "r") as tabfile:
        for row in csv.reader(tabfile, delimiter='\t'):
            segment_depths.setdefault(row[0], []).append(float(row[2]))
            domain_max += 1

    with matplotlib.pyplot.style.context(plot_style):
        fig = matplotlib.pyplot.gcf()
        DPI = plot_dpi or fig.get_dpi()
        fig.set_size_inches(float(plot_width) / float(DPI), float(plot_height) / float(DPI))

        font_size = (2.5 * plot_height) / float(DPI)

        ax = matplotlib.pyplot.subplot()    # Defines ax variable by creating an empty plot

        # Set the tick labels font
        for label in (ax.get_xticklabels() + ax.get_yticklabels()):
            label.set_fontsize(font_size)
            
        # Binning
        bin_size = 1
        if bin_large_plots:
            # Bin locations and take summary value (maximum or minimum) in each bin
            binning_action = eval(binning_summary_statistic)
            
            inner_plot_width_inches = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted()).width
            inner_plot_width_px = inner_plot_width_inches * fig.dpi # width of actual plot (sans whitespace and y axis text)
            bins_per_pixel = 1 # increase to make smaller (but less visible) bins
            bin_size = 1 + int(domain_max/(inner_plot_width_px * bins_per_pixel))
            
            binned_segment_depths = OrderedDict()
            for segment_num, (segment_name, position_depths) in enumerate(segment_depths.items()):
                summary_depths_in_bins = [binning_action(position_depths[i:i + bin_size]) for i in range(0, len(position_depths), bin_size)]
                binned_segment_depths[segment_name] = summary_depths_in_bins
            segment_depths = binned_segment_depths
        
        # Plotting
        domain_max = 0
        for segment_num, (segment_name, position_depths) in enumerate(segment_depths.items()):
            prior_domain_max = domain_max
            domain_max += len(position_depths)

            colors = list(matplotlib.pyplot.rcParams['axes.prop_cycle'].by_key()['color'])    # get the colors for this style
            segment_color = colors[segment_num % len(colors)]    # pick a color, offset by the segment index

            x_values = range(prior_domain_max, domain_max)
            x_values = [x * bin_size for x in x_values]

            if plot_data_style == "filled":
                matplotlib.pyplot.fill_between(
                    x_values,
                    position_depths, [0] * len(position_depths),
                    linewidth=0,
                    antialiased=True,
                    color=segment_color
                )
            elif plot_data_style == "line":
                matplotlib.pyplot.plot(
                    x_values,
                    position_depths,
                    antialiased=True,
                    color=segment_color
                )
            elif plot_data_style == "dots":
                matplotlib.pyplot.plot(
                    x_values,
                    position_depths,
                    'ro',
                    antialiased=True,
                    color=segment_color
                )

        matplotlib.pyplot.title(plot_title, fontsize=font_size * 1.2)
        matplotlib.pyplot.xlabel("bp", fontsize=font_size * 1.1)
        
        ylabel = "read depth"
        if(bin_size > 1):
        	ylabel = "read depth ({summary} in {size}-bp bin)".format(size=bin_size, summary=binning_summary_statistic)
        matplotlib.pyplot.ylabel(ylabel, fontsize=font_size * 1.1)

        if plot_x_limits is not None:
            x_min, x_max = plot_x_limits
            matplotlib.pyplot.xlim(x_min, x_max)
        if plot_y_limits is not None:
            y_min, y_max = plot_y_limits
            matplotlib.pyplot.ylim(y_min, y_max)

        # to squash a backend renderer error on OSX related to tight layout
        if matplotlib.pyplot.get_backend().lower() in ['agg', 'macosx']:
            fig.set_tight_layout(True)
        else:
            fig.tight_layout()

        matplotlib.pyplot.savefig(out_plot_file, format=plot_format, dpi=DPI)    #, bbox_inches='tight')
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
    plot_x_limits,
    plot_y_limits,
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
    bin_large_plots=False,
    binning_summary_statistic="max",
    JVMmemory=None,
    picardOptions=None,
    min_score_to_filter=None,
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
            aligner_options = '-1' # hidden option to work around kernel/cpu bug; disables multithreaded file read: https://github.com/lh3/bwa/issues/102

    samtools = tools.samtools.SamtoolsTool()

    ref_indexed = util.file.mkstempfname('.reference.fasta')
    shutil.copyfile(ref_fasta, ref_indexed)

    aln_bam = util.file.mkstempfname('.bam')
    if aligner=="bwa":
        bwa = tools.bwa.Bwa()
        

        bwa.index(ref_indexed)

        bwa_opts = aligner_options.split()
        if sensitive:
            bwa_opts += "-k 12 -A 1 -B 1 -O 1 -E 1".split()

        bwa.align_mem_bam(in_bam, ref_indexed, aln_bam, options=bwa_opts,
                          min_score_to_filter=min_score_to_filter)
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
        plot_x_limits, plot_y_limits, base_q_threshold, mapping_q_threshold, max_coverage_depth, read_length_threshold,
        excludeDuplicates, bin_large_plots, binning_summary_statistic, out_summary
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
        '--minScoreToFilter',
        dest="min_score_to_filter",
        type=int,
        help=("Filter bwa alignments using this value as the minimum allowed "
              "alignment score. Specifically, sum the alignment scores across "
              "all alignments for each query (including reads in a pair, "
              "supplementary and secondary alignments) and then only include, "
              "in the output, queries whose summed alignment score is at least "
              "this value. This is only applied when the aligner is 'bwa'. "
              "The filtering on a summed alignment score is sensible for reads "
              "in a pair and supplementary alignments, but may not be "
              "reasonable if bwa outputs secondary alignments (i.e., if '-a' "
              "is in the aligner options). (default: not set - i.e., do not "
              "filter bwa's output)")
    )
    parser.add_argument('--aligner', choices=['novoalign', 'bwa'], default='bwa', help='aligner (default: %(default)s)')
    parser.add_argument('--aligner_options', default=None, help='aligner options (default for novoalign: "-r Random -l 40 -g 40 -x 20 -t 100 -k", bwa: bwa defaults')
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



def parser_fastqc(parser=argparse.ArgumentParser()):
    parser.add_argument('inBam', help='Input reads, BAM format.')
    parser.add_argument('out_html',  help='Output report, HTML format.')
    parser.add_argument('--out_zip', help='Output data, zip archive.')
    parser.add_argument('--threads', type=int, help='Number of threads.')
    util.cmd.attach_main(parser, tools.fastqc.FastQC().execute, split_args=True)
    return parser
__commands__.append(('fastqc', parser_fastqc))

# =======================


def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)


if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)
