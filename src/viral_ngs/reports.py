#!/usr/bin/env python
''' Reports
'''

__author__ = "dpark@broadinstitute.org"
__commands__ = []

import argparse
import logging
import glob
import os
import time
import pysam
import Bio.SeqIO
import util.cmd
import util.file
import util.misc
import tools.samtools
import assembly
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
    out['reads_raw'] = sum(samtools.count(bam) for bam in glob.glob(os.path.join(raw_reads_dir, sample + "*.bam")))

    # pre-assembly stats
    out['assembled_trinity'] = os.path.isfile(os.path.join(assembly_tmp, sample +
                                                           '.assembly1-trinity.fasta')) and 1 or 0
    sub_bam = os.path.join(assembly_tmp, sample + '.subsamp.bam')
    if os.path.isfile(sub_bam):
        out['trinity_in_reads'] = samtools.count(sub_bam)

    # assembly stats
    assembly_fname = os.path.join(assembly_dir, sample + '.fasta')
    if not os.path.isfile(assembly_fname):
        assembly_fname = os.path.join(assembly_tmp, sample + '.assembly2-vfat.fasta')
        if not os.path.isfile(assembly_fname):
            out['n_contigs'] = 0
            return (header, out)
    with open(assembly_fname, 'rt') as inf:
        counts = [(len(s), assembly.unambig_count(s.seq)) for s in Bio.SeqIO.parse(inf, 'fasta') if len(s) > 0]
    out['n_contigs'] = len(counts)
    out['contig_len'] = ','.join(str(x) for x, y in counts)
    out['unambig_bases'] = ','.join(str(y) for x, y in counts)
    out['pct_unambig'] = ','.join(str(float(y) / x) for x, y in counts)

    # read counts from align-to-self
    bam_fname = os.path.join(align_dir, sample + '.bam')
    if not os.path.isfile(bam_fname):
        return (header, out)
    out['aln2self_reads_tot'] = samtools.count(bam_fname)
    out['aln2self_reads_aln'] = samtools.count(bam_fname, opts=['-F', '4'])
    out['aln2self_reads_rmdup'] = samtools.count(bam_fname, opts=['-F', '1028'])
    if out['aln2self_reads_aln']:
        out['aln2self_pct_nondup'] = float(out['aln2self_reads_rmdup']) / out['aln2self_reads_aln']

    # genome coverage stats
    bam_fname = os.path.join(align_dir, sample + '.mapped.bam')
    with pysam.AlignmentFile(bam_fname, 'rb') as bam:
        coverages = list([pcol.nsegments for pcol in bam.pileup()])
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


def get_refalign_stats(sample):
    pass


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
                        raise
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
                    if not line.startswith('Input bam'):
                        spike, count = line.strip().split('\t')
                        outf.write('\t'.join([s, spike, count]) + '\n')


def parser_consolidate_spike_count(parser=argparse.ArgumentParser()):
    parser.add_argument('inDir', help='Input spike count directory.')
    parser.add_argument('outFile', help='Output report file.')
    util.cmd.attach_main(parser, consolidate_spike_count, split_args=True)
    return parser


__commands__.append(('consolidate_spike_count', parser_consolidate_spike_count))

# =======================


def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)


if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)
