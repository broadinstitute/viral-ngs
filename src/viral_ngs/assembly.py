#!/usr/bin/env python3
''' This script contains a number of utilities for viral sequence assembly
    from NGS reads.  Primarily used for Lassa and Ebola virus analysis in
    the Sabeti Lab / Broad Institute Viral Genomics.
'''

__author__ = "dpark@broadinstitute.org, rsealfon@broadinstitute.org"
__commands__ = []

# built-ins
import argparse
import logging
import random
import os
import os.path
import shutil
import subprocess
import functools
import operator
import concurrent.futures
import csv
import itertools
from collections import OrderedDict, defaultdict

# from viral-core
import util.cmd
import util.file
import util.misc
import read_utils
import tools
import tools.picard
import tools.samtools
import tools.gatk
import tools.novoalign
import tools.trimmomatic

# intra-module
import assemble.vcf
import assemble.spades
import assemble.mafft
import assemble.mummer
import assemble.muscle
import assemble.gap2seq

# third-party
import numpy
import Bio.AlignIO
import Bio.SeqIO
import Bio.Data.IUPACData

log = logging.getLogger(__name__)


class DenovoAssemblyError(RuntimeError):

    '''Indicates a failure of the de novo assembly step.  Can indicate an internal error in the assembler, but also
    insufficient input data (for assemblers which cannot report that condition separately).'''

    def __init__(self, reason):
        super(DenovoAssemblyError, self).__init__(reason)


def trim_rmdup_subsamp_reads(inBam, clipDb, outBam, n_reads=100000, trim_opts=None):
    ''' Take reads through Trimmomatic, Prinseq, and subsampling.
        This should probably move over to read_utils.
    '''

    downsamplesam = tools.picard.DownsampleSamTool()
    samtools = tools.samtools.SamtoolsTool()

    if n_reads < 1:
        raise Exception()

    # BAM -> fastq
    infq = list(map(util.file.mkstempfname, ['.in.1.fastq', '.in.2.fastq']))
    tools.picard.SamToFastqTool().execute(inBam, infq[0], infq[1])
    n_input = util.file.count_fastq_reads(infq[0])

    # --- Trimmomatic ---
    trimfq = list(map(util.file.mkstempfname, ['.trim.1.fastq', '.trim.2.fastq']))
    trimfq_unpaired = list(map(util.file.mkstempfname, ['.trim.unpaired.1.fastq', '.trim.unpaired.2.fastq']))
    if n_input == 0:
        for i in range(2):
            shutil.copyfile(infq[i], trimfq[i])
    else:
        tools.trimmomatic.TrimmomaticTool().execute(
            infq[0],
            infq[1],
            trimfq[0],
            trimfq[1],
            clipDb,
            unpairedOutFastq1=trimfq_unpaired[0],
            unpairedOutFastq2=trimfq_unpaired[1],
            **(trim_opts or {})
        )

    n_trim = max(map(util.file.count_fastq_reads, trimfq))    # count is pairs
    n_trim_unpaired = sum(map(util.file.count_fastq_reads, trimfq_unpaired))    # count is individual reads

    # --- Prinseq duplicate removal ---
    # the paired reads from trim, de-duplicated
    rmdupfq = list(map(util.file.mkstempfname, ['.rmdup.1.fastq', '.rmdup.2.fastq']))

    # the unpaired reads from rmdup, de-duplicated with the other singleton files later on
    rmdupfq_unpaired_from_paired_rmdup = list(
        map(util.file.mkstempfname, ['.rmdup.unpaired.1.fastq', '.rmdup.unpaired.2.fastq'])
    )

    prinseq = tools.prinseq.PrinseqTool()
    prinseq.rmdup_fastq_paired(
        trimfq[0],
        trimfq[1],
        rmdupfq[0],
        rmdupfq[1],
        includeUnmated=False,
        unpairedOutFastq1=rmdupfq_unpaired_from_paired_rmdup[0],
        unpairedOutFastq2=rmdupfq_unpaired_from_paired_rmdup[1]
    )

    n_rmdup_paired = max(map(util.file.count_fastq_reads, rmdupfq))    # count is pairs
    n_rmdup = n_rmdup_paired    # count is pairs

    tmp_header = util.file.mkstempfname('.header.sam')
    tools.samtools.SamtoolsTool().dumpHeader(inBam, tmp_header)

    # convert paired reads to bam
    # stub out an empty file if the input fastqs are empty
    tmp_bam_paired = util.file.mkstempfname('.paired.bam')
    if all(os.path.getsize(x) > 0 for x in rmdupfq):
        tools.picard.FastqToSamTool().execute(rmdupfq[0], rmdupfq[1], 'Dummy', tmp_bam_paired)
    else:
        opts = ['INPUT=' + tmp_header, 'OUTPUT=' + tmp_bam_paired, 'VERBOSITY=ERROR']
        tools.picard.PicardTools().execute('SamFormatConverter', opts, JVMmemory='50m')

    n_paired_subsamp = 0
    n_unpaired_subsamp = 0
    n_rmdup_unpaired = 0

    # --- subsampling ---

    # if we have too few paired reads after trimming and de-duplication, we can incorporate unpaired reads to reach the desired count
    if n_rmdup_paired * 2 < n_reads:
        # the unpaired reads from the trim operation, and the singletons from Prinseq
        unpaired_concat = util.file.mkstempfname('.unpaired.fastq')

        # merge unpaired reads from trimmomatic and singletons left over from paired-mode prinseq
        util.file.cat(unpaired_concat, trimfq_unpaired + rmdupfq_unpaired_from_paired_rmdup)

        # remove the earlier singleton files
        for f in trimfq_unpaired + rmdupfq_unpaired_from_paired_rmdup:
            os.unlink(f)

        # the de-duplicated singletons
        unpaired_concat_rmdup = util.file.mkstempfname('.unpaired.rumdup.fastq')
        prinseq.rmdup_fastq_single(unpaired_concat, unpaired_concat_rmdup)
        os.unlink(unpaired_concat)

        n_rmdup_unpaired = util.file.count_fastq_reads(unpaired_concat_rmdup)

        did_include_subsampled_unpaired_reads = True
        # if there are no unpaired reads, simply output the paired reads
        if n_rmdup_unpaired == 0:
            shutil.copyfile(tmp_bam_paired, outBam)
            n_output = samtools.count(outBam)
        else:
            # take pooled unpaired reads and convert to bam
            tmp_bam_unpaired = util.file.mkstempfname('.unpaired.bam')
            tools.picard.FastqToSamTool().execute(unpaired_concat_rmdup, None, 'Dummy', tmp_bam_unpaired)

            tmp_bam_unpaired_subsamp = util.file.mkstempfname('.unpaired.subsamp.bam')
            reads_to_add = (n_reads - (n_rmdup_paired * 2))

            downsamplesam.downsample_to_approx_count(tmp_bam_unpaired, tmp_bam_unpaired_subsamp, reads_to_add)
            n_unpaired_subsamp = samtools.count(tmp_bam_unpaired_subsamp)
            os.unlink(tmp_bam_unpaired)

            # merge the subsampled unpaired reads into the bam to be used as
            tmp_bam_merged = util.file.mkstempfname('.merged.bam')
            tools.picard.MergeSamFilesTool().execute([tmp_bam_paired, tmp_bam_unpaired_subsamp], tmp_bam_merged)
            os.unlink(tmp_bam_unpaired_subsamp)

            tools.samtools.SamtoolsTool().reheader(tmp_bam_merged, tmp_header, outBam)
            os.unlink(tmp_bam_merged)

            n_paired_subsamp = n_rmdup_paired    # count is pairs
            n_output = n_rmdup_paired * 2 + n_unpaired_subsamp    # count is individual reads

    else:
        did_include_subsampled_unpaired_reads = False
        log.info("PRE-SUBSAMPLE COUNT: %s read pairs", n_rmdup_paired)

        tmp_bam_paired_subsamp = util.file.mkstempfname('.unpaired.subsamp.bam')
        downsamplesam.downsample_to_approx_count(tmp_bam_paired, tmp_bam_paired_subsamp, n_reads)
        n_paired_subsamp = samtools.count(tmp_bam_paired_subsamp) // 2    # count is pairs
        n_output = n_paired_subsamp * 2    # count is individual reads

        tools.samtools.SamtoolsTool().reheader(tmp_bam_paired_subsamp, tmp_header, outBam)
        os.unlink(tmp_bam_paired_subsamp)

    os.unlink(tmp_bam_paired)
    os.unlink(tmp_header)

    n_final_individual_reads = samtools.count(outBam)

    log.info("Pre-DeNovoAssembly read filters: ")
    log.info("    {} read pairs at start ".format(n_input))
    log.info(
        "    {} read pairs after Trimmomatic {}".format(
            n_trim, "(and {} unpaired)".format(n_trim_unpaired) if n_trim_unpaired > 0 else ""
        )
    )
    log.info(
        "    {} read pairs after Prinseq rmdup {} ".format(
            n_rmdup, "(and {} unpaired from Trimmomatic+Prinseq)".format(n_rmdup_unpaired)
            if n_rmdup_unpaired > 0 else ""
        )
    )
    if did_include_subsampled_unpaired_reads:
        log.info(
            "   Too few individual reads ({}*2={}) from paired reads to reach desired threshold ({}), so including subsampled unpaired reads".format(
                n_rmdup_paired, n_rmdup_paired * 2, n_reads
            )
        )
    else:
        log.info("  Paired read count sufficient to reach threshold ({})".format(n_reads))
    log.info(
        "    {} individual reads for de novo assembly ({}{})".format(
            n_output, "paired subsampled {} -> {}".format(n_rmdup_paired, n_paired_subsamp)
            if not did_include_subsampled_unpaired_reads else "{} read pairs".format(n_rmdup_paired),
            " + unpaired subsampled {} -> {}".format(n_rmdup_unpaired, n_unpaired_subsamp)
            if did_include_subsampled_unpaired_reads else ""
        )
    )
    log.info("    {} individual reads".format(n_final_individual_reads))

    if did_include_subsampled_unpaired_reads:
        if n_final_individual_reads < n_reads:
            log.warning(
                "NOTE: Even with unpaired reads included, there are fewer unique trimmed reads than requested for de novo assembly input."
            )

    # clean up temp files
    for i in range(2):
        for f in infq, trimfq, rmdupfq:
            if os.path.exists(f[i]):
                os.unlink(f[i])

    # multiply counts so all reflect individual reads
    return (n_input * 2, n_trim * 2, n_rmdup * 2, n_output, n_paired_subsamp * 2, n_unpaired_subsamp)


def parser_trim_rmdup_subsamp(parser=argparse.ArgumentParser()):
    parser.add_argument('inBam', help='Input reads, unaligned BAM format.')
    parser.add_argument('clipDb', help='Trimmomatic clip DB.')
    parser.add_argument(
        'outBam',
        help="""Output reads, unaligned BAM format (currently, read groups and other
                header information are destroyed in this process)."""
    )
    parser.add_argument(
        '--n_reads',
        default=100000,
        type=int,
        help='Subsample reads to no more than this many individual reads. Note that paired reads are given priority, and unpaired reads are included to reach the count if there are too few paired reads to reach n_reads. (default %(default)s)'
    )
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, trim_rmdup_subsamp_reads, split_args=True)
    return parser


__commands__.append(('trim_rmdup_subsamp', parser_trim_rmdup_subsamp))


def assemble_spades(
    in_bam,
    clip_db,
    out_fasta,
    spades_opts='',
    contigs_trusted=None, contigs_untrusted=None,
    filter_contigs=False,
    min_contig_len=0,
    kmer_sizes=(55,65),
    n_reads=10000000,
    outReads=None,
    always_succeed=False,
    max_kmer_sizes=1,
    mem_limit_gb=8,
    threads=None,
):
    '''De novo RNA-seq assembly with the SPAdes assembler.
    '''

    if outReads:
        trim_rmdup_bam = outReads
    else:
        trim_rmdup_bam = util.file.mkstempfname('.subsamp.bam')

    trim_rmdup_subsamp_reads(in_bam, clip_db, trim_rmdup_bam, n_reads=n_reads,
                             trim_opts=dict(maxinfo_target_length=35, maxinfo_strictness=.2))

    with tools.picard.SamToFastqTool().execute_tmp(trim_rmdup_bam, includeUnpaired=True, illuminaClipping=True
                                                   ) as (reads_fwd, reads_bwd, reads_unpaired):
        try:
            assemble.spades.SpadesTool().assemble(reads_fwd=reads_fwd, reads_bwd=reads_bwd, reads_unpaired=reads_unpaired,
                                               contigs_untrusted=contigs_untrusted, contigs_trusted=contigs_trusted,
                                               contigs_out=out_fasta, filter_contigs=filter_contigs,
                                               min_contig_len=min_contig_len,
                                               kmer_sizes=kmer_sizes, always_succeed=always_succeed, max_kmer_sizes=max_kmer_sizes,
                                               spades_opts=spades_opts, mem_limit_gb=mem_limit_gb,
                                               threads=threads)
        except subprocess.CalledProcessError as e:
            raise DenovoAssemblyError('SPAdes assembler failed: ' + str(e))

    if not outReads:
        os.unlink(trim_rmdup_bam)



def parser_assemble_spades(parser=argparse.ArgumentParser()):
    parser.add_argument('in_bam', help='Input unaligned reads, BAM format. May include both paired and unpaired reads.')
    parser.add_argument('clip_db', help='Trimmomatic clip db')
    parser.add_argument('out_fasta', help='Output assembled contigs. Note that, since this is RNA-seq assembly, for each assembled genomic region there may be several contigs representing different variants of that region.')
    parser.add_argument('--contigsTrusted', dest='contigs_trusted', 
                        help='Optional input contigs of high quality, previously assembled from the same sample')
    parser.add_argument('--contigsUntrusted', dest='contigs_untrusted', 
                        help='Optional input contigs of high medium quality, previously assembled from the same sample')
    parser.add_argument('--nReads', dest='n_reads', type=int, default=10000000, 
                        help='Before assembly, subsample the reads to at most this many')
    parser.add_argument('--outReads', default=None, help='Save the trimmomatic/prinseq/subsamp reads to a BAM file')
    parser.add_argument('--filterContigs', dest='filter_contigs', default=False, action='store_true', 
                        help='only output contigs SPAdes is sure of (drop lesser-quality contigs from output)')
    parser.add_argument('--alwaysSucceed', dest='always_succeed', default=False, action='store_true',
                        help='if assembly fails for any reason, output an empty contigs file, rather than failing with '
                        'an error code')
    parser.add_argument('--minContigLen', dest='min_contig_len', type=int, default=0,
                        help='only output contigs longer than this many bp')
    parser.add_argument('--spadesOpts', dest='spades_opts', default='', help='(advanced) Extra flags to pass to the SPAdes assembler')
    parser.add_argument('--memLimitGb', dest='mem_limit_gb', default=4, type=int, help='Max memory to use, in GB (default: %(default)s)')
    util.cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, assemble_spades, split_args=True)
    return parser


__commands__.append(('assemble_spades', parser_assemble_spades))

def gapfill_gap2seq(in_scaffold, in_bam, out_scaffold, threads=None, mem_limit_gb=4, time_soft_limit_minutes=60.0,
                    random_seed=0, gap2seq_opts='', always_succeed=False):
    ''' This step runs the Gap2Seq tool to close gaps between contigs in a scaffold.
    '''
    try:
        assemble.gap2seq.Gap2SeqTool().gapfill(in_scaffold, in_bam, out_scaffold, gap2seq_opts=gap2seq_opts, threads=threads,
                                            mem_limit_gb=mem_limit_gb, time_soft_limit_minutes=time_soft_limit_minutes, 
                                            random_seed=random_seed)
    except Exception as e:
        if not always_succeed:
            raise
        log.warning('Gap-filling failed (%s); ignoring error, emulating successful run where simply no gaps were filled.')
        shutil.copyfile(in_scaffold, out_scaffold)

def parser_gapfill_gap2seq(parser=argparse.ArgumentParser(description='Close gaps between contigs in a scaffold')):
    parser.add_argument('in_scaffold', help='FASTA file containing the scaffold.  Each FASTA record corresponds to one '
                        'segment (for multi-segment genomes).  Contigs within each segment are separated by Ns.')
    parser.add_argument('in_bam', help='Input unaligned reads, BAM format.')
    parser.add_argument('out_scaffold', help='Output assemble.')
    parser.add_argument('--memLimitGb', dest='mem_limit_gb', default=4.0, help='Max memory to use, in gigabytes %(default)s')
    parser.add_argument('--timeSoftLimitMinutes', dest='time_soft_limit_minutes', default=60.0,
                        help='Stop trying to close more gaps after this many minutes (default: %(default)s); this is a soft/advisory limit')
    parser.add_argument('--maskErrors', dest='always_succeed', default=False, action='store_true',
                        help='In case of any error, just copy in_scaffold to out_scaffold, emulating a successful run that simply could not '
                        'fill any gaps.')
    parser.add_argument('--gap2seqOpts', dest='gap2seq_opts', default='', help='(advanced) Extra command-line options to pass to Gap2Seq')
    parser.add_argument('--randomSeed', dest='random_seed', default=0, help='Random seed; 0 means use current time')
    util.cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, gapfill_gap2seq, split_args=True)
    return parser


__commands__.append(('gapfill_gap2seq', parser_gapfill_gap2seq))


def _order_and_orient_orig(inFasta, inReference, outFasta,
        outAlternateContigs=None,
        breaklen=None, # aligner='nucmer', circular=False, trimmed_contigs=None,
        maxgap=200, minmatch=10, mincluster=None,
        min_pct_id=0.6, min_contig_len=200, min_pct_contig_aligned=0.3):
    ''' This step cleans up the de novo assembly with a known reference genome.
        Uses MUMmer (nucmer or promer) to create a reference-based consensus
        sequence of aligned contigs (with runs of N's in between the de novo
        contigs).
    '''
    mummer = assemble.mummer.MummerTool()
    #if trimmed_contigs:
    #    trimmed = trimmed_contigs
    #else:
    #    trimmed = util.file.mkstempfname('.trimmed.contigs.fasta')
    #mummer.trim_contigs(inReference, inFasta, trimmed,
    #        aligner=aligner, circular=circular, extend=False, breaklen=breaklen,
    #        min_pct_id=min_pct_id, min_contig_len=min_contig_len,
    #        min_pct_contig_aligned=min_pct_contig_aligned)
    #mummer.scaffold_contigs(inReference, trimmed, outFasta,
    #        aligner=aligner, circular=circular, extend=True, breaklen=breaklen,
    #        min_pct_id=min_pct_id, min_contig_len=min_contig_len,
    #        min_pct_contig_aligned=min_pct_contig_aligned)
    mummer.scaffold_contigs_custom(
        inReference,
        inFasta,
        outFasta,
        outAlternateContigs=outAlternateContigs,
        extend=True,
        breaklen=breaklen,
        min_pct_id=min_pct_id,
        min_contig_len=min_contig_len,
        maxgap=maxgap,
        minmatch=minmatch,
        mincluster=mincluster,
        min_pct_contig_aligned=min_pct_contig_aligned
    )
    #if not trimmed_contigs:
    #    os.unlink(trimmed)
    return 0

def _call_order_and_orient_orig(inReference, outFasta, outAlternateContigs, **kwargs):
    return _order_and_orient_orig(inReference=inReference, outFasta=outFasta, outAlternateContigs=outAlternateContigs, **kwargs)

def order_and_orient(inFasta, inReference, outFasta,
        outAlternateContigs=None, outReference=None,
        breaklen=None, # aligner='nucmer', circular=False, trimmed_contigs=None,
        maxgap=200, minmatch=10, mincluster=None,
        min_pct_id=0.6, min_contig_len=200, min_pct_contig_aligned=0.3, n_genome_segments=0, 
        outStats=None, threads=None):
    ''' This step cleans up the de novo assembly with a known reference genome.
        Uses MUMmer (nucmer or promer) to create a reference-based consensus
        sequence of aligned contigs (with runs of N's in between the de novo
        contigs).
    '''

    chk = util.cmd.check_input

    ref_segments_all = [tuple(Bio.SeqIO.parse(inRef, 'fasta')) for inRef in util.misc.make_seq(inReference)]
    chk(ref_segments_all, 'No references given')
    chk(len(set(map(len, ref_segments_all))) == 1, 'All references must have the same number of segments')
    if n_genome_segments:
        if len(ref_segments_all) == 1:
            chk((len(ref_segments_all[0]) % n_genome_segments) == 0,
                'Number of reference sequences in must be a multiple of the number of genome segments')
        else:
            chk(len(ref_segments_all[0]) == n_genome_segments,
                'Number of reference segments must match the n_genome_segments parameter')
    else:
        n_genome_segments = len(ref_segments_all[0])
    ref_segments_all = functools.reduce(operator.concat, ref_segments_all, ())

    n_refs = len(ref_segments_all) // n_genome_segments
    log.info('n_genome_segments={} n_refs={}'.format(n_genome_segments, n_refs))
    ref_ids = []

    with util.file.tempfnames(suffixes=[ '.{}.ref.fasta'.format(ref_num) for ref_num in range(n_refs)]) as refs_fasta, \
         util.file.tempfnames(suffixes=[ '.{}.scaffold.fasta'.format(ref_num) for ref_num in range(n_refs)]) as scaffolds_fasta, \
         util.file.tempfnames(suffixes=[ '.{}.altcontig.fasta'.format(ref_num) for ref_num in range(n_refs)]) as alt_contigs_fasta:

        for ref_num in range(n_refs):
            this_ref_segs = ref_segments_all[ref_num*n_genome_segments : (ref_num+1)*n_genome_segments]
            ref_ids.append(this_ref_segs[0].id)
            Bio.SeqIO.write(this_ref_segs, refs_fasta[ref_num], 'fasta')

        with concurrent.futures.ProcessPoolExecutor(max_workers=util.misc.sanitize_thread_count(threads)) as executor:
            retvals = executor.map(functools.partial(_call_order_and_orient_orig, inFasta=inFasta,
                breaklen=breaklen, maxgap=maxgap, minmatch=minmatch, mincluster=mincluster, min_pct_id=min_pct_id,
                min_contig_len=min_contig_len, min_pct_contig_aligned=min_pct_contig_aligned),
                refs_fasta, scaffolds_fasta, alt_contigs_fasta)
        for r in retvals:
            # if an exception is raised by _call_order_and_orient_contig, the
            # concurrent.futures.Executor.map function documentations states
            # that the same exception will be raised when retrieving that entry
            # of the retval iterator. This for loop is intended to reveal any
            # CalledProcessErrors from mummer itself.
            assert r==0

        scaffolds = [tuple(Bio.SeqIO.parse(scaffolds_fasta[ref_num], 'fasta')) for ref_num in range(n_refs)]
        base_counts = [sum([len(seg.seq.ungap('N')) for seg in scaffold]) \
            if len(scaffold)==n_genome_segments else 0 for scaffold in scaffolds]
        best_ref_num = numpy.argmax(base_counts)
        if len(scaffolds[best_ref_num]) != n_genome_segments:
            raise IncompleteAssemblyError(len(scaffolds[best_ref_num]), n_genome_segments)
        log.info('base_counts={} best_ref_num={}'.format(base_counts, best_ref_num))
        shutil.copyfile(scaffolds_fasta[best_ref_num], outFasta)
        if outAlternateContigs:
            shutil.copyfile(alt_contigs_fasta[best_ref_num], outAlternateContigs)
        if outReference:
            shutil.copyfile(refs_fasta[best_ref_num], outReference)
        if outStats:
            ref_ranks = (-numpy.array(base_counts)).argsort().argsort()
            with open(outStats, 'w') as stats_f:
                stats_w = csv.DictWriter(stats_f, fieldnames='ref_num ref_name base_count rank'.split(), delimiter='\t')
                stats_w.writeheader()
                for ref_num, (ref_id, base_count, rank) in enumerate(zip(ref_ids, base_counts, ref_ranks)):
                    stats_w.writerow({'ref_num': ref_num, 'ref_name': ref_id, 'base_count': base_count, 'rank': rank})


def parser_order_and_orient(parser=argparse.ArgumentParser()):
    parser.add_argument('inFasta', help='Input de novo assembly/contigs, FASTA format.')
    parser.add_argument(
        'inReference', nargs='+',
        help=('Reference genome for ordering, orienting, and merging contigs, FASTA format.  Multiple filenames may be listed, each '
              'containing one reference genome. Alternatively, multiple references may be given by specifying a single filename, '
              'and giving the number of reference segments with the nGenomeSegments parameter.  If multiple references are given, '
              'they must all contain the same number of segments listed in the same order.')
    )
    parser.add_argument(
        'outFasta',
        help="""Output assembly, FASTA format, with the same number of
                chromosomes as inReference, and in the same order."""
    )
    parser.add_argument(
        '--outAlternateContigs',
        help="""Output sequences (FASTA format) from alternative contigs that mapped,
                but were not chosen for the final output.""",
        default=None
    )

    parser.add_argument('--nGenomeSegments', dest='n_genome_segments', type=int, default=0,
                        help="""Number of genome segments.  If 0 (the default), the `inReference` parameter is treated as one genome.
                        If positive, the `inReference` parameter is treated as a list of genomes of nGenomeSegments each.""")

    parser.add_argument('--outReference', help='Output the reference chosen for scaffolding to this file')
    parser.add_argument('--outStats', help='Output stats used in reference selection')
    #parser.add_argument('--aligner',
    #                    help='nucmer (nucleotide) or promer (six-frame translations) [default: %(default)s]',
    #                    choices=['nucmer', 'promer'],
    #                    default='nucmer')
    #parser.add_argument("--circular",
    #                    help="""Allow contigs to wrap around the ends of the chromosome.""",
    #                    default=False,
    #                    action="store_true",
    #                    dest="circular")
    parser.add_argument(
        "--breaklen",
        "-b",
        help="""Amount to extend alignment clusters by (if --extend).
                        nucmer default 200, promer default 60.""",
        type=int,
        default=None,
        dest="breaklen"
    )
    parser.add_argument(
        "--maxgap",
        "-g",
        help="""Maximum gap between two adjacent matches in a cluster.
                        Our default is %(default)s.
                        nucmer default 90, promer default 30. Manual suggests going to 1000.""",
        type=int,
        default=200,
        dest="maxgap"
    )
    parser.add_argument(
        "--minmatch",
        "-l",
        help="""Minimum length of an maximal exact match.
                        Our default is %(default)s.
                        nucmer default 20, promer default 6.""",
        type=int,
        default=10,
        dest="minmatch"
    )
    parser.add_argument(
        "--mincluster",
        "-c",
        help="""Minimum cluster length.
                        nucmer default 65, promer default 20.""",
        type=int,
        default=None,
        dest="mincluster"
    )
    parser.add_argument(
        "--min_pct_id",
        "-i",
        type=float,
        default=0.6,
        help="show-tiling: minimum percent identity for contig alignment (0.0 - 1.0, default: %(default)s)"
    )
    parser.add_argument(
        "--min_contig_len",
        type=int,
        default=200,
        help="show-tiling: reject contigs smaller than this (default: %(default)s)"
    )
    parser.add_argument(
        "--min_pct_contig_aligned",
        "-v",
        type=float,
        default=0.3,
        help="show-tiling: minimum percent of contig length in alignment (0.0 - 1.0, default: %(default)s)"
    )
    #parser.add_argument("--trimmed_contigs",
    #                    default=None,
    #                    help="optional output file for trimmed contigs")
    util.cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, order_and_orient, split_args=True)
    return parser


__commands__.append(('order_and_orient', parser_order_and_orient))


class IncompleteAssemblyError(Exception):
    def __init__(self, actual_n, expected_n):
        super(IncompleteAssemblyError, self).__init__(
            'All computed scaffolds are incomplete. Best assembly has {} contigs, expected {}'.format(
                actual_n, expected_n)
        )

class PoorAssemblyError(Exception):
    def __init__(self, chr_idx, seq_len, non_n_count, min_length, segment_length):
        super(PoorAssemblyError, self).__init__(
            'Error: poor assembly quality, chr {}: contig length {}, unambiguous bases {}; bases required of reference segment length: {}/{} ({:.0%})'.format(
                chr_idx, seq_len, non_n_count, min_length, int(segment_length), float(min_length)/float(segment_length)
            )
        )


def impute_from_reference(
    inFasta,
    inReference,
    outFasta,
    minLengthFraction,
    minUnambig,
    replaceLength,
    newName=None,
    aligner='muscle',
    index=False
):
    '''
        This takes a de novo assembly, aligns against a reference genome, and
        imputes all missing positions (plus some of the chromosome ends)
        with the reference genome. This provides an assembly with the proper
        structure (but potentially wrong sequences in areas) from which
        we can perform further read-based refinement.
        Two steps:
        filter_short_seqs: We then toss out all assemblies that come out to
            < 15kb or < 95% unambiguous and fail otherwise.
        modify_contig: Finally, we trim off anything at the end that exceeds
            the length of the known reference assemble.  We also replace all
            Ns and everything within 55bp of the chromosome ends with the
            reference sequence.  This is clearly incorrect consensus sequence,
            but it allows downstream steps to map reads in parts of the genome
            that would otherwise be Ns, and we will correct all of the inferred
            positions with two steps of read-based refinement (below), and
            revert positions back to Ns where read support is lacking.
        FASTA indexing: output assembly is indexed for Picard, Samtools, Novoalign.
    '''
    tempFastas = []

    pmc = parser_modify_contig(argparse.ArgumentParser())
    assert aligner in ('muscle', 'mafft', 'mummer')

    with open(inFasta, 'r') as asmFastaFile:
        with open(inReference, 'r') as refFastaFile:
            asmFasta = Bio.SeqIO.parse(asmFastaFile, 'fasta')
            refFasta = Bio.SeqIO.parse(refFastaFile, 'fasta')
            for idx, (refSeqObj, asmSeqObj) in enumerate(itertools.zip_longest(refFasta, asmFasta)):
                # our zip fails if one file has more seqs than the other
                if not refSeqObj or not asmSeqObj:
                    raise KeyError("inFasta and inReference do not have the same number of sequences. This could be because the de novo assembly process was unable to create contigs for all segments.")

                # error if PoorAssembly
                minLength = len(refSeqObj) * minLengthFraction
                non_n_count = util.misc.unambig_count(asmSeqObj.seq)
                seq_len = len(asmSeqObj)
                log.info(
                    "Assembly Quality - segment {idx} - name {segname} - contig len {len_actual} / {len_desired} ({min_frac}) - unambiguous bases {unamb_actual} / {unamb_desired} ({min_unamb})".format(
                        idx=idx + 1,
                        segname=refSeqObj.id,
                        len_actual=seq_len,
                        len_desired=len(refSeqObj),
                        min_frac=minLengthFraction,
                        unamb_actual=non_n_count,
                        unamb_desired=seq_len,
                        min_unamb=minUnambig
                    )
                )
                if seq_len < minLength or non_n_count < seq_len * minUnambig:
                    raise PoorAssemblyError(idx + 1, seq_len, non_n_count, minLength, len(refSeqObj))

                # prepare temp input and output files
                tmpOutputFile = util.file.mkstempfname(prefix='seq-out-{idx}-'.format(idx=idx), suffix=".fasta")
                concat_file = util.file.mkstempfname('.ref_and_actual.fasta')
                ref_file = util.file.mkstempfname('.ref.fasta')
                actual_file = util.file.mkstempfname('.actual.fasta')
                aligned_file = util.file.mkstempfname('.'+aligner+'.fasta')
                refName = refSeqObj.id
                with open(concat_file, 'wt') as outf:
                    Bio.SeqIO.write([refSeqObj, asmSeqObj], outf, "fasta")
                with open(ref_file, 'wt') as outf:
                    Bio.SeqIO.write([refSeqObj], outf, "fasta")
                with open(actual_file, 'wt') as outf:
                    Bio.SeqIO.write([asmSeqObj], outf, "fasta")

                # align scaffolded genome to reference (choose one of three aligners)
                if aligner == 'mafft':
                    assemble.mafft.MafftTool().execute(
                        [ref_file, actual_file], aligned_file, False, True, True, False, False, None
                    )
                elif aligner == 'muscle':
                    if len(refSeqObj) > 40000:
                        assemble.muscle.MuscleTool().execute(
                            concat_file, aligned_file, quiet=False,
                            maxiters=2, diags=True
                        )
                    else:
                        assemble.muscle.MuscleTool().execute(concat_file, aligned_file, quiet=False)
                elif aligner == 'mummer':
                    assemble.mummer.MummerTool().align_one_to_one(ref_file, actual_file, aligned_file)

                # run modify_contig
                args = [
                    aligned_file, tmpOutputFile, refName, '--call-reference-ns', '--trim-ends', '--replace-5ends',
                    '--replace-3ends', '--replace-length', str(replaceLength), '--replace-end-gaps'
                ]
                if newName:
                    # renames the segment name "sampleName-idx" where idx is the segment number
                    args.extend(['--name', newName + "-" + str(idx + 1)])
                args = pmc.parse_args(args)
                args.func_main(args)

                # clean up
                os.unlink(concat_file)
                os.unlink(ref_file)
                os.unlink(actual_file)
                os.unlink(aligned_file)
                tempFastas.append(tmpOutputFile)

    # merge outputs
    util.file.concat(tempFastas, outFasta)
    for tmpFile in tempFastas:
        os.unlink(tmpFile)

    # Index final output FASTA for Picard/GATK, Samtools, and Novoalign
    if index:
        tools.samtools.SamtoolsTool().faidx(outFasta, overwrite=True)
        tools.picard.CreateSequenceDictionaryTool().execute(outFasta, overwrite=True)
        tools.novoalign.NovoalignTool().index_fasta(outFasta)

    return 0


def parser_impute_from_reference(parser=argparse.ArgumentParser()):
    parser.add_argument(
        'inFasta',
        help='Input assembly/contigs, FASTA format, already ordered, oriented and merged with inReference.'
    )
    parser.add_argument('inReference', help='Reference genome to impute with, FASTA format.')
    parser.add_argument('outFasta', help='Output assembly, FASTA format.')
    parser.add_argument("--newName", default=None, help="rename output chromosome (default: do not rename)")
    parser.add_argument(
        "--minLengthFraction",
        type=float,
        default=0.5,
        help="minimum length for contig, as fraction of reference (default: %(default)s)"
    )
    parser.add_argument(
        "--minUnambig",
        type=float,
        default=0.5,
        help="minimum percentage unambiguous bases for contig (default: %(default)s)"
    )
    parser.add_argument(
        "--replaceLength",
        type=int,
        default=0,
        help="length of ends to be replaced with reference (default: %(default)s)"
    )
    parser.add_argument(
        '--aligner',
        help="""which method to use to align inFasta to
                        inReference. "muscle" = MUSCLE, "mafft" = MAFFT,
                        "mummer" = nucmer.  [default: %(default)s]""",
        choices=['muscle', 'mafft', 'mummer'],
        default='muscle'
    )
    parser.add_argument(
        "--index",
        help="""Index outFasta for Picard/GATK, Samtools, and Novoalign.""",
        default=False,
        action="store_true",
        dest="index"
    )
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, impute_from_reference, split_args=True)
    return parser


__commands__.append(('impute_from_reference', parser_impute_from_reference))


def refine_assembly(
    inFasta,
    inBam,
    outFasta,
    outVcf=None,
    outBam=None,
    novo_params='',
    min_coverage=2,
    major_cutoff=0.5,
    chr_names=None,
    keep_all_reads=False,
    already_realigned_bam=None,
    JVMmemory=None,
    threads=None,
    gatk_path=None,
    novoalign_license_path=None
):
    ''' This a refinement step where we take a crude assembly, align
        all reads back to it, and modify the assembly to the majority
        allele at each position based on read pileups.
        This step considers both SNPs as well as indels called by GATK
        and will correct the consensus based on GATK calls.
        Reads are aligned with Novoalign, then PCR duplicates are removed
        with Picard (in order to debias the allele counts in the pileups),
        and realigned with GATK's IndelRealigner (in order to call indels).
        Output FASTA file is indexed for Picard, Samtools, and Novoalign.
    '''
    chr_names = chr_names or []

    # if the input fasta is empty, create an empty output fasta and return
    if (os.path.getsize(inFasta) == 0):
        util.file.touch(outFasta)
        return 0

    # Get tools
    picard_index = tools.picard.CreateSequenceDictionaryTool()
    picard_mkdup = tools.picard.MarkDuplicatesTool()
    samtools = tools.samtools.SamtoolsTool()
    novoalign = tools.novoalign.NovoalignTool(license_path=novoalign_license_path)
    gatk = tools.gatk.GATKTool(path=gatk_path)

    # Create deambiguated genome for GATK
    deambigFasta = util.file.mkstempfname('.deambig.fasta')
    deambig_fasta(inFasta, deambigFasta)
    picard_index.execute(deambigFasta, overwrite=True)
    samtools.faidx(deambigFasta, overwrite=True)

    if already_realigned_bam:
        realignBam = already_realigned_bam
    else:
        # Novoalign reads to self
        novoBam = util.file.mkstempfname('.novoalign.bam')
        min_qual = 0 if keep_all_reads else 1
        novoalign.execute(inBam, inFasta, novoBam, options=novo_params.split(), min_qual=min_qual, JVMmemory=JVMmemory)
        rmdupBam = util.file.mkstempfname('.rmdup.bam')
        opts = ['CREATE_INDEX=true']
        if not keep_all_reads:
            opts.append('REMOVE_DUPLICATES=true')
        picard_mkdup.execute([novoBam], rmdupBam, picardOptions=opts, JVMmemory=JVMmemory)
        os.unlink(novoBam)
        realignBam = util.file.mkstempfname('.realign.bam')
        gatk.local_realign(rmdupBam, deambigFasta, realignBam, JVMmemory=JVMmemory, threads=threads)
        os.unlink(rmdupBam)
        if outBam:
            shutil.copyfile(realignBam, outBam)

    # Modify original assembly with VCF calls from GATK
    tmpVcf = util.file.mkstempfname('.vcf.gz')
    tmpFasta = util.file.mkstempfname('.fasta')
    gatk.ug(realignBam, deambigFasta, tmpVcf, JVMmemory=JVMmemory, threads=threads)
    if already_realigned_bam is None:
        os.unlink(realignBam)
    os.unlink(deambigFasta)
    name_opts = []
    if chr_names:
        name_opts = ['--name'] + chr_names
    main_vcf_to_fasta(
        parser_vcf_to_fasta(argparse.ArgumentParser(
        )).parse_args([
            tmpVcf,
            tmpFasta,
            '--trim_ends',
            '--min_coverage',
            str(min_coverage),
            '--major_cutoff',
            str(major_cutoff)
        ] + name_opts)
    )
    if outVcf:
        shutil.copyfile(tmpVcf, outVcf)
        if outVcf.endswith('.gz'):
            shutil.copyfile(tmpVcf + '.tbi', outVcf + '.tbi')
    os.unlink(tmpVcf)
    shutil.copyfile(tmpFasta, outFasta)
    os.unlink(tmpFasta)

    # Index final output FASTA for Picard/GATK, Samtools, and Novoalign
    picard_index.execute(outFasta, overwrite=True)
    # if the input bam is empty, an empty fasta will be created, however
    # faidx cannot index an empty fasta file, so only index if the fasta
    # has a non-zero size
    if (os.path.getsize(outFasta) > 0):
        samtools.faidx(outFasta, overwrite=True)
        novoalign.index_fasta(outFasta)
        
    return 0


def parser_refine_assembly(parser=argparse.ArgumentParser()):
    parser.add_argument(
        'inFasta', help='Input assembly, FASTA format, pre-indexed for Picard, Samtools, and Novoalign.'
    )
    parser.add_argument('inBam', help='Input reads, unaligned BAM format.')
    parser.add_argument(
        'outFasta',
         help='Output refined assembly, FASTA format, indexed for Picard, Samtools, and Novoalign.'
    )
    parser.add_argument(
        '--already_realigned_bam',
        default=None,
        help="""BAM with reads that are already aligned to inFasta.
            This bypasses the alignment process by novoalign and instead uses the given
            BAM to make an assemble. When set, outBam is ignored."""
    )
    parser.add_argument(
        '--outBam',
        default=None,
        help='Reads aligned to inFasta. Unaligned and duplicate reads have been removed. GATK indel realigned.'
    )
    parser.add_argument('--outVcf', default=None, help='GATK genotype calls for genome in inFasta coordinate space.')
    parser.add_argument(
        '--min_coverage',
        default=3,
        type=int,
        help='Minimum read coverage required to call a position unambiguous.'
    )
    parser.add_argument(
        '--major_cutoff',
        default=0.5,
        type=float,
        help="""If the major allele is present at a frequency higher than this cutoff,
            we will call an unambiguous base at that position.  If it is equal to or below
            this cutoff, we will call an ambiguous base representing all possible alleles at
            that position. [default: %(default)s]"""
    )
    parser.add_argument(
        '--novo_params',
        default='-r Random -l 40 -g 40 -x 20 -t 100',
        help='Alignment parameters for Novoalign.'
    )
    parser.add_argument(
        "--chr_names",
        dest="chr_names",
        nargs="*",
        help="Rename all output chromosomes (default: retain original chromosome names)",
        default=[]
    )
    parser.add_argument(
        "--keep_all_reads",
        help="""Retain all reads in BAM file? Default is to remove unaligned and duplicate reads.""",
        default=False,
        action="store_true",
        dest="keep_all_reads"
    )
    parser.add_argument(
        '--JVMmemory',
        default=tools.gatk.GATKTool.jvmMemDefault,
        help='JVM virtual memory size (default: %(default)s)'
    )
    parser.add_argument(
        '--GATK_PATH',
        default=None,
        dest="gatk_path",
        help='A path containing the GATK jar file. This overrides the GATK_ENV environment variable or the GATK conda package. (default: %(default)s)'
    )
    parser.add_argument(
        '--NOVOALIGN_LICENSE_PATH',
        default=None,
        dest="novoalign_license_path",
        help='A path to the novoalign.lic file. This overrides the NOVOALIGN_LICENSE_PATH environment variable. (default: %(default)s)'
    )
    util.cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, refine_assembly, split_args=True)
    return parser


__commands__.append(('refine_assembly', parser_refine_assembly))


def parser_filter_short_seqs(parser=argparse.ArgumentParser()):
    parser.add_argument("inFile", help="input sequence file")
    parser.add_argument("minLength", help="minimum length for contig", type=int)
    parser.add_argument("minUnambig", help="minimum percentage unambiguous bases for contig", type=float)
    parser.add_argument("outFile", help="output file")
    parser.add_argument("-f", "--format", help="Format for input sequence (default: %(default)s)", default="fasta")
    parser.add_argument(
        "-of", "--output-format",
        help="Format for output sequence (default: %(default)s)",
        default="fasta"
    )
    util.cmd.common_args(parser, (('loglevel', None), ('version', None)))
    util.cmd.attach_main(parser, main_filter_short_seqs)
    return parser


def main_filter_short_seqs(args):
    '''Check sequences in inFile, retaining only those that are at least minLength'''
    # orig by rsealfon, edited by dpark
    # TO DO: make this more generalized to accept multiple minLengths (for multiple chromosomes/segments)
    with util.file.open_or_gzopen(args.inFile) as inf:
        with util.file.open_or_gzopen(args.outFile, 'w') as outf:
            Bio.SeqIO.write(
                [
                    s for s in Bio.SeqIO.parse(inf, args.format)
                    if len(s) >= args.minLength and util.misc.unambig_count(s.seq) >= len(s) * args.minUnambig
                ], outf, args.output_format
            )
    return 0


__commands__.append(('filter_short_seqs', parser_filter_short_seqs))


def parser_modify_contig(parser=argparse.ArgumentParser()):
    parser.add_argument("input", help="input alignment of reference and contig (should contain exactly 2 sequences)")
    parser.add_argument("output", help="Destination file for modified contigs")
    parser.add_argument("ref", help="reference sequence name (exact match required)")
    parser.add_argument(
        "-n", "--name",
        type=str, help="fasta header output name (default: existing header)",
        default=None
    )
    parser.add_argument(
        "-cn",
        "--call-reference-ns",
        help="""should the reference sequence be called if there is an
        N in the contig and a more specific base in the reference (default: %(default)s)""",
        default=False,
        action="store_true",
        dest="call_reference_ns"
    )
    parser.add_argument(
        "-t",
        "--trim-ends",
        help="should ends of contig.fasta be trimmed to length of reference (default: %(default)s)",
        default=False,
        action="store_true",
        dest="trim_ends"
    )
    parser.add_argument(
        "-r5",
        "--replace-5ends",
        help="should the 5'-end of contig.fasta be replaced by reference (default: %(default)s)",
        default=False,
        action="store_true",
        dest="replace_5ends"
    )
    parser.add_argument(
        "-r3",
        "--replace-3ends",
        help="should the 3'-end of contig.fasta be replaced by reference (default: %(default)s)",
        default=False,
        action="store_true",
        dest="replace_3ends"
    )
    parser.add_argument(
        "-l",
        "--replace-length",
        help="length of ends to be replaced (if replace-ends is yes) (default: %(default)s)",
        default=10,
        type=int
    )
    parser.add_argument("-f", "--format", help="Format for input alignment (default: %(default)s)", default="fasta")
    parser.add_argument(
        "-r",
        "--replace-end-gaps",
        help="Replace gaps at the beginning and end of the sequence with reference sequence (default: %(default)s)",
        default=False,
        action="store_true",
        dest="replace_end_gaps"
    )
    parser.add_argument(
        "-rn",
        "--remove-end-ns",
        help="Remove leading and trailing N's in the contig (default: %(default)s)",
        default=False,
        action="store_true",
        dest="remove_end_ns"
    )
    parser.add_argument(
        "-ca",
        "--call-reference-ambiguous",
        help="""should the reference sequence be called if the contig seq is ambiguous and
        the reference sequence is more informative & consistant with the ambiguous base
        (ie Y->C) (default: %(default)s)""",
        default=False,
        action="store_true",
        dest="call_reference_ambiguous"
    )
    util.cmd.common_args(parser, (('tmp_dir', None), ('loglevel', None), ('version', None)))
    util.cmd.attach_main(parser, main_modify_contig)
    return parser


def main_modify_contig(args):
    ''' Modifies an input contig. Depending on the options
        selected, can replace N calls with reference calls, replace ambiguous
        calls with reference calls, trim to the length of the reference, replace
        contig ends with reference calls, and trim leading and trailing Ns.
        Author: rsealfon.
    '''
    aln = Bio.AlignIO.read(args.input, args.format)

    # TODO?: take list of alignments in, one per chromosome, rather than
    #       single alignment

    if len(aln) != 2:
        raise Exception("alignment does not contain exactly 2 sequences, %s found" % len(aln))
    elif aln[0].name == args.ref:
        ref_idx = 0
        consensus_idx = 1
    elif aln[1].name == args.ref:
        ref_idx = 1
        consensus_idx = 0
    else:
        raise NameError("reference name '%s' not in alignment" % args.ref)

    mc = ContigModifier(str(aln[ref_idx].seq), str(aln[consensus_idx].seq))
    if args.remove_end_ns:
        mc.remove_end_ns()
    if args.call_reference_ns:
        mc.call_reference_ns()
    if args.call_reference_ambiguous:
        mc.call_reference_ambiguous()
    if args.trim_ends:
        mc.trim_ends()
    if args.replace_end_gaps:
        mc.replace_end_gaps()
    if args.replace_5ends:
        mc.replace_5ends(args.replace_length)
    if args.replace_3ends:
        mc.replace_3ends(args.replace_length)

    with open(args.output, "wt") as f:
        if hasattr(args, "name"):
            name = args.name
        else:
            name = aln[consensus_idx].name
        for line in util.file.fastaMaker([(name, mc.get_stripped_consensus())]):
            f.write(line)
    return 0


__commands__.append(('modify_contig', parser_modify_contig))


class ContigModifier(object):
    ''' Initial modifications to spades+MUMmer assembly output based on
        MUSCLE alignment to known reference genome
        author: rsealfon
    '''

    def __init__(self, ref, consensus):
        if len(ref) != len(consensus):
            raise Exception("improper alignment")
        self.ref = list(ref)
        self.consensus = list(consensus)
        self.len = len(ref)

    def get_stripped_consensus(self):
        return ''.join(self.consensus).replace('-', '')

    def call_reference_ns(self):
        log.debug("populating N's from reference...")
        for i in range(self.len):
            if self.consensus[i].upper() == "N":
                self.consensus[i] = self.ref[i]

    def call_reference_ambiguous(self):
        ''' This is not normally used by default in our pipeline '''
        log.debug("populating ambiguous bases from reference...")
        for i in range(self.len):
            if self.ref[i].upper() in Bio.Data.IUPACData.ambiguous_dna_values.get(self.consensus[i].upper(), []):
                self.consensus[i] = self.ref[i]

    def trim_ends(self):
        ''' This trims down the consensus so it cannot go beyond the given reference genome '''
        log.debug("trimming ends...")
        for end_iterator in (range(self.len), reversed(range(self.len))):
            for i in end_iterator:
                if self.ref[i] != "-":
                    break
                else:
                    self.consensus[i] = "-"

    def replace_end_gaps(self):
        ''' This fills out the ends of the consensus with reference sequence '''
        log.debug("populating leading and trailing gaps from reference...")
        for end_iterator in (range(self.len), reversed(range(self.len))):
            for i in end_iterator:
                if self.consensus[i] != "-":
                    break
                self.consensus[i] = self.ref[i]

    def replace_5ends(self, replace_length):
        ''' This replaces everything within <replace_length> of the ends of the
            reference genome with the reference genome.
        '''
        log.debug("replacing 5' ends...")
        ct = 0
        for i in range(self.len):
            if self.ref[i] != "-":
                ct = ct + 1
            if ct == replace_length:
                for j in range(i + 1):
                    self.consensus[j] = self.ref[j]
                break

    def replace_3ends(self, replace_length):
        log.debug("replacing 3' ends...")
        ct = 0
        for i in reversed(range(self.len)):
            if self.ref[i] != "-":
                ct = ct + 1
            if ct == replace_length:
                for j in range(i, self.len):
                    self.consensus[j] = self.ref[j]
                break

    def remove_end_ns(self):
        ''' This clips off any N's that begin or end the consensus.
            Not normally used in our pipeline
        '''
        log.debug("removing leading and trailing N's...")
        for end_iterator in (range(self.len), reversed(range(self.len))):
            for i in end_iterator:
                if (self.consensus[i].upper() == "N" or self.consensus[i] == "-"):
                    self.consensus[i] = "-"
                else:
                    break


class MutableSequence(object):

    def __init__(self, name, start, stop, init_seq=None):
        if not (stop >= start >= 1):
            raise IndexError("coords out of bounds")
        if init_seq is None:
            self.seq = list('N' * (stop - start + 1))
        else:
            self.seq = list(init_seq)
        if stop - start + 1 != len(self.seq):
            raise Exception("wrong length")
        self.start = start
        self.stop = stop
        self.name = name
        self.deletions = []

    def modify(self, p, new_base):
        if not (self.start <= p <= self.stop):
            raise IndexError("position out of bounds")
        i = p - self.start
        self.seq[i] = new_base

    def replace(self, start, stop, new_seq):
        if stop > start:
            self.deletions.append((start, stop, new_seq))
        self.__change__(start, stop, new_seq)

    def __change__(self, start, stop, new_seq):
        if not (self.start <= start <= stop <= self.stop):
            raise IndexError("positions out of bounds")
        start -= self.start
        stop -= self.start
        if start == stop:
            self.seq[start] = new_seq
        for i in range(max(stop - start + 1, len(new_seq))):
            if start + i <= stop:
                if i < len(new_seq):
                    if start + i == stop:
                        # new allele is >= ref length, fill out the rest of the bases
                        self.seq[start + i] = new_seq[i:]
                    else:
                        self.seq[start + i] = new_seq[i]
                else:
                    # new allele is shorter than ref, so delete extra bases
                    self.seq[start + i] = ''

    def replay_deletions(self):
        for start, stop, new_seq in self.deletions:
            self.__change__(start, stop, new_seq)

    def emit(self):
        return (self.name, ''.join(self.seq))


def alleles_to_ambiguity(allelelist):
    ''' Convert a list of DNA bases to a single ambiguity base.
        All alleles must be one base long.  '''
    for a in allelelist:
        if len(a) != 1:
            raise Exception("all alleles must be one base long")
    if len(allelelist) == 1:
        return allelelist[0]
    else:
        convert = dict([(tuple(sorted(v)), k) for k, v in Bio.Data.IUPACData.ambiguous_dna_values.items() if k != 'X'])
        key = tuple(sorted(set(a.upper() for a in allelelist)))
        return convert[key]


def vcfrow_parse_and_call_snps(vcfrow, samples, min_dp=0, major_cutoff=0.5, min_dp_ratio=0.0):
    ''' Parse a single row of a VCF file, emit an iterator over each sample,
        call SNP genotypes using custom viral method based on read counts.
    '''
    # process this row
    c = vcfrow[0]
    alleles = [vcfrow[3]] + [a for a in vcfrow[4].split(',') if a not in '.']
    start = int(vcfrow[1])
    stop = start + len(vcfrow[3]) - 1
    format_col = vcfrow[8].split(':')
    format_col = dict((format_col[i], i) for i in range(len(format_col)))
    assert 'GT' in format_col and format_col['GT'] == 0    # required by VCF spec
    assert len(vcfrow) == 9 + len(samples)
    info = [x.split('=') for x in vcfrow[7].split(';') if x != '.']
    info = dict(x for x in info if len(x) == 2)
    info_dp = int(info.get('DP', 0))

    # process each sample
    for i in range(len(samples)):
        sample = samples[i]
        rec = vcfrow[i + 9].split(':')

        # require a minimum read coverage
        if len(alleles) == 1:
            # simple invariant case
            dp = ('DP' in format_col and len(rec) > format_col['DP']) and int(rec[format_col['DP']]) or 0
            if dp < min_dp:
                continue
            geno = alleles
            if info_dp and float(dp) / info_dp < min_dp_ratio:
                log.warning(
                    "dropping invariant call at %s:%s-%s %s (%s) due to low DP ratio (%s / %s = %s < %s)", c, start,
                    stop, sample, geno, dp, info_dp, float(dp) / info_dp, min_dp_ratio
                )
                continue
        else:
            # variant: manually call the highest read count allele if it exceeds a threshold
            assert ('AD' in format_col and len(rec) > format_col['AD'])
            allele_depths = list(map(int, rec[format_col['AD']].split(',')))
            assert len(allele_depths) == len(alleles)
            allele_depths = [(allele_depths[i], alleles[i]) for i in range(len(alleles)) if allele_depths[i] > 0]
            allele_depths = list(reversed(sorted((n, a) for n, a in allele_depths if n >= min_dp)))
            if not allele_depths:
                continue
            dp = sum(n for n, a in allele_depths)

            if allele_depths[0][0] > (dp * major_cutoff):
                # call a single allele at this position if it is a clear winner
                geno = [allele_depths[0][1]]
            else:
                # call multiple alleles at this position if there is no clear winner
                geno = [a for _, a in allele_depths]
        if geno:
            yield (c, start, stop, sample, geno)


def vcf_to_seqs(vcfIter, chrlens, samples, min_dp=0, major_cutoff=0.5, min_dp_ratio=0.0):
    ''' Take a VCF iterator and produce an iterator of chromosome x sample full sequences.'''
    seqs = {}
    cur_c = None
    for vcfrow in vcfIter:
        try:
            for c, start, stop, s, alleles in vcfrow_parse_and_call_snps(
                vcfrow, samples, min_dp=min_dp,
                major_cutoff=major_cutoff,
                min_dp_ratio=min_dp_ratio
            ):
                # changing chromosome?
                if c != cur_c:
                    if cur_c is not None:
                        # dump the previous chromosome before starting a new one
                        for s in samples:
                            seqs[s].replay_deletions()    # because of the order of VCF rows with indels
                            yield seqs[s].emit()

                    # prepare base sequences for this chromosome
                    cur_c = c
                    for s in samples:
                        name = len(samples) > 1 and ("%s-%s" % (c, s)) or c
                        seqs[s] = MutableSequence(name, 1, chrlens[c])

                # modify sequence for this chromosome/sample/position
                if len(alleles) == 1:
                    # call a single allele
                    seqs[s].replace(start, stop, alleles[0])
                elif all(len(a) == 1 for a in alleles):
                    # call an ambiguous SNP
                    seqs[s].replace(start, stop, alleles_to_ambiguity(alleles))
                else:
                    # mix of indels with no clear winner... force the most popular one
                    seqs[s].replace(start, stop, alleles[0])
        except:
            log.exception("Exception occurred while parsing VCF file.  Row: '%s'", vcfrow)
            raise

    # at the end, dump the last chromosome
    if cur_c is not None:
        for s in samples:
            seqs[s].replay_deletions()    # because of the order of VCF rows with indels
            yield seqs[s].emit()


def parser_vcf_to_fasta(parser=argparse.ArgumentParser()):
    parser.add_argument("inVcf", help="Input VCF file")
    parser.add_argument("outFasta", help="Output FASTA file")
    parser.add_argument(
        "--trim_ends",
        action="store_true",
        dest="trim_ends",
        default=False,
        help="""If specified, we will strip off continuous runs of N's from the beginning
        and end of the sequences before writing to output.  Interior N's will not be
        changed."""
    )
    parser.add_argument(
        "--min_coverage",
        dest="min_dp",
        type=int,
        help="""Specify minimum read coverage (with full agreement) to make a call.
        [default: %(default)s]""",
        default=3
    )
    parser.add_argument(
        "--major_cutoff",
        dest="major_cutoff",
        type=float,
        help="""If the major allele is present at a frequency higher than this cutoff,
        we will call an unambiguous base at that position.  If it is equal to or below
        this cutoff, we will call an ambiguous base representing all possible alleles at
        that position. [default: %(default)s]""",
        default=0.5
    )
    parser.add_argument(
        "--min_dp_ratio",
        dest="min_dp_ratio",
        type=float,
        help="""The input VCF file often reports two read depth values (DP)--one for
        the position as a whole, and one for the sample in question.  We can optionally
        reject calls in which the sample read count is below a specified fraction of the
        total read count.  This filter will not apply to any sites unless both DP values
        are reported.  [default: %(default)s]""",
        default=0.0
    )
    parser.add_argument(
        "--name",
        dest="name",
        nargs="*",
        help="output sequence names (default: reference names in VCF file)",
        default=[]
    )
    util.cmd.common_args(parser, (('loglevel', None), ('version', None)))
    util.cmd.attach_main(parser, main_vcf_to_fasta)
    return parser


def main_vcf_to_fasta(args):
    ''' Take input genotypes (VCF) and construct a consensus sequence
        (fasta) by using majority-read-count alleles in the VCF.
        Genotypes in the VCF will be ignored--we will use the allele
        with majority read support (or an ambiguity base if there is no clear majority).
        Uncalled positions will be emitted as N's.
        Author: dpark.
    '''
    assert args.min_dp >= 0
    assert 0.0 <= args.major_cutoff < 1.0

    with assemble.vcf.VcfReader(args.inVcf) as vcf:
        chrlens = dict(vcf.chrlens())
        samples = vcf.samples()

    assert len(
        samples
    ) == 1, """Multiple sample columns were found in the intermediary VCF file
        of the refine_assembly step, suggesting multiple sample names are present
        upstream in the BAM file. Please correct this so there is only one sample in the BAM file."""

    with open(args.outFasta, 'wt') as outf:
        chr_idx = 0
        for chr_idx, (header, seq) in enumerate(
            vcf_to_seqs(
                util.file.read_tabfile(args.inVcf),
                chrlens,
                samples,
                min_dp=args.min_dp,
                major_cutoff=args.major_cutoff,
                min_dp_ratio=args.min_dp_ratio
            )
        ):
            if args.trim_ends:
                seq = seq.strip('Nn')
            if args.name:
                header = args.name[chr_idx % len(args.name)]
            for line in util.file.fastaMaker([(header, seq)]):
                outf.write(line)

    # done
    log.info("done")
    return 0


__commands__.append(('vcf_to_fasta', parser_vcf_to_fasta))


def parser_trim_fasta(parser=argparse.ArgumentParser()):
    parser.add_argument("inFasta", help="Input fasta file")
    parser.add_argument("outFasta", help="Output (trimmed) fasta file")
    util.cmd.common_args(parser, (('loglevel', None), ('version', None)))
    util.cmd.attach_main(parser, trim_fasta, split_args=True)
    return parser


def trim_fasta(inFasta, outFasta):
    ''' Take input sequences (fasta) and trim any continuous sections of
        N's from the ends of them.  Write trimmed sequences to an output fasta file.
    '''
    with open(outFasta, 'wt') as outf:
        with open(inFasta, 'rt') as inf:
            for record in Bio.SeqIO.parse(inf, 'fasta'):
                for line in util.file.fastaMaker([(record.id, str(record.seq).strip('Nn'))]):
                    outf.write(line)
    log.info("done")
    return 0


__commands__.append(('trim_fasta', parser_trim_fasta))


def deambig_base(base):
    ''' Take a single base (possibly a IUPAC ambiguity code) and return a random
        non-ambiguous base from among the possibilities '''
    return random.choice(Bio.Data.IUPACData.ambiguous_dna_values[base.upper()])


def deambig_fasta(inFasta, outFasta):
    ''' Take input sequences (fasta) and replace any ambiguity bases with a
        random unambiguous base from among the possibilities described by the ambiguity
        code.  Write output to fasta file.
    '''
    with util.file.open_or_gzopen(outFasta, 'wt') as outf:
        with util.file.open_or_gzopen(inFasta, 'rt') as inf:
            for record in Bio.SeqIO.parse(inf, 'fasta'):
                for line in util.file.fastaMaker([(record.id, ''.join(map(deambig_base, str(record.seq))))]):
                    outf.write(line)
    return 0


def parser_deambig_fasta(parser=argparse.ArgumentParser()):
    parser.add_argument("inFasta", help="Input fasta file")
    parser.add_argument("outFasta", help="Output fasta file")
    util.cmd.common_args(parser, (('loglevel', None), ('version', None)))
    util.cmd.attach_main(parser, deambig_fasta, split_args=True)
    return parser


__commands__.append(('deambig_fasta', parser_deambig_fasta))


def vcf_dpdiff(vcfs):
    for vcf in vcfs:
        with assemble.vcf.VcfReader(vcf) as v:
            samples = v.samples()
        assert len(samples) == 1
        for row in util.file.read_tabfile(vcf):
            dp1 = int(dict(x.split('=') for x in row[7].split(';') if x != '.').get('DP', 0))
            dp2 = 0
            if 'DP' in row[8].split(':'):
                dpidx = row[8].split(':').index('DP')
                if len(row[9].split(':')) > dpidx:
                    dp2 = int(row[9].split(':')[dpidx])
            ratio = ''
            if dp1:
                ratio = float(dp2) / dp1
            yield (row[0], row[1], samples[0], dp1, dp2, dp1 - dp2, ratio)


def parser_dpdiff(parser=argparse.ArgumentParser()):
    parser.add_argument("inVcfs", help="Input VCF file", nargs='+')
    parser.add_argument("outFile", help="Output flat file")
    util.cmd.common_args(parser, (('loglevel', None), ('version', None)))
    util.cmd.attach_main(parser, dpdiff, split_args=True)
    return parser


def dpdiff(inVcfs, outFile):
    ''' Take input VCF files (all with only one sample each) and report
        on the discrepancies between the two DP fields (one in INFO and one in the
        sample's genotype column).
    '''
    header = ['chr', 'pos', 'sample', 'dp_info', 'dp_sample', 'diff', 'ratio']
    with open(outFile, 'wt') as outf:
        outf.write('#' + '\t'.join(header) + '\n')
        for row in vcf_dpdiff(inVcfs):
            outf.write('\t'.join(map(str, row)) + '\n')
    return 0


#__commands__.append(('dpdiff', parser_dpdiff))



def alignment_summary(inFastaFileOne, inFastaFileTwo, outfileName=None, printCounts=False):
    """ Write or print pairwise alignment summary information for sequences in two FASTA
        files, including SNPs, ambiguous bases, and indels.
    """
    gap = '-'
    ambiguous = 'N'
    unambiguous = set(list("GATCgatc"))
    aligner = assemble.muscle.MuscleTool()

    per_chr_fastas = util.file.transposeChromosomeFiles([inFastaFileOne, inFastaFileTwo])

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

                if (c1 in unambiguous
                   and c2 in unambiguous):
                    unambig_both += 1
                    if c1 == c2:
                        same_unambig += 1
                    else:
                        snp_unambig += 1

                if ((c1 == gap and
                    c2 in unambiguous) or
                   (c2 == gap and
                    c1 in unambiguous)):
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
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, alignment_summary, split_args=True)
    return parser
__commands__.append(('alignment_summary', parser_alignment_summary))




def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)


if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)
