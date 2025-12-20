'''
    The minimap2 aligner.

'''

import collections
import logging
import os
import os.path
import shutil
import subprocess

from Bio import SeqIO

import tools
import tools.samtools
import tools.picard
import util.file
import util.misc
from errors import *

TOOL_NAME = 'minimap2'

log = logging.getLogger(__name__)

class Minimap2(tools.Tool):

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = []
            install_methods.append(tools.PrexistingUnixCommand(shutil.which(TOOL_NAME), require_executability=True))
        super(Minimap2, self).__init__(install_methods=install_methods)

    def _get_tool_version(self):
        self.tool_version = subprocess.check_output([TOOL_NAME ,'--version']).decode('UTF-8').strip()

    def execute(self, args, stdout=None, stdin=None, background=False):    # pylint: disable=W0221
        tool_cmd = [self.install_and_get_path()] + args
        log.debug(' '.join(tool_cmd))
        if stdout:
            stdout = open(stdout, 'w')
        if background:
            subprocess.Popen(tool_cmd, stdout=stdout, stdin=stdin)
        else:
            subprocess.check_call(tool_cmd, stdout=stdout, stdin=stdin)
            if stdout:
                stdout.close()

    def align_bam(self, inBam, refDb, outBam, options=None,
                      threads=None, JVMmemory=None, should_index=True):
        options = options or []

        samtools = tools.samtools.SamtoolsTool()
        threads = util.misc.sanitize_thread_count(threads)

        # fetch list of RGs
        rgs = list(samtools.getReadGroups(inBam).keys())

        if len(rgs) == 0:
            # Can't do this
            raise InvalidBamHeaderError("{} lacks read groups".format(inBam))

        elif len(rgs) == 1:
            # Only one RG, keep it simple
            self.align_one_rg(inBam, refDb, outBam, options=options, threads=threads, should_index=should_index)

        else:
            # Multiple RGs, align one at a time and merge
            align_bams = []

            for rg in rgs:
                tmp_bam = util.file.mkstempfname('.{}.bam'.format(rg))
                self.align_one_rg(
                    inBam,
                    refDb,
                    tmp_bam,
                    rgid=rg,
                    options=options,
                    threads=threads,
                    should_index=False  # Don't index intermediate BAMs that will be merged
                )
                if not samtools.isEmpty(tmp_bam):
                    align_bams.append(tmp_bam)
                else:
                    log.warning("No alignment output for RG %s in file %s against %s", rg, inBam, refDb)

            if len(align_bams)==0:
                log.warning("All read groups in file %s appear to be empty.", inBam)
                with util.file.tempfname('.empty.sam') as empty_sam:
                    samtools.dumpHeader(inBam, empty_sam)
                    samtools.sort(empty_sam, outBam)
            else:
                # Merge BAMs, sort, and index
                picardOptions = ['SORT_ORDER=coordinate', 'USE_THREADING=true', 'CREATE_INDEX=true']
                tools.picard.MergeSamFilesTool().execute(
                    align_bams,
                    outBam,
                    picardOptions=picardOptions,
                    JVMmemory=JVMmemory
                )

                for bam in align_bams:
                    os.unlink(bam)

    def align_one_rg(self, inBam, refDb, outBam, rgid=None, preset=None, options=None,
                         threads=None, JVMmemory=None, should_index=True):
        """
            Performs an alignment of one read group in a bam file to a reference fasta file using minimap2.
            Emits alignments in sorted, index bam files.
            inBam may contain more read groups, but we will subset input to the specified rgid.
            preset may be specified as a valid value for "minimap2 -x" which depends on the type of
                data (short accurate reads vs long noisy reads). If preset is set to None, we will autodetect
                based on the PL (platform) tag in the read group header (e.g. illumina, ont, pacbio)
        """
        options = list(options).copy() or []

        samtools = tools.samtools.SamtoolsTool()

        # Require exactly one RG
        rgs = samtools.getReadGroups(inBam)
        if len(rgs) == 0:
            raise InvalidBamHeaderError("{} lacks read groups".format(inBam))
        elif len(rgs) == 1:
            if not rgid:
                rgid = list(rgs.keys())[0]
        elif not rgid:
            raise InvalidBamHeaderError("{} has {} read groups, but we require exactly one".format(inBam, len(rgs)))
        if rgid not in rgs:
            raise InvalidBamHeaderError("{} has read groups, but not {}".format(inBam, rgid))

        headerFile = util.file.mkstempfname('.{}.header.txt'.format(rgid))
        # Strip inBam to just one RG (if necessary)
        removeInput = False
        if len(rgs) == 1:
            one_rg_inBam = inBam
            tools.samtools.SamtoolsTool().dumpHeader(one_rg_inBam, headerFile)
        else:
            # strip inBam to one read group
            with util.file.tempfname('.onebam.bam') as tmp_bam:
                samtools.view(['-1', '-r', rgid], inBam, tmp_bam)
                # special exit if this file is empty
                if samtools.isEmpty(tmp_bam):
                    log.warning("No reads present for RG %s in file: %s", rgid, inBam)
                    shutil.copyfile(tmp_bam, outBam)
                    return
                # simplify BAM header otherwise Novoalign gets confused
                one_rg_inBam = util.file.mkstempfname('.{}.in.bam'.format(rgid))
                removeInput = True
                
                with open(headerFile, 'wt') as outf:
                    for row in samtools.getHeader(inBam):
                        if len(row) > 0 and row[0] == '@RG':
                            if rgid != list(x[3:] for x in row if x.startswith('ID:'))[0]:
                                # skip all read groups that are not rgid
                                continue
                        outf.write('\t'.join(row) + '\n')
                samtools.reheader(tmp_bam, headerFile, one_rg_inBam)

        # get the read group line to give to mm2
        readgroup_line = ""
        with open(headerFile) as inf:
            for line in inf:
                if line.startswith("@RG"):
                    readgroup_line = line.rstrip("\r\n")
        if not readgroup_line:
            raise Exception()
        # rather than reheader the alignment bam file later so it has the readgroup information
        # from the original bam file, we'll pass the RG line to minimap2 to write out
        options.extend(('-R', readgroup_line.replace('\t','\\t')))

        # dynamically determine the mode of operation
        if '-x' not in options:
            if preset is None:
                platform = list(x for x in readgroup_line.split('\t') if x.startswith('PL:'))
                if len(platform) != 1:
                    raise Exception("cannot autodetect minimap2 aligner mode when PL: tag is not set in the read group header for {}: {}".format(inBam, readgroup_line))
                else:
                    platform = platform[0][3:].lower()
                    if platform == 'illumina':
                        preset = 'sr'
                    elif platform == 'ont':
                        preset = 'map-ont'
                    elif platform == 'pacbio':
                        preset = 'map-pb'
                    else:
                        raise Exception("PL: tag {} for read group {} in bam {} refers to a data type we do not know how to map with minimap2".format(platform, rgid, inBam))
            options.extend(('-x', preset))

        # perform actual alignment
        if samtools.isEmpty(one_rg_inBam):
            log.warning("Input file %s appears to lack reads for RG '%s'", inBam, rgid)
            # minimap doesn't like empty inputs, so copy empty bam through
            # samtools.sort(one_rg_inBam, outBam)
            self.align_cmd(one_rg_inBam, refDb, outBam, options=options, threads=threads, should_index=should_index)
        else:
            self.align_cmd(one_rg_inBam, refDb, outBam, options=options, threads=threads, should_index=should_index)

        # if there was more than one RG in the input, we had to create a temporary file with the one RG specified
        # and we can safely delete it this file
        # if there was only one RG in the input, we used it directly and should not delete it
        if removeInput:
            os.unlink(one_rg_inBam)

    def align_cmd(self, inReads, refDb, outAlign, options=None, threads=None, should_index=True):
        options = [] if not options else options

        threads = util.misc.sanitize_thread_count(threads)
        if '-t' not in options:
            options.extend(('-t', str(threads)))
        if '-2' not in options:
            options.append('-2')

        samtools = tools.samtools.SamtoolsTool()

        with util.file.tempfname('.aligned.sam') as aln_sam:
            fastq_pipe = samtools.bam2fq_pipe(inReads)
            options.extend(('-a', refDb, '-', '-o', aln_sam))
            self.execute(options, stdin=fastq_pipe.stdout)
            if fastq_pipe.wait():
                raise subprocess.CalledProcessError(fastq_pipe.returncode, "samtools.bam2fq_pipe() for {}".format(inReads))
            samtools.sort(aln_sam, outAlign, threads=threads)

        # cannot index sam files; only do so if a bam/cram is desired
        if should_index and (outAlign.endswith(".bam") or outAlign.endswith(".cram")):
            samtools.index(outAlign, threads=threads)

    def scaffold(self, contigs_fasta, ref_fasta, outAlign, divergence=20, options=None, threads=None):
        options = [] if not options else options

        threads = util.misc.sanitize_thread_count(threads)
        if '-t' not in options:
            options.extend(('-t', str(threads)))
        if '-2' not in options:
            options.append('-2')

        if divergence >= 20:
            options.extend(('-x', 'asm20'))
        elif divergence >= 10:
            options.extend(('-x', 'asm10'))
        else:
            options.extend(('-x', 'asm5'))

        with util.file.tempfname('.aligned.sam') as aln_sam:
            options.extend(('-a', ref_fasta, contigs_fasta, '-o', aln_sam))
            self.execute(options)
            samtools.sort(aln_sam, outAlign, threads=threads)

        # cannot index sam files; only do so if a bam/cram is desired
        if (outAlign.endswith(".bam") or outAlign.endswith(".cram")):
            samtools.index(outAlign)

    def idxstats(self, inReads, refDb, outIdxstats, outReadlist=None, threads=None):
        """
        Align reads and produce idxstats-like counts and optional read list.

        This method is optimized for counting reads per reference by:
        - Using PAF output (no -a flag) for faster alignment
        - Streaming and parsing PAF directly without intermediate files
        - Only counting primary alignments
        - For paired-end reads: requiring both R1 and R2 to have primary alignments
          (properly paired behavior)
        - For single-end reads: counting any primary alignment

        Args:
            inReads: Input reads (BAM format)
            refDb: Reference database (FASTA format)
            outIdxstats: Output file in samtools idxstats format (tab-separated:
                         seq_name, seq_length, mapped_count, unmapped_count)
                         Note: unmapped_count will always be 0 in this implementation
            outReadlist: Optional output file listing all read IDs that mapped
                         as primary alignments. If None, this output is skipped.
            threads: Number of threads for minimap2 (default: auto-detect)
        """
        log.info("Starting idxstats: aligning %s to %s", inReads, refDb)

        threads = util.misc.sanitize_thread_count(threads)
        samtools = tools.samtools.SamtoolsTool()

        # Parse reference FASTA for sequence lengths
        ref_lengths = collections.OrderedDict()
        for record in SeqIO.parse(refDb, 'fasta'):
            ref_lengths[record.id] = len(record.seq)

        # Handle empty reference case
        if not ref_lengths:
            log.warning("Reference FASTA %s contains no sequences", refDb)
            with open(outIdxstats, 'wt') as outf:
                pass  # write empty file
            if outReadlist:
                with open(outReadlist, 'wt') as outf:
                    pass  # write empty file
            return

        # Initialize counter for reads per reference
        reads_per_ref = collections.Counter()

        # Build minimap2 command for PAF output (no -a flag)
        options = ['-t', str(threads), '-2', '-x', 'sr', refDb, '-']

        # Open readlist file for streaming if requested
        readlist_file = open(outReadlist, 'wt') if outReadlist else None

        try:
            # Start samtools bam2fq pipe for input
            fastq_pipe = samtools.bam2fq_pipe(inReads)

            # Start minimap2 process with PAF output to stdout
            tool_cmd = [self.install_and_get_path()] + options
            log.debug(' '.join(tool_cmd))
            mm2_proc = subprocess.Popen(
                tool_cmd,
                stdin=fastq_pipe.stdout,
                stdout=subprocess.PIPE,
                text=True
            )

            # Parse PAF output line by line
            # For paired-end reads: require both R1 and R2 to have primary alignments
            # (properly paired behavior). For single-end reads: count any primary alignment.
            # Track recent read names to deduplicate readlist output.
            recent_reads = collections.deque(maxlen=4)
            pending_alignment = None  # Tuple: (base_read_id, target_name, is_r1) for paired-end
            lines_processed = 0
            paired_count = 0
            single_count = 0
            for line in mm2_proc.stdout:
                lines_processed += 1
                if lines_processed % 100000 == 0:
                    log.info("Processed %d PAF lines, %d paired + %d single reads so far",
                             lines_processed, paired_count, single_count)

                fields = line.rstrip('\n').split('\t')
                if len(fields) < 12:
                    continue

                read_id = fields[0]
                target_name = fields[5]

                # Check for primary alignment tag (tp:A:P)
                is_primary = False
                for field in fields[12:]:
                    if field == 'tp:A:P':
                        is_primary = True
                        break

                if not is_primary:
                    continue

                # Detect if paired-end by checking for /1 or /2 suffix
                # bam2fq -n adds these suffixes for paired-end reads
                if read_id.endswith('/1'):
                    is_paired = True
                    is_r1 = True
                    base_read_id = read_id[:-2]
                elif read_id.endswith('/2'):
                    is_paired = True
                    is_r1 = False
                    base_read_id = read_id[:-2]
                else:
                    is_paired = False
                    base_read_id = read_id

                if not is_paired:
                    # SINGLE-END: count immediately
                    single_count += 1
                    reads_per_ref[target_name] += 1
                    if readlist_file and base_read_id not in recent_reads:
                        recent_reads.append(base_read_id)
                        readlist_file.write(base_read_id + '\n')
                else:
                    # PAIRED-END: require both R1 and R2 to have primary alignments
                    if pending_alignment is None:
                        # First of pair: store it
                        pending_alignment = (base_read_id, target_name, is_r1)
                    else:
                        pending_base, pending_target, pending_is_r1 = pending_alignment

                        if base_read_id == pending_base and is_r1 != pending_is_r1:
                            # Mate found! Both R1 and R2 aligned - count BOTH
                            paired_count += 2
                            reads_per_ref[pending_target] += 1
                            reads_per_ref[target_name] += 1
                            if readlist_file and base_read_id not in recent_reads:
                                recent_reads.append(base_read_id)
                                readlist_file.write(base_read_id + '\n')
                            pending_alignment = None
                        else:
                            # Not a matching mate - discard pending (singleton), store current
                            pending_alignment = (base_read_id, target_name, is_r1)

            # Wait for processes to complete
            mm2_proc.wait()
            if fastq_pipe.wait():
                raise subprocess.CalledProcessError(
                    fastq_pipe.returncode,
                    "samtools.bam2fq_pipe() for {}".format(inReads)
                )

        finally:
            if readlist_file:
                readlist_file.close()

        # Write idxstats output
        with open(outIdxstats, 'wt') as outf:
            for ref_name, ref_len in ref_lengths.items():
                mapped_count = reads_per_ref.get(ref_name, 0)
                outf.write('{}\t{}\t{}\t0\n'.format(ref_name, ref_len, mapped_count))

        log.info("Completed idxstats: %d paired + %d single-end reads counted from %d PAF lines",
                 paired_count, single_count, lines_processed)
