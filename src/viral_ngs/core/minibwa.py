'''
    The minibwa aligner.

    minibwa (https://github.com/lh3/minibwa) is a short-read (and accurate
    long-read) aligner that combines BWA-MEM variable-length seeding with
    minimap2's chaining and SIMD base alignment. It is exposed here as an
    opt-in alternative to bwa / minimap2 / novoalign for aligning reads to a
    reference FASTA.

    CLI shape (see the minibwa README):
        minibwa index ref.fa
        minibwa map -x sr -t <threads> ref.fa reads.fq > aln.sam

    Notes on the assumptions this wrapper makes about the minibwa CLI. minibwa
    is new and its help text does not (yet) document every flag; these mirror
    the behavior of its minimap2 lineage and should be re-verified against the
    installed binary:
      - `minibwa map` emits SAM to stdout by default (PAF only with `-f`), so
        this wrapper redirects stdout rather than passing an output-file flag.
      - `-x sr` selects the short-read preset and (like minimap2's `sr` preset)
        smart-pairs interleaved reads read from stdin (`-`).
      - `-R <line>` attaches a read-group header line, as in bwa/minimap2.
'''

import logging
import os
import os.path
import shutil
import subprocess

from . import samtools
from . import picard
from . import file as util_file, misc as util_misc
from .errors import *
from . import Tool, PrexistingUnixCommand

TOOL_NAME = 'minibwa'

log = logging.getLogger(__name__)


class Minibwa(Tool):

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [PrexistingUnixCommand(shutil.which(TOOL_NAME), require_executability=True)]
        super(Minibwa, self).__init__(install_methods=install_methods)

    def _get_tool_version(self):
        self.tool_version = subprocess.check_output([TOOL_NAME, '--version']).decode('UTF-8').strip()

    def execute(self, command, args, stdout=None, stdin=None, background=False):    # pylint: disable=W0221
        tool_cmd = [self.install_and_get_path(), command] + args
        log.debug(' '.join(tool_cmd))
        if stdout:
            stdout = open(stdout, 'w')
        if background:
            subprocess.Popen(tool_cmd, stdout=stdout, stdin=stdin)
        else:
            subprocess.check_call(tool_cmd, stdout=stdout, stdin=stdin)
            if stdout:
                stdout.close()

    def index(self, inFasta):
        self.execute('index', [inFasta])

    def align_bam(self, inBam, refDb, outBam, options=None,
                  threads=None, JVMmemory=None, should_index=True):
        options = options or []

        samtools_tool = samtools.SamtoolsTool()
        threads = util_misc.sanitize_thread_count(threads)

        # fetch list of RGs
        rgs = list(samtools_tool.getReadGroups(inBam).keys())

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
                tmp_bam = util_file.mkstempfname('.{}.bam'.format(rg))
                self.align_one_rg(
                    inBam,
                    refDb,
                    tmp_bam,
                    rgid=rg,
                    options=options,
                    threads=threads,
                    should_index=False  # Don't index intermediate BAMs that will be merged
                )
                if not samtools_tool.isEmpty(tmp_bam):
                    align_bams.append(tmp_bam)
                else:
                    log.warning("No alignment output for RG %s in file %s against %s", rg, inBam, refDb)

            if len(align_bams) == 0:
                log.warning("All read groups in file %s appear to be empty.", inBam)
                with util_file.tempfname('.empty.sam') as empty_sam:
                    samtools_tool.dumpHeader(inBam, empty_sam)
                    samtools_tool.sort(empty_sam, outBam)
            else:
                # Merge BAMs, sort, and index
                picardOptions = ['SORT_ORDER=coordinate', 'USE_THREADING=true', 'CREATE_INDEX=true']
                picard.MergeSamFilesTool().execute(
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
            Performs an alignment of one read group in a bam file to a reference fasta file using minibwa.
            Emits alignments in sorted, indexed bam files.
            inBam may contain more read groups, but we will subset input to the specified rgid.
            preset may be specified as a valid value for "minibwa map -x" (e.g. "sr" for short reads,
                "lr" for accurate long reads). If preset is None, we autodetect based on the PL
                (platform) tag in the read group header (e.g. illumina, ont, pacbio).
        """
        options = list(options).copy() if options else []

        samtools_tool = samtools.SamtoolsTool()

        # Require exactly one RG
        rgs = samtools_tool.getReadGroups(inBam)
        if len(rgs) == 0:
            raise InvalidBamHeaderError("{} lacks read groups".format(inBam))
        elif len(rgs) == 1:
            if not rgid:
                rgid = list(rgs.keys())[0]
        elif not rgid:
            raise InvalidBamHeaderError("{} has {} read groups, but we require exactly one".format(inBam, len(rgs)))
        if rgid not in rgs:
            raise InvalidBamHeaderError("{} has read groups, but not {}".format(inBam, rgid))

        headerFile = util_file.mkstempfname('.{}.header.txt'.format(rgid))
        # Strip inBam to just one RG (if necessary)
        removeInput = False
        if len(rgs) == 1:
            one_rg_inBam = inBam
            samtools_tool.dumpHeader(one_rg_inBam, headerFile)
        else:
            # strip inBam to one read group
            with util_file.tempfname('.onebam.bam') as tmp_bam:
                samtools_tool.view(['-1', '-r', rgid], inBam, tmp_bam)
                # special exit if this file is empty
                if samtools_tool.isEmpty(tmp_bam):
                    log.warning("No reads present for RG %s in file: %s", rgid, inBam)
                    shutil.copyfile(tmp_bam, outBam)
                    return
                one_rg_inBam = util_file.mkstempfname('.{}.in.bam'.format(rgid))
                removeInput = True

                with open(headerFile, 'wt') as outf:
                    for row in samtools_tool.getHeader(inBam):
                        if len(row) > 0 and row[0] == '@RG':
                            if rgid != list(x[3:] for x in row if x.startswith('ID:'))[0]:
                                # skip all read groups that are not rgid
                                continue
                        outf.write('\t'.join(row) + '\n')
                samtools_tool.reheader(tmp_bam, headerFile, one_rg_inBam)

        # get the read group line to give to minibwa
        readgroup_line = ""
        with open(headerFile) as inf:
            for line in inf:
                if line.startswith("@RG"):
                    readgroup_line = line.rstrip("\r\n")
        if not readgroup_line:
            raise Exception("{} lacks an @RG header line for read group {}".format(inBam, rgid))
        # rather than reheader the alignment bam file later, pass the RG line to minibwa to write out
        options.extend(('-R', readgroup_line.replace('\t', '\\t')))

        # dynamically determine the mode of operation
        if '-x' not in options:
            if preset is None:
                platform = list(x for x in readgroup_line.split('\t') if x.startswith('PL:'))
                if len(platform) != 1:
                    raise Exception("cannot autodetect minibwa aligner mode when PL: tag is not set in the read group header for {}: {}".format(inBam, readgroup_line))
                platform = platform[0][3:].lower()
                if platform == 'illumina':
                    preset = 'sr'
                elif platform in ('ont', 'pacbio'):
                    preset = 'lr'
                else:
                    raise Exception("PL: tag {} for read group {} in bam {} refers to a data type we do not know how to map with minibwa".format(platform, rgid, inBam))
            options.extend(('-x', preset))

        self.align_cmd(one_rg_inBam, refDb, outBam, options=options, threads=threads, should_index=should_index)

        # if there was more than one RG in the input, we created a temporary one-RG file and can delete it;
        # if there was only one RG, we used the input directly and should not delete it
        if removeInput:
            os.unlink(one_rg_inBam)

    def align_cmd(self, inReads, refDb, outAlign, options=None, threads=None, should_index=True):
        options = list(options) if options else []

        threads = util_misc.sanitize_thread_count(threads)
        if '-t' not in options:
            options.extend(('-t', str(threads)))

        samtools_tool = samtools.SamtoolsTool()

        with util_file.tempfname('.aligned.sam') as aln_sam:
            fastq_pipe = samtools_tool.bam2fq_pipe(inReads)
            # `minibwa map` emits SAM to stdout; positional args are the reference then the reads
            # ('-' = interleaved reads from stdin)
            self.execute('map', options + [refDb, '-'], stdout=aln_sam, stdin=fastq_pipe.stdout)
            if fastq_pipe.wait():
                raise subprocess.CalledProcessError(fastq_pipe.returncode, "samtools.bam2fq_pipe() for {}".format(inReads))
            samtools_tool.sort(aln_sam, outAlign, threads=threads)

        # cannot index sam files; only do so if a bam/cram is desired
        if should_index and (outAlign.endswith(".bam") or outAlign.endswith(".cram")):
            samtools_tool.index(outAlign, threads=threads)
