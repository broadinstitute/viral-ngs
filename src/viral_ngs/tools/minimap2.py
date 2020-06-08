'''
    The minimap2 aligner.

'''

from collections import defaultdict
import logging
import os
import os.path
import subprocess
import shutil
import concurrent.futures

import tools
import tools.samtools
import tools.picard
import util.file
import util.misc

TOOL_NAME = 'minimap2'

log = logging.getLogger(__name__)

class Bwa(tools.Tool):

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = []
            install_methods.append(tools.PrexistingUnixCommand(TOOL_NAME, verifycmd=False, require_executability=True))
        super(Bwa, self).__init__(install_methods=install_methods)

    def version(self):
        return subprocess.check_output([TOOL_NAME ,'--version']).decode('UTF-8').strip()

    def execute(self, args, stdout=None, stdin=None, background=False):    # pylint: disable=W0221
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

    def align_bam(self, inBam, refDb, outBam, options=None,
                      min_score_to_filter=None, threads=None, JVMmemory=None, invert_filter=False, should_index=True):
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
            self.align_one_rg(inBam, refDb, outBam, options=options,
                                  min_score_to_filter=min_score_to_filter,
                                  threads=threads, invert_filter=invert_filter)

        else:
            # Multiple RGs, align one at a time and merge
            align_bams = []

            threads_for_chunk = int(round(min(max(threads / len(rgs),1),threads),0))+1
            # worker count limited to 1 for now to reduce in-memory index size resulting from
            # running multiple copies of bwa in parallel
            workers = 1 #len(rgs) if len(rgs)<threads else threads
            with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
                futures = []# executor.submit(util.file.count_occurrences_in_tsv, filePath, include_noise=includeNoise) for rg in rgs]

                for rg in rgs:
                    tmp_bam = util.file.mkstempfname('.{}.bam'.format(rg))
                    futures.append(executor.submit(
                        self.align_one_rg,
                        inBam,
                        refDb,
                        tmp_bam,
                        rgid=rg,
                        options=options,
                        min_score_to_filter=min_score_to_filter,
                        threads=threads_for_chunk, 
                        invert_filter=invert_filter
                    ))

                for future in concurrent.futures.as_completed(futures):
                    if future.result():
                        rg, aln_bam = future.result()
                        if os.path.getsize(aln_bam) > 0:
                            align_bams.append(aln_bam)
                        else:
                            log.warning("No alignment output for RG %s in file %s against %s", rg, inBam, refDb)

            if len(align_bams) == 0:
                util.file.touch(outBam)
            else:
                # Merge BAMs, sort, and index
                picardOptions = ['SORT_ORDER=coordinate', 'USE_THREADING=true']
                if should_index:
                    picardOptions.append('CREATE_INDEX=true')
                tools.picard.MergeSamFilesTool().execute(
                    align_bams,
                    outBam,
                    picardOptions=picardOptions,
                    JVMmemory=JVMmemory
                )
                # no longer required since MergeSamFiles creates the index
                #if outBam.endswith(".bam") or outBam.endswith(".cram"):
                #    samtools.index(outBam)
                for bam in align_bams:
                    os.unlink(bam)

    def align_one_rg(self, inBam, refDb, outBam, rgid=None, options=None,
                         min_score_to_filter=None, threads=None, JVMmemory=None, invert_filter=False, should_index=True):
        """
            Performs an alignment of one read group in a bam file to a reference fasta file

            TODO: With the addition of a third aligner to viral-ngs, the functionality
            common to this method and to the comparable method in the Novoalign wrapper should
            be broken out as an "aligner" superclass, capable of aligning bam or fastq files with an arbitrary
            aligner, while preserving read groups. 
        """
        options = options or []

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
                samtools.view(['-b', '-r', rgid], inBam, tmp_bam)
                # special exit if this file is empty
                if samtools.count(tmp_bam) == 0:
                    log.warning("No reads present for RG %s in file: %s", rgid, inBam)
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

        # perform actual alignment

        # get the read group line to give to BWA
        readgroup_line = ""
        with open(headerFile) as inf:
            for line in inf:
                if line.startswith("@RG"):
                    readgroup_line = line

        assert len(readgroup_line) > 0

        #with util.file.tempfname('.aligned.bam') as tmp_bam_aligned:
        # rather than reheader the alignment bam file later so it has the readgroup information
        # from the original bam file, we'll pass the RG line to bwa to write out
        self.align_cmd(one_rg_inBam, refDb, outBam, options=options+['-R',
                 readgroup_line.rstrip("\r\n").replace('\t','\\t')],
                 threads=threads)

        return (rgid, outBam)
        # if there was more than one RG in the input, we had to create a temporary file with the one RG specified
        # and we can safely delete it this file
        # if there was only one RG in the input, we used it directly and should not delete it
        if removeInput:
            os.unlink(one_rg_inBam)

    def align_cmd(self, inReads, refDb, outAlign, options=None, threads=None):
        options = [] if not options else options

        threads = util.misc.sanitize_thread_count(threads)
        if '-t' not in options:
            options.extend(('-t', str(threads)))
        if '-2' not in options:
            options.append('-2')

        samtools = tools.samtools.SamtoolsTool()

        with util.file.tempfname('.aligned.sam') as aln_sam:
            fastq_pipe = samtools.bam2fq_pipe(inReads)
            options.extend((refDb, '-', '-o', aln_sam))
            self.execute(options, stdin=fastq_pipe)
            if fastq_pipe.poll():
                raise subprocess.CalledProcessError(fastq_pipe.returncode, "samtools.bam2fq_pipe() for {}".format(inReads))
            samtools.sort(aln_sam, outAlign, threads=threads)

        # cannot index sam files; only do so if a bam/cram is desired
        if (outAlign.endswith(".bam") or outAlign.endswith(".cram")):
            samtools.index(outAlign)
