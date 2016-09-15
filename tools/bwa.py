'''
    The BWA aligner.

'''

import logging
import os
import os.path
import subprocess

import tools
import tools.samtools
import tools.picard
import util.file
import util.misc

TOOL_NAME = 'bwa'
TOOL_VERSION = '0.7.15'

log = logging.getLogger(__name__)


class Bwa(tools.Tool):

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [tools.CondaPackage(TOOL_NAME, version=TOOL_VERSION)]
        tools.Tool.__init__(self, install_methods=install_methods)

    def version(self):
        return TOOL_VERSION

    def execute(self, command, args, stdout=None):    # pylint: disable=W0221
        tool_cmd = [self.install_and_get_path(), command] + args
        log.debug(' '.join(tool_cmd))
        if stdout:
            stdout = open(stdout, 'w')
        subprocess.check_call(tool_cmd, stdout=stdout)
        if stdout:
            stdout.close()

    def index(self, inFasta, algorithm=None):
        cmd = []
        if algorithm is not None:
            if algorithm not in ('is', 'bwtsw'):
                raise NameError(algorithm + " is not a recognized algorithm")
            cmd.extend(('-a', algorithm))
        cmd.append(inFasta)
        self.execute('index', cmd)


    def align_mem_bam(self, inBam, refDb, outBam, options=None, min_qual=30, threads=None, JVMmemory=None):
        options = options or []

        samtools = tools.samtools.SamtoolsTool()

        # fetch list of RGs
        rgs = list(samtools.getReadGroups(inBam).keys())

        if len(rgs) == 0:
            # Can't do this
            raise InvalidBamHeaderError("{} lacks read groups".format(inBam))

        elif len(rgs) == 1:
            # Only one RG, keep it simple
            self.align_mem_one_rg(inBam, refDb, outBam, options=options, min_qual=min_qual, threads=threads)

        else:
            # Multiple RGs, align one at a time and merge
            align_bams = []
            for rg in rgs:
                tmp_bam = util.file.mkstempfname('.{}.bam'.format(rg))
                self.align_mem_one_rg(
                    inBam,
                    refDb,
                    tmp_bam,
                    rgid=rg,
                    options=options,
                    min_qual=min_qual,
                    threads=threads
                )
                if os.path.getsize(tmp_bam) > 0:
                    align_bams.append(tmp_bam)

            # Merge BAMs, sort, and index
            tools.picard.MergeSamFilesTool().execute(
                align_bams,
                outBam,
                picardOptions=['SORT_ORDER=coordinate', 'USE_THREADING=true', 'CREATE_INDEX=true'],
                JVMmemory=JVMmemory
            )
            for bam in align_bams:
                os.unlink(bam)


    def align_mem_one_rg(self, inBam, refDb, outBam, rgid=None, options=None, min_qual=30, threads=None, JVMmemory=None):
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
            tmp_bam = util.file.mkstempfname('.onebam.bam')
            samtools.view(['-b', '-r', rgid], inBam, tmp_bam)
            # special exit if this file is empty
            if samtools.count(tmp_bam) == 0:
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
            os.unlink(tmp_bam)

        # perform actual alignment

        # get the read group line to give to BWA
        readgroup_line = ""
        with open(headerFile) as inf:
            for line in inf:
                if line.startswith("@RG"):
                    readgroup_line = line

        assert len(readgroup_line) > 0
        
        aln_bam_prefilter = util.file.mkstempfname('.prefiltered.bam')
        # rather than reheader the alignment bam file later so it has the readgroup information
        # from the original bam file, we'll pass the RG line to bwa to write out
        self.mem(one_rg_inBam, refDb, aln_bam_prefilter, opts=options+['-R', readgroup_line.rstrip("\n").rstrip("\r")], min_qual=min_qual, threads=threads)
        
        # if there was more than one RG in the input, we had to create a temporary file with the one RG specified
        # and we can safely delete it this file
        # if there was only one RG in the input, we used it directly and should not delete it
        if removeInput:
            os.unlink(one_rg_inBam)

        # @haydenm says: 
        # For some reason (particularly when the --sensitive option is on), bwa
        # doesn't listen to its '-T' flag and outputs alignments with score less
        # than the '-T 30' threshold. So filter these:
        if min_qual > 0:
            tmp_bam_aligned = util.file.mkstempfname('.aligned.bam')
            tools.samtools.SamtoolsTool().view(["-b", "-h", "-q", str(min_qual)], aln_bam_prefilter, tmp_bam_aligned)
            os.unlink(aln_bam_prefilter)
        else:
            shutil.move(aln_bam_prefilter, tmp_bam_aligned)

        # if the aligned bam file contains no reads after filtering
        # just create an empty file
        if tools.samtools.SamtoolsTool().count(tmp_bam_aligned) == 0:
             util.file.touch(outBam)
        else:
            # samtools reheader seems to segfault on some alignments created by bwa
            # so rather than reheader, BWA will write out the RG given to it via '-R'
            # reheadered_bam = util.file.mkstempfname('.reheadered.bam')
            # tools.samtools.SamtoolsTool().reheader(tmp_bam_aligned, headerFile, reheadered_bam)
            # os.unlink(tmp_bam_aligned)
            # os.unlink(headerFile)
            # os.system("samtools view -h {} > /Users/tomkinsc/Desktop/test_reheader.bam".format(reheadered_bam))

            # sort
            sorter = tools.picard.SortSamTool()
            sorter.execute(
                tmp_bam_aligned,
                outBam,
                sort_order='coordinate',
                picardOptions=['CREATE_INDEX=true', 'VALIDATION_STRINGENCY=SILENT'],
                JVMmemory=JVMmemory
            )
            #os.unlink(reheadered_bam)

    def mem(self, inReads, refDb, outAlign, opts=None, min_qual=None, threads=None):
        opts = [] if not opts else opts

        threads = threads or util.misc.available_cpu_count()
        samtools = tools.samtools.SamtoolsTool()
        fq1 = util.file.mkstempfname('.1.fastq')
        fq2 = util.file.mkstempfname('.2.fastq')
        aln_sam = util.file.mkstempfname('.sam')
        samtools.bam2fq(inReads, fq1, fq2)

        if '-t' not in opts:
            threads = threads or utils.misc.available_cpu_count()
            opts.extend( ('-t', str(threads)) )

        if '-T' not in opts:
            min_qual = min_qual or 30
            opts.extend( ('-T', str(min_qual)) )

        self.execute('mem', opts + [refDb, fq1, fq2], stdout=aln_sam)

        os.unlink(fq1)
        os.unlink(fq2)
        samtools.sort(aln_sam, outAlign)
        os.unlink(aln_sam)
        # cannot index sam files; only do so if a bam/cram is desired
        if outAlign.endswith(".bam") or outAlign.endswith(".cram"):
            samtools.index(outAlign)
