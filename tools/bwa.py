'''
    The BWA aligner.

'''

from collections import defaultdict
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

    def index(self, inFasta, output=None, algorithm=None):
        cmd = []
        if algorithm is not None:
            if algorithm not in ('is', 'bwtsw'):
                raise NameError(algorithm + " is not a recognized algorithm")
            cmd.extend(('-a', algorithm))
        if output:
            cmd.extend(('-p', output))
        cmd.append(inFasta)
        self.execute('index', cmd)

    def align_mem_bam(self, inBam, refDb, outBam, options=None,
                      min_score_to_filter=None, threads=None, JVMmemory=None):
        options = options or []

        samtools = tools.samtools.SamtoolsTool()

        # fetch list of RGs
        rgs = list(samtools.getReadGroups(inBam).keys())

        if len(rgs) == 0:
            # Can't do this
            raise InvalidBamHeaderError("{} lacks read groups".format(inBam))

        elif len(rgs) == 1:
            # Only one RG, keep it simple
            self.align_mem_one_rg(inBam, refDb, outBam, options=options,
                                  min_score_to_filter=min_score_to_filter,
                                  threads=threads)

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
                    min_score_to_filter=min_score_to_filter,
                    threads=threads
                )
                if os.path.getsize(tmp_bam) > 0:
                    align_bams.append(tmp_bam)
                else:
                    log.warning("No alignment output for RG %s in file %s against %s", rg, inBam, refDb)

            if len(align_bams) == 0:
                util.file.touch(outBam)
            else:
                # Merge BAMs, sort, and index
                tools.picard.MergeSamFilesTool().execute(
                    align_bams,
                    outBam,
                    picardOptions=['SORT_ORDER=coordinate', 'USE_THREADING=true', 'CREATE_INDEX=true'],
                    JVMmemory=JVMmemory
                )
                for bam in align_bams:
                    os.unlink(bam)


    def align_mem_one_rg(self, inBam, refDb, outBam, rgid=None, options=None,
                         min_score_to_filter=None, threads=None, JVMmemory=None):
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
            os.unlink(tmp_bam)

        # perform actual alignment

        # get the read group line to give to BWA
        readgroup_line = ""
        with open(headerFile) as inf:
            for line in inf:
                if line.startswith("@RG"):
                    readgroup_line = line

        assert len(readgroup_line) > 0

        tmp_bam_aligned = util.file.mkstempfname('.aligned.bam')
        # rather than reheader the alignment bam file later so it has the readgroup information
        # from the original bam file, we'll pass the RG line to bwa to write out
        self.mem(one_rg_inBam, refDb, tmp_bam_aligned, options=options+['-R',
                 readgroup_line.rstrip("\r\n")],
                 min_score_to_filter=min_score_to_filter, threads=threads)

        # if there was more than one RG in the input, we had to create a temporary file with the one RG specified
        # and we can safely delete it this file
        # if there was only one RG in the input, we used it directly and should not delete it
        if removeInput:
            os.unlink(one_rg_inBam)

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

    def mem(self, inReads, refDb, outAlign, options=None, min_score_to_filter=None,
            threads=None):
        options = [] if not options else options

        threads = threads or util.misc.available_cpu_count()
        samtools = tools.samtools.SamtoolsTool()
        fq1 = util.file.mkstempfname('.1.fastq')
        fq2 = util.file.mkstempfname('.2.fastq')
        aln_sam = util.file.mkstempfname('.sam')
        samtools.bam2fq(inReads, fq1, fq2)

        if '-t' not in options:
            options.extend(('-t', str(threads)))

        self.execute('mem', options + [refDb, fq1, fq2], stdout=aln_sam)
        os.unlink(fq1)
        os.unlink(fq2)

        if min_score_to_filter:
            # Filter reads in the alignment based on on their alignment score
            aln_sam_filtered = util.file.mkstempfname('.sam')
            self.filter_sam_on_alignment_score(aln_sam, aln_sam_filtered,
                                               min_score_to_filter, options)
        else:
            aln_sam_filtered = aln_sam

        samtools.sort(aln_sam_filtered, outAlign, threads=threads)
        os.unlink(aln_sam)
        # cannot index sam files; only do so if a bam/cram is desired
        if outAlign.endswith(".bam") or outAlign.endswith(".cram"):
            samtools.index(outAlign)

    def filter_sam_on_alignment_score(self, in_sam, out_sam, min_score_to_filter,
                                      bwa_options):
        """Filter reads in an alignment based on their alignment score.

        This is useful because bwa mem ignores its -T option on paired-end
        data (https://github.com/lh3/bwa/issues/133); -T normally specifies
        the threshold on score below which alignments are not output.

        This sums the alignment score across all alignments for each
        query name in the aligned SAM. This has the effect of summing
        the score across both reads in a pair (their query name is the same)
        as well as supplementary alignments (SAM flag 2048) and secondary
        alignments (SAM flag 256). The use of a summed score is reasonable for
        supplementary alignments, which bwa mem outputs. bwa mem normally
        does not output secondary alignments, but if it does then the
        use of a summed score might not be reasonable.
        """
        if '-a' in bwa_options:
            log.warning(("'bwa mem -a' will output secondary alignments, "
                         "and the filter on alignment score will use "
                         "a score that is summed across the primary and "
                         "secondary alignments for each read; this might "
                         "not be desired"))

        # Determine an alignment score for each query name
        qname_alignment_scores = defaultdict(int)
        with open(in_sam) as in_sam_f:
            for line in in_sam_f:
                line = line.rstrip()
                if line.startswith('@'):
                    # Skip headers
                    continue
                ls = line.split('\t')

                # bwa's output should have optional fields for all
                # alignments
                assert len(ls) >= 12

                qname = ls[0]
                aln_score = None
                for opt_field in ls[11:]:
                    opt_field_tag, opt_field_type, opt_field_val = opt_field.split(':')
                    if opt_field_tag == 'AS':
                        # The alignment score output by bwa should be a
                        # signed integer
                        assert opt_field_type == 'i'
                        aln_score = int(opt_field_val)
                        break
                if aln_score is None:
                    raise Exception(("Unknown alignment score for query "
                                    "name %s") % qname)
                qname_alignment_scores[qname] += aln_score

        # Write a SAM file that only includes reads whose query has a
        # sufficient alignment score summed across its alignments
        with open(out_sam, 'w') as out_sam_f:
            with open(in_sam) as in_sam_f:
                for line in in_sam_f:
                    line = line.rstrip()
                    if line.startswith('@'):
                        # Rewrite all headers
                        out_sam_f.write(line + '\n')
                        continue
                    ls = line.split('\t')

                    qname = ls[0]
                    if qname_alignment_scores[qname] >= min_score_to_filter:
                        # Write this query name
                        out_sam_f.write(line + '\n')
