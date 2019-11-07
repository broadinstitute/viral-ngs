'''
    Tool wrapper for the BBMap aligner and related tools.
'''

import logging
import os
import os.path
import subprocess

import util.file
import tools
import tools.samtools
import tools.picard

TOOL_NAME = 'bbmap'
TOOL_VERSION = '38.71'

_log = logging.getLogger(__name__)  # pylint: disable=invalid-name

class BBMapTool(tools.Tool):
    '''Tool wrapper for the BBMap aligner and related tools.'''

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [tools.CondaPackage(TOOL_NAME, version=TOOL_VERSION, executable='bbmap.sh')]
        tools.Tool.__init__(self, install_methods=install_methods)

    def version(self):
        return TOOL_VERSION

    def execute(self, tool, JVMmemory=None, **kwargs):  # pylint: disable=arguments-differ
        tool_dir = os.path.dirname(self.install_and_get_path())
        tool_cmd = [os.path.join(tool_dir, tool)] + \
                   ['{}={}'.format('in' if arg=='in_' else arg,
                                   (val is True and 't') or (val is False and 'f') or val)
                    for arg, val in kwargs.items()]
        if JVMmemory:
            tool_cmd.append('-Xmx'+JVMmemory)
        _log.debug('Running BBMap tool: %s', ' '.join(tool_cmd))
        subprocess.check_call(tool_cmd)

    def align(self, in_bam, ref_fasta, out_bam, min_qual=0, nodisk=True, JVMmemory=None, **kwargs):
        with tools.samtools.SamtoolsTool().bam2fq_tmp(in_bam) as (in1, in2), \
             util.file.tmp_dir('_bbmap_align') as t_dir:
            tmp_bam = os.path.join(t_dir, 'bbmap_out.bam')
            self.execute(tool='bbmap.sh', in1=in1, in2=in2, ref=ref_fasta, out=tmp_bam, nodisk=nodisk, JVMmemory=JVMmemory, **kwargs)
            
            # Samtools filter (optional)
            if min_qual:
                samtools = tools.samtools.SamtoolsTool()
                tmp_bam2 = os.path.join(tdir, 'bbmap.filtered.bam')
                samtools.view(['-b', '-S', '-1', '-q', str(min_qual)], tmp_bam, tmp_bam2)
                os.unlink(tmp_bam)
                tmp_bam = tmp_bam2

            # Picard SortSam
            sorter = tools.picard.SortSamTool()
            sorter.execute(
                tmp_bam,
                out_bam,
                sort_order='coordinate',
                picardOptions=['CREATE_INDEX=true', 'VALIDATION_STRINGENCY=SILENT'],
                JVMmemory=JVMmemory
            )

    def dedup_clumpify(self, in_bam, out_bam, optical=False, subs=3, passes=4, dupedist=40, kmer_size=31, spany=False, adjacent=False, treat_as_unpaired=False, containment=True, JVMmemory=None, **kwargs):
        '''
            clumpify-based deduplication
            see:
                https://www.biostars.org/p/225338/
                https://www.biostars.org/p/225338/#230178
            and also:
                https://www.biostars.org/p/229842/#229940
                https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/clumpify-guide/

            From clumpify.sh usage:
                optical=False   If true, *only* mark or remove *optical* duplicates.
                                Optical duplicate removal is limited by xy-position, 
                                and is intended for single-end sequencing; 
                                especially useful for NextSeq data.
                dedupe=True     Remove duplicate reads.  For pairs, both must match.
                dupedist=40     (dist) Max distance to consider for optical duplicates.
                                Higher removes more duplicates but is more likely to
                                remove PCR rather than optical duplicates.
                                This is platform-specific; recommendations:
                                   NextSeq      40  (and spany=t)
                                   HiSeq 1T     40
                                   HiSeq 2500   40
                                   HiSeq 3k/4k  2500
                                   Novaseq      12000
                k=31              Use kmers of this length (1-31).  Shorter kmers may
                                  increase compression, but 31 is recommended for error
                                  correction.
                containment=True  Allow containments (where one sequence is shorter).
        '''
        unpair = treat_as_unpaired
        repair = treat_as_unpaired

        with tools.samtools.SamtoolsTool().bam2fq_tmp(in_bam) as (in1, in2), \
            util.file.tmp_dir('_bbmap_clumpify') as t_dir:

            # We may want to merge overlapping paired reads via BBMerge first; per clumpify docs:
            #   Clumpify supports paired reads, in which case it will clump based on read 1 only.
            #   However, it's much more effective to treat reads as unpaired. For example, merge 
            #   the reads with BBMerge, then concatenate the merged reads with the unmerged pairs,
            #   and clump them all together as unpaired.

            self.execute(tool='clumpify.sh', 
                            in1=in1, in2=in2,
                            out=out_bam,
                            dedupe=True,
                            subs=subs,
                            passes=passes,
                            dupedist=dupedist,
                            k=kmer_size,
                            optical=optical,
                            spany=spany,
                            adjacent=adjacent,
                            usetmpdir=True,
                            tmpdir=t_dir,
                            # if reads should be treated as unpaired, both 'unpair','repair' should be set to True
                            unpair=treat_as_unpaired,
                            repair=treat_as_unpaired,
                            containment=containment,
                            JVMmemory=JVMmemory,
                            **kwargs)
