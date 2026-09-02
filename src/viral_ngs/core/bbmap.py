'''
    Tool wrapper for the BBMap aligner and related tools.
'''

import logging
import os
import os.path
import shutil
import subprocess

from . import file as util_file, misc as util_misc  # was: from viral_ngs import util
from . import samtools, picard  # was: from viral_ngs import tools
from . import Tool, PrexistingUnixCommand

TOOL_NAME = 'bbmap.sh'

_log = logging.getLogger(__name__)  # pylint: disable=invalid-name

class BBMapTool(Tool):
    '''Tool wrapper for the BBMap aligner and related tools.'''

    # BBTools' own heap autodetection (calcmem.sh) subtracts a fixed 500MB floor
    # from detected free RAM, so under container memory limits it can emit a
    # negative -Xmx and the JVM refuses to start. Callers that do not specify
    # memory should pass this instead of relying on autodetection. Modest on
    # purpose: clumpify's groups=auto bounds its own footprint by spilling to
    # disk, so a large heap buys little.
    memDefault = '2g'

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [PrexistingUnixCommand(shutil.which(TOOL_NAME), require_executability=True)]
        super(BBMapTool, self).__init__(install_methods=install_methods)

    def _get_tool_version(self):
        self.tool_version = subprocess.check_output([os.path.join(os.path.dirname(self.install_and_get_path()), 'bbversion.sh')]).decode('UTF-8').strip()

    def execute(self, tool, **kwargs):  # pylint: disable=arguments-differ
        tool_dir = os.path.dirname(self.install_and_get_path())
        tool_cmd = [os.path.join(tool_dir, tool)] + \
                   ['{}={}'.format('in' if arg=='in_' else arg,
                                   (val is True and 't') or (val is False and 'f') or val)
                    for arg, val in kwargs.items()]
        _log.debug('Running BBMap tool: %s', ' '.join(tool_cmd))
        subprocess.check_call(tool_cmd)

    def align(self, inBam, refFasta, outBam, min_qual=0, nodisk=True, JVMmemory=None, **kwargs):
        with samtools.SamtoolsTool().bam2fq_tmp(inBam) as (in1, in2), \
             util_file.tmp_dir('_bbmap_align') as t_dir:
            tmp_bam = os.path.join(t_dir, 'bbmap_out.bam')
            self.execute(tool='bbmap.sh', in1=in1, in2=in2, ref=refFasta, out=tmp_bam, nodisk=nodisk, **kwargs)
            
            # Samtools filter (optional)
            if min_qual:
                tmp_bam2 = os.path.join(tdir, 'bbmap.filtered.bam')
                cmd = [samtools.install_and_get_path(), 'view', '-b', '-S', '-1', '-q', str(min_qual), tmp_bam]
                _log.debug('%s > %s', ' '.join(cmd), tmp_bam2)
                with open(tmp_bam2, 'wb') as outf:
                    util_misc.run_and_save(cmd, outf=outf)
                os.unlink(tmp_bam)
                tmp_bam = tmp_bam2

            # Picard SortSam
            sorter = picard.SortSamTool()
            sorter.execute(
                tmp_bam,
                outBam,
                sort_order='coordinate',
                picardOptions=['CREATE_INDEX=true', 'VALIDATION_STRINGENCY=SILENT'],
                JVMmemory=JVMmemory
            )

    def bbnorm(self, inFastq, outFastq, tmpdir=None, target=None, k=None,
               passes=None, mindepth=0, threads=None, memory=None):
        """
        Run bbnorm for read normalization/deduplication.

        Args:
            inFastq: Input FASTQ file (interleaved for paired-end, or single-end)
            outFastq: Output FASTQ file (same format as input)
            tmpdir: Temporary directory for multipass processing
            target: Target normalization depth (default: bbnorm default of 100)
            k: Kmer length (default: bbnorm default of 31)
            passes: Number of passes (default: bbnorm default of 2)
            mindepth: Min kmer depth threshold (default: 0 to include all)
            threads: Number of threads (default: auto)
            memory: Java memory allocation string (e.g., "4g")

        Note: bbnorm auto-detects interleaved vs single-end format.
        """
        tool_dir = os.path.dirname(self.install_and_get_path())
        tool_cmd = [os.path.join(tool_dir, 'bbnorm.sh'), '-eoom']

        if memory:
            tool_cmd.append('-Xmx{}'.format(memory))

        # Build arguments
        tool_cmd.append('in={}'.format(inFastq))
        tool_cmd.append('out={}'.format(outFastq))
        tool_cmd.append('mindepth={}'.format(mindepth))

        if tmpdir is not None:
            tool_cmd.append('tmpdir={}'.format(tmpdir))
        if target is not None:
            tool_cmd.append('target={}'.format(target))
        if k is not None:
            tool_cmd.append('k={}'.format(k))
        if passes is not None:
            tool_cmd.append('passes={}'.format(passes))
        if threads is not None:
            tool_cmd.append('threads={}'.format(threads))

        _log.debug('Running bbnorm: %s', ' '.join(tool_cmd))
        subprocess.check_call(tool_cmd)

    def clumpify(self, inFastq, outFastq, inFastq2=None, outFastq2=None,
                 interleaved=None, tmpdir=None, subs=None, containment=False,
                 optical=False, dupedist=None, threads=None, memory=None):
        """
        Run clumpify in deduplication mode.

        Args:
            inFastq: Input FASTQ (read 1, or interleaved/single-end)
            outFastq: Output FASTQ corresponding to inFastq
            inFastq2: Optional read 2 input, for split (non-interleaved) pairs
            outFastq2: Optional read 2 output, required when inFastq2 is given
            interleaved: Force interleaved interpretation of a single input file
                         (clumpify auto-detects by default)
            tmpdir: Temporary directory for clumpify's intermediate group files
            subs: Max substitutions allowed between duplicates (clumpify default: 2)
            containment: Also treat a shorter sequence contained in a longer one
                         as a duplicate (clumpify default: off)
            optical: Restrict removal to optical duplicates only. Requires
                     Illumina-style read names carrying flowcell coordinates.
            dupedist: Max flowcell distance for optical duplicates (clumpify
                      default: 40). Only meaningful with optical=True.
            threads: Number of threads (default: auto)
            memory: Java memory allocation string (e.g., "4g")

        Note: clumpify's `passes` parameter is deliberately not exposed. Despite
        appearing in duplicate-removal recipes online, it controls the number of
        *error-correction* passes and has no effect on deduplication.

        Note: clumpify reorders its output into clumps. Callers must not rely on
        read order -- rmdup_clumpify_bam only harvests read IDs from the output
        and filters the original BAM.
        """
        tool_dir = os.path.dirname(self.install_and_get_path())
        tool_cmd = [os.path.join(tool_dir, 'clumpify.sh'), '-eoom']

        if memory:
            tool_cmd.append('-Xmx{}'.format(memory))

        # Build arguments
        if inFastq2 is not None:
            if outFastq2 is None:
                raise ValueError('outFastq2 is required when inFastq2 is given')
            tool_cmd.append('in1={}'.format(inFastq))
            tool_cmd.append('in2={}'.format(inFastq2))
            tool_cmd.append('out1={}'.format(outFastq))
            tool_cmd.append('out2={}'.format(outFastq2))
        else:
            tool_cmd.append('in={}'.format(inFastq))
            tool_cmd.append('out={}'.format(outFastq))
        tool_cmd.append('dedupe=t')

        if interleaved is not None:
            tool_cmd.append('interleaved={}'.format('t' if interleaved else 'f'))
        if tmpdir is not None:
            tool_cmd.append('usetmpdir=t')
            tool_cmd.append('tmpdir={}'.format(tmpdir))
        if subs is not None:
            tool_cmd.append('subs={}'.format(subs))
        if containment:
            tool_cmd.append('containment=t')
        if optical:
            tool_cmd.append('optical=t')
        if dupedist is not None:
            tool_cmd.append('dupedist={}'.format(dupedist))
        if threads is not None:
            tool_cmd.append('threads={}'.format(threads))

        _log.debug('Running clumpify: %s', ' '.join(tool_cmd))
        subprocess.check_call(tool_cmd)
