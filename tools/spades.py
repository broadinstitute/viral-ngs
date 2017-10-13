'''
    SPAdes - St. Petersburg Assembler

'''

import logging
import os
import os.path
import subprocess
import shutil
import random
import shlex
import tempfile

import tools
import tools.samtools
import tools.picard
import util.file
import util.misc

TOOL_NAME = 'spades'
TOOL_VERSION = '3.11.1'

log = logging.getLogger(__name__)

class SpadesTool(tools.Tool):
    '''Tool wrapper for SPAdes tool (St. Petersburg Assembler)'''

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [tools.CondaPackage(TOOL_NAME, version=TOOL_VERSION, executable='spades.py',
                                                  verifycmd='spades.py --version')]
        tools.Tool.__init__(self, install_methods=install_methods)

    def version(self):
        return TOOL_VERSION

    def execute(self, args):    # pylint: disable=W0221
        tool_cmd = [self.install_and_get_path()] + args
        log.debug(' '.join(tool_cmd))
        subprocess.check_call(tool_cmd)

    def assemble(self, reads_fwd, reads_bwd, contigs_out, reads_unpaired=None, contigs_trusted=None,
                 contigs_untrusted=None, mem_limit_gb=4, threads=0, spades_opts=''):
        '''Assemble contigs from reads and (optionally) pre-existing contigs.

        Inputs:
            reads_fwd, reads_bwd - paired reads in fasta format
            reads_unpaired - optionally, additional unpaired reads in fasta format
            contigs_trusted - optionally, already-assembled contigs of high quality
            contigs_untrusted - optionally, already-assembled contigs of average quality
        Params:
            mem_limit_gb - max memory to use, in gigabytes
            threads - number of threads to use (0 means use all available CPUs)
            spades_opts - additional options to pass to spades
        Outputs:
            contigs_out - assembled contigs in fasta format.  Note that, since we use the
                RNA-seq assembly mode, for some genome regions we may get several contigs
                representing alternative transcripts.  Fasta record name of each contig indicates
                its length, coverage, and the group of alternative transcripts to which it belongs.
                See details at 
                http://cab.spbu.ru/files/release3.10.1/rnaspades_manual.html#sec2.2 .
        '''

        if not threads: threads = util.misc.available_cpu_count()

        if (reads_fwd and reads_bwd
            and os.path.getsize(reads_fwd) > 0 and os.path.getsize(reads_bwd) > 0
            or reads_unpaired and os.path.getsize(reads_unpaired) > 0):

            with util.file.tmp_dir('_spades') as spades_dir:
                log.debug('spades_dir=' + spades_dir)
                args = []
                if reads_fwd and reads_bwd and os.path.getsize(reads_fwd) > 0 and os.path.getsize(reads_bwd) > 0:
                    args += ['-1', reads_fwd, '-2', reads_bwd ]
                if reads_unpaired and os.path.getsize(reads_unpaired) > 0:
                    args += [ '--s1', reads_unpaired ]
                if contigs_trusted: args += [ '--trusted-contigs', contigs_trusted ]
                if contigs_untrusted: args += [ '--untrusted-contigs', contigs_untrusted ]
                if spades_opts: args += shlex.split(spades_opts)
                args += [ '--rna', '-m' + str(mem_limit_gb), '-t', str(threads), '-o', spades_dir ]

                self.execute(args=args)

                shutil.copyfile(src=os.path.join(spades_dir, 'transcripts.fasta'), dst=contigs_out)

        else:
            # spades crashes on empty input, so just return empty output
            util.file.make_empty(contigs_out)
            return

