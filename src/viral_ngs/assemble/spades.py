'''
    Tool wrapper for SPAdes, St. Petersburg Assembler ( http://cab.spbu.ru/software/spades/ )
'''

import logging
import os
import os.path
import subprocess
import shutil
import random
import shlex
import tempfile

import Bio.SeqIO

import tools
import tools.samtools
import tools.picard
import util.file
import util.misc

TOOL_NAME = 'spades.py'

log = logging.getLogger(__name__)

class SpadesTool(tools.Tool):
    '''Tool wrapper for SPAdes tool (St. Petersburg Assembler)'''

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [tools.PrexistingUnixCommand(shutil.which(TOOL_NAME), require_executability=True)]
        super(SpadesTool, self).__init__(install_methods=install_methods)

    def version(self):
        if self.tool_version is None:
            self._get_tool_version()
        return self.tool_version

    def _get_tool_version(self):
        self.tool_version = subprocess.check_output([self.install_and_get_path(), '--version']).decode('UTF-8').strip().split()[1]

    def execute(self, args):    # pylint: disable=W0221
        tool_cmd = [self.install_and_get_path()] + list(map(str, args))
        log.debug(' '.join(tool_cmd))
        subprocess.check_call(tool_cmd)

    def assemble(self, reads_fwd, reads_bwd, contigs_out, reads_unpaired=None, contigs_trusted=None,
                 contigs_untrusted=None, kmer_sizes=None, always_succeed=True,
                 filter_contigs=False, min_contig_len=0, mem_limit_gb=8, threads=None, spades_opts=''):
        '''Assemble contigs from RNA-seq reads and (optionally) pre-existing contigs.

        Inputs:
            Required:
              reads_fwd, reads_bwd (fasta/q): paired reads
            Optional:
              reads_unpaired (fasta/q): optionally, additional unpaired reads
              contigs_trusted (fasta/q): optionally, already-assembled contigs of high quality
              contigs_untrusted (fasta/q): optionally, already-assembled contigs of average quality
        Params:
            kmer_sizes: if given, use these comma-separated ascending kmer sizes--all values must be odd and < 127.  kmer size of 0 or None means use size auto-selected
              by SPAdes based on read length.
            always_succeed: if True, if spades fails with an error, pretend it just produced no contigs
            filter_contigs: if True, outputs only "long and reliable transcripts with rather high expression" (per SPAdes docs)
            min_contig_len: drop contigs shorter than this many bp
            mem_limit_gb: max memory to use, in gigabytes
            threads: number of threads to use
            spades_opts: additional options to pass to spades
        Outputs:
            contigs_out: assembled contigs in fasta format.  Note that, since we use the
                RNA-seq assembly mode, for some genome regions we may get several contigs
                representing alternative transcripts.  Fasta record name of each contig indicates
                its length, coverage, and the group of alternative transcripts to which it belongs.
                See details at 
                http://cab.spbu.ru/files/release3.11.1/rnaspades_manual.html#sec2.4 .
        '''

        threads = util.misc.sanitize_thread_count(threads)

        util.file.make_empty(contigs_out)
        contigs_cumul_count = 0

        if ((reads_fwd and reads_bwd and os.path.getsize(reads_fwd) > 0 and os.path.getsize(reads_bwd) > 0) or
            (reads_unpaired and os.path.getsize(reads_unpaired) > 0)):

            with util.file.tmp_dir('_spades') as spades_dir:
                log.debug('spades_dir=' + spades_dir)
                args = []
                if reads_fwd and reads_bwd and os.path.getsize(reads_fwd) > 0 and os.path.getsize(reads_bwd) > 0:
                    args += ['-1', reads_fwd, '-2', reads_bwd ]
                if reads_unpaired and os.path.getsize(reads_unpaired) > 0:
                    args += [ '--s1', reads_unpaired ]
                if contigs_trusted: args += [ '--trusted-contigs', contigs_trusted ]
                if contigs_untrusted: args += [ '--untrusted-contigs', contigs_untrusted ]
                if kmer_sizes: args += [ '-k', kmer_sizes ]
                if spades_opts: args += shlex.split(spades_opts)
                args += [ '-m' + str(mem_limit_gb), '-t', str(threads), '-o', spades_dir ]

                transcripts_fname = os.path.join(spades_dir, ('hard_filtered_' if filter_contigs else '') + 'transcripts.fasta')

                try:
                    self.execute(args=args)
                except Exception as e:
                    if always_succeed:
                        log.warning('SPAdes failed: {}'.format(e))
                        util.file.make_empty(transcripts_fname)
                    else:
                        raise

                # work around the bug that spades may succeed yet not create the transcripts.fasta file
                if not os.path.isfile(transcripts_fname):
                    msg = 'SPAdes failed to make transcripts.fasta'
                    if always_succeed:
                        log.warning(msg)
                        util.file.make_empty(transcripts_fname)
                    else:
                        raise RuntimeError(msg)

                if min_contig_len:
                    transcripts = Bio.SeqIO.parse(transcripts_fname, 'fasta')
                    transcripts_sans_short = [r for r in transcripts if len(r.seq) >= min_contig_len]
                    transcripts_fname = os.path.join(spades_dir, 'transcripts_over_{}bp.fasta'.format(min_contig_len))
                    Bio.SeqIO.write(transcripts_sans_short, transcripts_fname, 'fasta')

                contigs_cumul = os.path.join(spades_dir, 'contigs_cumul.{}.fasta'.format(contigs_cumul_count))
                contigs_cumul_count += 1

                util.file.concat(inputFilePaths=(contigs_out, transcripts_fname), outputFilePath=contigs_cumul, append=True)
                shutil.copyfile(contigs_cumul, contigs_out)

            # end: with util.file.tmp_dir('_spades') as spades_dir
        # if input non-empty
    # end: def assemble(self, reads_fwd, reads_bwd, contigs_out, reads_unpaired=None, contigs_trusted=None, ...)
# end: class SpadesTool(tools.Tool)


