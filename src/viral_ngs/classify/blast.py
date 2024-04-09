"Tools in the blast+ suite."

import logging
import os
import shutil
import subprocess

import tools
import tools.samtools
import util.misc

TOOL_NAME = "blastn"
#Creating task.log
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("blast_py.log"),  
        logging.StreamHandler() 
    ]
)
_log = logging.getLogger(__name__)

class BlastTools(tools.Tool):
    """'Abstract' base class for tools in the blast+ suite.
       Subclasses must define class member subtool_name."""

    def __init__(self, install_methods=None):
        unwanted = [
            'blast_formatter', 'blastdb_aliastool', 'blastdbcheck', 'blastdbcmd', 'convert2blastmask', 'deltablast',
            'legacy_blast.pl', 'makembindex', 'makeprofiledb', 'psiblast', 'rpsblast', 'rpstblastn', 'segmasker',
            'tblastn', 'tblastx', 'update_blastdb.pl', 'windowmasker'
        ]
        self.subtool_name = self.subtool_name if hasattr(self, "subtool_name") else "blastn"
        if install_methods is None:
            install_methods = [tools.PrexistingUnixCommand(shutil.which(self.subtool_name), require_executability=False)]
        super(BlastTools, self).__init__(install_methods=install_methods)

    def execute(self, *args):
        cmd = [self.install_and_get_path()]
        cmd.extend(args)
        util.misc.run_and_print(cmd, buffered=True, check=True)


class BlastnTool(BlastTools):
    """ Tool wrapper for blastn """
    subtool_name = 'blastn'

    def get_hits_pipe(self, inPipe, db, threads=None, task=None, outfmt='6', max_target_seqs=1, output_type="read_id"):
        _log.debug("Executing get_hits_pipe function.")
        #Validates output type run during this function
        _log.debug(f"get_hits_pipe called with outfmt: {outfmt}")
        #toggle between extracting read IDs only or full blast output (all lines)
        if output_type not in ['read_id', 'full_line']:
            _log.warning(f"Invalid output_type '{output_type}' specified. Defaulting to 'read_id'.")
            output_type = 'read_id'
        # run blastn and emit list of read IDs
        threads = util.misc.sanitize_thread_count(threads)
        cmd = [self.install_and_get_path(),
            '-db', db,
            '-word_size', 16,
            '-num_threads', threads,
            '-evalue', '1e-6',
            '-outfmt', str(outfmt),
            '-max_target_seqs', str(max_target_seqs),
            '-task', str(task) if task else 'blastn',
        ]
        cmd = [str(x) for x in cmd]
        #Log BLAST command executed
        _log.debug('Running blastn command: {}'.format(' '.join(cmd)))
        _log.debug('| ' + ' '.join(cmd) + ' |')
        blast_pipe = subprocess.Popen(cmd, stdin=inPipe, stdout=subprocess.PIPE)
        output, error = blast_pipe.communicate()

        #Display error message if BLAST failed
        if blast_pipe.returncode!= 0:
            _log.error('Error running blastn command: {}'.format(error))
            raise subprocess.CalledProcessError(blast_pipe.returncode, cmd)
        
        # If read_id is defined, strip tab output to just query read ID names and emit (default)
        if output_type == 'read_id':
            last_read_id = None
            for line in output.decode('UTF-8').splitlines():
                line = line.strip()
                read_id = line.split('\t')[0]
            #only emit if it is not a duplicate of the previous read ID
                if read_id != last_read_id:
                    last_read_id = read_id
                    yield read_id
        #Else, writes entire line of BLAST output
        elif output_type == 'full_line':
            for line in output.decode('UTF-8').splitlines():
                yield line 
        #Display on CMD if BLAST fails
        if blast_pipe.returncode!= 0:
            _log.error('Error running blastn command: {}'.format(error))
            raise subprocess.CalledProcessError(blast_pipe.returncode, cmd)
        #Logging configuration written to blast_py.log if BLAST passes/fails
        if blast_pipe.returncode == 0:
            _log.info("Blastn process completed succesfully.")
        else:
            _log.error("Blastn process failed with exit code: %s", blast_pipe.returncode)
            raise subprocess.CalledProcessError(blast_pipe.returncode, cmd)
       
    def get_hits_bam(self, inBam, db, threads=None):
        return self.get_hits_pipe(
            tools.samtools.SamtoolsTool().bam2fa_pipe(inBam),
            db,
            threads=threads)

    def get_hits_fasta(self, inFasta, db, threads=None, task=None, outfmt='6', max_target_seqs=1, output_type='read_id'):
        _log.debug("Executing get_hits_fasta function.")
        _log.debug(f"get_hits_fasta called with outfmt: {outfmt}")
        with open(inFasta, 'rt') as inf:
            for hit in self.get_hits_pipe(inf, db, threads=threads, task=None, outfmt=outfmt, max_target_seqs=1, output_type=output_type):
                yield hit


class MakeblastdbTool(BlastTools):
    """ Tool wrapper for makeblastdb """
    subtool_name = 'makeblastdb'

    def build_database(self, fasta_files, database_prefix_path):
        """ builds a srprism database """

        input_fasta = ""

        # we can pass in a string containing a fasta file path
        # or a list of strings
        if 'basestring' not in globals():
           basestring = str
        if isinstance(fasta_files, basestring):
            fasta_files = [fasta_files]
        elif isinstance(fasta_files, list):
            pass
        else:
            raise TypeError("fasta_files was not a single fasta file, nor a list of fasta files") # or something along that line

        # if more than one fasta file is specified, join them
        # otherwise if only one is specified, just use it
        if len(fasta_files) > 1:
            input_fasta = util.file.mkstempfname("fasta")
            util.file.cat(input_fasta, fasta_files)
        elif len(fasta_files) == 1:
            input_fasta = fasta_files[0]
        else:
            raise IOError("No fasta file provided")

        args = ['-dbtype', 'nucl', '-in', input_fasta, '-out', database_prefix_path]
        self.execute(*args)

        return database_prefix_path
