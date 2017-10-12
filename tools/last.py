"Tools in the 'last' suite."

# built-ins
import os
import logging
import subprocess

# within this module
import util.file
import util.misc
import tools
import tools.samtools

_log = logging.getLogger(__name__)

TOOL_NAME = "last"
TOOL_VERSION = "876"


class LastTools(tools.Tool):
    """
    "Abstract" base class for tools in the 'last' suite.
    Subclasses must define class members subtool_name #and subtool_name_on_broad.
    """

    def __init__(self, install_methods=None):
        self.subtool_name = self.subtool_name if hasattr(self, "subtool_name") else None
        self.subtool_name_on_broad = self.subtool_name_on_broad if hasattr(self, "subtool_name_on_broad") else None
        if install_methods is None:
            install_methods = []
            install_methods.append(tools.CondaPackage(TOOL_NAME, executable=self.subtool_name, version=TOOL_VERSION))
        tools.Tool.__init__(self, install_methods=install_methods)


class Lastal(LastTools):
    """ wrapper for lastal subtool """
    subtool_name = 'lastal'
    subtool_name_on_broad = 'lastal'

    def get_hits(self, inBam, db,
            max_gapless_alignments_per_position=1,
            min_length_for_initial_matches=5,
            max_length_for_initial_matches=50,
            max_initial_matches_per_position=100
        ):

            # convert BAM to interleaved FASTQ
            fastq_pipe = tools.samtools.SamtoolsTool().bam2fq_pipe(inBam)

            # run lastal and emit list of read IDs
            # -P 0 = use threads = core count
            # -N 1 = report at most one alignment per query sequence
            # -i 1G = perform work in batches of at most 1GB of query sequence at a time
            # -f tab = write output in tab format instead of maf format
            cmd = [self.install_and_get_path(),
                '-n', max_gapless_alignments_per_position,
                '-l', min_length_for_initial_matches,
                '-L', max_length_for_initial_matches,
                '-m', max_initial_matches_per_position,
                '-Q', '1', '-P', '0', '-N', '1', '-i', '1G', '-f', 'tab',
                db,
            ]
            cmd = [str(x) for x in cmd]
            _log.debug('| ' + ' '.join(cmd) + ' |')
            lastal_pipe = subprocess.Popen(cmd, stdin=fastq_pipe, stdout=subprocess.PIPE)

            # strip tab output to just query read ID names and emit
            for line in lastal_pipe.stdout:
                line = line.decode('UTF-8').rstrip('\n\r')
                if not line.startswith('#'):
                    read_id = line.split('\t')[6]
                    if read_id.endswith('/1') or read_id.endswith('/2'):
                        read_id = read_id[:-2]
                    yield read_id


class Lastdb(LastTools):
    """ wrapper for lastdb subtool """
    subtool_name = 'lastdb'
    subtool_name_on_broad = 'lastdb'

    def build_database(self, fasta_files, database_prefix_path): # pylint: disable=W0221
        output_file_prefix = os.path.basename(database_prefix_path)
        output_directory = os.path.dirname(database_prefix_path)


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

        self.execute(input_fasta, output_directory, output_file_prefix)    

        return database_prefix_path


    def execute(self, inputFasta, outputDirectory, outputFilePrefix):    # pylint: disable=W0221
        # get the path to the binary
        tool_cmd = [self.install_and_get_path()]

        # if the output directory (and its parents) do not exist, create them
        if not os.path.exists(outputDirectory):
            os.makedirs(outputDirectory)

        # append the prefix given to files created by lastdb
        tool_cmd.append(outputFilePrefix)

        # append the input filepath
        tool_cmd.append(os.path.realpath(inputFasta))

        # execute the lastdb command
        # lastdb writes files to the current working directory, so we need to set
        # it to the desired output location
        with util.file.pushd_popd(os.path.realpath(outputDirectory)):
            _log.debug(" ".join(tool_cmd))
            subprocess.check_call(tool_cmd)


