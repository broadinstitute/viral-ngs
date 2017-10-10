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

    def lastal_get_hits(self, inBam, db, outList,
            max_gapless_alignments_per_position=1,
            min_length_for_initial_matches=5,
            max_length_for_initial_matches=50,
            max_initial_matches_per_position=100
        ):

        # convert BAM to paired FASTQ
        inFastq = util.file.mkstempfname('.all.fastq')
        tools.samtools.SamtoolsTool().bam2fq(inBam, inFastq)

        # run lastal and save output to lastal_out
        # -P 0 = use threads = core count
        # -N 1 = report at most one alignment per query sequence
        # -i 1G = perform work in batches of at most 1GB of query sequence at a time
        # -f tab = write output in tab format instead of maf format
        lastal_out = util.file.mkstempfname('.lastal')
        with open(lastal_out, 'wt') as outf:
            cmd = [self.install_and_get_path(),
                '-n', max_gapless_alignments_per_position,
                '-l', min_length_for_initial_matches,
                '-L', max_length_for_initial_matches,
                '-m', max_initial_matches_per_position,
                '-Q', '1', '-P', '0', '-N', '1', '-i', '1G', '-f', 'tab',
                db,
                inFastq,
            ]
            cmd = [str(x) for x in cmd]
            _log.debug(' '.join(cmd) + ' > ' + lastal_out)
            util.misc.run_and_save(cmd, outf=outf)

        # strip tab output to just query read ID names
        with open(outList, 'wt') as outf:
            with open(lastal_out, 'rt') as inf:
                for row in util.file.read_tabfile(lastal_out):
                    read_id = row[6]
                    if read_id.endswith('/1') or read_id.endswith('/2'):
                        read_id = read_id[:-2]
                    outf.write(read_id + '\n')
        os.unlink(lastal_out)


class MafSort(LastTools):
    """ wrapper for maf-sort subtool """
    subtool_name = 'maf-sort'
    subtool_name_on_broad = 'scripts/maf-sort.sh'


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

        # store the cwd because we will be changing it to the file destination
        cwd_before_lastdb = os.getcwd()

        # append the prefix given to files created by lastdb
        tool_cmd.append(outputFilePrefix)

        # append the input filepath
        tool_cmd.append(os.path.realpath(inputFasta))

        # lastdb writes files to the current working directory, so we need to set
        # it to the desired output location
        os.chdir(os.path.realpath(outputDirectory))

        # execute the lastdb command
        _log.debug(" ".join(tool_cmd))
        subprocess.check_call(tool_cmd)

        # restore cwd
        os.chdir(cwd_before_lastdb)


class MafConvert(LastTools):
    """ wrapper for maf-convert subtool """
    subtool_name = 'maf-convert'
    subtool_name_on_broad = 'scripts/maf-convert.py'
