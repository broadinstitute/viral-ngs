'''
DIAMOND - blastx replacement for large database protein sequence queries
'''
from builtins import super
import itertools
import logging
import os
import os.path
import shlex
import shutil
import subprocess
import tools
import util.file

TOOL_VERSION = '0.9.10'

log = logging.getLogger(__name__)

@tools.skip_install_test(condition=tools.is_osx())
class Diamond(tools.Tool):

    SUBCOMMANDS = ['makedb', 'blastx', 'blastp', 'view']

    def __init__(self, install_methods=None):
        if not install_methods:
            install_methods = [
                tools.CondaPackage("diamond", version=TOOL_VERSION)
            ]
        super(Diamond, self).__init__(install_methods=install_methods)

    def version(self):
        return TOOL_VERSION

    def build(self, db, protein_fastas, options=None, option_string=None):
        '''Create a diamond database.

        Args:
          db: Diamond database file to create.
          protein_fastas: List of input fasta files to process.
        '''
        assert len(protein_fastas), ('Diamond requires input files to create a database.')
        options = options or {}
        temp_file = util.file.temp_catted_files(protein_fastas, prefix='diamond_', suffix='.fasta')
        with temp_file as input_fasta:
            options['--in'] = input_fasta
            options['--db'] = db

            self.execute('makedb', options=options, option_string=option_string)

    def view(self, diamond_alignment, output_file, output_format='tab', options=None, option_string=None):
        '''Perform translation between diamond output and blast tab/sam output.
        '''

        assert output_format in ('tab', 'sam'), 'Invalid diamond view format'
        options = options or {}
        options['--out'] = output_file
        options['--daa'] = diamond_alignment
        options['--outfmt'] = output_format
        self.execute('view', options=options, option_string=option_string)

    def execute(self, command, options=None, option_string=None, return_stdout=False):
        '''Run a diamond command

        Args:
          options: Dict of command line options to values. Set value to None
            for an option with no value.
          return_stdout: Whether to return stdout as well as in
            (exitcode, stdout).
        '''
        assert command in Diamond.SUBCOMMANDS, 'Diamond command is unknown'

        cmd = [self.install_and_get_path(), command]
        if options:
            # We need some way to allow empty options args like --log, hence
            # we filter out on 'x is None'.
            cmd.extend([str(x) for x in itertools.chain(*options.items()) if x is not None])
        if option_string:
            cmd.extend(shlex.split(option_string))
        log.debug("Calling {}: {}".format(command, " ".join(cmd)))
        subprocess.check_call(cmd)
