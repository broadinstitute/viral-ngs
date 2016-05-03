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
import util.misc

URL = 'https://github.com/bbuchfink/diamond/archive/b576a5c03177603554f4627ae367f7bcbc6b8dcb.zip'
TOOL_VERSION = '0.7.9'
CONDA_VERSION = tools.CondaPackageVersion('0.7.10', 'boost1.60_1')
DIAMOND_COMMIT_DIR = 'diamond-b576a5c03177603554f4627ae367f7bcbc6b8dcb'
DIAMOND_DIR = 'diamond-{}'.format(TOOL_VERSION)

log = logging.getLogger(__name__)


class Diamond(tools.Tool):

    SUBCOMMANDS = ['makedb', 'blastx', 'blastp', 'view']

    def __init__(self, install_methods=None):
        if not install_methods:
            install_methods = [
                tools.CondaPackage("diamond", version=CONDA_VERSION),
                DownloadAndBuildDiamond(URL, os.path.join(DIAMOND_DIR, 'bin', 'diamond'))]
        super().__init__(install_methods=install_methods)

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

            return self.execute('makedb', options=options, option_string=option_string)

    def blastx(self, db, query_files, diamond_alignment, options=None, option_string=None):
        '''Perform a blastx-like search from query file to database.

        Args:
          db: Diamond database file.
          query_files: List of input fastq files.
          diamond_alignment: Diamond alignment output file. Must end in .daa
        '''
        assert diamond_alignment.endswith('.daa'), 'Output must end in .daa'
        options = options or {}
        temp_file = util.file.temp_catted_files(query_files, prefix='diamond_', suffix='.fasta')
        with temp_file as query:
            options['--db'] = db
            options['--query'] = query
            options['--daa'] = diamond_alignment
            return self.execute('blastx', options=options, option_string=option_string)

    def view(self, diamond_alignment, output_file, output_format='tab', options=None, option_string=None):
        '''Perform translation between diamond output and blast tab/sam output.
        '''

        assert output_format in ('tab', 'sam'), 'Invalid diamond view format'
        options = options or {}
        options['--out'] = output_file
        options['--daa'] = diamond_alignment
        options['--outfmt'] = output_format
        return self.execute('view', options=options, option_string=option_string)

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
        return util.misc.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)


class DownloadAndBuildDiamond(tools.DownloadPackage):

    # We need to refactor to have a generic cmake installer.
    def post_download(self):
        diamond_dir = os.path.join(self.destination_dir, DIAMOND_DIR)
        # We should rather have a way to rename self.download_file in
        # DownloadPackage generically.
        if not os.path.exists(diamond_dir):
            shutil.move(os.path.join(self.destination_dir, DIAMOND_COMMIT_DIR), diamond_dir)
        build_dir = os.path.join(diamond_dir, 'src')
        #util.file.mkdir_p(build_dir)
        env = os.environ.copy()
        # The default travis gcc version is 4.6, which is too old to build
        # diamond properly.
        if os.environ.get('TRAVIS') == 'true':
            env['CC'] = 'gcc-4.9'
            env['CXX'] = 'g++-4.9'
        #util.misc.run_and_print(['cmake', '..'], env=env, cwd=build_dir)
        util.misc.run_and_print(['make'], env=env, cwd=build_dir)
