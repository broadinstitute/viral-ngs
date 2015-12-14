'''
KRAKEN metagenomics classifier
'''
import itertools
import logging
import os
import os.path
import shlex
import shutil
import stat
import tools
import util.file
import util.misc
from builtins import super

URL = 'https://github.com/yesimon/kraken/archive/75154106773b41b1d0e55b3274178134eb14723d.zip'
TOOL_VERSION = '0.10.5-beta'
KRAKEN_COMMIT_DIR = 'kraken-75154106773b41b1d0e55b3274178134eb14723d'
KRAKEN_DIR = 'kraken-{}'.format(TOOL_VERSION)

JELLYFISH_URL = 'https://github.com/gmarcais/Jellyfish/archive/43fc99e4d44d11f115dc6741ff705cf7e113f251.zip'
JELLYFISH_VERSION = '1.1.11'
JELLYFISH_COMMIT_DIR = 'Jellyfish-43fc99e4d44d11f115dc6741ff705cf7e113f251'
JELLYFISH_DIR = 'jellyfish-{}'.format(JELLYFISH_VERSION)

YAGGO_URL = 'https://github.com/gmarcais/yaggo/releases/download/v1.5.9/yaggo'
YAGGO_VERSION = '1.5.9'

log = logging.getLogger(__name__)


class Yaggo(tools.Tool):

    def __init__(self, install_methods=None):
        if not install_methods:
            install_methods = [DownloadAndInstallYaggo(YAGGO_URL, 'yaggo')]
        super().__init__(install_methods=install_methods)


class DownloadAndInstallYaggo(tools.DownloadPackage):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.verifycmd = '{}/yaggo -v > /dev/null 2>& 1'.format(util.file.get_build_path())

    def post_download(self):
        yaggo_path = os.path.join(self.destination_dir, 'yaggo')
        os.chmod(yaggo_path, 0o755)


class Jellyfish(tools.Tool):

    def __init__(self, install_methods=None):
        if not install_methods:
            install_methods = [
                DownloadAndInstallJellyfish(
                    JELLYFISH_URL, os.path.join(JELLYFISH_DIR, 'bin', 'jellyfish'))
            ]
        super().__init__(install_methods=install_methods)


class DownloadAndInstallJellyfish(tools.DownloadPackage):

    def post_download(self):
        yaggo_path = Yaggo().install_and_get_path()
        env = os.environ.copy()
        env['PATH'] = '{}:{}'.format(os.path.dirname(yaggo_path), env['PATH'])
        jellyfish_dir = os.path.join(self.destination_dir, JELLYFISH_DIR)

        shutil.move(os.path.join(self.destination_dir, JELLYFISH_COMMIT_DIR), jellyfish_dir)

        install_dir = os.path.join(jellyfish_dir, 'local')
        util.file.replace_in_file(
            os.path.join(jellyfish_dir, 'Makefile.am'), 'AM_CXXFLAGS = -g -O3',
            'AM_CXXFLAGS = -g -O3 -Wno-maybe-uninitialized')
        util.misc.run_and_print(['autoreconf', '-i'], cwd=jellyfish_dir, env=env)
        util.misc.run_and_print(['./configure', '--prefix={}'.format(install_dir)], cwd=jellyfish_dir, env=env)
        util.misc.run_and_print(['make', 'install'], cwd=jellyfish_dir, env=env)


class Kraken(tools.Tool):

    BINS = ['kraken', 'kraken-build', 'kraken-filter', 'kraken-mpa-report', 'kraken-report', 'kraken-translate']

    def __init__(self, install_methods=None):
        if not install_methods:
            install_methods = [DownloadAndInstallKraken(URL, os.path.join(KRAKEN_DIR, 'bin', 'kraken'))]
        super().__init__(install_methods=install_methods)

    def version(self):
        return TOOL_VERSION

    @property
    def libexec(self):
        return os.path.dirname(self.executable_path())

    def build(self, db, options=None, option_string=None):
        '''Create a kraken database.

        Args:
          db: Kraken database directory to build. Must have library/ and
            taxonomy/ subdirectories to build from.
          *args: List of input filenames to process.
        '''
        return self.execute('kraken-build', db, options=options, option_string=option_string)

    def classify(self, db, args=None, options=None, option_string=None):
        """Classify input fasta/fastq

        Args:
          db: Kraken built database directory.
          args: List of input filenames to process.
        """
        assert len(args), 'Kraken requires input filenames.'
        return self.execute('kraken', db, args=args, options=options, option_string=option_string)

    def execute(self, command, db, args=None, options=None, option_string=None):
        '''Run a kraken-* command.

        Args:
          db: Kraken database directory.
          args: List of positional args.
          options: List of keyword options.
          option_string: Raw strip command line options.
        '''
        assert command in Kraken.BINS, 'Kraken command is unknown'
        options = options or {}
        option_string = option_string or ''
        args = args or []

        jellyfish_path = Jellyfish().install_and_get_path()
        env = os.environ.copy()
        env['PATH'] = '{}:{}'.format(os.path.dirname(jellyfish_path), env['PATH'])
        cmd = [os.path.join(self.libexec, command), '--db', db]
        # We need some way to allow empty options args like --build, hence
        # we filter out on 'x is None'.
        cmd.extend([str(x) for x in itertools.chain(*options.items()) if x is not None])
        cmd.extend(shlex.split(option_string))
        cmd.extend(args)
        log.debug('Calling %s: %s', command, ' '.join(cmd))
        return util.misc.run_and_print(cmd, env=env)


class DownloadAndInstallKraken(tools.DownloadPackage):

    def post_download(self):
        jellyfish_path = Jellyfish().install_and_get_path()
        env = os.environ.copy()
        env['PATH'] = '{}:{}'.format(os.path.dirname(jellyfish_path), env['PATH'])
        kraken_dir = os.path.join(self.destination_dir, KRAKEN_DIR)

        shutil.move(os.path.join(self.destination_dir, KRAKEN_COMMIT_DIR), kraken_dir)
        libexec_dir = os.path.join(kraken_dir, 'libexec')
        bin_dir = os.path.join(kraken_dir, 'bin')
        util.misc.run_and_print(['./install_kraken.sh', 'libexec'], cwd=kraken_dir, env=env)
        util.file.mkdir_p(bin_dir)
        for bin_name in Kraken.BINS:
            libexec_bin = os.path.join(libexec_dir, bin_name)
            os.symlink(libexec_bin, os.path.join(bin_dir, bin_name))
