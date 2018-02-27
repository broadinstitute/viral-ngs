'''
Prodigal ORF finder
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

log = logging.getLogger(__name__)

TOOL_VERSION = '2.6.3'

class Prodigal(tools.Tool):

    def __init__(self, install_methods=None):
        if not install_methods:
            install_methods = [
                tools.CondaPackage('prodigal', version=TOOL_VERSION)
            ]
        super(Prodigal, self).__init__(install_methods=install_methods)

    def execute(self, input_fn, output_translated_fn=None, options=None, option_string=None):
        '''Find ORFs

        Args:
          input_fn: Input fasta filename
          output_translated_fn: Output fasta filename
        '''
        cmd = [self.install_and_get_path()]
        cmd.extend(['-i', input_fn])
        if output_translated_fn:
            cmd.extend(['-a', output_translated_fn])

        options = options or {}
        if options:
            # We need some way to allow empty options args like --log, hence
            # we filter out on 'x is None'.
            cmd.extend([str(x) for x in itertools.chain(*options.items()) if x is not None])
        if option_string:
            cmd.extend(shlex.split(option_string))
        log.debug('Calling: {}'.format(' '.join(cmd)))
        util.misc.run(cmd)
