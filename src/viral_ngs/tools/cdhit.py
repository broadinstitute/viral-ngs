'''
CD-HIT
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

log = logging.getLogger(__name__)

class CdHit(tools.Tool):

    COMMANDS = [
        'cd-hit',
        'cd-hit-est',
        'cd-hit-2d',
        'cd-hit-est-2d',
        'cd-hit-dup',
        ]

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [tools.PrexistingUnixCommand(shutil.which('cd-hit'), require_executability=True)]
        super(CdHit, self).__init__(install_methods=install_methods)

    def execute(self, command, input_fn, output_fn, options=None, option_string=None, background=None):
        '''Perform a clustering on DNA/RNA sequences

        Args:
          input_fn: Input fasta filename
          output_fn: Output fasta filename
        '''
        assert command in self.COMMANDS
        cmd = [command]
        cmd.extend(['-i', input_fn])
        cmd.extend(['-o', output_fn])

        options = options or {}
        if options:
            # We need some way to allow empty options args like --log, hence
            # we filter out on 'x is None'.
            cmd.extend([str(x) for x in itertools.chain(*options.items()) if x is not None])
        if option_string:
            cmd.extend(shlex.split(option_string))
        log.debug("Calling {}: {}".format(command, " ".join(cmd)))
        if background:
            return subprocess.Popen(cmd)
        else:
            return subprocess.check_call(cmd)
