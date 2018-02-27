'''
Infernal
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

CONDA_VERSION='1.1.2'

class Infernal(tools.Tool):

    COMMANDS = [
        'cmscan',
        'cmbuild',
        'cmpress',
        ]

    def __init__(self, install_methods=None):
        if not install_methods:
            install_methods = [
                tools.CondaPackage('infernal', version=CONDA_VERSION, executable='cmscan')
            ]
        super(Infernal, self).__init__(install_methods=install_methods)

    def execute(self, command, options=None, option_string=None):
        '''Execute an infernal command

        '''
        assert command in self.COMMANDS
        cmd = [os.path.join(os.path.dirname(self.install_and_get_path()), command)]

        options = options or {}
        if options:
            # We need some way to allow empty options args like --log, hence
            # we filter out on 'x is None'.
            cmd.extend([str(x) for x in itertools.chain(*options.items()) if x is not None])
        if option_string:
            cmd.extend(shlex.split(option_string))
        log.debug("Calling {}: {}".format(command, " ".join(cmd)))
        subprocess.check_call(cmd)

    def cmscan(self, db, input_fn, output_fn, num_threads=None, options=None, option_string=None):
        '''Perform an RNA structure search

        '''
        cmd = [os.path.join(os.path.dirname(self.install_and_get_path()), 'cmscan'),
               '--tblout', output_fn, db, input_fn]

        if num_threads:
            cmd.insert(1, str(num_threads))
            cmd.insert(1, '--cpu')

        options = options or {}
        if options:
            # We need some way to allow empty options args like --log, hence
            # we filter out on 'x is None'.
            cmd.extend([str(x) for x in itertools.chain(*options.items()) if x is not None])
        if option_string:
            cmd.extend(shlex.split(option_string))
        log.debug("Calling 'cmscan {}'".format(" ".join(cmd)))
        subprocess.check_call(cmd)
