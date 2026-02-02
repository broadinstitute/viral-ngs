'''class Tool, class InstallMethod, and related subclasses and methods'''

__author__ = "dpark@broadinstitute.org,irwin@broadinstitute.org"

import collections
import json
import operator
import os
import re
import logging
import tempfile
import shutil
import shlex
import subprocess
import util.file
import util.misc

from urllib.request import urlretrieve    # pylint: disable=E0611
from urllib.parse import urlparse    # pylint: disable=E0611

# Put all tool files in __all__
# allows "from tools import *" to import all tools for testtools
__all__ = sorted(
    [
        filename[:-3]    # Remove .py
        for filename in os.listdir(os.path.dirname(__file__))    # tools directory
        if filename.endswith(
            '.py') and filename != '__init__.py' and filename not in [    # Add any files to exclude here:
    # e.g. 'sometool.py',
            ]
    ]
)
installed_tools = {}

_log = logging.getLogger(__name__)


def iter_leaf_subclasses(a_class):
    "Iterate over subclasses at all levels that don't themselves have a subclass"
    is_leaf = True
    for subclass in sorted(a_class.__subclasses__(), key=operator.attrgetter("__name__")):
        is_leaf = False
        for leaf_class in iter_leaf_subclasses(subclass):
            if not getattr(leaf_class, '_skiptest', False):
                yield leaf_class
    if is_leaf:
        yield a_class

def all_tool_classes():
    return iter_leaf_subclasses(Tool)


def get_tool_by_name(name):
    if name not in installed_tools:
        raise NotImplementedError
    return installed_tools[name]


def skip_install_test(condition=None):
    '''Decorate the Tool class to skip the installation test.'''

    def decorator(klass):
        if callable(condition) and not condition():
            return klass
        klass._skiptest = True
        return klass

    return decorator


def is_osx():
    return os.uname()[0] == 'Darwin'


class Tool(object):
    ''' Base tool class that includes install machinery.
    '''

    def __init__(self, install_methods=None):
        install_methods = install_methods or []

        self.install_methods = install_methods
        self.installed_method = None
        self.exec_path = None
        self.tool_version = None

    def is_installed(self):
        return (self.installed_method is not None)

    def install(self):
        if not self.is_installed():
            for m in self.install_methods:
                if not m.is_attempted():
                    m.attempt_install()
                if m.is_installed():
                    self.installed_method = m
                    self.exec_path = m.executable_path()
                    installed_tools[self.__class__.__name__] = self
                    break

    def get_install_methods(self):
        return self.install_methods

    def set_install_methods(self, methods):
        self.install_methods = methods

    def version(self):
        if self.tool_version is None:
            self._get_tool_version()
        return self.tool_version

    def _get_tool_version():
        raise NotImplementedError

    def executable_path(self):
        return self.exec_path

    def execute(self, *args):
        assert not os.system(self.exec_path + ' ' + args)

    def install_and_get_path(self):
        self.install()
        if self.executable_path() == None:
            raise NameError("unsuccessful in installing " + type(self).__name__)
        return self.executable_path()


class InstallMethod(object):
    ''' Base class for installation methods for a given tool.
        None of these methods should ever fail/error. attempt_install should
        return silently regardless of the outcome (is_installed must be
        called to verify success or failure).
    '''

    def __init__(self):
        self.attempts = 0

    def is_attempted(self):
        return self.attempts

    def attempt_install(self):    # Override _attempt_install, not this.
        self.attempts += 1
        self._attempt_install()

    def _attempt_install(self):
        raise NotImplementedError

    def is_installed(self):
        raise NotImplementedError

    def executable_path(self):
        raise NotImplementedError


class PrexistingUnixCommand(InstallMethod):
    ''' This is an install method that tries to find whether an executable
        binary already exists for free on the unix file system--it doesn't
        actually try to install anything.
    '''

    def __init__(self, path, verifycmd=None, verifycode=0, require_executability=True):
        self.path = path
        self.verifycmd = verifycmd
        self.verifycode = verifycode
        self.installed = False
        self.require_executability = require_executability
        InstallMethod.__init__(self)

    def _attempt_install(self):
        if os.access(self.path, (os.X_OK | os.R_OK) if self.require_executability else os.R_OK):
            if self.verifycmd:
                self.installed = (os.system(self.verifycmd) == self.verifycode)
            else:
                self.installed = True
        else:
            self.installed = False

    def is_installed(self):
        if not self.is_attempted():
            self.attempt_install()
        return self.installed

    def executable_path(self):
        return self.installed and self.path or None

