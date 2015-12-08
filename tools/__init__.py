'''class Tool, class InstallMethod, and related subclasses and methods'''

__author__ = "dpark@broadinstitute.org,irwin@broadinstitute.org"

import os
import logging
import tempfile
import shutil
import util.file
import util.misc
import json

try:
    # Python 3.x
    from urllib.request import urlretrieve # pylint: disable=E0611
    from urllib.parse import urlparse # pylint: disable=E0611
except ImportError:
    # Python 2.x
    from urllib import urlretrieve
    from urlparse import urlparse

# Put all tool files in __all__
# allows "from tools import *" to import all tooles for testtools
__all__ = [filename[:-3]  # Remove .py
           for filename in os.listdir(os.path.dirname(__file__))  # tools directory
           if filename.endswith('.py') and filename != '__init__.py' and filename not in [  # Add any files to exclude here:
               # e.g. 'sometool.py',
           ]]
installed_tools = {}

log = logging.getLogger(__name__)


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

        TO DO: add something about dependencies..
    '''

    def __init__(self, install_methods=None):
        install_methods = install_methods or []

        self.install_methods = install_methods
        self.installed_method = None
        self.exec_path = None

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
        return None

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

    def attempt_install(self):  # Override _attempt_install, not this.
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

class CondaPackage(InstallMethod):
    ''' This is an install method for tools that can be installed via 
        conda.
    '''

    def __init__(self, 
                 package, 
                 channel="bioconda", 
                 executable=None, 
                 version="", 
                 verifycmd=None, verifycode=0, require_executability=True,
                 env_path=None):
        # if the executable name is specifed, use it; otherwise use the package name
        self.executable = executable or package
        self.package    = package        
        self.channel    = channel
        self.version    = version

        self.verifycmd = verifycmd
        self.verifycode = verifycode
        self.require_executability = require_executability

        env_path = env_path or os.path.join(util.file.get_build_path(), 'conda-tools')
        self.env_path   = os.path.realpath(os.path.expanduser(env_path))

        self.installed = False

        # conda must be installed
        try:
            util.misc.run_and_print(["conda", "-h"], silent=True)
            log.debug("conda is installed")
        except:  
            log.error("conda must be installed")
            raise

        try:
            util.misc.run_and_print(["conda", "build", "-h"], silent=True)
        except:  
            log.warning("conda-build must be installed; installing...")
            util.misc.run_and_print(["conda", "install", "-y", "conda-build"])

        InstallMethod.__init__(self)

    @property
    def _package_str(self):
        if len(self.version) and float(self.version) > 0.0:
            ver_str = "{pkg}={ver}".format(pkg=self.package, ver=self.version)
        else:
            ver_str = self.package
        return ver_str

    @property
    def bin_path(self):
        return os.path.join(self.env_path, "bin")

    def executable_path(self):
        return os.path.join(self.bin_path, self.executable)

    def is_installed(self):
        result = util.misc.run_and_print(["conda", "list", "-p", self.env_path, "--json", self.package], silent=True)
        if result.returncode == 0:
            data = json.loads(result.stdout.decode("UTF-8"))
            if len(data) > 0:
                return True
        return False

    def verify_install(self):
        if os.access(self.executable_path(), (os.X_OK | os.R_OK) if self.require_executability else os.R_OK):
            if self.verifycmd:
                log.debug("validating")
                self.installed = (os.system(self.verifycmd) == self.verifycode)
            else:
                self.installed = True
        else:
            self.installed = False
        return self.installed

    def _attempt_install(self):
        self.installed = self.is_installed()

        if not self.installed:
            # try to create the environment and install the package
            result = util.misc.run_and_print(["conda", "create", "-q", "-y", "--json", "-c", self.channel, "-p", self.env_path, self._package_str], silent=True)
            data = json.loads(result.stdout.decode("UTF-8"))
            if "error" in data.keys() and "prefix already exists" in data["error"]:
                # the environment already exists
                # the package may not be installed...
                log.debug("Conda environment already exists...")                
                result = util.misc.run_and_print(["conda", "install", "--json", "-c", self.channel, "-y", "-q", "-p", self.env_path, self._package_str], silent=True)
                if result.returncode == 0:
                    data = json.loads(result.stdout.decode("UTF-8"))
                    if data["success"] == True:
                        # set self.installed = True
                        self.verify_install()
            else:
                if "success" in data.keys() and data["success"]:
                    # we were able to create the environment and install the package
                    log.debug("Conda environment created.")
                    if self.is_installed():
                        # set self.installed = True
                        self.verify_install()

class DownloadPackage(InstallMethod):
    ''' This is an install method for downloading, unpacking, and post-
            processing straight from the source.
        target_rel_path is the executable's path relative to destination_dir
        destination_dir defaults to the project build directory
        post_download_command will be executed if it isn't None, in
            destination_dir.
        if post_download_ret != None, assert it is returned by
            post_download_command
    '''

    def __init__(self, url, target_rel_path, destination_dir=None,
                 verifycmd=None, verifycode=0, require_executability=True,
                 post_download_command=None, post_download_ret=0):
        if destination_dir is None:
            destination_dir = util.file.get_build_path()
        self.url = url
        self.targetpath = os.path.join(destination_dir, target_rel_path)
        self.destination_dir = destination_dir
        self.verifycmd = verifycmd
        self.verifycode = verifycode
        self.installed = False
        self.require_executability = require_executability
        self.post_download_command = post_download_command
        self.post_download_ret = post_download_ret
        self.download_file = None
        InstallMethod.__init__(self)

    def is_installed(self):
        return self.installed

    def executable_path(self):
        return self.installed and self.targetpath or None

    def verify_install(self):
        if os.access(self.targetpath, (os.X_OK | os.R_OK) if self.require_executability else os.R_OK):
            if self.verifycmd:
                log.debug("validating")
                self.installed = (os.system(self.verifycmd) == self.verifycode)
            else:
                self.installed = True
        else:
            self.installed = False
        return self.installed

    def _attempt_install(self):
        if not self.verify_install():
            self.pre_download()
            self.download()
            self.post_download()
            self.verify_install()

    def pre_download(self):
        pass

    def download(self):
        download_dir = tempfile.gettempdir()
        util.file.mkdir_p(download_dir)
        filepath = urlparse(self.url).path
        file_basename = filepath.split('/')[-1]
        log.info("Downloading from %s to %s/%s ...", self.url, download_dir, file_basename)
        urlretrieve(self.url, os.path.join(download_dir, file_basename))
        self.download_file = file_basename
        self.unpack(download_dir)

    def post_download(self):
        if self.post_download_command:
            return_code = os.system('cd "{}" && {}'.format(self.destination_dir, self.post_download_command))
            if self.post_download_ret is not None:
                assert return_code == self.post_download_ret

    def unpack(self, download_dir):
        log.debug("unpacking")
        util.file.mkdir_p(self.destination_dir)
        if self.download_file.endswith('.zip'):
            if os.system("unzip -o %s/%s -d %s > /dev/null" % (download_dir, self.download_file, self.destination_dir
                                                              )):

                return
            else:
                os.unlink(os.path.join(download_dir, self.download_file))
        elif (self.download_file.endswith('.tar.gz') or self.download_file.endswith('.tgz') or
              self.download_file.endswith('.tar.bz2') or self.download_file.endswith('.tar')):
            if self.download_file.endswith('.tar'):
                compression_option = ''
            elif self.download_file.endswith('.tar.bz2'):
                compression_option = 'j'
            else:
                compression_option = 'z'
            untar_cmd = "tar -C {} -x{}pf {}/{}".format(self.destination_dir, compression_option, download_dir,
                                                        self.download_file)
            log.debug("Untaring with command: %s", untar_cmd)
            exitCode = os.system(untar_cmd)
            if exitCode:
                log.info("tar returned non-zero exitcode %s", exitCode)
                return
            else:
                log.debug("tar returned with exit code 0")
                os.unlink(os.path.join(download_dir, self.download_file))
        else:
            shutil.move(os.path.join(download_dir, self.download_file),
                        os.path.join(self.destination_dir, self.download_file))
