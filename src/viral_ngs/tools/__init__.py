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


class CondaPackageVersion(object):

    def __init__(self, version, build_type=None):
        self.version = version
        self.build_type = build_type

    def __cmp__(self, other):
        return self.__dict__ == other.__dict__

    def __repr__(self):
        if self.build_type:
            return '{}-{}'.format(self.version, self.build_type)
        else:
            return self.version

    def satisfies(self, spec):
        if spec.build_type:
            if self.build_type != spec.build_type:
                return False
        # if the specified version is blank/None (not specified), any version will work
        return (spec.version=="") or (spec.version==None) or (self.version == spec.version)

    @property
    def version_spec(self):
        if self.build_type:
            return '{}={}'.format(self.version, self.build_type)
        else:
            return self.version


class CondaPackage(InstallMethod):
    ''' This is an install method for tools that can be installed via
        conda.
    '''

    QUIET_COMMANDS = [
        'install',
        'create',
        'remove',
        ]

    def execute(self, cmd, loglevel=logging.DEBUG, buffered=None, check=None, silent=None, stderr=None):
        run_cmd = ['conda']
        if cmd[0] in self.QUIET_COMMANDS:
            run_cmd.extend(['-q', '-y'])
        run_cmd.extend(cmd)
        result = util.misc.run_and_print(
            run_cmd, loglevel=loglevel, env=self.conda_env, buffered=buffered, check=check, silent=silent, stderr=stderr)

        if result.returncode == 0:
            try:
                command_output = result.stdout.decode("UTF-8")
                if cmd[0] == '-V':
                    return command_output
                return json.loads(self._string_from_start_of_json(command_output.strip()))
            except Exception as e:
                _log.warning("Failed to decode JSON during conda command '%s' emitting output: %s", " ".join(run_cmd),result.stdout.decode("UTF-8"))
                return  # return rather than raise so we can fall back to the next install method

    def __init__(
        self,
        package,
        channel="bioconda",
        executable=None,
        version="",
        verifycmd=None,
        verifycode=0,
        require_executability=True,
        env=None,
        env_root_path=None,
        conda_cache_path=None,
        patches=None,
        post_install_command=None,
        post_install_ret=0,
        post_verify_command=None,
        post_verify_ret=0
    ):
        # if the executable name is specifed, use it; otherwise use the package name
        self.executable = executable or package
        self.package = package
        self.channel = channel
        if type(version) == CondaPackageVersion:
            self.version = version
        else:
            self.version = CondaPackageVersion(version)

        self.post_install_command = post_install_command
        self.post_install_ret = post_install_ret

        # call the post-verification command.
        # Useful for copying in license files, building databases, etc.
        # The post-verify command is executed relative to the conda environment bin/
        # And has the active conda environment on the PATH
        self.post_verify_command = post_verify_command
        self.post_verify_ret = post_verify_ret
        self.post_verify_cmd_executed = False

        self.verifycmd = verifycmd
        self.verifycode = verifycode
        self.require_executability = require_executability
        self.patches = patches or []

        # if we have specified a conda env/root, use it
        # alternatively, use cms tools dir as default env root if os.environ["CONDA_PREFIX"] is not defined
        #
        # as it is always the path
        # CONDA_PREFIX is always a full path, but may not be present in older versions
        # CONDA_ENV_PATH may be used instead
        # in even older versions of conda, CONDA_DEFAULT_ENV can be used as another fallback
        self.env_path = None
        if "CONDA_PREFIX" in os.environ and len(os.environ["CONDA_PREFIX"]):
            #_log.debug('CONDA_PREFIX found')
            last_path_component = os.path.basename(os.path.normpath(os.environ["CONDA_PREFIX"]))
            self.env_path = os.path.dirname(os.environ["CONDA_PREFIX"]) if last_path_component == "bin" else os.environ["CONDA_PREFIX"]
        elif "CONDA_ENV_PATH" in os.environ and len(os.environ["CONDA_ENV_PATH"]):
            #_log.debug('CONDA_ENV_PATH found')
            last_path_component = os.path.basename(os.path.normpath(os.environ["CONDA_ENV_PATH"]))
            self.env_path = os.path.dirname(os.environ["CONDA_ENV_PATH"]) if last_path_component == "bin" else os.environ["CONDA_ENV_PATH"]
        elif "CONDA_DEFAULT_ENV" in os.environ and len(os.environ["CONDA_DEFAULT_ENV"]):
            #_log.debug('CONDA_PREFIX not found, using CONDA_DEFAULT_ENV')
            conda_env_path = os.environ.get('CONDA_DEFAULT_ENV')  # path to current conda environment
            if conda_env_path:
                if os.path.isdir(conda_env_path):
                    #_log.debug('Conda env found is specified as dir: %s' % conda_env_path)
                    conda_env_path = os.path.abspath(conda_env_path)
                    last_path_component = os.path.basename(os.path.normpath(conda_env_path))
                    self.env_path = os.path.dirname(last_path_component) if last_path_component == "bin" else conda_env_path
                else: # if conda env is an environment name, infer the path
                    #_log.debug('Conda env found is specified by name: %s' % conda_env_path)
                    result = util.misc.run_and_print(["conda", "env", "list", "--json"], silent=True, env=os.environ)
                    if result.returncode == 0:
                        command_output = result.stdout.decode("UTF-8")
                        data = json.loads(self._string_from_start_of_json(command_output))
                        if len(data) > 0:
                            if "envs" in data and len(data["envs"]):
                                for item in data["envs"]:
                                    if os.path.basename(os.path.realpath(item)) == conda_env_path:
                                        self.env_path = os.path.realpath(item)
                                        break

        # if the env is being overridden, or if we could not find an active conda env
        if env_root_path or env or not self.env_path:
            env_root_path = env_root_path or os.path.join(util.file.get_project_path(), 'tools', 'conda-tools')
            env = env or 'default'
            self.env_path = os.path.realpath(os.path.expanduser(
                os.path.join(env_root_path, env)))

        # set an env variable to the conda cache path. this env gets passed to the
        # the subprocess, and the variable instructs conda where to place its cache files
        conda_cache_path = conda_cache_path or os.path.join(util.file.get_project_path(), 'tools', 'conda-cache')
        self.conda_cache_path = os.path.realpath(os.path.expanduser(conda_cache_path))
        self.conda_env = os.environ
        old_envs_path = os.environ.get('CONDA_DEFAULT_ENV')
        self.conda_env["CONDA_ENVS_PATH"] = conda_cache_path+":"+os.path.dirname(self.env_path)

        #_log.info("Tool install conda env path: %s", self.env_path)
        self.installed = False
        self._is_attempted = False

        super(CondaPackage, self).__init__()

    @staticmethod
    def _string_from_start_of_json(string_with_json):
        # JSON can start with "{" or "["
        # via http://www.json.org/
        try:
            matches = re.compile("\{|\[").search(string_with_json)
            if matches:
                return string_with_json[matches.start():]
            else:
                _log.warning("Does not look like json: %s" % string_with_json)
                return None
        except:
            _log.warning("Does not look like json: %s" % string_with_json)
            return None

    @property
    def _package_str(self):
        if self.version:
            ver_str = "{pkg}={ver}".format(
                pkg=self.package, ver=self.version.version_spec)
        else:
            ver_str = self.package
        return ver_str

    def is_attempted(self):
        return self._is_attempted

    def is_installed(self):
        return self.installed

    @property
    def bin_path(self):
        return os.path.join(self.env_path, "bin")

    def executable_path(self):
        return os.path.join(self.bin_path, self.executable)

    def apply_patches(self):
        for path, patch in self.patches:
            self._patch(path, patch)


    def _patch(self, path, patch):
        """Patch a path relative to conda root.

        The patch is relative to the tools/patches directory.
        """
        file_path = os.path.join(self.env_path, path)
        patch_path = os.path.join(
            util.file.get_project_path(), 'tools', 'patches', patch)
        subprocess.check_call(['patch', file_path, patch_path])

    @property
    def _package_installed(self):
        data = self.execute(['list', "-f", "-c", "-p", self.env_path, "--json", self.package], silent=True)
        if len(data) > 0:
            _log.debug('Conda package found: {}'.format(data))
            return True
        return False

    def verify_install(self):
        # if conda does not the package as installed, no need to check further
        installed_version = self.get_installed_version()
        if not installed_version:
            # report the package as not installed
            return False

        # if the package is installed, check the binary for executability
        if os.access(self.executable_path(), (os.X_OK | os.R_OK) if self.require_executability else os.R_OK):
            # optionally use the verify command, if specified
            if self.verifycmd:
                if os.system(self.verifycmd) == self.verifycode:
                    _log.debug("Validating with cmd: {}".format(self.verifycmd))
                    self.installed = installed_version
            else:
                self.installed = installed_version
        else:
            self.installed = False
        if self.installed:
            # call the post-verification command.
            # Useful for copying in license files, building databases, etc.
            # This is executed relative to the conda environment bin/
            # And has the active conda environment on the PATH
            if (not self.post_verify_cmd_executed) and self.post_verify_command:
                post_verify_command = shlex.split(self.post_verify_command)
                _log.debug("Running post-verification cmd: {}".format(self.post_verify_command))

                result = util.misc.run_and_print(post_verify_command, silent=False, check=False, env=self.conda_env, cwd=self.bin_path)
                post_verify_cmd_return_code = result.returncode
                if post_verify_cmd_return_code == self.post_verify_ret:
                    self.post_verify_cmd_executed = True
                else:
                    raise subprocess.CalledProcessError(post_verify_cmd_return_code, "Post verification command failed with exit %s: %s" % (post_verify_cmd_return_code, self.post_verify_command))

            return installed_version
        return False

    def _attempt_install(self):
        try:
            # check for presence of conda command
            data = self.execute(['-V'], check=True, silent=True)
        except:
            _log.warning("conda NOT installed")
            self._is_attempted = True
            self.installed = False
            return

        # conda-build is not needed for pre-built binaries from conda channels
        # though we may will need it in the future for custom local builds
        # try:
        #     util.misc.run_and_print(["conda", "build", "-V"], silent=True, check=True, env=self.conda_env)
        # except:
        #     _log.warning("conda-build must be installed; installing...")
        #     util.misc.run_and_print(["conda", "install", "-y", "conda-build"], check=True)

        # if the package is already installed, we need to check if the version is correct
        pkg_version = self.verify_install()
        if pkg_version:
            _log.debug("Currently installed version of {package}: {version}".format(
                package=self.package, version=pkg_version))
            # if the installed version is not the one specified
            if not pkg_version.satisfies(self.version):
                _log.debug("Expected version of {package}:            {version}".format(package=self.package, version=self.version))
                _log.debug("Incorrect version of {package} installed. Removing it...".format(package=self.package) )

                # uninstall the current (incorrect) version
                self.uninstall_package()
                # and continue to install...
            else:
                # if the package is installed and is the correct version
                # return so we don't bother installing
                return

        # install the package and verify
        _log.debug("Attempting install...")
        self.install_package()
        self.verify_install()
        self.post_install()

    def get_installed_version(self):
        # If we ever use conda to install pip packages as tools, "-c" needs to be removed

        data = self.execute(["list", "-c", "--json", "-f", "-p", self.env_path, self.package], check=True, silent=True)
        if data is None or not len(data):
            return
        if isinstance(data[0], dict):
            installed_package_string = data[0]["dist_name"]
        else:
            installed_package_string = data[0]
        # regex to match package specs in the format bioconda::biopython-1.68-py35_0
        package_info_re = re.compile(r"(?:(?P<channel>.*)::)?(?P<package_name>.*)-(?P<version>.*)-(?P<build_type>.*)")
        matches = package_info_re.match(installed_package_string)
        if matches:
            installed_version = matches.group("version")
            installed_package = matches.group("package_name")
            installed_build_type = matches.group("build_type")
            return CondaPackageVersion(installed_version, installed_build_type)

    def package_available(self):
        # If we ever use conda to install pip packages as tools, "-c" needs to be removed
        data = self.execute(["search", "--json", "-c", self.channel, self.package], check=False, loglevel=logging.INFO)
        if data and len(data) and self.package in data and "error" not in data:
            for sub_pkg in data[self.package]:
                if sub_pkg.get("version", "") == self.version.version_spec:
                    return True
        _log.info("Conda package for %s is not available on this platform.", self.package)
        return False

    def uninstall_package(self):
        data = self.execute(["remove", "-q", "-y", "--json", "-p", self.env_path, self.package],
                            loglevel=logging.INFO)
        if not data:
            return
        if data["success"] == True:
            _log.debug("Package removed.")
            self.installed = False

    def install_package(self):
        data = self.execute(["list", "--json", "-p", self.env_path], silent=True, check=True)
        if not data:
            return
        for d in data:
            if "error" in d and "Not a conda environment" in d["message"]:
                _log.warning("Conda environment doesn't exist")
                return

        # the environment already exists
        # the package may not be installed...
        _log.debug("Conda environment already exists. Installing package...")

        data = self.execute(["install", "--json", "-c", self.channel, "-y", "-q", "--no-update-deps", "-p", self.env_path, self._package_str])
        if not data:
            return

        if data["success"] == True:
            _log.debug("Package installed.")
        self.apply_patches()

    def post_install(self):
        """
            Runs a shell command after package installation,
            relative to the directory containing the executable
        """
        if self.post_install_command:
            destination_dir = os.path.dirname(os.path.realpath(self.executable_path()))
            return_code = os.system('cd "{}" && {}'.format(destination_dir, self.post_install_command))
            if self.post_install_ret is not None:
                assert return_code == self.post_install_ret


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

    def __init__(
        self,
        url,
        target_rel_path,
        destination_dir=None,
        verifycmd=None,
        verifycode=0,
        require_executability=True,
        post_download_command=None,
        post_download_ret=0
    ):
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
                _log.debug("validating")
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
        _log.info("Downloading from %s to %s/%s ...", self.url, download_dir, file_basename)
        urlretrieve(self.url, os.path.join(download_dir, file_basename))
        self.download_file = file_basename
        self.unpack(download_dir)

    def post_download(self):
        if self.post_download_command:
            return_code = os.system('cd "{}" && {}'.format(self.destination_dir, self.post_download_command))
            if self.post_download_ret is not None:
                assert return_code == self.post_download_ret

    def unpack(self, download_dir):
        _log.debug("unpacking")
        util.file.mkdir_p(self.destination_dir)
        if self.download_file.endswith('.zip'):
            if os.system("unzip -o %s/%s -d %s > /dev/null" % (download_dir, self.download_file, self.destination_dir)):

                return
            else:
                os.unlink(os.path.join(download_dir, self.download_file))
        elif (
            self.download_file.endswith('.tar.gz') or self.download_file.endswith('.tgz') or
            self.download_file.endswith('.tar.bz2') or self.download_file.endswith('.tar')
        ):
            if self.download_file.endswith('.tar'):
                compression_option = ''
            elif self.download_file.endswith('.tar.bz2'):
                compression_option = 'j'
            else:
                compression_option = 'z'
            untar_cmd = "tar -C {} -x{}pf {}/{}".format(
                self.destination_dir, compression_option, download_dir, self.download_file
            )
            _log.debug("Untaring with command: %s", untar_cmd)
            exitCode = os.system(untar_cmd)
            if exitCode:
                _log.info("tar returned non-zero exitcode %s", exitCode)
                return
            else:
                _log.debug("tar returned with exit code 0")
                os.unlink(os.path.join(download_dir, self.download_file))
        else:
            shutil.move(
                os.path.join(download_dir, self.download_file), os.path.join(self.destination_dir, self.download_file)
            )
