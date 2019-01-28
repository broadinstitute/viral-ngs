import logging
import os
import pytest
import tools


class StubCondaPackage(tools.Tool):

    # Skip gathering in all_tool_classes
    _skiptest = True

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
        self.executable = executable or package
        if type(version) == tools.CondaPackageVersion:
            self.version = version
        else:
            self.version = tools.CondaPackageVersion(version)

        self.env_path = None
        if 'CONDA_PREFIX' in os.environ and len(os.environ['CONDA_PREFIX']):
            last_path_component = os.path.basename(os.path.normpath(os.environ['CONDA_PREFIX']))
            self.env_path = os.path.dirname(os.environ['CONDA_PREFIX']) if last_path_component == "bin" else os.environ['CONDA_PREFIX']

        else:
            raise Exception

    def is_attempted(self):
        return True

    def is_installed(self):
        return True

    @property
    def bin_path(self):
        return os.path.join(self.env_path, 'bin')

    def executable_path(self):
        return os.path.join(self.bin_path, self.executable)

    def apply_patches(self):
        pass

    def verify_install(self):
        return

    def package_available(self):
        return True

    def uninstall_package(self):
        return True

    def install_package(self):
        return True

    def get_installed_version(self):
        return self.version

    def post_install(self):
        pass

    def execute(self, cmd, loglevel=logging.DEBUG):
        return {}
