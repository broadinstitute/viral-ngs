import tools
import os.path
import subprocess
from os.path import join
from builtins import super

TOOL_NAME = 'krona'
CONDA_TOOL_VERSION = '2.7'


class Krona(tools.Tool):
    def __init__(self, install_methods=None):
        if not install_methods:
            install_methods = []
            install_methods.append(
                tools.CondaPackage(
                    TOOL_NAME,
                    version=CONDA_TOOL_VERSION,
                    executable='ktImportTaxonomy'))
        super(Krona, self).__init__(install_methods=install_methods)

    @property
    def opt(self):
        if not self.executable_path():
            self.install_and_get_path()
        bin_path = os.path.dirname(self.executable_path())
        # Get at the opt directory from the conda env root
        opt = os.path.abspath(join(bin_path, '..', 'opt', 'krona'))
        return opt

    def import_taxonomy(self,
                        db,
                        input_tsvs,
                        output,
                        query_column=None,
                        taxid_column=None,
                        score_column=None,
                        magnitude_column=None,
                        root_name=None,
                        no_hits=None,
                        no_rank=None):
        self.install()
        bin_path = os.path.dirname(self.executable_path())
        env = os.environ.copy()
        env['PATH'] = '{}:{}'.format(bin_path, env['PATH'])
        cmd = ['ktImportTaxonomy', '-tax', db, '-o', output]
        if query_column is not None:
            cmd.extend(['-q', str(query_column)])
        if taxid_column is not None:
            cmd.extend(['-t', str(taxid_column)])
        if score_column is not None:
            cmd.extend(['-s', str(score_column)])
        if magnitude_column is not None:
            cmd.extend(['-m', str(magnitude_column)])
        if root_name is not None:
            cmd.extend(['-n', root_name])
        if no_hits is not None:
            cmd.append('-i')
        if no_rank is not None:
            cmd.append('-k')
        cmd.extend(input_tsvs)

        subprocess.check_call(cmd, env=env)

    def create_db(self, db_dir):
        """Caution - this deletes the original .dmp files."""
        bin_path = os.path.dirname(self.executable_path())
        env = os.environ.copy()
        env['PATH'] = '{}:{}'.format(bin_path, env['PATH'])

        sh = join(self.opt, 'updateTaxonomy.sh')
        cmd = [sh, '--only-build', os.path.abspath(db_dir)]
        subprocess.check_call(cmd, env=env)
