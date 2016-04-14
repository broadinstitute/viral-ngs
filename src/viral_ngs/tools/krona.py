import tools
import os.path
from os.path import join
import util.misc
from builtins import super

TOOL_NAME = 'krona'
CONDA_TOOL_VERSION = '2.6'

class Krona(tools.Tool):

    BINS = ['kraken', 'kraken-build', 'kraken-filter', 'kraken-mpa-report', 'kraken-report', 'kraken-translate']

    def __init__(self, install_methods=None):
        if not install_methods:
            install_methods = []
            install_methods.append(tools.CondaPackage(TOOL_NAME, version=CONDA_TOOL_VERSION,
                                                      executable='ktImportTaxonomy'))
        super().__init__(install_methods=install_methods)


    def import_taxonomy(self, input_tsvs, output, query_column=None, taxid_column=None,
                        score_column=None, no_hits=None, no_rank=None):
        self.install()
        bin_path = os.path.dirname(self.executable_path())
        cmd = [join(bin_path, 'ktImportTaxonomy'), '-o', output]
        if query_column is not None:
            cmd.extend(['-q', str(query_column)])
        if taxid_column is not None:
            cmd.extend(['-t', str(taxid_column)])
        if score_column is not None:
            cmd.extend(['-s', str(score_column)])
        if no_hits is not None:
            cmd.append('-i')
        if no_rank is not None:
            cmd.append('-k')
        cmd.extend(input_tsvs)
        util.misc.run_and_print(cmd)
