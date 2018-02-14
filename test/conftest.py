import functools
import copy
import os

import util.cmd
import util.misc
import util.file

import pytest
import coverage

def add_measure_cmd_cov(cmd_impl, the_cov):
    
    @util.misc.wraps(cmd_impl)
    def run_with_cov(args):
        with util.file.tempfname('.coverage') as cov_fname:

            the_cov.get_data().write_file(cov_fname)
            the_cov.get_data().erase()

            try:
                cmd_impl(args)
            finally:
                cmd_cov = the_cov.get_data()
                #print('COV data for cmd', util.misc.unwrap(cmd_impl))
                #cmd_cov.write_fileobj(sys.stdout)

                # next: merge cmd_cov with other cov data for the same cmd 

                #print('---------end cov data--------')
                old_cov = coverage.CoverageData()
                old_cov.read_file(cov_fname)
                the_cov.get_data().update(old_cov)

    return run_with_cov

def pytest_configure(config):
    if 'VIRAL_NGS_GATHER_CMD_COVERAGE' not in os.environ: return
    plugin = config.pluginmanager.getplugin('_cov')
    if plugin and plugin.cov_controller:
        util.cmd.cmd_decorators.append(functools.partial(add_measure_cmd_cov, the_cov=plugin.cov_controller.cov))
