import operator
import os
import shutil
import sys
import tempfile
import time

import pytest

def timer():
    if sys.version_info < (3, 3):
        return time.time()
    return time.perf_counter()


def pytest_addoption(parser):
    parser.addoption("--runslow", action="store_true", help="run slow tests")

    group = parser.getgroup("terminal reporting", "reporting", after="general")
    group.addoption(
        '--fixture-durations',
        action="store",
        type=int,
        default=None,
        metavar="N",
        help="show N slowest fixture durations (N=0 for all)."
    ),


def pytest_configure(config):
    reporter = FixtureReporter(config)
    config.pluginmanager.register(reporter, 'fixturereporter')


@pytest.fixture(scope='session')
def tmpdir_session(request, tmpdir_factory):
    basetemp = str(tmpdir_factory.getbasetemp())
    tmpdir = tempfile.mkdtemp(dir=str(basetemp), prefix='test-session-'+str(id(request.session)+'-'))
    yield tmpdir
    if os.path.isdir(tmpdir): shutil.rmtree(tmpdir)

@pytest.fixture(scope='module')
def tmpdir_module(request, tmpdir_factory):
    basetemp = str(tmpdir_factory.getbasetemp())
    tmpdir = tempfile.mkdtemp(dir=str(basetemp), prefix='test-module-'+request.module.__name__+'-')
    yield tmpdir
    if os.path.isdir(tmpdir): shutil.rmtree(tmpdir)

@pytest.fixture(autouse=True)
def tmpdir_function(request, tmpdir, monkeypatch):
    tmpdir_str = str(tmpdir)
    monkeypatch.setattr(tempfile, 'tempdir', tmpdir_str)
    monkeypatch.setenv('TMPDIR', tmpdir_str)
    return tmpdir_str

class FixtureReporter:

    def __init__(self, config):
        import _pytest.config
        self.config = config
        self.stats = {}
        self.writer = _pytest.config.create_terminal_writer(config)
        self.durations = config.option.fixture_durations

    @pytest.hookimpl(hookwrapper=True)
    def pytest_fixture_setup(self, fixturedef, request):
        funcname = request._pyfuncitem.name
        fixname = fixturedef.argname

        fixturedef._timer_start = timer()
        yield
        duration = timer() - fixturedef._timer_start
        fixturedef._timer_duration = duration

        self.stats[(fixname, funcname)] = duration

    def pytest_terminal_summary(self, terminalreporter, exitstatus):
        if self.durations is None:
            return

        writer = terminalreporter.writer

        slowest = sorted(self.stats.items(), key=operator.itemgetter(1), reverse=True)
        if not self.durations:
            writer.sep("=", "slowest fixture durations")
        else:
            writer.sep("=", "slowest %s fixture durations" % self.durations)
            slowest = slowest[:self.durations]


        rows = []
        for (fixname, funcname), duration in slowest:
            row = ['{:.2f}s'.format(duration), fixname, funcname]
            rows.append(row)

        widths = [max(map(len, col)) for col in zip(*rows)]
        for row in rows:
            writer.write(" ".join((val.ljust(width) for val, width in zip(row, widths))))
            writer.line()


            
