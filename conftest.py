import operator
import os
import shutil
import sys
import tempfile
import time
import contextlib
import string

import util.file

import pytest

def timer():
    if sys.version_info < (3, 3):
        return time.time()
    return time.perf_counter()


def pytest_addoption(parser):
    group = parser.getgroup("reporting", description="terminal reporting", after="general")
    group.addoption(
        '--fixture-durations',
        action="store",
        type=int,
        default=None,
        metavar="N",
        help="show N slowest fixture durations (N=0 for all)."
    )

    group = parser.getgroup("test_selection", description="test selection", after="general")
    group.addoption(
        '--runslow',
        action="store_true",
        default=False,
        help="run slow tests."
    ),

def pytest_configure(config):
    reporter = FixtureReporter(config)
    config.pluginmanager.register(reporter, 'fixturereporter')

#
# Fixtures for creating a temp dir at session/module/class/function scope.
# Unlike pytest's tmpdir fixture, they use tempfile.mkdtemp to create a 
# tempdir in the most secure/race-condition-free manner possible.
# Also, since util.file.tmp_dir() is used, the tempdir contens can be
# preserved for debugging by setting the environment variable VIRAL_NGS_TMP_DIRKEEP.
#

@contextlib.contextmanager
def _tmpdir_aux(base_dir, scope, name):
    """Create and return a temporary directory; remove it and its contents on context exit."""
    with util.file.tmp_dir(dir=base_dir,
                           prefix='test-{}-{}-'.format(scope, name)) as tmpdir:
        yield tmpdir

@pytest.fixture(scope='session')
def tmpdir_session(request, tmpdir_factory):
    """Create a session-scope temporary directory."""
    with _tmpdir_aux(str(tmpdir_factory.getbasetemp()),
                     'session', id(request.session)) as tmpdir:
        yield tmpdir

@pytest.fixture(scope='module')
def tmpdir_module(request, tmpdir_session):
    """Create a module-scope temporary directory."""
    with _tmpdir_aux(tmpdir_session, 'module', request.module.__name__) as tmpdir:
        yield tmpdir

@pytest.fixture(scope='class')
def tmpdir_class(request, tmpdir_module):
    """Create a class-scope temporary directory."""
    with _tmpdir_aux(tmpdir_module, 'class',
                     request.cls.__name__ if request.cls else '__noclass__') as tmpdir:
        yield tmpdir

@pytest.fixture(autouse=True)
def tmpdir_function(request, tmpdir_class, monkeypatch):
    """Create a temporary directory and set it to be used by the tempfile module and as the TMPDIR environment variable."""
    with _tmpdir_aux(tmpdir_class, 'node', request.node.name) as tmpdir:
        monkeypatch.setattr(tempfile, 'tempdir', tmpdir)
        monkeypatch.setenv('TMPDIR', tmpdir)
        yield tmpdir

#############################################################################################

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