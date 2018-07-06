import operator
import os
import shutil
import sys
import tempfile
import time
import contextlib
import string

import pytest

def timer():
    if sys.version_info < (3, 3):
        return time.time()
    return time.perf_counter()


def pytest_addoption(parser):
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

def _max_fname_len(file_system_path, default_max_len=80):
    name_max_str = [s for s in os.pathconf_names if s.endswith('_NAME_MAX')]
    if len(name_max_str) == 1:
        try:
            return os.pathconf(file_system_path, name_max_str[0])
        except OSError:
            pass
    return default_max_len

def _make_fname_valid(file_system_path, fname, len_margin):
    valid_fname_chars = set(string.ascii_letters+string.digits+'-.')
    max_len = _max_fname_len(file_system_path)-len_margin
    return ''.join((c if c in valid_fname_chars else '_') for c in fname[:max_len])

@contextlib.contextmanager
def _tmpdir_aux(request, tmpdir_factory, scope, name):
    basetemp = str(tmpdir_factory.getbasetemp())
    name = _make_fname_valid(file_system_path=basetemp, fname=name, len_margin=50)
    tmpdir = tempfile.mkdtemp(dir=basetemp, prefix='test-{}-{}-'.format(scope, name))
    yield tmpdir
    if os.path.isdir(tmpdir): shutil.rmtree(tmpdir)

@pytest.fixture(scope='session')
def tmpdir_session(request, tmpdir_factory):
    with _tmpdir_aux(request, tmpdir_factory, 'session', id(request.session)) as tmpdir:
        yield tmpdir

@pytest.fixture(scope='module')
def tmpdir_module(request, tmpdir_factory):
    with _tmpdir_aux(request, tmpdir_factory, 'module', request.module.__name__) as tmpdir:
        yield tmpdir

@pytest.fixture(autouse=True)
def tmpdir_function(request, tmpdir_factory, monkeypatch):
    with _tmpdir_aux(request, tmpdir_factory, 'node', request.node.name) as tmpdir:
        monkeypatch.setattr(tempfile, 'tempdir', tmpdir)
        monkeypatch.setenv('TMPDIR', tmpdir)
        yield tmpdir

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
