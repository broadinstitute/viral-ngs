import operator
import os
import pytest
import shutil
import sys
import tempfile
import time


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


@pytest.fixture(scope='session')
def tmpdir_session(request, tmpdir_factory):
    tmpdir = str(tmpdir_factory.mktemp('test-session'))

    def reset():
        shutil.rmtree(tmpdir)

    request.addfinalizer(reset)
    return tmpdir


@pytest.fixture(scope='module')
def tmpdir_module(request, tmpdir_factory):
    tmpdir = str(tmpdir_factory.mktemp('test-module'))

    def reset():
        shutil.rmtree(tmpdir)

    request.addfinalizer(reset)
    return tmpdir


@pytest.fixture(autouse=True)
def tmpdir_function(request, tmpdir_factory):
    old_tempdir = tempfile.tempdir
    old_env_tmpdir = os.environ.get('TMPDIR')
    new_tempdir = str(tmpdir_factory.mktemp('test-function'))
    tempfile.tempdir = new_tempdir
    os.environ['TMPDIR'] = new_tempdir

    def reset():
        shutil.rmtree(new_tempdir)
        tempfile.tmpdir = old_tempdir
        if old_env_tmpdir:
            os.environ['TMPDIR'] = old_env_tmpdir

    request.addfinalizer(reset)
    return new_tempdir


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
