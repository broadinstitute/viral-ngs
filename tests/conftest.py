import operator
import os
import shutil
import sys
import tempfile
import time
import contextlib
import string
import inspect
import copy
import functools

import util.misc
import util.file

import pytest

import tools

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
    parser.addoption(
        "--runslow", action="store_true", default=False, help="run slow tests"
    )


def pytest_configure(config):
    reporter = FixtureReporter(config)
    config.pluginmanager.register(reporter, 'fixturereporter')


def pytest_collection_modifyitems(config, items):
    if config.getoption("--runslow"):
        # --runslow given in cli: do not skip slow tests
        return
    skip_slow = pytest.mark.skip(reason="need --runslow option to run")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)


# Fixtures for creating a temp dir at session/module/class/function scope.
# Unlike pytest's tmpdir fixture, they use tempfile.mkdtemp to create a
# tempdir in the most secure/race-condition-free manner possible.
# Also, since util.file.tmp_dir() is used, the tempdir contens can be
# preserved for debugging by setting the environment variable VIRAL_NGS_TMP_DIRKEEP.


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


@pytest.fixture
def monkeypatch_function_result(monkeypatch):
    """Patches result of a function for specified args"""

    @contextlib.contextmanager
    def _set_function_result(f, *patch_args, **patch_kwargs):
        """Within the context, calls to function `f` with `patch_args` and `patch_kwargs` will return the
        `patch_result` (a keyword-only arg), or raise the `patch_exception` (a keyword-only arg) if that is given.
        Note that for the patched args, f will not be called, so any side effects of it won't happen.
        The function will be patched in the module `patch_module` (a keyword-only arg), or in `f` 's module
        if `patch_module` is None.  The keyword-only args `patch_result`, `patch_exception` and `patch_module`
        will be removed from `patch_kwargs` before calling f.  This context manager monkeypatches f for exactly
        one argument combination; multiple argument combos can be patched by nesting calls to this context manager.

        For a usage example, see test.unit.test_util_misc.test_available_cpu_count() .
        """

        patch_kwargs = copy.copy(patch_kwargs)
        patch_result = patch_kwargs.pop('patch_result', None)
        patch_exception = patch_kwargs.pop('patch_exception', None)
        util.misc.chk(patch_exception is None or patch_result is None)
        patch_module = patch_kwargs.pop('patch_module', inspect.getmodule(f))
        get_call_args = functools.partial(inspect.getcallargs, util.misc.unwrap(f))
        patch_call_args = get_call_args(*patch_args, **patch_kwargs)

        @util.misc.wraps(f)
        def patched_f(*args, **kwargs):
            if get_call_args(*args, **kwargs) != patch_call_args:
                return f(*args, **kwargs)

            if patch_exception:
                raise patch_exception

            return patch_result

        util.misc.chk(inspect.ismodule(patch_module) and getattr(patch_module, f.__name__) == f)
        with monkeypatch.context() as m:
            m.setattr(patch_module, f.__name__, patched_f)
            yield

    return _set_function_result



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

        writer = terminalreporter

        slowest = sorted(self.stats.items(), key=operator.itemgetter(1), reverse=True)
        if not self.durations:
            writer.write_sep("=", "slowest fixture durations")
        else:
            writer.write_sep("=", "slowest %s fixture durations" % self.durations)
            slowest = slowest[:self.durations]


        rows = []
        for (fixname, funcname), duration in slowest:
            row = ['{:.2f}s'.format(duration), fixname, funcname]
            rows.append(row)

        widths = [max(map(len, col)) for col in zip(*rows)]
        for row in rows:
            writer.write(" ".join((val.ljust(width) for val, width in zip(row, widths))))
            writer.line()
