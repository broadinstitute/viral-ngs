import pytest

import util.metadata

def pytest_runtest_setup(item):
    util.metadata.record_test_start(item.nodeid)

def pytest_unconfigure(config):
    util.metadata.tests_ended()


