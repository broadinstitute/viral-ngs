import warnings

import pytest

pytest_plugins = ["test.plugins.metadata_tester"]

@pytest.fixture(scope='session', autouse='true')
def warnings_as_errors():
    """Turn warnings into errors, so we see them during testing"""
    with warnings.catch_warnings():
        warnings.simplefilter('error', UserWarning)
        yield
