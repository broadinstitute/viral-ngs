"""Unit tests for pytest customizations"""

import sys
import os
import os.path
import functools
import pytest

import util.misc
from util.file import slurp_file, get_test_input_path


def test_monkeypatch_function_result(monkeypatch_function_result):
    assert list(util.misc.unique([])) == []
    with monkeypatch_function_result(util.misc.unique, [], patch_result=[1]):
        assert list(util.misc.unique([])) == [1]
        assert list(util.misc.unique([2,1,2])) == [2, 1]
    assert list(util.misc.unique([])) == []

    inp = functools.partial(os.path.join, get_test_input_path())

    with pytest.raises(Exception):
        slurp_file('/some/file')

    with monkeypatch_function_result(slurp_file, '/some/file', patch_result='some_content',
                                     patch_module=sys.modules[__name__]):
        assert slurp_file('/some/file') == 'some_content'
        assert slurp_file(inp('ebov-makona.fasta')).startswith('>KJ660346.2')
        with monkeypatch_function_result(slurp_file, inp('ebov-makona.fasta'), patch_result='something_else',
                                         patch_module=sys.modules[__name__]):
            assert slurp_file('/some/file') == 'some_content'
            assert slurp_file(inp('ebov-makona.fasta')) == 'something_else'
        assert slurp_file(inp('ebov-makona.fasta')).startswith('>KJ660346.2')

        with monkeypatch_function_result(slurp_file, inp('ebov-makona.fasta'), patch_exception=RuntimeError(),
                                         patch_module=sys.modules[__name__]):
            with pytest.raises(Exception):
                slurp_file(inp('ebov-makona.fasta'))
            assert slurp_file('/some/file') == 'some_content'

    with pytest.raises(Exception):
        slurp_file('/some/file')
