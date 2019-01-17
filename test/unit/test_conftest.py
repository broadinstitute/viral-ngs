"""Unit tests for pytest customizations"""

import pytest

def func_A(x, y=0):
    return x+y

def test_monkeypatch_function_result(monkeypatch_function_result):
    assert func_A(1) == 1
    with monkeypatch_function_result(func_A, 1, patch_result=2):
        assert func_A(1) == 2
    assert func_A(1) == 1

    with monkeypatch_function_result(func_A, x=1, patch_result=2):
        assert func_A(1) == 2
    assert func_A(1) == 1
    
    assert func_A(2, y=2) == 4
    with monkeypatch_function_result(func_A, 2, y=2, patch_result=100):
        assert func_A(2, y=2) == 100
        assert func_A(2) == 2
        assert func_A(2, 5) == 7

    with pytest.raises(RuntimeError):
        with monkeypatch_function_result(func_A, 2, y=2, patch_exception=RuntimeError()):
            func_A(2, y=2)

    with monkeypatch_function_result(func_A, 2, y=2, patch_result=100), \
         monkeypatch_function_result(func_A, x=3, patch_result=200):
        assert func_A(2, y=2) == 100
        assert func_A(3, y=0) == 200
        assert func_A(4, y=0) == 4

