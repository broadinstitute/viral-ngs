# Unit tests for util.misc.py

__author__ = "dpark@broadinstitute.org"

import os, random, collections
import unittest
import subprocess

import util.misc
import util.file

import pytest

class TestRunAndPrint(unittest.TestCase):
    
    def testBasicRunSuccess(self):
        util.misc.run_and_print(['cat', '/dev/null'],
            silent=False, buffered=False, check=False)
        util.misc.run_and_print(['cat', '/dev/null'],
            silent=False, buffered=False, check=True)
        util.misc.run_and_print(['cat', '/dev/null'],
            silent=True, buffered=False, check=False)
        util.misc.run_and_print(['cat', '/dev/null'],
            silent=True, buffered=False, check=True)
        util.misc.run_and_print(['cat', '/dev/null'],
            silent=False, buffered=True, check=False)
        util.misc.run_and_print(['cat', '/dev/null'],
            silent=False, buffered=True, check=True)
        util.misc.run_and_print(['cat', '/dev/null'],
            silent=True, buffered=True, check=False)
        util.misc.run_and_print(['cat', '/dev/null'],
            silent=True, buffered=True, check=True)
            
    def testBasicRunFailDontCare(self):
        util.misc.run_and_print(['cat', '/notdev/notnull'],
            silent=False, buffered=False, check=False)
        util.misc.run_and_print(['cat', '/notdev/notnull'],
            silent=True, buffered=False, check=False)
        util.misc.run_and_print(['cat', '/notdev/notnull'],
            silent=False, buffered=True, check=False)
        util.misc.run_and_print(['cat', '/notdev/notnull'],
            silent=True, buffered=True, check=False)

    def testBasicRunFailAndCatch(self):
        self.assertRaises(subprocess.CalledProcessError,
            util.misc.run_and_print, ['cat', '/notdev/notnull'],
            silent=False, buffered=False, check=True)
        self.assertRaises(subprocess.CalledProcessError,
            util.misc.run_and_print, ['cat', '/notdev/notnull'],
            silent=False, buffered=True, check=True)
        self.assertRaises(subprocess.CalledProcessError,
            util.misc.run_and_print, ['cat', '/notdev/notnull'],
            silent=True, buffered=False, check=True)
        self.assertRaises(subprocess.CalledProcessError,
            util.misc.run_and_print, ['cat', '/notdev/notnull'],
            silent=True, buffered=True, check=True)


class TestFeatureSorter(unittest.TestCase):

    def testBasicSortingWithOverlap(self):
        fs = util.misc.FeatureSorter((
            ('abca', 10, 20),
            ('abca', 25, 35),
            ('abca', 15, 30),
        ))
        self.assertEqual(
            list(fs.get_features()),
            [
                ('abca', 10, 20, '+', None),
                ('abca', 15, 30, '+', None),
                ('abca', 25, 35, '+', None),
            ]
        )

    def testBasicIntervalsWithOverlap(self):
        fs = util.misc.FeatureSorter((
            ('abca', 10, 20),
            ('abca', 25, 35),
            ('abca', 15, 30),
        ))
        self.assertEqual(
            list(fs.get_intervals()),
            [
                ('abca', 10, 14, 1, [('abca', 10, 20, '+', None),]),
                ('abca', 15, 20, 2, [('abca', 10, 20, '+', None),('abca', 15, 30, '+', None),]),
                ('abca', 21, 24, 1, [('abca', 15, 30, '+', None),]),
                ('abca', 25, 30, 2, [('abca', 15, 30, '+', None),('abca', 25, 35, '+', None),]),
                ('abca', 31, 35, 1, [('abca', 25, 35, '+', None),]),
            ]
        )

    def testDisjointAndOverlappingIntervals(self):
        fs = util.misc.FeatureSorter((
            ('abca', 10, 20),
            ('abca', 80, 90, '+', None),
            ('abca', 25, 35, '-'),
            ('abca', 15, 30),
        ))
        self.assertEqual(
            list(fs.get_intervals()),
            [
                ('abca', 10, 14, 1, [('abca', 10, 20, '+', None),]),
                ('abca', 15, 20, 2, [('abca', 10, 20, '+', None),('abca', 15, 30, '+', None),]),
                ('abca', 21, 24, 1, [('abca', 15, 30, '+', None),]),
                ('abca', 25, 30, 2, [('abca', 15, 30, '+', None),('abca', 25, 35, '-', None),]),
                ('abca', 31, 35, 1, [('abca', 25, 35, '-', None),]),
                ('abca', 36, 79, 0, []),
                ('abca', 80, 90, 1, [('abca', 80, 90, '+', None),]),
            ]
        )

    def testMultiChrWindowedFeatures(self):
        fs = util.misc.FeatureSorter((
            ('abca', 10, 20),
            ('aaaa', 17, 33),
            ('abca', 80, 90),
            ('abca', 25, 35),
            ('abca', 15, 30),
        ))
        self.assertEqual(
            list(fs.get_features('abca', 11, 22)),
            [
                ('abca', 10, 20, '+', None),
                ('abca', 15, 30, '+', None),
            ]
        )

    def testOpenWindowRight(self):
        fs = util.misc.FeatureSorter((
            ('abca', 10, 20),
            ('aaaa', 17, 33),
            ('abca', 80, 90),
            ('abca', 25, 35),
            ('abca', 15, 30),
        ))
        self.assertEqual(
            list(fs.get_features('abca', 22)),
            [
                ('abca', 15, 30, '+', None),
                ('abca', 25, 35, '+', None),
                ('abca', 80, 90, '+', None),
            ]
        )

    def testOpenWindowLeft(self):
        fs = util.misc.FeatureSorter((
            ('abca', 10, 20),
            ('aaaa', 17, 33),
            ('abca', 80, 90),
            ('abca', 25, 35),
            ('abca', 15, 30),
        ))
        self.assertEqual(
            list(fs.get_features('abca', right=18)),
            [
                ('abca', 10, 20, '+', None),
                ('abca', 15, 30, '+', None),
            ]
        )

    def testMultiChrWithPayloadIntervals(self):
        fs = util.misc.FeatureSorter((
            ('abca', 10, 20),
            ('aaaa', 17, 33, '-', [100, 'name', []]),
            ('abca', 80, 90, '+', ['other info']),
            ('abca', 25, 35, '-'),
            ('abca', 15, 30),
        ))
        self.assertEqual(
            list(fs.get_intervals()),
            [
                ('abca', 10, 14, 1, [('abca', 10, 20, '+', None),]),
                ('abca', 15, 20, 2, [('abca', 10, 20, '+', None),('abca', 15, 30, '+', None),]),
                ('abca', 21, 24, 1, [('abca', 15, 30, '+', None),]),
                ('abca', 25, 30, 2, [('abca', 15, 30, '+', None),('abca', 25, 35, '-', None),]),
                ('abca', 31, 35, 1, [('abca', 25, 35, '-', None),]),
                ('abca', 36, 79, 0, []),
                ('abca', 80, 90, 1, [('abca', 80, 90, '+', ['other info']),]),
                ('aaaa', 17, 33, 1, [('aaaa', 17, 33, '-', [100, 'name', []]),]),
            ]
        )


class TestConfigIncludes(unittest.TestCase):

    def testConfigIncludes(self):

        def test_fn(f): return os.path.join(util.file.get_test_input_path(), 'TestUtilMisc', f)
        cfg1 = util.misc.load_config(test_fn('cfg1.yaml'))
        cfg2 = util.misc.load_config(test_fn('cfg2.yaml'), std_includes=[test_fn('cfg_std.yaml')],
                                     param_renamings={'std_param_A_old': 'std_param_A_new'})
        
        self.assertIn('paramA', cfg2)
        self.assertEqual(cfg2["env_vars"]["var_A"],1)
        self.assertEqual(cfg2["env_vars"]["var_B"],3)
        self.assertEqual(cfg2["env_vars"]["var_C"],4)
        with self.assertRaises(KeyError): cfg2["env_vars"]["var_E"]
        with self.assertRaises(KeyError): cfg2["var_A"]
        self.assertFalse(cfg2["paramZ"])
        self.assertTrue(cfg1["paramZ"])
        self.assertEqual(cfg2["empty_subtree"]["x"], 1)

        self.assertEqual(cfg2["std_methods"], ['a','b','c'])
        self.assertEqual(cfg2["stage1"]["stage2"]["step_num"],3)
        self.assertEqual(cfg2["stage1"]["stage2"]["step_list"],[5,10,15])
        self.assertEqual(cfg2["stage1"]["stage3"]["step_list"],[3,33])
        self.assertEqual(cfg1["stage1"]["stage3"]["step_list"],[51,101,151])

        self.assertEqual(cfg2["std_param_A_new"], 111)  # specified as std_param_A_old in cfg1.yaml

        self.assertEqual(util.misc.load_config(test_fn('empty.yaml')), {})

def test_tmp_set_env():
    """Test tmp_set_env()"""

    var = 'VIRAL_NGS_TEST_VAR'
    tmp_set_env = util.misc.tmp_set_env

    save_environ = dict(os.environ)

    assert var not in os.environ

    with tmp_set_env(var, 239) as old_val:
        assert not old_val
        assert os.environ[var] == '239'

        with tmp_set_env(var, 300) as old_val1:
            assert old_val1 == '239'
            assert os.environ[var] == '300'
            
        assert os.environ[var] == '239'
        with tmp_set_env(var, None):
            assert var not in os.environ
        assert os.environ[var] == '239'

        with pytest.raises(AssertionError):
            with tmp_set_env(var, None, sep=':'):
                pass

        with tmp_set_env(var, 300, sep=':') as old_val1:
            assert old_val1 == '239'
            assert os.environ[var] == '239:300'

            with tmp_set_env(var, 400, sep=':', append=False) as old_val2:
                assert old_val2 == '239:300'
                assert os.environ[var] == '400:239:300'
            assert os.environ[var] == '239:300'
        assert os.environ[var] == '239'
        
    assert var not in os.environ

    with tmp_set_env(var, 300, sep=':') as old_val:
        assert old_val is None
        assert os.environ[var] == '300'
    assert var not in os.environ

    with tmp_set_env(var, 300, sep=':', append=False) as old_val:
        assert old_val is None
        assert os.environ[var] == '300'
    assert var not in os.environ

    with tmp_set_env(var, None) as old_val:
        assert old_val is None
        assert var not in os.environ

    assert var not in os.environ

    with pytest.raises(AssertionError):
        with tmp_set_env(var, 1):
            os.environ[var] = '2'

    with pytest.raises(AssertionError):
        with tmp_set_env(var, None):
            os.environ[var] = '2'

    with pytest.raises(AssertionError):
        with tmp_set_env(var, 1):
            del os.environ[var]

    assert os.environ == save_environ
# end: def test_tmp_set_env()

def test_func_wrap():
    """Test wrap/unrap decorator helpers"""

    def f(): 
        print('in f')

    @util.misc.wraps(f)
    def g():
        print('in g')

    assert util.misc.unwrap(g) == f

    @util.misc.wraps(g)
    def h(): pass

    assert util.misc.unwrap(h) == f

def test_flatten_dict():
    """Test dict flattening"""

    flatten_dict = util.misc.flatten_dict
    
    inp_out = (
        ([{}], {}),
        ([{1:2, 3: {4:5}}],
         {(1,): 2, (3, 4): 5}),
        ([{1:2, 3: {4:5, 6:7}}],
         {(1,): 2, (3, 4): 5, (3,6):7}),
        ([{1:2, 3: {4:5, 6:7}, 8: {}}],
         {(1,): 2, (3, 4): 5, (3,6):7}),
        ([{1:2, 3: {4:5, 6:7}, 8: {9: {}}}],
         {(1,): 2, (3, 4): 5, (3,6):7}),
        ([{1:2, 3: {4:5, 6:7}, 8: {9: {10: (11, 12)}}}],
         {(1,): 2, (3, 4): 5, (3,6):7, (8,9,10): (11,12)}),
        ([(1,2), (tuple,)],
         {(0,): 1, (1,): 2}),
        ([{1:2, 3: {4:[5,8], 6:(7,10)}, 8: {9: []}}, (tuple,list)],
         {(1,): 2, (3, 4, 0): 5, (3, 4, 1): 8, (3, 6, 0): 7, (3, 6, 1): 10}),
    )
    for inp, out in inp_out:
        assert flatten_dict(*inp) == out

def test_flatten():
    """Test flatten()"""

    flatten = util.misc.flatten
    assert flatten([]) == []
    assert flatten(()) == []
    assert flatten([[]]) == []
    assert flatten([(),[],[]]) == []
    assert flatten([1]) == [1]
    assert flatten((1,)) == [1]
    assert flatten([1, [2]]) == [1, 2]
    assert flatten([1, [2, [3]]]) == [1, 2, 3]
    assert flatten([1, [], [2, [3, []]]]) == [1, 2, 3]
    assert flatten((1, [], (2, (3,)))) == [1, 2, 3]
    assert flatten([1, [2, [3]]], types_to_flatten=(tuple,)) == [[1, [2, [3]]]]
    assert flatten([1, [2, [3]]], types_to_flatten=(tuple,list)) == [1, 2, 3]

def test_get_func_arg_names():
    """Test getting function arg names"""

    gnfa = util.misc.get_named_func_args

    assert gnfa(lambda : None) == []
    assert gnfa(lambda a: None) == ['a']
    assert gnfa(lambda a, b: None) == ['a', 'b']
    assert gnfa(lambda a, b, c=None: None) == ['a', 'b', 'c']
    assert gnfa(lambda a, b, c=None, *args, **kwargs: None) == ['a', 'b', 'c']
    def f(a, b, c=None, *args, **kwargs):
        return 239
    assert gnfa(lambda a, b, c=None, *args, **kwargs: None) == ['a', 'b', 'c']

def test_dict_subset():
    """Test dict_subset()"""
    dict_subset = util.misc.dict_subset
    assert dict_subset({}, {}) == {}
    assert dict_subset({1:2}, {}) == {}
    assert dict_subset({1:2}, {1,2}) == {1:2}
    assert dict_subset({1:2,3:4}, {1,2,3}) == {1:2,3:4}
    assert dict_subset({1:2,3:4}, {2,3,5}) == {3:4}

def test_make_list():
    """Test make_list()"""
    assert util.misc.make_list() == []
    assert util.misc.make_list(1) == [1]
    assert util.misc.make_list(1, 2) == [1, 2]
    assert util.misc.make_list(1, 1) == [1, 1]
