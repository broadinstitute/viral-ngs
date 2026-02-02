# Unit tests for util.misc.py

__author__ = "dpark@broadinstitute.org"

import os, random, collections
import unittest
import subprocess
import multiprocessing
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

def test_as_type():
    """Test util.misc.as_type()"""

    as_type = util.misc.as_type

    test_data = (
        ('1', int, 1, int),
        ('1', (int,), 1, int),
        ('1', (int, float), 1, int),
        ('1', (float, int), 1., float),
        ('1.2', (float, int), 1.2, float),
        ('1.2', (int, float), 1.2, float),
        (1, int, 1, int),
        ('1.', (int, float), 1., float),
        (1., int, 1, int),
        (1.2, (int, float), 1, int),
        ('1e3', (int, float), 1000., float),
        ('1e3', (float, int), 1000., float),
        ('1.1e3', (int, float), 1100., float),
        ('-1.1e3', (float, int), -1100., float),
    )

    for val, types, out_val, out_type in test_data:
        result = as_type(val, types)
        assert result == out_val
        assert type(result) == out_type

    err_data = (
        ('1.', int), ('1.', (int,)), ('1e3', int), ('1e3', (int,)), ([1], int),
        ('', int), ('', (int, float)), ('', (float, int))
    )

    for val, types in err_data:
        with pytest.raises(TypeError):
            as_type(val, types)

@pytest.mark.parametrize("iter_d", [False, True])
@pytest.mark.parametrize("iter_subset", [False, True])
def test_subdict(iter_d, iter_subset):
    """Test util.misc.subdict()"""
    def subdict(d, subset): 
        return util.misc.subdict(iter(d.items()) if iter_d else d, 
                                 iter(subset) if iter_subset else subset)

    test_data = (
        ({}, {}, {}),
        ({1:2}, {}, {}),
        ({1:2}, {1}, {1:2}),
        ({1:2}, {2}, {}),
        ({1:2}, {1,2}, {1:2}),
        ({1:2,3:4}, {1,2,3}, {1:2,3:4}),
        ({1:2,3:4}, {2,3,5}, {3:4}),
    )

    for d, subset, expected in test_data:
        assert subdict(d, subset) == expected

        assert set(subdict(d, subset).keys()) == (set(d.keys()) & set(subset))
        assert subdict(d, {}) == {}
        assert subdict(d, d.keys()) == d
        assert subdict(d, list(d.keys())*2) == d

def test_chk():
    chk = util.misc.chk
    chk(True, 'no error')
    with pytest.raises(RuntimeError):
        chk(False)
    with pytest.raises(RuntimeError):
        chk(2 == 3, 'Something wrong')
    with pytest.raises(TypeError):
        chk(isinstance(None, int), 'Expected an int', TypeError)

def test_available_cpu_count(monkeypatch_function_result):
    reported_cpu_count = util.misc.available_cpu_count()
    assert reported_cpu_count >= int(os.environ.get('PYTEST_XDIST_WORKER_COUNT', '1'))
    assert util.misc.available_cpu_count() == reported_cpu_count

    # cgroup v2 limited to 1 cpu
    with monkeypatch_function_result(os.path.exists, "/sys/fs/cgroup/cgroup.controllers", patch_result=True, patch_module=os.path), \
         monkeypatch_function_result(util.file.slurp_file, '/sys/fs/cgroup/cpu.max', patch_result="100000 100000"):
        assert util.misc.available_cpu_count() == 1

    # cgroup v2 limited to 2 cpu
    with monkeypatch_function_result(os.path.exists, "/sys/fs/cgroup/cgroup.controllers", patch_result=True, patch_module=os.path), \
         monkeypatch_function_result(util.file.slurp_file, '/sys/fs/cgroup/cpu.max', patch_result="200000 100000"):
        assert util.misc.available_cpu_count() == 2

    # cgroup v2 with no CPU limit imposed on cgroup
    # (fall back to /proc/self/status method, with limit imposed there):
    #   'Cpus_allowed:  d' = 0b1101 bitmask (meaning execution allowed on 3 CPUs)
    with monkeypatch_function_result(os.path.exists, "/sys/fs/cgroup/cgroup.controllers", patch_result=True, patch_module=os.path), \
         monkeypatch_function_result(util.file.slurp_file, '/sys/fs/cgroup/cpu.max', patch_result="max 100000"), \
         monkeypatch_function_result(util.file.slurp_file, '/proc/self/status', patch_result='Cpus_allowed:  d'):
        assert util.misc.available_cpu_count() == 3

    # cgroup v1 limited to 2 CPUs
    with monkeypatch_function_result(os.path.exists, "/sys/fs/cgroup/cgroup.controllers", patch_result=False, patch_module=os.path), \
         monkeypatch_function_result(util.file.slurp_file, '/sys/fs/cgroup/cpu/cpu.cfs_quota_us', patch_result='200000'), \
         monkeypatch_function_result(util.file.slurp_file, '/sys/fs/cgroup/cpu/cpu.cfs_period_us', patch_result='100000'):
         
        assert util.misc.available_cpu_count() == 2

    # cgroup v1 limited to 1 CPU
    with monkeypatch_function_result(os.path.exists, "/sys/fs/cgroup/cgroup.controllers", patch_result=False, patch_module=os.path), \
         monkeypatch_function_result(util.file.slurp_file, '/sys/fs/cgroup/cpu/cpu.cfs_quota_us', patch_result='1'), \
         monkeypatch_function_result(util.file.slurp_file, '/sys/fs/cgroup/cpu/cpu.cfs_period_us', patch_result='1'):
         
        assert util.misc.available_cpu_count() == 1

    # cgroup v1 with no limit imposed on the cgroup
    # (fall back to /proc/self/status method, with limit imposed there):
    #   'Cpus_allowed:  c' = 0b1100 bitmask (meaning execution allowed on 2 CPUs)
    with monkeypatch_function_result(os.path.exists, "/sys/fs/cgroup/cgroup.controllers", patch_result=False, patch_module=os.path), \
         monkeypatch_function_result(util.file.slurp_file, '/sys/fs/cgroup/cpu/cpu.cfs_quota_us', patch_result='-1'), \
         monkeypatch_function_result(util.file.slurp_file, '/sys/fs/cgroup/cpu/cpu.cfs_period_us', patch_result='1'), \
         monkeypatch_function_result(util.file.slurp_file, '/proc/self/status', patch_result='Cpus_allowed:  c'):
         
        assert util.misc.available_cpu_count() == 2

    # cgroup v1 with no limit imposed on the cgoup or via /proc/self/status
    # (fall back to /proc/self/status method, with no limit imposed there)
    with monkeypatch_function_result(os.path.exists, "/sys/fs/cgroup/cgroup.controllers", patch_result=False, patch_module=os.path), \
         monkeypatch_function_result(util.file.slurp_file, '/sys/fs/cgroup/cpu/cpu.cfs_quota_us', patch_result='-1'), \
         monkeypatch_function_result(util.file.slurp_file, '/sys/fs/cgroup/cpu/cpu.cfs_period_us', patch_result='1'):
        
        assert util.misc.available_cpu_count() == reported_cpu_count

    # cgroup v1 with no limit imposed on the cgoup
    # with 'Cpus_allowed' not present in /proc/self/status
    # (fall back to multiprocessing.cpu_count() method)
    with monkeypatch_function_result(os.path.exists, "/sys/fs/cgroup/cgroup.controllers", patch_result=False, patch_module=os.path), \
         monkeypatch_function_result(util.file.slurp_file, '/sys/fs/cgroup/cpu/cpu.cfs_quota_us', patch_result='-1'), \
         monkeypatch_function_result(util.file.slurp_file, '/sys/fs/cgroup/cpu/cpu.cfs_period_us', patch_result='1'), \
         monkeypatch_function_result(util.file.slurp_file, '/proc/self/status', patch_result='unexpected_key:  1'):
        
        assert util.misc.available_cpu_count() == reported_cpu_count

    # cgroup v1 with no limit imposed on the cgoup
    # with 'Cpus_allowed' not present in /proc/self/status
    # (fall back to multiprocessing.cpu_count() method with CPU count of 2 reported)
    with monkeypatch_function_result(os.path.exists, "/sys/fs/cgroup/cgroup.controllers", patch_result=False, patch_module=os.path), \
         monkeypatch_function_result(util.file.slurp_file, '/sys/fs/cgroup/cpu/cpu.cfs_quota_us', patch_result='-1'), \
         monkeypatch_function_result(util.file.slurp_file, '/sys/fs/cgroup/cpu/cpu.cfs_period_us', patch_result='1'), \
         monkeypatch_function_result(util.file.slurp_file, '/proc/self/status', patch_result='unexpected_key:  1'), \
         monkeypatch_function_result(multiprocessing.cpu_count, patch_result=2, patch_module=multiprocessing):
        
        assert util.misc.available_cpu_count() == 2