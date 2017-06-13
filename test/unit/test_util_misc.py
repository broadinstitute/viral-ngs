# Unit tests for util.misc.py

__author__ = "dpark@broadinstitute.org"

import os, random, collections
import unittest
import subprocess
import util.misc
import util.file


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


class TestDictFlatten(unittest.TestCase):

    def testDictFlatten(self):
        z=dict(a=1,b=dict(c=2,e=3))
        self.assertEqual(util.misc.unflatten_dict(util.misc.flatten_dict(z)), z)

        self.assertEqual(util.misc.flatten_dict({}),{})
        self.assertEqual(util.misc.flatten_dict({'a':1}),{('a',): 1})
        self.assertEqual(util.misc.flatten_dict({'a':1,'b':2}),{('a',): 1, ('b',):2})
        self.assertEqual(util.misc.flatten_dict({'a':1,'b':{'c':2}}),{('a',): 1, ('b','c'):2})
        self.assertEqual(util.misc.flatten_dict({'a':1,'b':{'c':2,'d':{'e':3}}}),
                         {('a',): 1, ('b','c'):2, ('b','d','e'):3})

        # test random dicts

        def rand_dict(max_depth=5):
            z=dict()
            for i in range(random.randint(1,3)):
                z[str(random.randint(0,5))] = random.randint(0,100) if max_depth==0 or random.randint(0,1) \
                                              else rand_dict(max_depth-1)
            return z

        for randTestNum in range(1000):
            z = rand_dict()
            z_flat=util.misc.flatten_dict(z)
            assert all( isinstance(k,tuple) and isinstance(v,int) for k,v in z_flat.items())
            z_unflat=util.misc.unflatten_dict(z_flat)
            z_reflat=util.misc.flatten_dict(z_unflat)
            
            self.assertEqual(z_unflat, z)
            self.assertEqual(z_reflat, z_flat)

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

        self.assertEqual(cfg2["std_methods"], ['a','b','c'])
        self.assertEqual(cfg2["stage1"]["stage2"]["step_num"],3)
        self.assertEqual(cfg2["stage1"]["stage2"]["step_list"],[5,10,15])
        self.assertEqual(cfg2["stage1"]["stage3"]["step_list"],[3,33])
        self.assertEqual(cfg1["stage1"]["stage3"]["step_list"],[51,101,151])

        self.assertEqual(cfg2["std_param_A_new"], 111)  # specified as std_param_A_old in cfg1.yaml




        


