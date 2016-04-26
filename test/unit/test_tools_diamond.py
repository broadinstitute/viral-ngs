# Unit tests for diamond

from builtins import super
import collections
import glob
import os.path
from os.path import join
import tempfile
import unittest
import six
from mock import patch
from util import file
import util.misc
import tools.diamond
from test import TestCaseWithTmp


CompletedProcess = collections.namedtuple(
    'CompletedProcess', ['args', 'returncode', 'stdout', 'stderr'])


class TestToolDiamond(TestCaseWithTmp):

    def setUp(self):
        super().setUp()
        self.diamond = tools.diamond.Diamond()
        self.diamond.install()

        patcher = patch('util.misc.run_and_print', autospec=True)
        self.addCleanup(patcher.stop)
        self.mock_run = patcher.start()
        self.mock_run.return_value = CompletedProcess([], 0, '', '')

        self.data_dir = util.file.get_test_input_path(self)
        self.db_dir = join(self.data_dir, 'db')

    def test_create_and_blastx(self):
        db = join(tempfile.tempdir, 'fake.dmnd')

        # To be replaced with recursive glob in Python 3.5.
        protein_fastas = glob.glob('{}/library/*/*/*.ffn'.format(self.db_dir))

        self.assertEqual(0, self.diamond.build(db, protein_fastas).returncode)
        args = self.mock_run.call_args[0][0]
        self.assertIn('makedb', args)
        self.assertIn('--db', args)
        self.assertIn(db, args)

        inputs = [join(self.data_dir, f)
                  for f in ['zaire_ebola.1.fastq', 'zaire_ebola.2.fastq']]

        output = join(tempfile.tempdir, 'zaire_ebola.daa')
        output_tab = join(tempfile.tempdir, 'zaire_ebola.m8')

        tmpdir = tempfile.gettempdir()
        self.assertEqual(
            0, self.diamond.blastx(db, inputs, output,
                                   options={'--tmpdir': tmpdir}).returncode)
        args = self.mock_run.call_args[0][0]
        self.assertIn('blastx', args)
        self.assertTrue(util.misc.list_contains(['--daa', output], args))
        self.assertTrue(util.misc.list_contains(['--db', db], args))
        self.assertTrue(util.misc.list_contains(['--tmpdir', tmpdir], args))

        self.assertEqual(0, self.diamond.view(output, output_tab).returncode)
        args = self.mock_run.call_args[0][0]
        self.assertIn('view', args)
        self.assertTrue(util.misc.list_contains(['--outfmt', 'tab'], args))
        self.assertTrue(util.misc.list_contains(['--daa', output], args))
        self.assertTrue(util.misc.list_contains(['--out', output_tab], args))


if __name__ == '__main__':
    unittest.main()
