# Unit tests for diamond

from builtins import super
import collections
import glob
import mock
import os.path
import tempfile
import unittest
import util.file
import util.misc
import tools.diamond
import test


CompletedProcess = collections.namedtuple(
    'CompletedProcess', ['args', 'returncode', 'stdout', 'stderr'])


class TestToolDiamond(test.TestCaseWithTmp):

    def setUp(self):
        super().setUp()

        patcher = mock.patch('subprocess.check_call', autospec=True)
        self.addCleanup(patcher.stop)
        self.mock_run = patcher.start()
        self.mock_run.return_value = CompletedProcess([], 0, '', '')

        patcher = mock.patch('tools.CondaPackage', autospec=True)
        self.addCleanup(patcher.stop)
        self.mock_conda = patcher.start()
        self.mock_conda.return_value.verify_install.return_value = "mock"
        self.mock_conda.return_value.is_attempted.return_value = True
        self.mock_conda.return_value.is_installed.return_value = True
        self.mock_conda.return_value.require_executability = False
        self.mock_conda.return_value.executable_path.return_value = "/dev/null"

        self.diamond = tools.diamond.Diamond()
        self.diamond.install()
        self.data_dir = os.path.join(util.file.get_test_input_path(), 'TestMetagenomicsSimple')
        self.db_dir = os.path.join(self.data_dir, 'db')

    def test_create(self):
        db = os.path.join(tempfile.tempdir, 'fake.dmnd')

        # To be replaced with recursive glob in Python 3.5.
        protein_fastas = glob.glob('{}/library/*/*/*.faa'.format(self.db_dir))

        self.diamond.build(db, protein_fastas)
        args = self.mock_run.call_args[0][0]
        self.assertIn('makedb', args)
        self.assertIn('--db', args)
        self.assertIn(db, args)
