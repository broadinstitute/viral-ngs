# Unit tests for kraken
from builtins import super
import os.path
import tempfile
import unittest
import mock
import util.file
import util.misc
import tools.kraken
from test import TestCaseWithTmp


class TestToolKraken(TestCaseWithTmp):

    def setUp(self):
        super().setUp()
        patcher = mock.patch('util.misc.run_and_print', autospec=True)
        self.addCleanup(patcher.stop)
        self.mock_run = patcher.start()

        patcher = mock.patch('tools.CondaPackage', autospec=True)
        self.addCleanup(patcher.stop)
        self.mock_conda = patcher.start()
        self.mock_conda.return_value.verify_install.return_value = "mock"
        self.mock_conda.return_value.is_attempted.return_value = True
        self.mock_conda.return_value.is_installed.return_value = True
        self.mock_conda.return_value.require_executability = False
        self.mock_conda.return_value.executable_path.return_value = "/dev/null"

        self.kraken = tools.kraken.Kraken()
        self.inBam = util.file.mkstempfname('.bam')
        self.db = tempfile.mkdtemp('db')

    def test_kraken_classify(self):
        in_bam = util.file.mkstempfname('.bam')
        out_reads = util.file.mkstempfname('.reads')
        self.kraken.classify(in_bam, self.db, out_reads)
        args = self.mock_run.call_args[0][0]
        self.assertEqual('kraken', os.path.basename(args[0]))
        self.assertTrue(util.misc.list_contains(['--db', self.db], args), args)
        self.assertTrue(util.misc.list_contains(['--output', out_reads], args), args)
        self.assertTrue(util.misc.list_contains(['--threads', str(util.misc.available_cpu_count())], args), args)
        self.assertIn('--paired', args)

