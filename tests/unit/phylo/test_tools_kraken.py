# Unit tests for kraken
from builtins import super
import os.path
import tempfile
import unittest
from mock import patch
import util.file
import util.misc
import tools.kraken
from test import TestCaseWithTmp


class TestToolKraken(TestCaseWithTmp):

    def setUp(self):
        super().setUp()
        self.kraken = tools.kraken.Kraken()
        self.kraken.install()

        patcher = patch('util.misc.run_and_print', autospec=True)
        self.addCleanup(patcher.stop)
        self.mock_run = patcher.start()

        patcher = patch('tools.kraken.Jellyfish', autospec=True)
        self.addCleanup(patcher.stop)
        self.mock_jelly = patcher.start()
        self.mock_jelly().install_and_get_path.return_value = tempfile.tempdir

        self.inBam = util.file.mkstempfname('.bam')
        self.db = tempfile.mkdtemp('db')

    def test_kraken_classify(self):
        in_fastq1 = util.file.mkstempfname('.1.fastq')
        in_fastq2 = util.file.mkstempfname('.2.fastq')
        out_reads = util.file.mkstempfname('.reads')
        self.kraken.classify(self.db, [in_fastq1, in_fastq2], out_reads, options={
            '--paired': None
        })
        args = self.mock_run.call_args[0][0]
        self.assertEqual('kraken', os.path.basename(args[0]))
        self.assertTrue(util.misc.list_contains(['--db', self.db], args))
        self.assertTrue(util.misc.list_contains(['--output', out_reads], args))
        self.assertIn('--paired', args)
        self.assertIn(in_fastq1, args)
        self.assertIn(in_fastq2, args)
