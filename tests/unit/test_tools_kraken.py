# Unit tests for kraken
from builtins import super
import os.path
import tempfile
import unittest
from mock import patch
import util.file
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

        patcher = patch('util.misc.run', new=self.mock_run)

        self.inBam = util.file.mkstempfname('.bam')
        self.db = tempfile.mkdtemp('db')

    def test_kraken_classify(self):
        in_fastq1 = util.file.mkstempfname('.1.fastq')
        in_fastq2 = util.file.mkstempfname('.2.fastq')
        out_reads = util.file.mkstempfname('.reads')
        self.kraken.classify(self.db, [in_fastq1, in_fastq2], out_reads, options={
            '--paired': None
        })
        called_cmd = self.mock_run.call_args[0][0]
        self.assertIn('--db', called_cmd[1:])
        self.assertIn(self.db, called_cmd[1:])
        self.assertIn('--output', called_cmd[1:])
        self.assertIn(out_reads, called_cmd[1:])
        self.assertIn(in_fastq1, called_cmd[1:])
        self.assertIn(in_fastq2, called_cmd[1:])
        self.assertIn('--paired', called_cmd[1:])
        self.assertEqual('kraken', os.path.basename(called_cmd[0]))


if __name__ == '__main__':
    unittest.main()
