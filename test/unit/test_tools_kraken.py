# Unit tests for kraken
from builtins import super
import os.path
import tempfile
import unittest
import mock
import subprocess
import util.file
import util.misc
import tools.kraken
from test import TestCaseWithTmp


class TestToolKrakenMocked(TestCaseWithTmp):

    def setUp(self):
        super().setUp()
        patcher = mock.patch('util.misc.run_and_print', autospec=True)
        self.addCleanup(patcher.stop)
        self.mock_run_and_print = patcher.start()
        patcher = mock.patch('util.misc.run', autospec=True)
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
        out_reads = util.file.mkstempfname('.reads')
        self.kraken.classify(self.inBam, self.db, out_reads)
        args = self.mock_run_and_print.call_args[0][0]
        self.assertEqual('kraken', os.path.basename(args[0]))
        self.assertTrue(util.misc.list_contains(['--db', self.db], args), args)
        self.assertTrue(util.misc.list_contains(['--output', out_reads], args), args)
        self.assertTrue(util.misc.list_contains(['--threads', str(util.misc.available_cpu_count())], args), args)
        self.assertIn('--paired', args)

    def test_kraken_filter(self):
        in_reads = util.file.mkstempfname('.kraken_reads.unfilt.txt')
        out_reads = util.file.mkstempfname('.kraken_reads.filt.txt')
        for thresh in (0.05, 0.3, 0.81):
            self.kraken.filter(in_reads, self.db, out_reads, thresh)
            args = self.mock_run.call_args[0][0]
            self.assertEqual('kraken-filter', os.path.basename(args[0]))
            self.assertIn(in_reads, args)
            self.assertTrue(util.misc.list_contains(['--db', self.db], args), args)
            self.assertTrue(util.misc.list_contains(['--threshold', str(thresh)], args), args)

    def test_kraken_report(self):
        in_reads = util.file.mkstempfname('.kraken_reads.txt')
        out_report = util.file.mkstempfname('.kraken_report.txt')
        self.kraken.report(in_reads, self.db, out_report)
        args = self.mock_run.call_args[0][0]
        self.assertEqual('kraken-report', os.path.basename(args[0]))
        self.assertIn(in_reads, args)
        self.assertTrue(util.misc.list_contains(['--db', self.db], args), args)

    def test_classify_num_threads(self):
        out_reads = util.file.mkstempfname('.reads')
        
        self.kraken.classify(self.inBam, self.db, out_reads)
        args = self.mock_run_and_print.call_args[0][0]
        self.assertEqual('kraken', os.path.basename(args[0]))
        self.assertIn('--threads', args)
        actual = args[args.index('--threads')+1]
        self.assertEqual(actual, str(util.misc.available_cpu_count()), args)

        for requested in (1,2,3,8,11,20):
            expected = min(util.misc.available_cpu_count(), requested)
            self.kraken.classify(self.inBam, self.db, out_reads, numThreads=requested)
            args = self.mock_run_and_print.call_args[0][0]
            self.assertEqual('kraken', os.path.basename(args[0]))
            self.assertIn('--threads', args)
            actual = args[args.index('--threads')+1]
            self.assertEqual(actual, str(expected), "failure for requested %s, expected %s, actual %s" % (requested, expected, actual))


class TestToolKrakenExecute(TestCaseWithTmp):

    def setUp(self):
        super().setUp()
        input_dir = os.path.join(util.file.get_test_input_path(), 'TestMetagenomicsViralMix')
        self.kraken_db_viral_mix = tempfile.mkdtemp('-kraken_db_viral_mix')
        cmd = ['tar', '-C', self.kraken_db_viral_mix, '-xzf', os.path.join(input_dir, 'kraken_db.tar.gz')
        subprocess.check_call(cmd)

    def test_kraken_execution(self):
        input_bam = os.path.join(util.file.get_test_input_path(), 'TestMetagenomicsViralMix', 'test-reads.bam')
        outdir = tempfile.mkdtemp('-kraken_execute')
        out = os.path.join(outdir, 'viral_mix.kraken')
        out_filtered = os.path.join(outdir, 'viral_mix.filtered-kraken')
        out_report = os.path.join(outdir, 'viral_mix.kraken-report')
        self.assertEqual(0, kraken.classify(input_bam, self.kraken_db_viral_mix, out))
        self.assertEqual(0, kraken.filter(out, self.kraken_db_viral_mix, out_filtered, 0.05))
        self.assertEqual(0, kraken.report(out_filtered, self.kraken_db_viral_mix, out_report))
        self.assertGreater(os.path.getsize(out_report), 0)
        self.assertGreater(os.path.getsize(out_filtered), 0)

    def test_kraken_on_empty(self):
        input_bam = os.path.join(util.file.get_test_input_path(), 'empty.bam')
        outdir = tempfile.mkdtemp('-kraken_empty')
        out = os.path.join(outdir, 'empty.kraken')
        out_filtered = os.path.join(outdir, 'empty.filtered-kraken')
        out_report = os.path.join(outdir, 'empty.kraken-report')
        self.assertEqual(0, kraken.classify(input_bam, self.kraken_db_viral_mix, out))
        self.assertEqual(0, kraken.filter(out, self.kraken_db_viral_mix, out_filtered, 0.05))
        self.assertEqual(0, kraken.report(out_filtered, self.kraken_db_viral_mix, out_report))
        self.assertTrue(os.path.isfile(out_report))
        self.assertTrue(os.path.isfile(out_filtered))

    def test_kraken_on_almost_empty(self):
        input_bam = os.path.join(util.file.get_test_input_path(), 'almost-empty.bam')
        outdir = tempfile.mkdtemp('-kraken_almost_empty')
        out = os.path.join(outdir, 'almost_empty.kraken')
        out_filtered = os.path.join(outdir, 'almost_empty.filtered-kraken')
        out_report = os.path.join(outdir, 'almost_empty.kraken-report')
        self.assertEqual(0, kraken.classify(input_bam, self.kraken_db_viral_mix, out))
        self.assertEqual(0, kraken.filter(out, self.kraken_db_viral_mix, out_filtered, 0.05))
        self.assertEqual(0, kraken.report(out_filtered, self.kraken_db_viral_mix, out_report))
        self.assertTrue(os.path.isfile(out_report))
        self.assertTrue(os.path.isfile(out_filtered))
