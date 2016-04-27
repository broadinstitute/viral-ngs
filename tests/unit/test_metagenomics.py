# Unit tests for metagenomics.py
from builtins import super
import six
import argparse
import os.path
import tempfile
import unittest

import mock
from mock import patch

import tools.picard
import metagenomics
import util.file
from test import TestCaseWithTmp

class TestCommandHelp(unittest.TestCase):

    def test_help_parser_for_each_command(self):
        for cmd_name, parser_fun in metagenomics.__commands__:
            print(cmd_name)
            parser = parser_fun(argparse.ArgumentParser())
            helpstring = parser.format_help()


class TestKrakenCalls(TestCaseWithTmp):

    def setUp(self):
        super().setUp()
        patcher = patch('tools.picard.SamToFastqTool', autospec=True)
        self.addCleanup(patcher.stop)
        self.mock_samtofastq = patcher.start()

        patcher = patch('tools.kraken.Kraken', autospec=True)
        self.addCleanup(patcher.stop)
        self.mock_kraken = patcher.start()

        self.inBam = util.file.mkstempfname('.bam')
        self.db = tempfile.mkdtemp('db')

    def test_assert(self):
        with self.assertRaises(AssertionError):
            metagenomics.kraken(self.inBam, self.db)

    def test_out_reads(self):
        out_reads = util.file.mkstempfname('.kraken_reads.gz')
        metagenomics.kraken(self.inBam, self.db, outReads=out_reads)
        self.mock_samtofastq().execute.assert_called_once_with(
            self.inBam, mock.ANY, mock.ANY, picardOptions=mock.ANY,
            JVMmemory=mock.ANY)
        picard_opts = self.mock_samtofastq().execute.call_args[1]['picardOptions']
        six.assertCountEqual(
            self, picard_opts,
            ['CLIPPING_ACTION=X', 'CLIPPING_ATTRIBUTE=%s'
             % tools.picard.SamToFastqTool.illumina_clipping_attribute])

        self.mock_kraken().classify.assert_called_once_with(
            self.db, mock.ANY, mock.ANY, options={
                '--paired': None,
                '--threads': 1,
        })

        # Exists due to test. Need a better tmpfile patcher.
        self.assertTrue(os.path.isfile(out_reads))

    def test_filter_threshold(self):
        out_reads = util.file.mkstempfname('.kraken_reads.gz')
        metagenomics.kraken(self.inBam, self.db, outReads=out_reads, filterThreshold=0.05)
        self.mock_kraken().execute.assert_called_with(
            'kraken-filter', self.db, mock.ANY, args=[mock.ANY], options={
                '--threshold': 0.05
        })
        self.mock_kraken().classify.assert_called_once_with(
            self.db, mock.ANY, mock.ANY, options={
                '--paired': None,
                '--threads': 1,
        })

    def test_out_report(self):
        out_report = util.file.mkstempfname('.kraken_report.txt')
        metagenomics.kraken(self.inBam, self.db, outReport=out_report)

        self.mock_kraken().classify.assert_called_once_with(
            self.db, mock.ANY, mock.ANY, options={
                '--paired': None,
                '--threads': 1,
        })

        self.mock_kraken().execute.assert_called_with(
            'kraken-report', self.db, out_report, args=[mock.ANY])

    def test_out_reads_and_report(self):
        out_reads = util.file.mkstempfname('.kraken_reads.gz')
        out_report = util.file.mkstempfname('.kraken_report.txt')
        metagenomics.kraken(self.inBam, self.db, outReads=out_reads, outReport=out_report)
        self.mock_kraken().classify.assert_called_once_with(
            self.db, mock.ANY, mock.ANY, options={
                '--paired': None,
                '--threads': 1,
        })
        self.mock_kraken().execute.assert_called_with(
            'kraken-report', self.db, out_report, args=[mock.ANY])

    def test_num_threads(self):
        out_reads = util.file.mkstempfname('.kraken_reads.gz')
        metagenomics.kraken(self.inBam, self.db, outReads=out_reads, numThreads=11)
        self.mock_kraken().classify.assert_called_once_with(
            self.db, mock.ANY, mock.ANY, options={
                '--paired': None,
                '--threads': 11,
        })


class TestDiamondCalls(TestCaseWithTmp):
    def setUp(self):
        super().setUp()
        patcher = patch('tools.picard.SamToFastqTool', autospec=True)
        self.addCleanup(patcher.stop)
        self.mock_samtofastq = patcher.start()

        patcher = patch('tools.diamond.Diamond', autospec=True)
        self.addCleanup(patcher.stop)
        self.mock_diamond = patcher.start()

        self.inBam = util.file.mkstempfname('.bam')
        self.db = tempfile.mkdtemp('db')

    def test_output_m8(self):
        out_m8 = util.file.mkstempfname('.')
        metagenomics.diamond(self.inBam, self.db, out_m8)
        self.mock_samtofastq().execute.assert_called_once_with(
            self.inBam, mock.ANY, mock.ANY, picardOptions=mock.ANY,
            JVMmemory=mock.ANY)
        picard_opts = self.mock_samtofastq().execute.call_args[1]['picardOptions']
        six.assertCountEqual(
            self, picard_opts,
            ['CLIPPING_ACTION=X', 'CLIPPING_ATTRIBUTE=%s'
             % tools.picard.SamToFastqTool.illumina_clipping_attribute])
        self.mock_diamond().install.assert_called_once_with()
        self.mock_diamond().blastx.assert_called_once_with(
            self.db, mock.ANY, mock.ANY, options={'threads': 1})
        self.assertTrue(self.mock_diamond().view.called)

    def test_num_threads(self):
        out_m8 = util.file.mkstempfname('.')
        metagenomics.diamond(self.inBam, self.db, out_m8, numThreads=11)
        self.mock_diamond().blastx.assert_called_once_with(
            self.db, mock.ANY, mock.ANY, options={'threads': 11})
        self.assertEqual(self.mock_diamond().view.call_args[1]['options']['threads'], 11)

class TestKronaCalls(TestCaseWithTmp):

    def setUp(self):
        super().setUp()
        patcher = patch('tools.krona.Krona', autospec=True)
        self.addCleanup(patcher.stop)
        self.mock_krona = patcher.start()

        self.inTsv = util.file.mkstempfname('.tsv')
        self.db = tempfile.mkdtemp('db')

    def test_krona_import_taxonomy(self):
        out_html = util.file.mkstempfname('.html')
        metagenomics.krona(self.inTsv, out_html, queryColumn=3, taxidColumn=5, scoreColumn=7,
                           noHits=True, noRank=True, db=self.db)
        self.mock_krona().import_taxonomy.assert_called_once_with(
            [self.inTsv], out_html, query_column=3, taxid_column=5, score_column=7,
            no_hits=True, no_rank=True, db=self.db)
