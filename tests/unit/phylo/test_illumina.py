# Unit tests for illumina.py

__author__ = "dpark@broadinstitute.org"

import unittest
import os
import os.path
import tempfile
import argparse
import filecmp
import util
import util.file
import illumina
from test import TestCaseWithTmp


class TestCommandHelp(unittest.TestCase):

    def test_help_parser_for_each_command(self):
        for cmd_name, parser_fun in illumina.__commands__:
            parser = parser_fun(argparse.ArgumentParser())
            helpstring = parser.format_help()




class TestSampleSheet(TestCaseWithTmp):

    def test_miseq(self):
        inDir = util.file.get_test_input_path(self)
        samples = illumina.SampleSheet(os.path.join(inDir, 'SampleSheet-miseq-1.csv'))
        self.assertEqual(samples.num_indexes(), 2)
        self.assertEqual(len(samples.get_rows()), 20)
    
    def test_broad_platform(self):
        inDir = util.file.get_test_input_path(self)
        samples = illumina.SampleSheet(os.path.join(inDir, 'SampleSheet-hiseq-1.csv'), only_lane=2)
        self.assertEqual(samples.num_indexes(), 2)
        self.assertEqual(len(samples.get_rows()), 24)
        samples = illumina.SampleSheet(os.path.join(inDir, 'SampleSheet-hiseq-1.csv'), allow_non_unique=True)
        self.assertEqual(len(samples.get_rows()), 48)

    def test_walkup_submission_no_header_no_lf(self):
        inDir = util.file.get_test_input_path(self)
        samples = illumina.SampleSheet(os.path.join(inDir, 'SampleSheet-submit-1.csv'))
        self.assertEqual(samples.num_indexes(), 2)
        self.assertEqual(len(samples.get_rows()), 24)

    def test_walkup_submission(self):
        inDir = util.file.get_test_input_path(self)
        samples = illumina.SampleSheet(os.path.join(inDir, 'SampleSheet-submit-2.csv'))
        self.assertEqual(samples.num_indexes(), 2)
        self.assertEqual(len(samples.get_rows()), 96)

    def test_walkup_submission_no_lf(self):
        inDir = util.file.get_test_input_path(self)
        samples = illumina.SampleSheet(os.path.join(inDir, 'SampleSheet-submit-3.csv'))
        self.assertEqual(samples.num_indexes(), 2)
        self.assertEqual(len(samples.get_rows()), 24)
    
    def test_tabfile(self):
        inDir = util.file.get_test_input_path(self)
        samples = illumina.SampleSheet(os.path.join(inDir, 'SampleSheet-custom-1.txt'))
        self.assertEqual(samples.num_indexes(), 2)
        self.assertEqual(len(samples.get_rows()), 24)


class TestRunInfo(TestCaseWithTmp):

    def test_miseq(self):
        inDir = util.file.get_test_input_path(self)
        runinfo = illumina.RunInfo(os.path.join(inDir, 'RunInfo-miseq.xml'))
        self.assertEqual(runinfo.get_flowcell(), 'AEF96')
        self.assertEqual(runinfo.get_rundate_american(), '08/21/2015')
        self.assertEqual(runinfo.get_rundate_iso(), '2015-08-21')
        self.assertEqual(runinfo.get_machine(), 'M04004')
        self.assertEqual(runinfo.get_read_structure(), '101T8B8B101T')
    
    def test_hiseq(self):
        inDir = util.file.get_test_input_path(self)
        runinfo = illumina.RunInfo(os.path.join(inDir, 'RunInfo-hiseq.xml'))
        self.assertEqual(runinfo.get_flowcell(), 'HVFF2ADXX')
        self.assertEqual(runinfo.get_rundate_american(), '08/21/2015')
        self.assertEqual(runinfo.get_rundate_iso(), '2015-08-21')
        self.assertEqual(runinfo.get_machine(), 'SL-HDF')
        self.assertEqual(runinfo.get_read_structure(), '101T8B8B101T')


class TestIlluminaDir(TestCaseWithTmp):

    def test_directory(self):
        inDir = util.file.get_test_input_path(self)
        test_in = os.path.join(inDir, 'empty_dir')
        with illumina.IlluminaDirectory(test_in) as idir:
            self.assertTrue(os.path.isdir(idir.get_BCLdir()))
    
    def test_tarball_normal(self):
        inDir = util.file.get_test_input_path(self)
        with illumina.IlluminaDirectory(os.path.join(inDir, 'bcl-plain.tgz')) as idir:
            self.assertTrue(os.path.isdir(idir.get_BCLdir()))
        with illumina.IlluminaDirectory(os.path.join(inDir, 'bcl-plain.tar.bz2')) as idir:
            self.assertTrue(os.path.isdir(idir.get_BCLdir()))
    
    def test_tarball_indented(self):
        inDir = util.file.get_test_input_path(self)
        with illumina.IlluminaDirectory(os.path.join(inDir, 'bcl-indented.tgz')) as idir:
            self.assertTrue(os.path.isdir(idir.get_BCLdir()))

    def test_tarball_sample_sheet(self):
        inDir = util.file.get_test_input_path(self)
        with illumina.IlluminaDirectory(os.path.join(inDir, 'bcl-samplesheet.tar.gz')) as idir:
            self.assertTrue(os.path.isdir(idir.get_BCLdir()))
            self.assertEqual(len(idir.get_SampleSheet().get_rows()), 15)
        with illumina.IlluminaDirectory(os.path.join(inDir, 'bcl-both.tar.gz')) as idir:
            self.assertTrue(os.path.isdir(idir.get_BCLdir()))
            self.assertEqual(len(idir.get_SampleSheet().get_rows()), 15)

    def test_tarball_run_info(self):
        inDir = util.file.get_test_input_path(self)
        with illumina.IlluminaDirectory(os.path.join(inDir, 'bcl-runinfo.tar.gz')) as idir:
            self.assertTrue(os.path.isdir(idir.get_BCLdir()))
            self.assertEqual(idir.get_RunInfo().get_flowcell(), 'AHVPA')
            self.assertEqual(idir.get_RunInfo().get_machine(), 'M04004')
            self.assertEqual(idir.get_RunInfo().get_rundate_iso(), '2015-08-27')
        with illumina.IlluminaDirectory(os.path.join(inDir, 'bcl-both.tar.gz')) as idir:
            self.assertTrue(os.path.isdir(idir.get_BCLdir()))
            self.assertEqual(idir.get_RunInfo().get_flowcell(), 'AHVPA')
            self.assertEqual(idir.get_RunInfo().get_read_structure(), '101T8B8B101T')

    def test_tarball_fail_missing_data(self):
        inDir = util.file.get_test_input_path(self)
        with illumina.IlluminaDirectory(os.path.join(inDir, 'bcl-runinfo.tar.gz')) as idir:
            self.assertRaises(Exception, idir.get_SampleSheet())
        with illumina.IlluminaDirectory(os.path.join(inDir, 'bcl-samplesheet.tar.gz')) as idir:
            self.assertRaises(Exception, idir.get_RunInfo())


class TestMiseqToBam(TestCaseWithTmp):
    
    def test_paired(self):
        pass
        



if __name__ == '__main__':
    unittest.main()
