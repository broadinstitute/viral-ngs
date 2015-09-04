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
import tools.samtools
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
        self.assertEqual(runinfo.num_reads(), 2)
    
    def test_hiseq(self):
        inDir = util.file.get_test_input_path(self)
        runinfo = illumina.RunInfo(os.path.join(inDir, 'RunInfo-hiseq.xml'))
        self.assertEqual(runinfo.get_flowcell(), 'HVFF2ADXX')
        self.assertEqual(runinfo.get_rundate_american(), '08/21/2015')
        self.assertEqual(runinfo.get_rundate_iso(), '2015-08-21')
        self.assertEqual(runinfo.get_machine(), 'SL-HDF')
        self.assertEqual(runinfo.get_read_structure(), '101T8B8B101T')
        self.assertEqual(runinfo.num_reads(), 2)


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
    
    def test_paired_1(self):
        inDir = util.file.get_test_input_path(self)
        outBam = util.file.mkstempfname('.bam')
        outHeader = util.file.mkstempfname('.txt')
        sampleSheet = os.path.join(inDir, 'SampleSheet.csv')
        runInfo = os.path.join(inDir, 'RunInfo.xml')
        fastq = (os.path.join(inDir, 'mebv-0-1_S5_L001_R1_001.fastq.gz'),
                 os.path.join(inDir, 'mebv-0-1_S5_L001_R2_001.fastq.gz'))
        illumina.miseq_fastq_to_bam(outBam, sampleSheet, fastq[0], inFastq2=fastq[1], runInfo=runInfo)
        tools.samtools.SamtoolsTool().dumpHeader(outBam, outHeader)
        self.assertEqualContents(outHeader, os.path.join(inDir, 'mebv.0.1.bam.header.txt'))
        
    def test_paired_2(self):
        inDir = util.file.get_test_input_path(self)
        outBam = util.file.mkstempfname('.bam')
        outHeader = util.file.mkstempfname('.txt')
        sampleSheet = os.path.join(inDir, 'SampleSheet.csv')
        runInfo = os.path.join(inDir, 'RunInfo.xml')
        fastq = (os.path.join(inDir, 'mebv-48-5_S17_L001_R1_001.fastq.gz'),
                 os.path.join(inDir, 'mebv-48-5_S17_L001_R2_001.fastq.gz'))
        illumina.miseq_fastq_to_bam(outBam, sampleSheet, fastq[0], inFastq2=fastq[1], runInfo=runInfo)
        tools.samtools.SamtoolsTool().dumpHeader(outBam, outHeader)
        self.assertEqualContents(outHeader, os.path.join(inDir, 'mebv.48.5.bam.header.txt'))

    def test_paired_custom_seq_center(self):
        inDir = util.file.get_test_input_path(self)
        outBam = util.file.mkstempfname('.bam')
        outHeader = util.file.mkstempfname('.txt')
        sampleSheet = os.path.join(inDir, 'SampleSheet.csv')
        runInfo = os.path.join(inDir, 'RunInfo.xml')
        fastq = (os.path.join(inDir, 'mebv-48-5_S17_L001_R1_001.fastq.gz'),
                 os.path.join(inDir, 'mebv-48-5_S17_L001_R2_001.fastq.gz'))
        illumina.miseq_fastq_to_bam(outBam, sampleSheet, fastq[0], inFastq2=fastq[1], runInfo=runInfo, sequencing_center='CustomSeqCenter')
        tools.samtools.SamtoolsTool().dumpHeader(outBam, outHeader)
        self.assertEqualContents(outHeader, os.path.join(inDir, 'mebv.48.5.custom.bam.header.txt'))

    def test_fail_missing_pair(self):
        inDir = util.file.get_test_input_path(self)
        outBam = util.file.mkstempfname('.bam')
        sampleSheet = os.path.join(inDir, 'SampleSheet.csv')
        runInfo = os.path.join(inDir, 'RunInfo.xml')
        fastq = (os.path.join(inDir, 'mebv-48-5_S17_L001_R1_001.fastq.gz'),)
        self.assertRaises(Exception, illumina.miseq_fastq_to_bam, outBam, sampleSheet, fastq[0], runInfo=runInfo)

    def test_fail_backwards_pair(self):
        inDir = util.file.get_test_input_path(self)
        outBam = util.file.mkstempfname('.bam')
        sampleSheet = os.path.join(inDir, 'SampleSheet.csv')
        runInfo = os.path.join(inDir, 'RunInfo.xml')
        fastq = (os.path.join(inDir, 'mebv-48-5_S17_L001_R2_001.fastq.gz'),
                 os.path.join(inDir, 'mebv-48-5_S17_L001_R1_001.fastq.gz'))
        self.assertRaises(Exception, illumina.miseq_fastq_to_bam, outBam, sampleSheet, fastq[0], fastq2=fastq[1], runInfo=runInfo)

    def test_fail_mismatched_pair(self):
        inDir = util.file.get_test_input_path(self)
        outBam = util.file.mkstempfname('.bam')
        sampleSheet = os.path.join(inDir, 'SampleSheet.csv')
        runInfo = os.path.join(inDir, 'RunInfo.xml')
        fastq = (os.path.join(inDir, 'mebv-48-5_S16_L001_R1_001.fastq.gz'),
                 os.path.join(inDir, 'mebv-48-5_S17_L001_R2_001.fastq.gz'))
        self.assertRaises(Exception, illumina.miseq_fastq_to_bam, outBam, sampleSheet, fastq[0], fastq2=fastq[1], runInfo=runInfo)

    def test_fail_oob_index(self):
        inDir = util.file.get_test_input_path(self)
        outBam = util.file.mkstempfname('.bam')
        sampleSheet = os.path.join(inDir, 'SampleSheet.csv')
        runInfo = os.path.join(inDir, 'RunInfo.xml')
        fastq = (os.path.join(inDir, 'mebv-48-5_S33_L001_R1_001.fastq.gz'),
                 os.path.join(inDir, 'mebv-48-5_S33_L001_R2_001.fastq.gz'))
        self.assertRaises(Exception, illumina.miseq_fastq_to_bam, outBam, sampleSheet, fastq[0], fastq2=fastq[1], runInfo=runInfo)

    def test_fail_bad_format(self):
        inDir = util.file.get_test_input_path(self)
        outBam = util.file.mkstempfname('.bam')
        sampleSheet = os.path.join(inDir, 'SampleSheet.csv')
        runInfo = os.path.join(inDir, 'RunInfo.xml')
        fastq = (os.path.join(inDir, 'mebv-48-5_S17_L001_R1_002.fastq.gz'),
                 os.path.join(inDir, 'mebv-48-5_S17_L001_R2_002.fastq.gz'))
        self.assertRaises(Exception, illumina.miseq_fastq_to_bam, outBam, sampleSheet, fastq[0], fastq2=fastq[1], runInfo=runInfo)
        fastq = (os.path.join(inDir, 'mebv-48-5_S17_L002_R1_001.fastq.gz'),
                 os.path.join(inDir, 'mebv-48-5_S17_L002_R2_001.fastq.gz'))
        self.assertRaises(Exception, illumina.miseq_fastq_to_bam, outBam, sampleSheet, fastq[0], fastq2=fastq[1], runInfo=runInfo)
        fastq = (os.path.join(inDir, 'mebv-48-5_17_L001_R1_001.fastq.gz'),
                 os.path.join(inDir, 'mebv-48-5_17_L001_R2_001.fastq.gz'))
        self.assertRaises(Exception, illumina.miseq_fastq_to_bam, outBam, sampleSheet, fastq[0], fastq2=fastq[1], runInfo=runInfo)



if __name__ == '__main__':
    unittest.main()
