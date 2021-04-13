# Unit tests for illumina.py
# -*- coding: utf-8 -*-

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

    def test_tabfile_win_endings(self):
        inDir = util.file.get_test_input_path(self)
        samples = illumina.SampleSheet(os.path.join(inDir, 'SampleSheet-custom-1_win-endings.txt'))
        self.assertEqual(samples.num_indexes(), 2)
        self.assertEqual(len(samples.get_rows()), 24)

    def test_gz_tabfile_win_endings(self):
        inDir = util.file.get_test_input_path(self)
        samples = illumina.SampleSheet(os.path.join(inDir, 'SampleSheet-custom-1_win-endings.txt.gz'))
        self.assertEqual(samples.num_indexes(), 2)
        self.assertEqual(len(samples.get_rows()), 24)

    def test_tabfile_macos9_endings(self):
        inDir = util.file.get_test_input_path(self)
        samples = illumina.SampleSheet(os.path.join(inDir, 'SampleSheet-custom-1_macos9-endings.txt'))
        self.assertEqual(samples.num_indexes(), 2)
        self.assertEqual(len(samples.get_rows()), 24)

    def test_tabfile_extras(self):
        inDir = util.file.get_test_input_path(self)
        samples = illumina.SampleSheet(os.path.join(inDir, 'SampleSheet-custom-2.txt'))
        self.assertEqual(samples.num_indexes(), 2)
        self.assertEqual(len(samples.get_rows()), 96)

    def test_tabfile_extras_win(self):
        inDir = util.file.get_test_input_path(self)
        samples = illumina.SampleSheet(os.path.join(inDir, 'SampleSheet-custom-2_win-endings.tsv'))
        self.assertEqual(samples.num_indexes(), 2)
        self.assertEqual(len(samples.get_rows()), 96)

    def test_blank_line_in_tabular_section(self):
        inDir = util.file.get_test_input_path(self)
        samples = illumina.SampleSheet(os.path.join(inDir, 'SampleSheet-with-blanklines.csv'))
        self.assertEqual(samples.num_indexes(), 2)
        self.assertEqual(len(samples.get_rows()), 12)

    def test_picard_block(self):
        inDir = util.file.get_test_input_path(self)
        samples = illumina.SampleSheet(os.path.join(inDir, 'SampleSheet-in-Broad-MiSeq-Format_with_Picard_Block.csv'))
        self.assertEqual(samples.num_indexes(), 2)
        self.assertEqual(len(samples.get_rows()), 11)


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
        self.assertEqual(runinfo.get_machine_model(), "Illumina MiSeq")
        self.assertEqual(runinfo.get_flowcell_lane_count(), 1)

    def test_hiseq(self):
        inDir = util.file.get_test_input_path(self)
        runinfo = illumina.RunInfo(os.path.join(inDir, 'RunInfo-hiseq.xml'))
        self.assertEqual(runinfo.get_flowcell(), 'HVFF2ADXX')
        self.assertEqual(runinfo.get_rundate_american(), '08/21/2015')
        self.assertEqual(runinfo.get_rundate_iso(), '2015-08-21')
        self.assertEqual(runinfo.get_machine(), 'SL-HDF')
        self.assertEqual(runinfo.get_read_structure(), '101T8B8B101T')
        self.assertEqual(runinfo.num_reads(), 2)
        self.assertEqual(runinfo.get_machine_model(), "Illumina HiSeq 2500")
        self.assertEqual(runinfo.get_flowcell_lane_count(), 2)

    def test_novaseq(self):
        inDir = util.file.get_test_input_path(self)
        runinfo = illumina.RunInfo(os.path.join(inDir, 'RunInfo-novaseq.xml'))
        self.assertEqual(runinfo.get_flowcell(), 'HCYTJDMXX')
        self.assertEqual(runinfo.get_rundate_american(), '06/27/2018')
        self.assertEqual(runinfo.get_rundate_iso(), '2018-06-27')
        self.assertEqual(runinfo.get_machine(), 'A00198')
        self.assertEqual(runinfo.get_read_structure(), '101T8B8B101T')
        self.assertEqual(runinfo.num_reads(), 2)
        self.assertEqual(runinfo.get_machine_model(), "Illumina NovaSeq 6000")
        self.assertEqual(runinfo.get_flowcell_lane_count(), 2)


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
        with illumina.IlluminaDirectory(os.path.join(inDir, 'bcl-plain.tar.lz4')) as idir:
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

    def test_tarball_uncompressed(self):
        inDir = util.file.get_test_input_path(self)
        with illumina.IlluminaDirectory(os.path.join(inDir, 'bcl-both-uncompressed.tar')) as idir:
            self.assertTrue(os.path.isdir(idir.get_BCLdir()))
            self.assertEqual(len(idir.get_SampleSheet().get_rows()), 15)

    def test_tarball_deep_dir_tree(self):
        inDir = util.file.get_test_input_path(self)
        with illumina.IlluminaDirectory(os.path.join(inDir, 'bcl-both-broad_full_path.tar.gz')) as idir:
            self.assertTrue(os.path.isdir(idir.get_BCLdir()))
            self.assertEqual(len(idir.get_SampleSheet().get_rows()), 15)

    def test_zip_archive(self):
        inDir = util.file.get_test_input_path(self)
        with illumina.IlluminaDirectory(os.path.join(inDir, 'bcl-both-2.zip')) as idir:
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
            self.assertRaises(Exception, idir.get_SampleSheet)
        with illumina.IlluminaDirectory(os.path.join(inDir, 'bcl-samplesheet.tar.gz')) as idir:
            self.assertRaises(Exception, idir.get_RunInfo)


class TestDifficultSampleNames(TestCaseWithTmp):

    def test_paired_1(self):
        inDir = util.file.get_test_input_path(self)
        outBam = util.file.mkstempfname('.bam')
        outHeader = util.file.mkstempfname('.txt')
        sampleSheet = os.path.join(inDir, 'SampleSheet.csv')
        runInfo = os.path.join(inDir, 'RunInfo.xml')
        fastq = (os.path.join(inDir, 'mebv-0-1_S5_L001_R1_001.fastq.gz'),
                 os.path.join(inDir, 'mebv-0-1_S5_L001_R2_001.fastq.gz'))
        illumina.miseq_fastq_to_bam(outBam, sampleSheet, fastq[0], inFastq2=fastq[1], runInfo=runInfo)
        rgs = list(tools.samtools.SamtoolsTool().getReadGroups(outBam).values())
        self.assertEqual(len(rgs), 1)
        rgs = rgs[0]
        self.assertEqual(rgs.get('ID'), 'AEF96')
        self.assertEqual(rgs.get('PL'), 'illumina')
        self.assertEqual(rgs.get('PU'), 'AEF96.1.CGTACTAG-CTAAGCCT')
        self.assertEqual(rgs.get('LB'), u'GID-14-E021.ldifficult-value+for_-Sénégalsample_name0.1')
        self.assertEqual(rgs.get('SM'), u'GID-14-E021')
        self.assertEqual(rgs.get('CN'), 'M04004')
        self.assertTrue(rgs.get('DT','').startswith('2015-08-2'))

    def test_inline_commas_strings(self):
        inDir = util.file.get_test_input_path(self)
        samples = illumina.SampleSheet(os.path.join(inDir, 'SampleSheet-inline-commas-strings.csv'))
        self.assertEqual(samples.num_indexes(), 2)
        self.assertEqual(len(samples.get_rows()), 18)

        sample_names = [r["sample"] for r in samples.get_rows()]
        names_to_validate = [
            'Zika "seedstock_1 (in K562, 5ng)',
            "Zika 'seedstock_3 (in K562, 5ng)",
            "Zika seedstock_7 (in K562, 10pg)"
        ]
        for sample_name in names_to_validate:
            assert util.file.string_to_file_name(sample_name) in sample_names


class TestIlluminaBarcodeHelper(TestCaseWithTmp):
    def test_one_correction(self):
        dir_prefix = "one_correction"
        in_dir = util.file.get_test_input_path(self)
        in_barcodes = os.path.join(in_dir,dir_prefix,"barcodes.txt")
        in_metrics = os.path.join(in_dir,dir_prefix,"metrics.txt")
        out_report = util.file.mkstempfname('.txt')
        expected = os.path.join(in_dir,dir_prefix,"expected.txt")

        args = [in_barcodes, in_metrics, out_report]
        args = illumina.parser_guess_barcodes(argparse.ArgumentParser()).parse_args(args)
        args.func_main(args)

        self.assertEqualContents(out_report, expected)

    def test_ambiguous(self):
        dir_prefix = "ambiguous"
        in_dir = util.file.get_test_input_path(self)
        in_barcodes = os.path.join(in_dir,dir_prefix,"barcodes.txt")
        in_metrics = os.path.join(in_dir,dir_prefix,"metrics.txt")
        out_report = util.file.mkstempfname('.txt')
        expected = os.path.join(in_dir,dir_prefix,"expected.txt")

        args = [in_barcodes, in_metrics, out_report]
        args = illumina.parser_guess_barcodes(argparse.ArgumentParser()).parse_args(args)
        args.func_main(args)

        self.assertEqualContents(out_report, expected)

    def test_single_index_run(self):
        dir_prefix = "single_index"
        in_dir = util.file.get_test_input_path(self)
        in_barcodes = os.path.join(in_dir,dir_prefix,"barcodes.txt")
        in_metrics = os.path.join(in_dir,dir_prefix,"metrics.txt")
        out_report = util.file.mkstempfname('.txt')
        expected = os.path.join(in_dir,dir_prefix,"expected.txt")

        args = [in_barcodes, in_metrics, out_report]
        args = illumina.parser_guess_barcodes(argparse.ArgumentParser()).parse_args(args)
        args.func_main(args)

        self.assertEqualContents(out_report, expected)

    def test_few_assigned(self):
        dir_prefix = "few_assigned"
        in_dir = util.file.get_test_input_path(self)
        in_barcodes = os.path.join(in_dir,dir_prefix,"barcodes.txt")
        in_metrics = os.path.join(in_dir,dir_prefix,"metrics.txt")
        out_report = util.file.mkstempfname('.txt')

        self.assertRaises(Exception, illumina.main_guess_barcodes, in_barcodes, in_metrics, out_report)

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
        rgs = list(tools.samtools.SamtoolsTool().getReadGroups(outBam).values())
        self.assertEqual(len(rgs), 1)
        rgs = rgs[0]
        self.assertEqual(rgs.get('ID'), 'AEF96')
        self.assertEqual(rgs.get('PL'), 'illumina')
        self.assertEqual(rgs.get('PU'), 'AEF96.1.CGTACTAG-CTAAGCCT')
        self.assertEqual(rgs.get('LB'), 'mebv.0.1')
        self.assertEqual(rgs.get('SM'), 'mebv.0.1')
        self.assertEqual(rgs.get('CN'), 'M04004')
        self.assertTrue(rgs.get('DT','').startswith('2015-08-2'))

    def test_paired_2(self):
        inDir = util.file.get_test_input_path(self)
        outBam = util.file.mkstempfname('.bam')
        outHeader = util.file.mkstempfname('.txt')
        sampleSheet = os.path.join(inDir, 'SampleSheet.csv')
        runInfo = os.path.join(inDir, 'RunInfo.xml')
        fastq = (os.path.join(inDir, 'mebv-48-5_S17_L001_R1_001.fastq.gz'),
                 os.path.join(inDir, 'mebv-48-5_S17_L001_R2_001.fastq.gz'))
        illumina.miseq_fastq_to_bam(outBam, sampleSheet, fastq[0], inFastq2=fastq[1], runInfo=runInfo)
        rgs = list(tools.samtools.SamtoolsTool().getReadGroups(outBam).values())
        self.assertEqual(len(rgs), 1)
        rgs = rgs[0]
        self.assertEqual(rgs.get('ID'), 'AEF96')
        self.assertEqual(rgs.get('PL'), 'illumina')
        self.assertEqual(rgs.get('PU'), 'AEF96.1.GGACTCCT-TATCCTCT')
        self.assertEqual(rgs.get('LB'), 'mebv.48.5')
        self.assertEqual(rgs.get('SM'), 'mebv.48.5')
        self.assertEqual(rgs.get('CN'), 'M04004')
        self.assertTrue(rgs.get('DT','').startswith('2015-08-2'))

    def test_paired_custom_seq_center(self):
        inDir = util.file.get_test_input_path(self)
        outBam = util.file.mkstempfname('.bam')
        outHeader = util.file.mkstempfname('.txt')
        sampleSheet = os.path.join(inDir, 'SampleSheet.csv')
        runInfo = os.path.join(inDir, 'RunInfo.xml')
        fastq = (os.path.join(inDir, 'mebv-48-5_S17_L001_R1_001.fastq.gz'),
                 os.path.join(inDir, 'mebv-48-5_S17_L001_R2_001.fastq.gz'))
        illumina.miseq_fastq_to_bam(outBam, sampleSheet, fastq[0], inFastq2=fastq[1], runInfo=runInfo, sequencing_center='CustomSeqCenter')
        rgs = list(tools.samtools.SamtoolsTool().getReadGroups(outBam).values())
        self.assertEqual(len(rgs), 1)
        rgs = rgs[0]
        self.assertEqual(rgs.get('ID'), 'AEF96')
        self.assertEqual(rgs.get('PL'), 'illumina')
        self.assertEqual(rgs.get('PU'), 'AEF96.1.GGACTCCT-TATCCTCT')
        self.assertEqual(rgs.get('LB'), 'mebv.48.5')
        self.assertEqual(rgs.get('SM'), 'mebv.48.5')
        self.assertEqual(rgs.get('CN'), 'CustomSeqCenter')
        self.assertTrue(rgs.get('DT','').startswith('2015-08-2'))

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
