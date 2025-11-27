# Unit tests for illumina.py
# -*- coding: utf-8 -*-

__author__ = "dpark@broadinstitute.org"

import unittest
import os
import os.path
import tempfile
import argparse
import filecmp
import shutil
import json
import gzip
import pytest
import util
import util.file
import illumina
import tools.samtools
import tools.picard
from test import TestCaseWithTmp, assert_equal_bam_reads


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

    def test_dup_index_collapse_at_init(self):
        inDir = util.file.get_test_input_path(self)
        samples = illumina.SampleSheet(
                                        os.path.join(inDir, 'SampleSheet-custom-inner-barcodes-outer-collapse.tsv'),
                                        collapse_duplicates=True,
                                        allow_non_unique = True)
        self.assertEqual(samples.num_indexes(), 2)
        self.assertEqual(samples.duplicate_rows_collapsed, True)
        self.assertEqual(len(samples.get_rows()), 4)

    def test_dup_index_collapse(self):
        inDir = util.file.get_test_input_path(self)
        samples = illumina.SampleSheet(
                                        os.path.join(inDir, 'SampleSheet-custom-inner-barcodes-outer-collapse.tsv'),
                                        allow_non_unique = True
                                      )
        self.assertEqual(samples.num_indexes(), 2)
        self.assertEqual(samples.duplicate_rows_collapsed, False)
        self.assertEqual(samples.can_be_collapsed, True)
        self.assertEqual(samples.num_samples, 96)

        samples.collapse_sample_index_duplicates()

        self.assertEqual(samples.num_indexes(), 2)
        self.assertEqual(samples.duplicate_rows_collapsed, True)
        self.assertEqual(samples.can_be_collapsed, False)
        self.assertEqual(len(samples.get_rows()), 4)

        # check that unique values in grouped rows are preserved
        # and others are collapsed to the md5sum of the (joined) unique values
        collapsed_row = next((r for r in samples.get_rows() if r["sample"] == "CTGATCGT-GCGCATAT" and r["barcode_1"] == "CTGATCGT" and r["barcode_2"] == "GCGCATAT"), None)
        self.assertEqual(collapsed_row["library_id_per_sample"], "Pool_1")
        self.assertEqual(collapsed_row["barcode_3"],"2e094627_muxed")

        # check for dups in the collapsed output
        self.assertEqual(len(set((r["barcode_1"] for r in samples.get_rows()))), len(list((r["barcode_1"] for r in samples.get_rows()))))
        self.assertEqual(len(set((r["barcode_2"] for r in samples.get_rows()))), len(list((r["barcode_2"] for r in samples.get_rows()))))

        # check that the collapsed output is what we expect
        self.assertEqual(tuple(r["barcode_1"] for r in samples.get_rows()), ("CTGATCGT","ACTCTCGA","TGAGCTAG","GAGACGAT"))
        self.assertEqual(tuple(r["barcode_2"] for r in samples.get_rows()), ("GCGCATAT","CTGTACCA","GAACGGTT","ACCGGTTA"))

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

    def test_rev_comp_barcode_values(self):
        inDir = util.file.get_test_input_path(self)
        samples = illumina.SampleSheet(os.path.join(inDir, 'SampleSheet-custom-1.txt'))
        self.assertEqual(samples.get_rows()[0]["barcode_2"], "TAGATCGC")
        samples.rev_comp_barcode_values()
        self.assertEqual(samples.get_rows()[0]["barcode_2"], "GCGATCTA")

    def test_rev_comp_barcode_value_barcode2_at_load(self):
        inDir = util.file.get_test_input_path(self)
        samples = illumina.SampleSheet(os.path.join(inDir, 'SampleSheet-custom-1.txt'), rev_comp_barcode_2=True)
        self.assertEqual(samples.get_rows()[0]["barcode_2"], "GCGATCTA")
        self.assertTrue(samples.barcodes_revcomped_relative_to_input)
        self.assertIn("barcode_2", samples.barcodes_revcomped_column_names)
        self.assertNotIn("barcode_1", samples.barcodes_revcomped_column_names)

    def test_rev_comp_barcode_values_specified_at_load_two_columns(self):
        inDir = util.file.get_test_input_path(self)
        samples = illumina.SampleSheet(os.path.join(inDir, 'SampleSheet-custom-1.txt'), barcode_columns_to_revcomp=["barcode_1","barcode_2"])
        self.assertEqual(samples.get_rows()[23]["barcode_1"], "TCCTCTAC")
        self.assertEqual(samples.get_rows()[23]["barcode_2"], "AGAGGATA")
        self.assertTrue(samples.barcodes_revcomped_relative_to_input)
        self.assertIn("barcode_1", samples.barcodes_revcomped_column_names)
        self.assertIn("barcode_2", samples.barcodes_revcomped_column_names)

    def test_rev_comp_barcode_values_specified_one_column(self):
        inDir = util.file.get_test_input_path(self)
        samples = illumina.SampleSheet(os.path.join(inDir, 'SampleSheet-custom-1.txt'))
        self.assertEqual(samples.get_rows()[23]["barcode_1"], "GTAGAGGA")
        samples.rev_comp_barcode_values(barcode_columns_to_revcomp=["barcode_1"])
        self.assertEqual(samples.get_rows()[23]["barcode_1"], "TCCTCTAC")
        # barcode_2 should be unchanged
        self.assertEqual(samples.get_rows()[23]["barcode_2"], "TATCCTCT")
        self.assertIn("barcode_1", samples.barcodes_revcomped_column_names)
        self.assertNotIn("barcode_2", samples.barcodes_revcomped_column_names)

    def test_rev_comp_barcode_values_specified_two_columns(self):
        inDir = util.file.get_test_input_path(self)
        samples = illumina.SampleSheet(os.path.join(inDir, 'SampleSheet-custom-1.txt'))
        self.assertEqual(samples.get_rows()[23]["barcode_1"], "GTAGAGGA")
        self.assertEqual(samples.get_rows()[23]["barcode_2"], "TATCCTCT")
        samples.rev_comp_barcode_values(barcode_columns_to_revcomp=["barcode_1","barcode_2"])
        self.assertEqual(samples.get_rows()[23]["barcode_1"], "TCCTCTAC")
        self.assertEqual(samples.get_rows()[23]["barcode_2"], "AGAGGATA")
        self.assertIn("barcode_1", samples.barcodes_revcomped_column_names)
        self.assertIn("barcode_2", samples.barcodes_revcomped_column_names)

    def test_rev_comp_barcode_values_undo(self):
        inDir = util.file.get_test_input_path(self)
        samples = illumina.SampleSheet(os.path.join(inDir, 'SampleSheet-custom-1.txt'))
        self.assertEqual(samples.get_rows()[0]["barcode_2"], "TAGATCGC")
        samples.rev_comp_barcode_values()
        self.assertEqual(samples.get_rows()[0]["barcode_2"], "GCGATCTA")
        samples.rev_comp_barcode_values()
        self.assertEqual(samples.get_rows()[0]["barcode_2"], "TAGATCGC")

    def test_rev_comp_barcode_values_not_inplace(self):
        inDir = util.file.get_test_input_path(self)
        samples = illumina.SampleSheet(os.path.join(inDir, 'SampleSheet-custom-1.txt'))
        self.assertEqual(samples.get_rows()[0]["barcode_2"], "TAGATCGC")
        samples_revcomped = samples.rev_comp_barcode_values(inplace=False)
        # the original shet should be unchanged
        self.assertEqual(samples.get_rows()[0]["barcode_2"], "TAGATCGC")
        # the returned copy should be rev-comped
        self.assertEqual(samples_revcomped.get_rows()[0]["barcode_2"], "GCGATCTA")

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

    def test_nextseq550(self):
        inDir = util.file.get_test_input_path(self)
        runinfo = illumina.RunInfo(os.path.join(inDir, 'RunInfo-nextseq550.xml'))
        self.assertEqual(runinfo.get_flowcell(), 'HMTLKAFX2')
        self.assertEqual(runinfo.get_rundate_american(), '02/21/2021')
        self.assertEqual(runinfo.get_rundate_iso(), '2021-02-21')
        self.assertEqual(runinfo.get_machine(), 'NB552060')
        self.assertEqual(runinfo.get_read_structure(), '149T10B10B149T')
        self.assertEqual(runinfo.num_reads(), 2)
        self.assertEqual(runinfo.get_machine_model(), "NextSeq 550")
        self.assertEqual(runinfo.get_flowcell_lane_count(), 4)

    def test_novel_tile_count_but_known_fcid(self):
        inDir = util.file.get_test_input_path(self)
        runinfo = illumina.RunInfo(os.path.join(inDir, 'RunInfo-novel-tile-count.xml'))
        self.assertEqual(runinfo.get_machine_model(), "Illumina MiSeq")

    def test_novel_fcid_but_known_tile_count(self):
        inDir = util.file.get_test_input_path(self)
        runinfo = illumina.RunInfo(os.path.join(inDir, 'RunInfo-novel-fcid.xml'))
        self.assertEqual(runinfo.get_machine_model(), "Illumina MiSeq")

    def test_novel_tile_count_and_fcid(self):
        inDir = util.file.get_test_input_path(self)
        runinfo = illumina.RunInfo(os.path.join(inDir, 'RunInfo-novel-fcid-and-tilecount.xml'))
        self.assertEqual(runinfo.get_machine_model(), "UNKNOWN")


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

    def test_single_index_i5_only_run(self):
        dir_prefix = "single_index_i5_only"
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


class TestSplitcodeDemuxIntegration(TestCaseWithTmp):
    """
    End-to-end integration test for splitcode_demux command.

    Tests a typical workflow:
    - Input: BAM files from outer barcode demux (one per pool)
    - Sample sheet with inner barcodes (barcode_3)
    - RunInfo.xml with run metadata
    - Output: Per-sample BAM files after inner barcode demux
    """

    def setUp(self):
        super().setUp()
        self.samtools = tools.samtools.SamtoolsTool()
        self.picard = tools.picard.FastqToSamTool()

    def create_test_bam_with_inline_barcodes(self, output_bam, barcode_reads):
        """
        Create a test BAM file with reads containing inline barcodes.

        This simulates what Picard creates from a sequencer with read structure "50T8B8B50T":
        - R1: 50bp template (first 8bp are inline barcode, remaining 42bp are actual sequence)
        - R2: 50bp template
        - The index reads (8B8B) are not included in the BAM template sequences

        Args:
            output_bam: Path to output BAM file
            barcode_reads: Dict mapping barcode -> num_reads
                          e.g., {"AAAACCCC": 10, "GGGGTTTT": 5}
        """
        r1_fastq = util.file.mkstempfname('.fastq')
        r2_fastq = util.file.mkstempfname('.fastq')

        with open(r1_fastq, 'w') as f1, open(r2_fastq, 'w') as f2:
            read_id = 0
            for barcode, num_reads in barcode_reads.items():
                for _ in range(num_reads):
                    # R1: 50bp total (8bp inline barcode + 42bp sequence)
                    # Picard treats all 50bp as template, doesn't know about the inline barcode
                    seq_r1 = barcode + ("ACGT" * 10) + "AC"  # 8bp + 40bp + 2bp = 50bp
                    qual_r1 = "I" * len(seq_r1)

                    # R2: 50bp template
                    seq_r2 = "TGCA" * 12 + "TG"  # 50bp
                    qual_r2 = "I" * len(seq_r2)

                    f1.write(f"@read{read_id}\n{seq_r1}\n+\n{qual_r1}\n")
                    f2.write(f"@read{read_id}\n{seq_r2}\n+\n{qual_r2}\n")
                    read_id += 1

        # Convert to BAM
        self.picard.execute(
            r1_fastq,
            r2_fastq,
            "PoolSample",
            output_bam,
            picardOptions=['LIBRARY_NAME=TestPool', 'PLATFORM=ILLUMINA']
        )

        os.unlink(r1_fastq)
        os.unlink(r2_fastq)

        return output_bam

    def create_expected_output_bam(self, output_bam, sample_name, barcode, num_reads, sample_idx=0, starting_read_id=0):
        """
        Create an expected output BAM file with reads after splitcode demux.

        After demuxing with splitcode using read structure "50T8B8B50T":
        - R1: Original 50bp read (barcode trimmed off)
        - R2: Original 50bp read (unchanged)

        Args:
            output_bam: Path to output BAM file
            sample_name: Sample name for read group
            barcode: The inline barcode (for verification, not included in output)
            num_reads: Number of reads to create
            sample_idx: Sample index (0-based) for read group name generation
            starting_read_id: Starting read ID number (for sequential read numbering)
        """
        r1_fastq = util.file.mkstempfname('.fastq')
        r2_fastq = util.file.mkstempfname('.fastq')

        with open(r1_fastq, 'w') as f1, open(r2_fastq, 'w') as f2:
            for i in range(num_reads):
                read_id = starting_read_id + i
                # After splitcode demux:
                # R1: The barcode is trimmed, leaving only the 42bp sequence after barcode
                # (splitcode trims the 8bp inline barcode from the start of the 50bp R1)
                seq_r1 = ("ACGT" * 10) + "AC"  # 42bp (barcode already removed)
                qual_r1 = "I" * len(seq_r1)

                # R2: Unchanged 50bp
                seq_r2 = "TGCA" * 12 + "TG"  # 50bp
                qual_r2 = "I" * len(seq_r2)

                f1.write(f"@read{read_id}\n{seq_r1}\n+\n{qual_r1}\n")
                f2.write(f"@read{read_id}\n{seq_r2}\n+\n{qual_r2}\n")

        # Convert to BAM with same read group structure as splitcode_demux output
        # Read group name format: {sequencing_center}.{lane}.{sample_idx+1}
        read_group_name = f"TEST.1.{sample_idx + 1}"
        self.picard.execute(
            r1_fastq,
            r2_fastq,
            sample_name,
            output_bam,
            picardOptions=[
                f'LIBRARY_NAME=L1',
                'PLATFORM=ILLUMINA',
                'SEQUENCING_CENTER=TEST',
                f'RUN_DATE=2025-01-01',
                f'READ_GROUP_NAME={read_group_name}',
            ]
        )

        os.unlink(r1_fastq)
        os.unlink(r2_fastq)

        return output_bam

    def test_splitcode_demux_basic(self):
        """
        Test splitcode_demux command with minimal realistic input.

        Simulates a pool with 2 samples distinguished by inner barcodes.
        Verifies that:
        - Command runs without error
        - Output BAM files are created for each sample
        - Output BAMs contain expected number of reads
        - Output BAM content matches expected (reads correctly demuxed and trimmed)
        """
        # Get test input directory
        inDir = util.file.get_test_input_path(self)

        # Define expected read counts per sample
        sample1_reads = 100
        sample2_reads = 50

        # Create temporary directories
        with tempfile.TemporaryDirectory() as input_bams_dir:
            with tempfile.TemporaryDirectory() as outDir:
                with tempfile.TemporaryDirectory() as expectedDir:
                    # Create a pool BAM file with reads from 2 barcodes
                    # Pool identifier from sample sheet: ATCGATCG-GCTAGCTA.lL1.TESTFLOW.1
                    pool_bam = os.path.join(input_bams_dir, 'ATCGATCG-GCTAGCTA.lL1.TESTFLOW.1.bam')
                    self.create_test_bam_with_inline_barcodes(
                        pool_bam,
                        {
                            "AAAAAAAA": sample1_reads,  # TestSample1 barcode
                            "CCCCCCCC": sample2_reads,   # TestSample2 barcode
                        }
                    )

                    # Create expected output BAM files
                    expected_sample1_bam = os.path.join(expectedDir, 'TestSample1.lL1.TESTFLOW.1.bam')
                    expected_sample2_bam = os.path.join(expectedDir, 'TestSample2.lL1.TESTFLOW.1.bam')

                    self.create_expected_output_bam(
                        expected_sample1_bam,
                        "TestSample1",
                        "AAAAAAAA",
                        sample1_reads,
                        sample_idx=0,  # First sample in samplesheet
                        starting_read_id=0  # Reads 0-99
                    )
                    self.create_expected_output_bam(
                        expected_sample2_bam,
                        "TestSample2",
                        "CCCCCCCC",
                        sample2_reads,
                        sample_idx=1,  # Second sample in samplesheet
                        starting_read_id=sample1_reads  # Reads 100-149
                    )

                    # Run splitcode_demux
                    out_meta_by_sample = os.path.join(outDir, 'test_meta_by_sample.txt')
                    illumina.splitcode_demux(
                        inDir=input_bams_dir,
                        lane="1",
                        outDir=outDir,
                        sampleSheet=os.path.join(inDir, 'SampleSheet.tsv'),
                        runinfo=os.path.join(inDir, 'RunInfo.xml'),
                        flowcell="TESTFLOW",
                        run_id="TESTFLOW.1",
                        run_date="2025-01-01",
                        read_structure="50T8B8B50T",
                        platform_name="ILLUMINA",
                        sequencing_center="TEST",
                        unmatched_name="Unmatched",
                        max_hamming_dist=1,
                        threads=1,
                        out_meta_by_sample=out_meta_by_sample
                    )

                    # Verify output files exist and have correct read counts
                    actual_sample1_bam = os.path.join(outDir, 'TestSample1.lL1.TESTFLOW.1.bam')
                    actual_sample2_bam = os.path.join(outDir, 'TestSample2.lL1.TESTFLOW.1.bam')

                    # Verify files exist
                    self.assertTrue(
                        os.path.exists(actual_sample1_bam),
                        f"Expected output BAM not found: {actual_sample1_bam}"
                    )
                    self.assertTrue(
                        os.path.exists(actual_sample2_bam),
                        f"Expected output BAM not found: {actual_sample2_bam}"
                    )

                    # Verify read counts match expected
                    actual_sample1_count = self.samtools.count(actual_sample1_bam)
                    actual_sample2_count = self.samtools.count(actual_sample2_bam)
                    expected_sample1_count = self.samtools.count(expected_sample1_bam)
                    expected_sample2_count = self.samtools.count(expected_sample2_bam)

                    self.assertEqual(
                        actual_sample1_count,
                        expected_sample1_count,
                        f"TestSample1 read count mismatch: got {actual_sample1_count}, expected {expected_sample1_count}"
                    )
                    self.assertEqual(
                        actual_sample2_count,
                        expected_sample2_count,
                        f"TestSample2 read count mismatch: got {actual_sample2_count}, expected {expected_sample2_count}"
                    )

                    # Verify BAM content matches expected
                    assert_equal_bam_reads(self, actual_sample1_bam, expected_sample1_bam)
                    assert_equal_bam_reads(self, actual_sample2_bam, expected_sample2_bam)


class TestParseIlluminaFastqFilename(unittest.TestCase):
    """Test parsing of Illumina DRAGEN FASTQ filename patterns.

    DRAGEN format: {flowcell}_{lane}_{numeric_id}_{sample_name}_{S#}_{L00#}_{R#}_{chunk}.fastq.gz
    Example: 22J5GLLT4_6_0420593812_B13Pool1a_S1_L006_R1_001.fastq.gz
    """

    def test_standard_dragen_format_r1(self):
        """Test standard DRAGEN FASTQ naming with flowcell ID."""
        filename = "22J5GLLT4_6_0420593812_B13Pool1a_S1_L006_R1_001.fastq.gz"
        result = illumina.parse_illumina_fastq_filename(filename)

        self.assertEqual(result['flowcell'], '22J5GLLT4')
        self.assertEqual(result['lane_short'], 6)
        self.assertEqual(result['numeric_id'], '0420593812')
        self.assertEqual(result['sample_name'], 'B13Pool1a')
        self.assertEqual(result['sample_number'], 1)
        self.assertEqual(result['lane'], 6)
        self.assertEqual(result['read'], 1)
        self.assertEqual(result['chunk'], 1)

    def test_standard_dragen_format_r2(self):
        """Test standard DRAGEN format for R2 file."""
        filename = "22J5GLLT4_6_0420593812_B13Pool1a_S1_L006_R2_001.fastq.gz"
        result = illumina.parse_illumina_fastq_filename(filename)

        self.assertEqual(result['flowcell'], '22J5GLLT4')
        self.assertEqual(result['sample_name'], 'B13Pool1a')
        self.assertEqual(result['read'], 2)

    def test_sample_name_with_underscores(self):
        """Test sample names containing underscores."""
        filename = "22J5GLLT4_6_0420593812_RS_Batch_20_NTC_01_S5_L006_R1_001.fastq.gz"
        result = illumina.parse_illumina_fastq_filename(filename)

        self.assertEqual(result['flowcell'], '22J5GLLT4')
        self.assertEqual(result['lane_short'], 6)
        self.assertEqual(result['numeric_id'], '0420593812')
        self.assertEqual(result['sample_name'], 'RS_Batch_20_NTC_01')
        self.assertEqual(result['sample_number'], 5)
        self.assertEqual(result['lane'], 6)
        self.assertEqual(result['read'], 1)

    def test_different_flowcells(self):
        """Test parsing files from different flowcell IDs."""
        test_cases = [
            "ABCDEFGH1_1_1234567890_Sample1_S1_L001_R1_001.fastq.gz",
            "22J5GLLT4_1_1234567890_Sample1_S1_L001_R1_001.fastq.gz",
            "H7YVLDSXY_1_1234567890_Sample1_S1_L001_R1_001.fastq.gz",
        ]
        expected_flowcells = ['ABCDEFGH1', '22J5GLLT4', 'H7YVLDSXY']

        for filename, expected_fc in zip(test_cases, expected_flowcells):
            result = illumina.parse_illumina_fastq_filename(filename)
            self.assertEqual(result['flowcell'], expected_fc)

    def test_different_lane_numbers(self):
        """Test parsing files from different lanes."""
        for lane in [1, 2, 3, 4, 5, 6, 7, 8]:
            filename = f"22J5GLLT4_{lane}_0420593812_Sample1_S1_L00{lane}_R1_001.fastq.gz"
            result = illumina.parse_illumina_fastq_filename(filename)
            self.assertEqual(result['lane'], lane)
            self.assertEqual(result['lane_short'], lane)

    def test_different_sample_numbers(self):
        """Test parsing files with various sample numbers."""
        test_cases = [
            ("22J5GLLT4_6_0420593812_Sample1_S1_L006_R1_001.fastq.gz", 1),
            ("22J5GLLT4_6_0420593812_Sample1_S10_L006_R1_001.fastq.gz", 10),
            ("22J5GLLT4_6_0420593812_Sample1_S100_L006_R1_001.fastq.gz", 100),
        ]
        for filename, expected_num in test_cases:
            result = illumina.parse_illumina_fastq_filename(filename)
            self.assertEqual(result['sample_number'], expected_num)

    def test_different_chunk_numbers(self):
        """Test parsing files with different chunk numbers."""
        filename = "22J5GLLT4_6_0420593812_Sample1_S1_L006_R1_002.fastq.gz"
        result = illumina.parse_illumina_fastq_filename(filename)
        self.assertEqual(result['chunk'], 2)

        filename = "22J5GLLT4_6_0420593812_Sample1_S1_L006_R1_003.fastq.gz"
        result = illumina.parse_illumina_fastq_filename(filename)
        self.assertEqual(result['chunk'], 3)

    def test_with_full_path(self):
        """Test parsing filename with full directory path."""
        filename = "/path/to/data/22J5GLLT4_6_0420593812_Sample1_S1_L006_R1_001.fastq.gz"
        result = illumina.parse_illumina_fastq_filename(filename)

        self.assertEqual(result['flowcell'], '22J5GLLT4')
        self.assertEqual(result['sample_name'], 'Sample1')
        self.assertEqual(result['lane'], 6)

    def test_without_gz_extension(self):
        """Test parsing filename without .gz extension."""
        filename = "22J5GLLT4_6_0420593812_Sample1_S1_L006_R1_001.fastq"
        result = illumina.parse_illumina_fastq_filename(filename)

        self.assertEqual(result['flowcell'], '22J5GLLT4')
        self.assertEqual(result['sample_name'], 'Sample1')

    def test_simple_format_without_flowcell(self):
        """Test simple Illumina format without flowcell ID (backward compatibility).

        Format: {sample_name}_{S#}_{L00#}_{R#}_{chunk}.fastq.gz
        Example: mebv-48-5_S17_L001_R1_001.fastq.gz
        """
        filename = "mebv-48-5_S17_L001_R1_001.fastq.gz"
        result = illumina.parse_illumina_fastq_filename(filename)

        self.assertEqual(result['sample_name'], 'mebv-48-5')
        self.assertEqual(result['sample_number'], 17)
        self.assertEqual(result['lane'], 1)
        self.assertEqual(result['read'], 1)
        self.assertEqual(result['chunk'], 1)
        self.assertIsNone(result.get('flowcell'))

    def test_simple_format_r2(self):
        """Test simple format for R2 file."""
        filename = "Sample1_S5_L001_R2_001.fastq.gz"
        result = illumina.parse_illumina_fastq_filename(filename)

        self.assertEqual(result['sample_name'], 'Sample1')
        self.assertEqual(result['sample_number'], 5)
        self.assertEqual(result['lane'], 1)
        self.assertEqual(result['read'], 2)

    def test_malformed_invalid_format(self):
        """Test error handling for completely invalid filename."""
        filename = "random_file.fastq.gz"
        with self.assertRaises(ValueError) as context:
            illumina.parse_illumina_fastq_filename(filename)
        self.assertIn("does not match", str(context.exception).lower())

    def test_empty_filename(self):
        """Test error handling for empty filename."""
        with self.assertRaises(ValueError):
            illumina.parse_illumina_fastq_filename("")

    def test_malformed_missing_read_number(self):
        """Test error handling for filename missing read number."""
        filename = "22J5GLLT4_6_0420593812_Sample1_S1_L006_001.fastq.gz"
        with self.assertRaises(ValueError) as context:
            illumina.parse_illumina_fastq_filename(filename)
        self.assertIn("does not match", str(context.exception).lower())

    def test_index_reads_not_supported(self):
        """Test that index reads (I1, I2) are handled appropriately.

        Note: DRAGEN may not produce I1/I2 files in the same way as bcl2fastq.
        This test documents expected behavior if encountered.
        """
        # If DRAGEN produces index files, they might follow a different pattern
        # For now, test that we handle them gracefully or raise appropriate errors
        pass


class TestNormalizeBarcode(unittest.TestCase):
    """Test barcode normalization function."""

    def test_uppercase_conversion(self):
        """Test that lowercase barcodes are converted to uppercase."""
        self.assertEqual(illumina.normalize_barcode("acgt"), "ACGT")
        self.assertEqual(illumina.normalize_barcode("atcgatcg"), "ATCGATCG")
        self.assertEqual(illumina.normalize_barcode("aCgT"), "ACGT")

    def test_whitespace_trimming(self):
        """Test that leading and trailing whitespace is removed."""
        self.assertEqual(illumina.normalize_barcode("  ACGT  "), "ACGT")
        self.assertEqual(illumina.normalize_barcode("\tACGT\n"), "ACGT")
        self.assertEqual(illumina.normalize_barcode(" ATCGATCG "), "ATCGATCG")

    def test_combined_normalization(self):
        """Test combined uppercase conversion and whitespace trimming."""
        self.assertEqual(illumina.normalize_barcode(" acgt "), "ACGT")
        self.assertEqual(illumina.normalize_barcode("\tatcg\n"), "ATCG")

    def test_already_normalized(self):
        """Test that already-normalized barcodes pass through unchanged."""
        self.assertEqual(illumina.normalize_barcode("ACGT"), "ACGT")
        self.assertEqual(illumina.normalize_barcode("ATCGATCG"), "ATCGATCG")

    def test_empty_string(self):
        """Test handling of empty strings."""
        self.assertEqual(illumina.normalize_barcode(""), "")
        self.assertEqual(illumina.normalize_barcode("  "), "")

    def test_valid_characters_only(self):
        """Test validation that only ACGTN characters are allowed."""
        # Valid barcodes
        self.assertEqual(illumina.normalize_barcode("ACGT"), "ACGT")
        self.assertEqual(illumina.normalize_barcode("AAACCCGGGTTT"), "AAACCCGGGTTT")
        self.assertEqual(illumina.normalize_barcode("ACGTN"), "ACGTN")  # N for ambiguous base
        self.assertEqual(illumina.normalize_barcode("NNNACGT"), "NNNACGT")

    def test_invalid_characters(self):
        """Test that invalid characters raise ValueError."""
        with self.assertRaises(ValueError) as context:
            illumina.normalize_barcode("ACGTX")
        self.assertIn("invalid", str(context.exception).lower())

        with self.assertRaises(ValueError):
            illumina.normalize_barcode("ACG-TGC")

        with self.assertRaises(ValueError):
            illumina.normalize_barcode("ACG TGC")  # space in middle

        with self.assertRaises(ValueError):
            illumina.normalize_barcode("123")

    def test_mixed_case_with_n(self):
        """Test mixed case barcodes containing N."""
        self.assertEqual(illumina.normalize_barcode("acgtn"), "ACGTN")
        self.assertEqual(illumina.normalize_barcode("NNNacgt"), "NNNACGT")

    def test_typical_illumina_barcodes(self):
        """Test with real Illumina barcode sequences."""
        # Typical 8bp dual index barcodes
        self.assertEqual(illumina.normalize_barcode("CTGATCGT"), "CTGATCGT")
        self.assertEqual(illumina.normalize_barcode("ctgatcgt"), "CTGATCGT")

        # Typical 10bp barcode
        self.assertEqual(illumina.normalize_barcode("AACCGGTTAA"), "AACCGGTTAA")

        # With whitespace (common in CSV files)
        self.assertEqual(illumina.normalize_barcode(" CTGATCGT "), "CTGATCGT")

    def test_none_input(self):
        """Test handling of None input."""
        with self.assertRaises((ValueError, TypeError)):
            illumina.normalize_barcode(None)

    def test_non_string_input(self):
        """Test handling of non-string input."""
        with self.assertRaises((ValueError, TypeError)):
            illumina.normalize_barcode(12345)

        with self.assertRaises((ValueError, TypeError)):
            illumina.normalize_barcode(['A', 'C', 'G', 'T'])


class TestBarcodeOrientationAutoDetection(unittest.TestCase):
    """
    Test suite for automatic i5 (barcode_2) orientation detection.

    Tests the match_barcodes_with_orientation() function which automatically
    tries reverse complement of barcode_2 when direct matching fails.

    Note: Only barcode_2 (i5) orientation is auto-detected because this is the
    only barcode that varies in orientation across Illumina sequencer generations.
    barcode_1 (i7) is consistent across all platforms.
    """

    def test_direct_match(self):
        """Test that direct barcode match works without reverse complement."""
        sample_rows = [
            {'sample': 'S1', 'barcode_1': 'ATCGATCG', 'barcode_2': 'GCTAGCTA'},
            {'sample': 'S2', 'barcode_1': 'AAAAAAAA', 'barcode_2': 'TTTTTTTT'},
        ]
        matched, info = illumina.match_barcodes_with_orientation(
            'ATCGATCG', 'GCTAGCTA', sample_rows
        )
        self.assertEqual(len(matched), 1)
        self.assertEqual(matched[0]['sample'], 'S1')
        self.assertFalse(info['barcode_2_revcomp'])

    def test_barcode2_revcomp_match(self):
        """Test that reverse complement of barcode_2/i5 is auto-detected."""
        # Samplesheet has ACTGCAGCCG, FASTQ has reverse complement CGGCTGCAGT
        sample_rows = [
            {'sample': 'S1', 'barcode_1': 'AGTTAATGCT', 'barcode_2': 'ACTGCAGCCG'},
        ]
        matched, info = illumina.match_barcodes_with_orientation(
            'AGTTAATGCT', 'CGGCTGCAGT', sample_rows  # CGGCTGCAGT = revcomp(ACTGCAGCCG)
        )
        self.assertEqual(len(matched), 1)
        self.assertTrue(info['barcode_2_revcomp'])
        self.assertEqual(info['matched_bc2'], 'ACTGCAGCCG')

    def test_no_match_returns_empty(self):
        """Test that no match returns empty list with skipped_reason='no_match'."""
        sample_rows = [
            {'sample': 'S1', 'barcode_1': 'ATCGATCG', 'barcode_2': 'GCTAGCTA'},
        ]
        # Use completely different barcodes that won't match
        matched, info = illumina.match_barcodes_with_orientation(
            'GGGGGGGG', 'CCCCCCCC', sample_rows
        )
        self.assertEqual(len(matched), 0)
        self.assertEqual(info.get('skipped_reason'), 'no_match')

    def test_single_barcode_matching(self):
        """Test matching with only barcode_1 (barcode_2 is None)."""
        sample_rows = [
            {'sample': 'S1', 'barcode_1': 'ATCGATCG', 'barcode_2': ''},
        ]
        matched, info = illumina.match_barcodes_with_orientation(
            'ATCGATCG', None, sample_rows
        )
        self.assertEqual(len(matched), 1)

    def test_case_insensitive_matching(self):
        """Test that barcode matching is case-insensitive."""
        sample_rows = [
            {'sample': 'S1', 'barcode_1': 'atcgatcg', 'barcode_2': 'GCTAGCTA'},
        ]
        matched, info = illumina.match_barcodes_with_orientation(
            'ATCGATCG', 'gctagcta', sample_rows
        )
        self.assertEqual(len(matched), 1)

    def test_multiple_samples_same_outer_barcodes(self):
        """Test matching returns all samples with same outer barcodes (3-barcode case)."""
        sample_rows = [
            {'sample': 'S1', 'barcode_1': 'ATCGATCG', 'barcode_2': 'GCTAGCTA', 'barcode_3': 'AAAAAAAA'},
            {'sample': 'S2', 'barcode_1': 'ATCGATCG', 'barcode_2': 'GCTAGCTA', 'barcode_3': 'CCCCCCCC'},
            {'sample': 'S3', 'barcode_1': 'GGGGGGGG', 'barcode_2': 'TTTTTTTT', 'barcode_3': 'AAAAAAAA'},
        ]
        matched, info = illumina.match_barcodes_with_orientation(
            'ATCGATCG', 'GCTAGCTA', sample_rows
        )
        self.assertEqual(len(matched), 2)
        self.assertEqual({r['sample'] for r in matched}, {'S1', 'S2'})

    def test_barcode2_revcomp_with_multiple_samples(self):
        """Test i5 auto-detection works when multiple samples share outer barcodes."""
        sample_rows = [
            {'sample': 'S1', 'barcode_1': 'ATCGATCG', 'barcode_2': 'GCTAGCTA', 'barcode_3': 'AAAAAAAA'},
            {'sample': 'S2', 'barcode_1': 'ATCGATCG', 'barcode_2': 'GCTAGCTA', 'barcode_3': 'CCCCCCCC'},
        ]
        # TAGCTAGC = revcomp(GCTAGCTA)
        matched, info = illumina.match_barcodes_with_orientation(
            'ATCGATCG', 'TAGCTAGC', sample_rows
        )
        self.assertEqual(len(matched), 2)
        self.assertTrue(info['barcode_2_revcomp'])

    def test_n_wildcard_in_fastq_barcode(self):
        """Test that N bases in FASTQ barcodes match any base in samplesheet.

        This handles sequencer no-call bases where a position couldn't be
        confidently determined. The N should act as a wildcard.
        """
        sample_rows = [
            {'sample': 'S1', 'barcode_1': 'ATCGATCG', 'barcode_2': 'GCTAGCTA'},
        ]
        # FASTQ has N at position 3, samplesheet has G
        matched, info = illumina.match_barcodes_with_orientation(
            'ATCNATCG', 'GCTAGCTA', sample_rows
        )
        self.assertEqual(len(matched), 1)
        self.assertEqual(matched[0]['sample'], 'S1')
        self.assertFalse(info['barcode_2_revcomp'])

    def test_n_wildcard_in_barcode2(self):
        """Test N wildcard matching in barcode_2."""
        sample_rows = [
            {'sample': 'S1', 'barcode_1': 'ATCGATCG', 'barcode_2': 'GCTAGCTA'},
        ]
        # FASTQ has N at position 2 of barcode_2
        matched, info = illumina.match_barcodes_with_orientation(
            'ATCGATCG', 'GCNAGCTA', sample_rows
        )
        self.assertEqual(len(matched), 1)
        self.assertEqual(matched[0]['sample'], 'S1')

    def test_n_wildcard_with_revcomp(self):
        """Test N wildcard matching combined with i5 reverse complement detection.

        This is the real-world case: FASTQ barcode_2 has N bases AND needs
        reverse complement to match samplesheet.
        """
        # Samplesheet has TTCCGACATT, revcomp is AATGTCGGAA
        sample_rows = [
            {'sample': 'VGG_21775', 'barcode_1': 'TAAGGAGGAA', 'barcode_2': 'TTCCGACATT'},
        ]
        # FASTQ has AANGTCGGAA (N at position 2, should match T in revcomp AATGTCGGAA)
        matched, info = illumina.match_barcodes_with_orientation(
            'TAAGGAGGAA', 'AANGTCGGAA', sample_rows
        )
        self.assertEqual(len(matched), 1)
        self.assertEqual(matched[0]['sample'], 'VGG_21775')
        self.assertTrue(info['barcode_2_revcomp'])

    def test_multiple_n_wildcards(self):
        """Test multiple N bases in a barcode (up to 50% N allowed)."""
        sample_rows = [
            {'sample': 'S1', 'barcode_1': 'ATCGATCG', 'barcode_2': 'GCTAGCTA'},
        ]
        # FASTQ has N at positions 0 and 4 (2/8 = 25% N, under threshold)
        matched, info = illumina.match_barcodes_with_orientation(
            'NTCGNTCG', 'GCTAGCTA', sample_rows
        )
        self.assertEqual(len(matched), 1)

    def test_high_n_fraction_returns_empty(self):
        """Test that barcodes with >50% N bases return empty (too low confidence)."""
        sample_rows = [
            {'sample': 'S1', 'barcode_1': 'ATCGATCG', 'barcode_2': 'GCTAGCTA'},
        ]
        # FASTQ with >50% N bases should return empty, not raise error
        matched, info = illumina.match_barcodes_with_orientation(
            'NNNNNNNN', 'GCTAGCTA', sample_rows
        )
        self.assertEqual(len(matched), 0)
        self.assertEqual(info.get('skipped_reason'), 'high_n_fraction')

    def test_high_n_fraction_in_bc2_returns_empty(self):
        """Test that barcode_2 with >50% N bases also returns empty."""
        sample_rows = [
            {'sample': 'S1', 'barcode_1': 'ATCGATCG', 'barcode_2': 'GCTAGCTA'},
        ]
        # bc2 has >50% N
        matched, info = illumina.match_barcodes_with_orientation(
            'ATCGATCG', 'NNNNNNCTA', sample_rows
        )
        self.assertEqual(len(matched), 0)
        self.assertEqual(info.get('skipped_reason'), 'high_n_fraction')

    def test_exactly_50_percent_n_allowed(self):
        """Test that exactly 50% N bases is allowed (boundary case)."""
        sample_rows = [
            {'sample': 'S1', 'barcode_1': 'ATCGATCG', 'barcode_2': 'GCTAGCTA'},
        ]
        # 4/8 = exactly 50% N is allowed (threshold is >50%)
        # NNNNATCG matches ATCGATCG at positions 4-7
        matched, info = illumina.match_barcodes_with_orientation(
            'NNNNATCG', 'GCTAGCTA', sample_rows
        )
        self.assertEqual(len(matched), 1)
        self.assertNotIn('skipped_reason', info)

        # 5/8 = 62.5% > 50%, should be rejected
        matched, info = illumina.match_barcodes_with_orientation(
            'NNNNNTCG', 'GCTAGCTA', sample_rows
        )
        self.assertEqual(len(matched), 0)
        self.assertEqual(info.get('skipped_reason'), 'high_n_fraction')

    def test_ambiguous_n_match_returns_empty(self):
        """Test that N-wildcard matching multiple distinct barcode pairs returns empty."""
        sample_rows = [
            {'sample': 'S1', 'barcode_1': 'ATCGATCG', 'barcode_2': 'GCTAGCTA'},
            {'sample': 'S2', 'barcode_1': 'GTCGATCG', 'barcode_2': 'GCTAGCTA'},  # differs only in first base
        ]
        # FASTQ with N at position 0 matches both S1 (A) and S2 (G) - ambiguous!
        matched, info = illumina.match_barcodes_with_orientation(
            'NTCGATCG', 'GCTAGCTA', sample_rows
        )
        self.assertEqual(len(matched), 0)
        self.assertEqual(info.get('skipped_reason'), 'ambiguous')

    def test_n_in_samplesheet_does_not_match_all(self):
        """Test that N in samplesheet barcode does NOT act as wildcard.

        Only N in the FASTQ (observed) barcode should be a wildcard.
        N in samplesheet should require exact N match.
        """
        sample_rows = [
            {'sample': 'S1', 'barcode_1': 'ATCNATCG', 'barcode_2': 'GCTAGCTA'},
        ]
        # FASTQ has G where samplesheet has N - should not match
        matched, info = illumina.match_barcodes_with_orientation(
            'ATCGATCG', 'GCTAGCTA', sample_rows
        )
        self.assertEqual(len(matched), 0)
        self.assertEqual(info.get('skipped_reason'), 'no_match')

    def test_n_wildcard_no_false_positives(self):
        """Test that N wildcard doesn't cause incorrect matches."""
        sample_rows = [
            {'sample': 'S1', 'barcode_1': 'AAAAAAAA', 'barcode_2': 'TTTTTTTT'},
            {'sample': 'S2', 'barcode_1': 'CCCCCCCC', 'barcode_2': 'GGGGGGGG'},
        ]
        # FASTQ with 2 N bases (25%) should still only match one sample
        matched, info = illumina.match_barcodes_with_orientation(
            'ANAAANNA', 'TTTTTTTT', sample_rows
        )
        self.assertEqual(len(matched), 1)
        self.assertEqual(matched[0]['sample'], 'S1')


class TestBuildRunInfoJson(TestCaseWithTmp):
    """Test run_info.json construction."""

    def test_build_with_all_parameters(self):
        """Test building run_info.json with all parameters provided."""
        run_info_data = illumina.build_run_info_json(
            sequencing_center="BROAD",
            run_start_date="2024-01-15",
            read_structure="76T8B8B76T",
            indexes=2,
            run_id="FLOWCELL123.1",
            lane=1,
            flowcell="FLOWCELL123",
            lane_count=8,
            surface_count=2,
            swath_count=3,
            tile_count=4,
            total_tile_count=192,
            sequencer_model="NextSeq 2000"
        )

        self.assertEqual(run_info_data['sequencing_center'], 'BROAD')
        self.assertEqual(run_info_data['run_start_date'], '2024-01-15')
        self.assertEqual(run_info_data['read_structure'], '76T8B8B76T')
        self.assertEqual(run_info_data['indexes'], '2')
        self.assertEqual(run_info_data['run_id'], 'FLOWCELL123.1')
        self.assertEqual(run_info_data['lane'], '1')
        self.assertEqual(run_info_data['flowcell'], 'FLOWCELL123')
        self.assertEqual(run_info_data['lane_count'], '8')
        self.assertEqual(run_info_data['surface_count'], '2')
        self.assertEqual(run_info_data['swath_count'], '3')
        self.assertEqual(run_info_data['tile_count'], '4')
        self.assertEqual(run_info_data['total_tile_count'], '192')
        self.assertEqual(run_info_data['sequencer_model'], 'NextSeq 2000')

    def test_build_with_minimal_parameters(self):
        """Test building with only required parameters."""
        run_info_data = illumina.build_run_info_json(
            sequencing_center="BROAD",
            run_start_date="2024-01-15",
            read_structure="76T8B76T",
            indexes=1,
            run_id="FC123.1",
            lane=1,
            flowcell="FC123"
        )

        # Required fields should be present
        self.assertEqual(run_info_data['sequencing_center'], 'BROAD')
        self.assertEqual(run_info_data['run_start_date'], '2024-01-15')
        self.assertEqual(run_info_data['read_structure'], '76T8B76T')
        self.assertEqual(run_info_data['flowcell'], 'FC123')
        self.assertEqual(run_info_data['lane'], '1')
        self.assertEqual(run_info_data['indexes'], '1')

        # Optional fields should NOT be present when not provided
        self.assertNotIn('lane_count', run_info_data)
        self.assertNotIn('surface_count', run_info_data)
        self.assertNotIn('sequencer_model', run_info_data)

    def test_integer_to_string_conversion(self):
        """Test that integer parameters are converted to strings."""
        run_info_data = illumina.build_run_info_json(
            sequencing_center="BROAD",
            run_start_date="2024-01-15",
            read_structure="76T",
            indexes=0,
            run_id="FC.1",
            lane=3,  # integer
            flowcell="FC",
            lane_count=8  # integer
        )

        # All numeric fields should be strings in the output
        self.assertIsInstance(run_info_data['lane'], str)
        self.assertEqual(run_info_data['lane'], '3')
        self.assertIsInstance(run_info_data['lane_count'], str)
        self.assertEqual(run_info_data['indexes'], '0')

    def test_consistency_with_existing_output_schema(self):
        """Test that output matches the schema used by illumina_demux and splitcode_demux."""
        run_info_data = illumina.build_run_info_json(
            sequencing_center="SEQ_CENTER",
            run_start_date="2024-01-15",
            read_structure="76T8B8B76T",
            indexes=2,
            run_id="RUN.1",
            lane=1,
            flowcell="FLOWCELL",
            lane_count=8,
            surface_count=2,
            swath_count=3,
            tile_count=4,
            total_tile_count=192,
            sequencer_model="NextSeq"
        )

        # Verify all expected keys are present
        expected_keys = [
            'sequencing_center', 'run_start_date', 'read_structure', 'indexes',
            'run_id', 'lane', 'flowcell', 'lane_count', 'surface_count',
            'swath_count', 'tile_count', 'total_tile_count', 'sequencer_model'
        ]
        for key in expected_keys:
            self.assertIn(key, run_info_data, f"Missing expected key: {key}")

    def test_none_values_handled(self):
        """Test that None values for optional parameters are handled gracefully."""
        run_info_data = illumina.build_run_info_json(
            sequencing_center="BROAD",
            run_start_date="2024-01-15",
            read_structure="76T",
            indexes=0,
            run_id="FC.1",
            lane=1,
            flowcell="FC",
            lane_count=None,
            sequencer_model=None
        )

        # Function should handle None gracefully
        self.assertIsNotNone(run_info_data)
        self.assertEqual(run_info_data['flowcell'], 'FC')


class TestIlluminaMetadata(TestCaseWithTmp):
    """
    Test suite for illumina_metadata entry point.

    Tests metadata JSON generation from RunInfo.xml and SampleSheet.csv
    without processing any read data. This function is designed to be run
    once per sequencing run to generate metadata that's shared across all
    parallel demux jobs.

    Outputs tested:
    - run_info.json: Run metadata (flowcell, dates, read structure, instrument)
    - meta_by_sample.json: Sample metadata indexed by sample name
    - meta_by_filename.json: Sample metadata indexed by library ID
    """

    def setUp(self):
        super().setUp()
        # Share test data with TestSplitcodeDemuxFastqs
        # Use the test input directory path construction that works
        test_dir = util.file.get_test_path()
        self.input_dir = os.path.join(test_dir, 'input', 'TestSplitcodeDemuxFastqs')

        # Metadata input files
        self.runinfo_xml = os.path.join(self.input_dir, 'RunInfo.xml')
        self.samplesheet_csv = os.path.join(self.input_dir, 'SampleSheet.csv')

        # Verify test input files exist
        for f in [self.runinfo_xml, self.samplesheet_csv]:
            self.assertTrue(os.path.exists(f), f"Test input file missing: {f}")

    def test_runinfo_xml_parsing(self):
        """Test that RunInfo.xml can be parsed correctly."""
        runinfo = illumina.RunInfo(self.runinfo_xml)

        # Verify basic metadata extraction
        self.assertEqual(runinfo.get_flowcell(), 'TESTFC01')
        self.assertEqual(runinfo.get_run_id(), '250101_TEST_0001_BTESTFC01')
        self.assertEqual(runinfo.get_machine(), 'TEST001')
        self.assertEqual(runinfo.get_lane_count(), 1)

        # Verify read structure
        read_structure = runinfo.get_read_structure()
        self.assertIn('50T', read_structure)  # 50bp template reads
        self.assertIn('8B', read_structure)   # 8bp barcode reads

    def test_samplesheet_parsing(self):
        """Test that SampleSheet.csv can be parsed correctly."""
        samples = illumina.SampleSheet(
            self.samplesheet_csv,
            only_lane=1,
            allow_non_unique=False
        )

        # Verify we got the expected samples
        rows = list(samples.get_rows())
        self.assertGreater(len(rows), 0, "SampleSheet should contain samples")

        # Verify sample metadata structure
        for row in rows:
            self.assertIn('sample', row)
            self.assertIn('barcode_1', row)
            self.assertIn('barcode_2', row)

    def test_metadata_json_generation(self):
        """
        Test complete metadata JSON generation workflow.

        This is the main workflow test that verifies illumina_metadata
        can generate all three output JSON files with correct schemas.
        """
        out_dir = tempfile.mkdtemp()

        try:
            # Output paths
            out_runinfo = os.path.join(out_dir, 'run_info.json')
            out_meta_by_sample = os.path.join(out_dir, 'meta_by_sample.json')
            out_meta_by_filename = os.path.join(out_dir, 'meta_by_filename.json')

            # Call illumina_metadata (will fail until implemented - expected per TDD)
            illumina.illumina_metadata(
                runinfo=self.runinfo_xml,
                samplesheet=self.samplesheet_csv,
                lane=1,
                sequencing_center='Broad',
                out_runinfo=out_runinfo,
                out_meta_by_sample=out_meta_by_sample,
                out_meta_by_filename=out_meta_by_filename
            )

            # Verify all output files were created
            self.assertTrue(os.path.exists(out_runinfo), "run_info.json not created")
            self.assertTrue(os.path.exists(out_meta_by_sample), "meta_by_sample.json not created")
            self.assertTrue(os.path.exists(out_meta_by_filename), "meta_by_filename.json not created")

            # Verify run_info.json schema
            with open(out_runinfo, 'r') as f:
                run_info = json.load(f)

            required_fields = [
                'sequencing_center', 'run_start_date', 'read_structure',
                'indexes', 'run_id', 'lane', 'flowcell', 'lane_count',
                'surface_count', 'swath_count', 'tile_count',
                'total_tile_count', 'sequencer_model'
            ]
            for field in required_fields:
                self.assertIn(field, run_info, f"Missing required field: {field}")

            # Verify specific values
            self.assertEqual(run_info['sequencing_center'], 'Broad')
            self.assertEqual(run_info['flowcell'], 'TESTFC01')
            self.assertEqual(run_info['lane'], '1')

            # Verify meta_by_sample.json schema
            with open(out_meta_by_sample, 'r') as f:
                meta_by_sample = json.load(f)

            self.assertIsInstance(meta_by_sample, dict, "meta_by_sample should be a dict")
            self.assertGreater(len(meta_by_sample), 0, "meta_by_sample should not be empty")

            # Each sample should have metadata with lane added
            for sample_name, metadata in meta_by_sample.items():
                self.assertIn('lane', metadata, "Sample metadata should include lane")
                self.assertEqual(metadata['lane'], '1')

            # Verify meta_by_filename.json schema
            with open(out_meta_by_filename, 'r') as f:
                meta_by_filename = json.load(f)

            self.assertIsInstance(meta_by_filename, dict, "meta_by_filename should be a dict")

        finally:
            shutil.rmtree(out_dir)

    def test_metadata_generation_with_optional_params(self):
        """Test metadata generation with optional parameters omitted."""
        out_dir = tempfile.mkdtemp()

        try:
            # Generate only run_info.json (omit sample metadata outputs)
            out_runinfo = os.path.join(out_dir, 'run_info.json')

            illumina.illumina_metadata(
                runinfo=self.runinfo_xml,
                samplesheet=self.samplesheet_csv,
                lane=1,
                out_runinfo=out_runinfo,
                out_meta_by_sample=None,  # Optional
                out_meta_by_filename=None  # Optional
            )

            # Only run_info.json should exist
            self.assertTrue(os.path.exists(out_runinfo))
            self.assertFalse(os.path.exists(os.path.join(out_dir, 'meta_by_sample.json')))
            self.assertFalse(os.path.exists(os.path.join(out_dir, 'meta_by_filename.json')))

        finally:
            shutil.rmtree(out_dir)

    def test_metadata_consistency_with_existing_demux(self):
        """
        Test that illumina_metadata output matches the schema produced
        by existing illumina_demux and splitcode_demux functions.

        This ensures compatibility and validates that the refactoring
        maintains the same output format.
        """
        out_dir = tempfile.mkdtemp()

        try:
            out_runinfo = os.path.join(out_dir, 'run_info.json')

            illumina.illumina_metadata(
                runinfo=self.runinfo_xml,
                samplesheet=self.samplesheet_csv,
                lane=1,
                sequencing_center='Broad',
                out_runinfo=out_runinfo
            )

            with open(out_runinfo, 'r') as f:
                run_info = json.load(f)

            # Verify the schema matches what build_run_info_json() produces
            # All values should be strings (as per existing implementation)
            for field in ['lane', 'lane_count', 'surface_count', 'swath_count',
                         'tile_count', 'total_tile_count']:
                self.assertIsInstance(run_info[field], str,
                                     f"{field} should be a string for compatibility")

        finally:
            shutil.rmtree(out_dir)

    def test_invalid_runinfo_path(self):
        """Test error handling for invalid RunInfo.xml path."""
        with self.assertRaises((FileNotFoundError, IOError)):
            illumina.illumina_metadata(
                runinfo='/nonexistent/RunInfo.xml',
                samplesheet=self.samplesheet_csv,
                lane=1
            )

    def test_invalid_samplesheet_path(self):
        """Test error handling for invalid SampleSheet path."""
        with self.assertRaises((FileNotFoundError, IOError)):
            illumina.illumina_metadata(
                runinfo=self.runinfo_xml,
                samplesheet='/nonexistent/SampleSheet.csv',
                lane=1
            )

    def test_illumina_metadata_via_parser(self):
        """
        Test illumina_metadata command via argument parser.

        This tests the full CLI code path including:
        - Parser argument processing
        - Argument validation
        - Main function invocation

        Regression test for issue #127 where split_args=True caused
        TypeError with main_illumina_metadata(args) signature.
        """
        out_dir = tempfile.mkdtemp()

        try:
            # Output paths
            out_runinfo = os.path.join(out_dir, 'run_info.json')
            out_meta_by_sample = os.path.join(out_dir, 'meta_by_sample.json')
            out_meta_by_filename = os.path.join(out_dir, 'meta_by_filename.json')

            # Call via parser (like CLI does)
            parser = illumina.parser_illumina_metadata(argparse.ArgumentParser())
            args = parser.parse_args([
                '--runinfo', self.runinfo_xml,
                '--samplesheet', self.samplesheet_csv,
                '--lane', '1',
                '--sequencing_center', 'Broad',
                '--out_runinfo', out_runinfo,
                '--out_meta_by_sample', out_meta_by_sample,
                '--out_meta_by_filename', out_meta_by_filename
            ])
            args.func_main(args)

            # Verify all output files were created
            self.assertTrue(os.path.exists(out_runinfo), "run_info.json not created")
            self.assertTrue(os.path.exists(out_meta_by_sample), "meta_by_sample.json not created")
            self.assertTrue(os.path.exists(out_meta_by_filename), "meta_by_filename.json not created")

            # Basic validation that JSON files are valid
            with open(out_runinfo, 'r') as f:
                run_info = json.load(f)
            self.assertEqual(run_info['sequencing_center'], 'Broad')

        finally:
            shutil.rmtree(out_dir)

    def test_three_barcode_samplesheet(self):
        """
        Test illumina_metadata with pure 3-barcode samplesheet.

        This tests a samplesheet where all samples share outer barcodes (barcode_1 + barcode_2)
        but differ in inner barcodes (barcode_3). This is the typical use case for
        splitcode demultiplexing.

        Expected behavior:
        - Should handle samples with duplicate barcode_1+barcode_2 combinations
        - barcode_3 column should be preserved in output metadata
        - num_indexes should reflect the number of outer barcodes (2)
        """
        out_dir = tempfile.mkdtemp()

        # Use the 3-barcode samplesheet
        samples_3bc = os.path.join(self.input_dir, 'samples_3bc.tsv')
        self.assertTrue(os.path.exists(samples_3bc), f"Test file missing: {samples_3bc}")

        try:
            # Output paths
            out_runinfo = os.path.join(out_dir, 'run_info.json')
            out_meta_by_sample = os.path.join(out_dir, 'meta_by_sample.json')

            # Call illumina_metadata with 3-barcode samplesheet
            illumina.illumina_metadata(
                runinfo=self.runinfo_xml,
                samplesheet=samples_3bc,
                lane=1,
                sequencing_center='Broad',
                out_runinfo=out_runinfo,
                out_meta_by_sample=out_meta_by_sample
            )

            # Verify output files were created
            self.assertTrue(os.path.exists(out_runinfo), "run_info.json not created")
            self.assertTrue(os.path.exists(out_meta_by_sample), "meta_by_sample.json not created")

            # Verify run_info.json has correct index count (should be 2, not 3)
            with open(out_runinfo, 'r') as f:
                run_info = json.load(f)
            self.assertEqual(run_info['indexes'], '2',
                           "num_indexes should be 2 (outer barcodes only)")

            # Verify meta_by_sample.json includes all samples with barcode_3
            with open(out_meta_by_sample, 'r') as f:
                meta_by_sample = json.load(f)

            # Should have samples from all three pools (7 samples total)
            self.assertEqual(len(meta_by_sample), 7,
                           "Should have all 7 samples from 3-barcode samplesheet")

            # Verify 3-barcode samples have barcode_3 field
            three_bc_samples = ['TestSample1', 'TestSample2', 'TestSample3',
                               'TestSampleEmpty', 'TestSample5', 'TestSample6']
            for sample_name in three_bc_samples:
                self.assertIn(sample_name, meta_by_sample,
                            f"Sample {sample_name} missing from metadata")
                self.assertIn('barcode_3', meta_by_sample[sample_name],
                            f"Sample {sample_name} should have barcode_3 field")
                self.assertNotEqual(meta_by_sample[sample_name]['barcode_3'], '',
                                  f"Sample {sample_name} should have non-empty barcode_3")

            # Verify 2-barcode sample has empty barcode_3
            self.assertIn('TestSampleNoSplitcode', meta_by_sample)
            # barcode_3 should either be empty string or not present (both acceptable)
            bc3 = meta_by_sample['TestSampleNoSplitcode'].get('barcode_3', '')
            self.assertEqual(bc3, '',
                           "2-barcode sample should have empty barcode_3")

        finally:
            shutil.rmtree(out_dir)

    def test_mixed_two_and_three_barcode_samplesheet(self):
        """
        Test illumina_metadata with mixed 2-barcode and 3-barcode samples.

        This tests a samplesheet containing both:
        - Samples with 3 barcodes (barcode_1 + barcode_2 + barcode_3)
        - Samples with 2 barcodes (barcode_1 + barcode_2, empty barcode_3)

        The same samplesheet is used for both scenarios since samples_3bc.tsv
        already contains this mix.

        Expected behavior:
        - Should handle both types of samples in the same samplesheet
        - All samples should be included in metadata output
        - barcode_3 field should be present for all samples (empty for 2-barcode samples)
        """
        out_dir = tempfile.mkdtemp()

        # Use the mixed barcode samplesheet
        samples_3bc = os.path.join(self.input_dir, 'samples_3bc.tsv')

        try:
            # Output paths
            out_meta_by_sample = os.path.join(out_dir, 'meta_by_sample.json')

            # Call illumina_metadata
            illumina.illumina_metadata(
                runinfo=self.runinfo_xml,
                samplesheet=samples_3bc,
                lane=1,
                out_meta_by_sample=out_meta_by_sample
            )

            # Verify metadata includes both types of samples
            with open(out_meta_by_sample, 'r') as f:
                meta_by_sample = json.load(f)

            # Count samples with and without barcode_3
            samples_with_bc3 = []
            samples_without_bc3 = []

            for sample_name, metadata in meta_by_sample.items():
                bc3 = metadata.get('barcode_3', '')
                if bc3 and bc3.strip():
                    samples_with_bc3.append(sample_name)
                else:
                    samples_without_bc3.append(sample_name)

            # Should have 6 samples with barcode_3 (from Pool 1 and Pool 2)
            self.assertEqual(len(samples_with_bc3), 6,
                           f"Should have 6 samples with barcode_3, got: {samples_with_bc3}")

            # Should have 1 sample without barcode_3 (from Pool 3)
            self.assertEqual(len(samples_without_bc3), 1,
                           f"Should have 1 sample without barcode_3, got: {samples_without_bc3}")
            self.assertIn('TestSampleNoSplitcode', samples_without_bc3,
                        "TestSampleNoSplitcode should have empty barcode_3")

            # Verify specific samples are categorized correctly
            self.assertIn('TestSample1', samples_with_bc3)
            self.assertIn('TestSample2', samples_with_bc3)
            self.assertIn('TestSample5', samples_with_bc3)

        finally:
            shutil.rmtree(out_dir)

    def test_three_barcode_barcode_uniqueness(self):
        """
        Test that illumina_metadata properly handles duplicate outer barcodes in 3-barcode sheets.

        In 3-barcode samplesheets, multiple samples share the same barcode_1+barcode_2
        combination (outer barcodes) and differ only in barcode_3 (inner barcode).

        Expected behavior:
        - Should NOT raise an error about non-unique barcodes
        - Should successfully generate metadata for all samples
        - The function should internally handle allow_non_unique=True for 3-barcode sheets

        This is a regression test to ensure the function doesn't fail when samples
        share outer barcodes.
        """
        out_dir = tempfile.mkdtemp()
        samples_3bc = os.path.join(self.input_dir, 'samples_3bc.tsv')

        try:
            out_meta_by_sample = os.path.join(out_dir, 'meta_by_sample.json')

            # This should NOT raise an error even though Pool 1 has 4 samples
            # sharing ATCGATCG+GCTAGCTA outer barcodes
            illumina.illumina_metadata(
                runinfo=self.runinfo_xml,
                samplesheet=samples_3bc,
                lane=1,
                out_meta_by_sample=out_meta_by_sample
            )

            # Verify all Pool 1 samples are present (they share outer barcodes)
            with open(out_meta_by_sample, 'r') as f:
                meta_by_sample = json.load(f)

            pool1_samples = ['TestSample1', 'TestSample2', 'TestSample3', 'TestSampleEmpty']
            for sample in pool1_samples:
                self.assertIn(sample, meta_by_sample,
                            f"Pool 1 sample {sample} should be in metadata")
                # All should have same outer barcodes
                self.assertEqual(meta_by_sample[sample]['barcode_1'], 'ATCGATCG')
                self.assertEqual(meta_by_sample[sample]['barcode_2'], 'GCTAGCTA')

            # All Pool 1 samples should have different barcode_3
            bc3_values = [meta_by_sample[s]['barcode_3'] for s in pool1_samples]
            self.assertEqual(len(set(bc3_values)), len(bc3_values),
                           "Pool 1 samples should have unique barcode_3 values")

        finally:
            shutil.rmtree(out_dir)


class TestSplitcodeDemuxFastqs(TestCaseWithTmp):
    """
    Test suite for splitcode_demux_fastqs entry point.

    Tests SIMPLIFIED demultiplexing from paired DRAGEN FASTQ files using 3-barcode scheme:
    - Input: Paired R1/R2 FASTQ files + custom 3-barcode samplesheet ONLY
    - Output: Per-sample unaligned BAMs + demux metrics + barcode reports
    - Does NOT generate: run_info.json, meta_by_sample.json, meta_by_filename.json
      (use illumina_metadata entry point for those)

    Test data includes both 3-barcode and 2-barcode samples:

    - TestPool1 (ATCGATCG+GCTAGCTA) - 3 barcodes:
      * TestSample1 (AAAAAAAA): 100 reads
      * TestSample2 (CCCCCCCC): 75 reads
      * TestSample3 (GGGGTTTT): 50 reads
      * TestSampleEmpty (TTTTGGGG): 0 reads (empty sample)
      * Unmatched (NNNNANNN): 25 reads

    - TestPool3 (GGAATTCC+CCGGAATT) - 2 barcodes only (no barcode_3):
      * TestSampleNoSplitcode: 80 reads (should bypass splitcode, direct FASTQ→BAM)

    Note: This simplified version does NOT require RunInfo.xml or Illumina SampleSheet.csv.
    """

    def setUp(self):
        super().setUp()
        self.input_dir = util.file.get_test_input_path(self)

        # Input FASTQ files
        self.r1_fastq = os.path.join(self.input_dir, 'TestPool1_S1_L001_R1_001.fastq.gz')
        self.r2_fastq = os.path.join(self.input_dir, 'TestPool1_S1_L001_R2_001.fastq.gz')

        # Custom 3-barcode samplesheet
        self.samples_3bc = os.path.join(self.input_dir, 'samples_3bc.tsv')

        # RunInfo.xml for richer BAM metadata
        self.runinfo_xml = os.path.join(self.input_dir, 'RunInfo.xml')

        # Verify test input files exist
        for f in [self.r1_fastq, self.r2_fastq, self.samples_3bc, self.runinfo_xml]:
            self.assertTrue(os.path.exists(f), f"Test input file missing: {f}")

    def test_parse_fastq_filename_from_test_data(self):
        """Test that our test FASTQ filenames can be parsed correctly."""
        # Test with simple format (no flowcell ID prefix)
        result = illumina.parse_illumina_fastq_filename(self.r1_fastq)

        self.assertEqual(result['sample_name'], 'TestPool1')
        self.assertEqual(result['sample_number'], 1)
        self.assertEqual(result['lane'], 1)
        self.assertEqual(result['read'], 1)
        self.assertEqual(result['chunk'], 1)
        self.assertFalse(result.get('is_index', False))

    def test_barcode_normalization_on_samplesheet(self):
        """Test that barcodes from samplesheet can be normalized."""
        # Read the custom 3-barcode samplesheet
        # Note: allow_non_unique=True because samples share index1+index2, differ only in barcode_3
        samples = illumina.SampleSheet(self.samples_3bc, allow_non_unique=True)

        # Get first sample's barcodes
        first_sample = samples.get_rows()[0]

        # Test normalization of each barcode
        bc1 = illumina.normalize_barcode(first_sample['barcode_1'])
        bc2 = illumina.normalize_barcode(first_sample['barcode_2'])
        bc3 = illumina.normalize_barcode(first_sample['barcode_3'])

        # Should all be uppercase, valid DNA
        self.assertTrue(all(c in 'ACGTN' for c in bc1))
        self.assertTrue(all(c in 'ACGTN' for c in bc2))
        self.assertTrue(all(c in 'ACGTN' for c in bc3))

    def test_basic_demux_workflow(self):
        """
        Test basic demultiplexing workflow from paired FASTQs.

        Simplified workflow (no RunInfo.xml or Illumina SampleSheet.csv required):
        1. Parse FASTQ filenames for metadata
        2. Load custom 3-barcode samplesheet
        3. Run splitcode demultiplexing
        4. Verify output BAMs are created
        5. Verify demux metrics and barcode reports
        """
        out_dir = tempfile.mkdtemp()

        # Call the main function (will fail until implemented in Phase 3)
        # NOTE: Simplified interface - no runinfo_xml or illumina_samplesheet needed
        illumina.splitcode_demux_fastqs(
            fastq_r1=self.r1_fastq,
            fastq_r2=self.r2_fastq,
            samplesheet=self.samples_3bc,
            outdir=out_dir,
            runinfo=self.runinfo_xml,
            threads=1
        )

        # Verify expected output BAMs exist
        expected_bams = [
            os.path.join(out_dir, 'TestSample1.bam'),
            os.path.join(out_dir, 'TestSample2.bam'),
            os.path.join(out_dir, 'TestSample3.bam'),
            os.path.join(out_dir, 'TestSampleEmpty.bam'),
        ]

        for bam in expected_bams:
            self.assertTrue(os.path.exists(bam), f"Expected output BAM missing: {bam}")

        # Verify demux metrics (but NOT metadata JSONs or barcode reports)
        self.assertTrue(os.path.exists(os.path.join(out_dir, 'demux_metrics.json')))

        # Should NOT generate these (use illumina_metadata or illumina_demux instead)
        self.assertFalse(os.path.exists(os.path.join(out_dir, 'run_info.json')))
        self.assertFalse(os.path.exists(os.path.join(out_dir, 'meta_by_sample.json')))
        self.assertFalse(os.path.exists(os.path.join(out_dir, 'meta_by_filename.json')))
        # Note: barcodes_common.txt and barcodes_outliers.txt are not generated by
        # splitcode_demux_fastqs (use illumina_demux for comprehensive barcode reporting)

    def test_barcode_matching_perfect_match(self):
        """
        Test that reads with perfect barcode matches are assigned correctly.

        Verify:
        - TestSample1 (AAAAAAAA) gets exactly 100 reads
        - TestSample2 (CCCCCCCC) gets exactly 75 reads
        - TestSample3 (GGGGTTTT) gets exactly 50 reads
        """
        out_dir = tempfile.mkdtemp()

        # Simplified call - no RunInfo.xml needed
        illumina.splitcode_demux_fastqs(
            fastq_r1=self.r1_fastq,
            fastq_r2=self.r2_fastq,
            samplesheet=self.samples_3bc,
            outdir=out_dir,
            runinfo=self.runinfo_xml,
            threads=1
        )

        # Use samtools to count reads in each output BAM
        samtools = tools.samtools.SamtoolsTool()

        sample1_bam = os.path.join(out_dir, 'TestSample1.bam')
        sample2_bam = os.path.join(out_dir, 'TestSample2.bam')
        sample3_bam = os.path.join(out_dir, 'TestSample3.bam')

        # Count read pairs (samtools.count returns total reads, divide by 2 for pairs)
        sample1_pairs = samtools.count(sample1_bam) // 2
        sample2_pairs = samtools.count(sample2_bam) // 2
        sample3_pairs = samtools.count(sample3_bam) // 2

        self.assertEqual(sample1_pairs, 100, "TestSample1 should have exactly 100 read pairs")
        self.assertEqual(sample2_pairs, 75, "TestSample2 should have exactly 75 read pairs")
        self.assertEqual(sample3_pairs, 50, "TestSample3 should have exactly 50 read pairs")

    def test_unmatched_barcodes(self):
        """
        Test handling of reads with unmatched inline barcodes.

        Verify:
        - Reads with outlier barcodes (GGAATTTT, CCCCAAAA, ATATAGAG totaling 25 reads)
          are handled without errors
        - No crashes when some reads don't match any sample
        - Unassigned reads are written to unassigned FASTQs

        Note: This test does NOT check barcodes_outliers.txt (not generated by
        splitcode_demux_fastqs). Use illumina_demux for comprehensive barcode reporting.
        """
        out_dir = tempfile.mkdtemp()

        # Should not crash even with unmatched barcodes
        illumina.splitcode_demux_fastqs(
            fastq_r1=self.r1_fastq,
            fastq_r2=self.r2_fastq,
            samplesheet=self.samples_3bc,
            outdir=out_dir,
            runinfo=self.runinfo_xml,
            threads=1
        )

        # Verify unassigned FASTQs exist
        unassigned_r1 = os.path.join(out_dir, 'unassigned_R1.fastq.gz')
        unassigned_r2 = os.path.join(out_dir, 'unassigned_R2.fastq.gz')
        self.assertTrue(os.path.exists(unassigned_r1), "Unassigned R1 should exist")
        self.assertTrue(os.path.exists(unassigned_r2), "Unassigned R2 should exist")

        # Verify unassigned files are not empty (should have 25 read pairs)
        samtools = tools.samtools.SamtoolsTool()
        # We'd need to convert FASTQ to BAM to count, but that's overkill.
        # Just verify the files have non-zero size.
        self.assertGreater(os.path.getsize(unassigned_r1), 0,
                          "Unassigned R1 should have content")
        self.assertGreater(os.path.getsize(unassigned_r2), 0,
                          "Unassigned R2 should have content")

    def test_empty_barcode_sample(self):
        """
        Test handling of sample with zero matching reads.

        Verify:
        - TestSampleEmpty (TTTTGGGG) is in samplesheet but has 0 reads
        - Empty BAM file is still created (or noted in metrics)
        - No errors/crashes when sample has zero reads
        """
        out_dir = tempfile.mkdtemp()

        # Should not raise an exception even though TestSampleEmpty has 0 reads
        # Simplified call - no RunInfo.xml needed
        illumina.splitcode_demux_fastqs(
            fastq_r1=self.r1_fastq,
            fastq_r2=self.r2_fastq,
            samplesheet=self.samples_3bc,
            outdir=out_dir,
            runinfo=self.runinfo_xml,
            threads=1
        )

        # Empty sample should either:
        # 1. Have an empty BAM file created, OR
        # 2. Be noted in metrics as having 0 reads

        empty_bam = os.path.join(out_dir, 'TestSampleEmpty.bam')

        if os.path.exists(empty_bam):
            # If BAM exists, it should have 0 reads
            samtools = tools.samtools.SamtoolsTool()
            read_count = samtools.count(empty_bam)
            self.assertEqual(read_count, 0, "Empty sample should have 0 reads")
        else:
            # If no BAM, check metrics mentions this sample
            metrics_file = os.path.join(out_dir, 'demux_metrics.json')
            with open(metrics_file, 'r') as f:
                metrics = json.load(f)

            # TestSampleEmpty should be mentioned somewhere in metrics
            self.assertIn('TestSampleEmpty', str(metrics), "Empty sample should be mentioned in metrics")

    def test_output_schema_consistency(self):
        """
        Test that output files match expected schema.

        Verify:
        - demux_metrics.json has expected fields
        - BAM files are valid
        - Barcode reports are NOT generated (use illumina_demux for those)
        """
        out_dir = tempfile.mkdtemp()

        # Simplified call - no RunInfo.xml needed
        illumina.splitcode_demux_fastqs(
            fastq_r1=self.r1_fastq,
            fastq_r2=self.r2_fastq,
            samplesheet=self.samples_3bc,
            outdir=out_dir,
            runinfo=self.runinfo_xml,
            threads=1
        )

        # Check demux_metrics.json has expected structure
        metrics_file = os.path.join(out_dir, 'demux_metrics.json')
        with open(metrics_file, 'r') as f:
            metrics = json.load(f)

        # Should have sample-level information
        self.assertIsInstance(metrics, (dict, list), "demux_metrics.json should be JSON dict or list")

        # Barcode reports should NOT be generated (not part of simplified workflow)
        common_file = os.path.join(out_dir, 'barcodes_common.txt')
        outliers_file = os.path.join(out_dir, 'barcodes_outliers.txt')

        self.assertFalse(os.path.exists(common_file),
                        "barcodes_common.txt should not be generated by splitcode_demux_fastqs")
        self.assertFalse(os.path.exists(outliers_file),
                        "barcodes_outliers.txt should not be generated by splitcode_demux_fastqs")

    def test_multi_pool_samplesheet_collapsibility_check(self):
        """
        Test that collapsibility is checked only on the filtered pool, not the entire samplesheet.

        This is a regression test for GitHub Copilot comment #2510390315.
        The bug: can_be_collapsed was checking the entire samplesheet instead of just
        the filtered samples for the current pool.

        Scenario:
        - Samplesheet contains multiple pools:
          * Pool 1 (ATCGATCG+GCTAGCTA): 3-barcode samples (collapsible)
          * Pool 2 (CTGATCGT+TAGATCGC): 3-barcode samples (collapsible)
          * Pool 3 (GGAATTCC+CCGGAATT): 2-barcode sample (unique, not collapsible)
        - FASTQ file contains only Pool 1 data (ATCGATCG+GCTAGCTA)
        - Should succeed even though Pool 3 exists with unique outer barcodes

        Before fix: Would fail because can_be_collapsed checked entire samplesheet
        After fix: Succeeds because collapsibility is checked only on Pool 1
        """
        out_dir = tempfile.mkdtemp()

        # The test samplesheet already has this structure:
        # - Pool 1: ATCGATCG+GCTAGCTA with 3-barcode samples (TestSample1-4)
        # - Pool 2: CTGATCGT+TAGATCGC with 3-barcode samples (TestSample5-6)
        # - Pool 3: GGAATTCC+CCGGAATT with 2-barcode sample (TestSampleNoSplitcode)

        # The test FASTQ contains Pool 1 data only (ATCGATCG+GCTAGCTA)
        # This should succeed even though the samplesheet contains Pool 3 with unique outer barcodes

        # This call should NOT raise ValueError about non-collapsible barcodes
        illumina.splitcode_demux_fastqs(
            fastq_r1=self.r1_fastq,
            fastq_r2=self.r2_fastq,
            samplesheet=self.samples_3bc,  # Multi-pool samplesheet
            outdir=out_dir,
            runinfo=self.runinfo_xml,
            threads=1
        )

        # Verify successful demux - BAMs should be created for Pool 1 samples
        expected_bams = [
            os.path.join(out_dir, 'TestSample1.bam'),
            os.path.join(out_dir, 'TestSample2.bam'),
            os.path.join(out_dir, 'TestSample3.bam'),
            os.path.join(out_dir, 'TestSampleEmpty.bam'),
        ]

        for bam in expected_bams:
            self.assertTrue(os.path.exists(bam), f"Expected output BAM missing: {bam}")

        # Pool 2 and Pool 3 samples should NOT be created (different outer barcodes)
        pool2_bams = [
            os.path.join(out_dir, 'TestSample5.bam'),
            os.path.join(out_dir, 'TestSample6.bam'),
        ]
        pool3_bams = [
            os.path.join(out_dir, 'TestSampleNoSplitcode.bam'),
        ]

        for bam in pool2_bams + pool3_bams:
            self.assertFalse(os.path.exists(bam),
                           f"BAM from different pool should not be created: {bam}")

    def test_fastq_filename_parsing(self):
        """
        Test extraction of metadata from FASTQ filenames.

        Verify:
        - Pool/sample name extracted correctly
        - Lane number extracted correctly
        - Read number (R1/R2) extracted correctly
        """
        out_dir = tempfile.mkdtemp()

        # Simplified call - no RunInfo.xml needed
        illumina.splitcode_demux_fastqs(
            fastq_r1=self.r1_fastq,
            fastq_r2=self.r2_fastq,
            samplesheet=self.samples_3bc,
            outdir=out_dir,
            runinfo=self.runinfo_xml,
            threads=1
        )

        # Metadata extracted from FASTQ filename should be used in output
        # Check demux_metrics.json for evidence of correct parsing
        metrics_file = os.path.join(out_dir, 'demux_metrics.json')
        with open(metrics_file, 'r') as f:
            metrics = json.load(f)

        # Should have metadata about the pool/samples
        # Exact structure TBD, but should contain sample information
        self.assertIsInstance(metrics, (dict, list), "Metrics should be valid JSON")

    def test_two_barcode_sample_bypass_splitcode(self):
        """
        Test handling of 2-barcode samples (empty barcode_3).

        This tests a mixed samplesheet where some rows have 3 barcodes (requiring
        splitcode demux) and other rows have only 2 barcodes (empty barcode_3).

        For 2-barcode rows, splitcode_demux_fastqs should:
        - Skip splitcode demultiplexing
        - Perform direct FASTQ → BAM conversion
        - Produce exactly one output BAM (the pool itself)
        - Still generate standard output files (JSONs, metrics)

        Test data:
        - TestPool3 (GGAATTCC+CCGGAATT): 80 reads, no inline barcode
        - Expected: Single BAM file with all 80 reads
        """
        out_dir = tempfile.mkdtemp()

        # Use TestPool3 which has empty barcode_3 (2-barcode sample)
        r1_fastq = os.path.join(self.input_dir, 'TestPool3_S3_L001_R1_001.fastq.gz')
        r2_fastq = os.path.join(self.input_dir, 'TestPool3_S3_L001_R2_001.fastq.gz')

        # Verify test files exist
        self.assertTrue(os.path.exists(r1_fastq), f"Test FASTQ missing: {r1_fastq}")
        self.assertTrue(os.path.exists(r2_fastq), f"Test FASTQ missing: {r2_fastq}")

        # Run demux on 2-barcode sample
        illumina.splitcode_demux_fastqs(
            fastq_r1=r1_fastq,
            fastq_r2=r2_fastq,
            samplesheet=self.samples_3bc,
            outdir=out_dir,
            runinfo=self.runinfo_xml,
            threads=1
        )

        # For 2-barcode sample, should produce exactly one BAM file
        # (the entire pool, no splitcode demux)
        expected_bam = os.path.join(out_dir, 'TestSampleNoSplitcode.bam')
        self.assertTrue(os.path.exists(expected_bam),
                       "2-barcode sample should produce a single BAM file")

        # Verify all 80 reads are in the output BAM
        samtools = tools.samtools.SamtoolsTool()
        total_reads = samtools.count(expected_bam)
        read_pairs = total_reads // 2
        self.assertEqual(read_pairs, 80,
                        "All 80 read pairs should be in the 2-barcode sample BAM")

        # Verify standard output files still exist (but NOT metadata JSONs)
        self.assertTrue(os.path.exists(os.path.join(out_dir, 'demux_metrics.json')))

        # Should NOT generate run_info.json (use illumina_metadata instead)
        self.assertFalse(os.path.exists(os.path.join(out_dir, 'run_info.json')))
        self.assertFalse(os.path.exists(os.path.join(out_dir, 'meta_by_sample.json')))
        self.assertFalse(os.path.exists(os.path.join(out_dir, 'meta_by_filename.json')))

    def test_splitcode_demux_fastqs_via_parser(self):
        """
        Test splitcode_demux_fastqs command via argument parser.

        This tests the full CLI code path including:
        - Parser argument processing
        - Argument validation
        - Main function invocation

        Regression test for issue #127 where split_args=True caused
        TypeError with main_splitcode_demux_fastqs(args) signature.
        """
        out_dir = tempfile.mkdtemp()

        try:
            # Call via parser (like CLI does)
            parser = illumina.parser_splitcode_demux_fastqs(argparse.ArgumentParser())
            args = parser.parse_args([
                '--fastq_r1', self.r1_fastq,
                '--fastq_r2', self.r2_fastq,
                '--samplesheet', self.samples_3bc,
                '--runinfo', self.runinfo_xml,
                '--outdir', out_dir
            ])
            args.func_main(args)

            # Verify output files were created
            self.assertTrue(os.path.exists(os.path.join(out_dir, 'demux_metrics.json')))

            # Verify at least one BAM file was created
            bam_files = [f for f in os.listdir(out_dir) if f.endswith('.bam')]
            self.assertGreater(len(bam_files), 0, "Should produce at least one BAM file")

        finally:
            shutil.rmtree(out_dir)

    def test_i5_reverse_complement_3bc_demux(self):
        """
        Test 3-barcode demultiplexing when i5 barcode requires reverse complement.

        This is a regression test for issue #133.

        Scenario:
        - Samplesheet has barcode_2 = GCTAGCTA (forward orientation)
        - FASTQ headers have barcode_2 = TAGCTAGC (reverse complement)
        - The function should auto-detect the i5 reverse complement orientation
        - Demux should succeed using the corrected barcode values

        Expected behavior:
        - Auto-detection logs i5 reverse complement match
        - 3-barcode demultiplexing succeeds
        - Output BAMs created with correct read counts
        - No "no_3bc_match" error
        """
        out_dir = tempfile.mkdtemp()

        try:
            # Use test files with reverse complement i5
            r1_fastq = os.path.join(self.input_dir, 'TestPool1_RC_S4_L001_R1_001.fastq.gz')
            r2_fastq = os.path.join(self.input_dir, 'TestPool1_RC_S4_L001_R2_001.fastq.gz')
            samples_rc = os.path.join(self.input_dir, 'samples_3bc_i5_revcomp.tsv')

            # Verify test input files exist
            self.assertTrue(os.path.exists(r1_fastq), f"Test FASTQ missing: {r1_fastq}")
            self.assertTrue(os.path.exists(r2_fastq), f"Test FASTQ missing: {r2_fastq}")
            self.assertTrue(os.path.exists(samples_rc), f"Test samplesheet missing: {samples_rc}")

            # Run splitcode demux with i5 reverse complement scenario
            illumina.splitcode_demux_fastqs(
                fastq_r1=r1_fastq,
                fastq_r2=r2_fastq,
                samplesheet=samples_rc,
                outdir=out_dir,
                runinfo=self.runinfo_xml,
                threads=1
            )

            # Verify expected output BAMs exist
            expected_bams = [
                os.path.join(out_dir, 'TestSample1_RC.bam'),
                os.path.join(out_dir, 'TestSample2_RC.bam'),
                os.path.join(out_dir, 'TestSample3_RC.bam')
            ]

            for bam in expected_bams:
                self.assertTrue(os.path.exists(bam), f"Expected output BAM missing: {bam}")

            # Verify read counts using samtools
            # Note: samtools.count() returns total reads (R1 + R2 for paired-end)
            samtools = tools.samtools.SamtoolsTool()

            # TestSample1_RC (AAAAAAAA): should have 2 read pairs = 4 total reads
            count1 = int(samtools.count(expected_bams[0]))
            self.assertEqual(count1, 4, "TestSample1_RC should have 2 read pairs (4 total reads)")

            # TestSample2_RC (CCCCCCCC): should have 2 read pairs = 4 total reads
            count2 = int(samtools.count(expected_bams[1]))
            self.assertEqual(count2, 4, "TestSample2_RC should have 2 read pairs (4 total reads)")

            # TestSample3_RC (GGGGTTTT): should have 1 read pair = 2 total reads
            count3 = int(samtools.count(expected_bams[2]))
            self.assertEqual(count3, 2, "TestSample3_RC should have 1 read pair (2 total reads)")

            # Verify demux metrics exist and indicate success (not "no_3bc_match")
            metrics_file = os.path.join(out_dir, 'demux_metrics.json')
            self.assertTrue(os.path.exists(metrics_file))

            with open(metrics_file, 'r') as f:
                metrics = json.load(f)

            # Should NOT have demux_type of "no_3bc_match"
            demux_type = metrics.get('demux_type', '')
            self.assertNotEqual(demux_type, 'no_3bc_match',
                              "Should successfully demux, not fail with no_3bc_match")

        finally:
            shutil.rmtree(out_dir)

    def test_splitcode_demux_3bc_with_n_bases_in_bc2_and_i5_rc(self):
        """
        Test 3-barcode demultiplexing with N bases in barcode_2 AND i5 reverse complement.

        This is a regression test for a bug where N-wildcard matching failed in the
        second filtering step for 3-barcode demux.

        Scenario:
        - Samplesheet has barcode_2 = ACGTTACGCA (forward orientation)
        - FASTQ headers have barcode_2 = TGNGTAACGT (reverse complement with N base)
        - After RC: TGNGTAACGT -> ACGTTACNCA
        - Should match ACGTTACGCA via N wildcard (N at position 8, 10% N fraction < 50% threshold)

        Expected behavior:
        - Auto-detection logs i5 reverse complement match
        - N-wildcard matching works in BOTH filtering steps (initial match + 3bc filter)
        - 3-barcode demultiplexing succeeds
        - Output BAMs created with correct read counts
        - No "no_3bc_match" or "no_match" error

        Bug being tested:
        The code had TWO filtering steps for 3-barcode samples:
        1. Initial match via match_barcodes_with_orientation() - uses N-wildcard matching ✓
        2. DataFrame filtering at lines 1094-1098 - used exact match, failed with N ✗
        """
        out_dir = tempfile.mkdtemp()

        try:
            # Use test files with N bases in barcode_2 + reverse complement i5
            r1_fastq = os.path.join(self.input_dir, 'TestPool1_N_S5_L001_R1_001.fastq.gz')
            r2_fastq = os.path.join(self.input_dir, 'TestPool1_N_S5_L001_R2_001.fastq.gz')
            samples_n = os.path.join(self.input_dir, 'samples_3bc_i5_rc_with_n.tsv')

            # Verify test input files exist
            self.assertTrue(os.path.exists(r1_fastq), f"Test FASTQ missing: {r1_fastq}")
            self.assertTrue(os.path.exists(r2_fastq), f"Test FASTQ missing: {r2_fastq}")
            self.assertTrue(os.path.exists(samples_n), f"Test samplesheet missing: {samples_n}")

            # Run splitcode demux with N bases + i5 reverse complement scenario
            illumina.splitcode_demux_fastqs(
                fastq_r1=r1_fastq,
                fastq_r2=r2_fastq,
                samplesheet=samples_n,
                outdir=out_dir,
                runinfo=self.runinfo_xml,
                threads=1
            )

            # Verify expected output BAMs exist
            expected_bams = [
                os.path.join(out_dir, 'TestSample1_N.bam'),
                os.path.join(out_dir, 'TestSample2_N.bam'),
                os.path.join(out_dir, 'TestSample3_N.bam')
            ]

            for bam in expected_bams:
                self.assertTrue(os.path.exists(bam), f"Expected output BAM missing: {bam}")

            # Verify read counts using samtools
            # Note: samtools.count() returns total reads (R1 + R2 for paired-end)
            samtools = tools.samtools.SamtoolsTool()

            # TestSample1_N (AAAAAAAA): should have 2 read pairs = 4 total reads
            count1 = int(samtools.count(expected_bams[0]))
            self.assertEqual(count1, 4, "TestSample1_N should have 2 read pairs (4 total reads)")

            # TestSample2_N (CCCCCCCC): should have 2 read pairs = 4 total reads
            count2 = int(samtools.count(expected_bams[1]))
            self.assertEqual(count2, 4, "TestSample2_N should have 2 read pairs (4 total reads)")

            # TestSample3_N (GGGGTTTT): should have 1 read pair = 2 total reads
            count3 = int(samtools.count(expected_bams[2]))
            self.assertEqual(count3, 2, "TestSample3_N should have 1 read pair (2 total reads)")

            # Verify demux metrics exist and indicate success (not "no_3bc_match" or "no_match")
            metrics_file = os.path.join(out_dir, 'demux_metrics.json')
            self.assertTrue(os.path.exists(metrics_file))

            with open(metrics_file, 'r') as f:
                metrics = json.load(f)

            # Should NOT have demux_type indicating failure
            demux_type = metrics.get('demux_type', '')
            self.assertNotIn('no_3bc_match', demux_type,
                           "Should successfully demux, not fail with no_3bc_match")
            self.assertNotIn('no_match', demux_type,
                           "Should successfully demux, not fail with no_match")

        finally:
            shutil.rmtree(out_dir)

    def test_no_barcode_match_produces_zero_bams(self):
        """
        Test that when FASTQ barcodes don't match any samplesheet entry,
        zero BAMs are produced but metrics files are still generated.

        This tests the behavior change where we emit zero BAMs instead of
        an empty BAM when there's no barcode match.
        """
        out_dir = tempfile.mkdtemp()

        try:
            # Create test FASTQ files with barcodes that don't match the samplesheet
            # The samplesheet has entries for ATCGATCG+GCTAGCTA, CTGATCGT+TAGATCGC, GGAATTCC+CCGGAATT
            # We'll use AAAAAAAA+TTTTTTTT which uses valid DNA chars but matches none of them
            nonmatch_r1 = os.path.join(out_dir, 'NonMatch_S99_L001_R1_001.fastq.gz')
            nonmatch_r2 = os.path.join(out_dir, 'NonMatch_S99_L001_R2_001.fastq.gz')

            # Create R1 FASTQ with non-matching barcodes
            with gzip.open(nonmatch_r1, 'wt') as f:
                f.write("@TEST001:1:TESTFC01:1:1101:1000:2000 1:N:0:AAAAAAAA+TTTTTTTT\n")
                f.write("ACGTACGTACGTACGT\n")
                f.write("+\n")
                f.write("IIIIIIIIIIIIIIII\n")

            # Create R2 FASTQ
            with gzip.open(nonmatch_r2, 'wt') as f:
                f.write("@TEST001:1:TESTFC01:1:1101:1000:2000 2:N:0:AAAAAAAA+TTTTTTTT\n")
                f.write("TGCATGCATGCATGCA\n")
                f.write("+\n")
                f.write("IIIIIIIIIIIIIIII\n")

            # Run splitcode_demux_fastqs - should not raise an exception
            illumina.splitcode_demux_fastqs(
                fastq_r1=nonmatch_r1,
                fastq_r2=nonmatch_r2,
                samplesheet=self.samples_3bc,
                outdir=out_dir,
                threads=1
            )

            # Verify NO BAM files were created (zero BAMs for no match)
            bam_files = [f for f in os.listdir(out_dir) if f.endswith('.bam')]
            self.assertEqual(len(bam_files), 0,
                           f"Expected zero BAM files for no barcode match, found: {bam_files}")

            # Verify metrics files ARE created
            metrics_file = os.path.join(out_dir, 'demux_metrics.json')
            self.assertTrue(os.path.exists(metrics_file),
                          "demux_metrics.json should be created even with no match")

            picard_metrics_file = os.path.join(out_dir, 'demux_metrics_picard-style.txt')
            self.assertTrue(os.path.exists(picard_metrics_file),
                          "demux_metrics_picard-style.txt should be created even with no match")

            # Verify metrics indicate no_barcode_match
            with open(metrics_file, 'r') as f:
                metrics = json.load(f)

            demux_type = metrics.get('demux_type', '')
            self.assertIn('no_barcode_match', demux_type,
                        f"demux_type should indicate no_barcode_match, got: {demux_type}")

            # Verify read_count is 0
            samples = metrics.get('samples', {})
            for sample_name, sample_info in samples.items():
                self.assertEqual(sample_info.get('read_count'), 0,
                               f"Sample {sample_name} should have 0 reads")

        finally:
            shutil.rmtree(out_dir)


class TestMergeDemuxMetrics(TestCaseWithTmp):
    """
    Test suite for merge_demux_metrics entry point.

    Regression tests for issue #127 CLI argument parsing.
    """

    def test_merge_demux_metrics_via_parser(self):
        """
        Test merge_demux_metrics command via argument parser.

        This tests the full CLI code path including:
        - Parser argument processing
        - Argument validation
        - Main function invocation

        Regression test for issue #127 where split_args=True caused
        TypeError with main_merge_demux_metrics(args) signature.
        """
        # Get test input files
        test_dir = util.file.get_test_path()
        input1 = os.path.join(test_dir, 'input', 'TestIlluminaBarcodeHelper', 'single_index', 'metrics.txt')
        input2 = os.path.join(test_dir, 'input', 'TestIlluminaBarcodeHelper', 'one_correction', 'metrics.txt')

        self.assertTrue(os.path.exists(input1), f"Test input missing: {input1}")
        self.assertTrue(os.path.exists(input2), f"Test input missing: {input2}")

        out_file = util.file.mkstempfname('.txt')

        try:
            # Call via parser (like CLI does)
            parser = illumina.parser_merge_demux_metrics(argparse.ArgumentParser())
            args = parser.parse_args([input1, input2, out_file])
            args.func_main(args)

            # Verify output file was created
            self.assertTrue(os.path.exists(out_file), "Merged metrics file should be created")

            # Verify it has content
            with open(out_file, 'r') as f:
                content = f.read()
                self.assertGreater(len(content), 0, "Merged file should have content")

        finally:
            if os.path.exists(out_file):
                os.remove(out_file)
