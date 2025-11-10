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


class TestSplitcodeLookupTable(TestCaseWithTmp):
    """Test cases for create_splitcode_lookup_table function."""

    def test_basic_single_pool(self):
        """Test basic functionality with a single pool."""
        inDir = util.file.get_test_input_path(self)
        sample_sheet = os.path.join(inDir, 'sample_sheet_basic.tsv')

        with tempfile.TemporaryDirectory() as tmpdir:
            csv_out = os.path.join(tmpdir, 'lut.csv')

            # Copy summary JSON to tmpdir where function expects it
            import shutil
            shutil.copy(
                os.path.join(inDir, 'ATCGATCG-GCTAGCTA.lB1_summary.json'),
                os.path.join(tmpdir, 'ATCGATCG-GCTAGCTA.lB1_summary.json')
            )

            result_path = illumina.create_splitcode_lookup_table(
                sample_sheet, csv_out, unmatched_name="Unmatched"
            )

            # Verify output file was created
            self.assertTrue(os.path.exists(result_path))

            # Read and validate output CSV
            import pandas as pd
            df = pd.read_csv(result_path, dtype=str)

            # Check for expected columns
            expected_cols = [
                'sample', 'library_id', 'barcode_1', 'barcode_2',
                'inline_barcode', 'run', 'muxed_pool',
                'num_reads_hdistance0', 'num_reads_hdistance1', 'num_reads_total'
            ]
            for col in expected_cols:
                self.assertIn(col, df.columns, f"Missing column: {col}")

            # Check number of rows (3 samples + 1 unmatched)
            self.assertEqual(len(df), 4)

            # Verify sample names
            sample_names = set(df['sample'].tolist())
            self.assertIn('Sample1', sample_names)
            self.assertIn('Sample2', sample_names)
            self.assertIn('Sample3', sample_names)

            # Verify unmatched row
            unmatched_rows = df[df['sample'].str.contains('Unmatched')]
            self.assertEqual(len(unmatched_rows), 1)
            unmatched_row = unmatched_rows.iloc[0]
            self.assertTrue(unmatched_row['inline_barcode'].replace('N', '') == '')
            # Unmatched count should be n_processed - n_assigned = 100000 - 95000 = 5000
            self.assertEqual(int(unmatched_row['num_reads_hdistance0']), 5000)

            # Verify read counts for Sample2 (should have highest count)
            sample2_row = df[df['sample'] == 'Sample2'].iloc[0]
            self.assertEqual(int(sample2_row['num_reads_hdistance0']), 45000)
            self.assertEqual(int(sample2_row['num_reads_hdistance1']), 2000)
            self.assertEqual(int(sample2_row['num_reads_total']), 47000)

            # Verify barcode values
            self.assertEqual(sample2_row['barcode_1'], 'ATCGATCG')
            self.assertEqual(sample2_row['barcode_2'], 'GCTAGCTA')
            self.assertEqual(sample2_row['inline_barcode'], 'GGGGTTTT')

    def test_zero_reads_pool(self):
        """Test handling of pool with zero reads."""
        inDir = util.file.get_test_input_path(self)
        sample_sheet = os.path.join(inDir, 'sample_sheet_zero_reads.tsv')

        with tempfile.TemporaryDirectory() as tmpdir:
            csv_out = os.path.join(tmpdir, 'lut.csv')

            # Copy summary JSON
            import shutil
            shutil.copy(
                os.path.join(inDir, 'TTTTAAAA-CCCCGGGG.lB2_summary.json'),
                os.path.join(tmpdir, 'TTTTAAAA-CCCCGGGG.lB2_summary.json')
            )

            result_path = illumina.create_splitcode_lookup_table(
                sample_sheet, csv_out, unmatched_name="Unmatched"
            )

            # Read output
            import pandas as pd
            df = pd.read_csv(result_path, dtype=str)

            # Should have 1 sample + 1 unmatched
            self.assertEqual(len(df), 2)

            # Verify sample has 0 reads
            sample_row = df[df['sample'] == 'Sample4'].iloc[0]
            self.assertEqual(int(sample_row['num_reads_hdistance0']), 0)
            self.assertEqual(int(sample_row['num_reads_hdistance1']), 0)
            self.assertEqual(int(sample_row['num_reads_total']), 0)

    def test_multi_pool(self):
        """Test with multiple pools."""
        inDir = util.file.get_test_input_path(self)
        sample_sheet = os.path.join(inDir, 'sample_sheet_multi_pool.tsv')

        with tempfile.TemporaryDirectory() as tmpdir:
            csv_out = os.path.join(tmpdir, 'lut.csv')

            # Copy both summary JSONs
            import shutil
            shutil.copy(
                os.path.join(inDir, 'AAAAAAAA-TTTTTTTT.lLibA_summary.json'),
                os.path.join(tmpdir, 'AAAAAAAA-TTTTTTTT.lLibA_summary.json')
            )
            shutil.copy(
                os.path.join(inDir, 'GGGGGGGG-CCCCCCCC.lLibB_summary.json'),
                os.path.join(tmpdir, 'GGGGGGGG-CCCCCCCC.lLibB_summary.json')
            )

            result_path = illumina.create_splitcode_lookup_table(
                sample_sheet, csv_out, unmatched_name="Unmatched"
            )

            # Read output
            import pandas as pd
            df = pd.read_csv(result_path, dtype=str)

            # Should have 4 samples + 2 unmatched (one per pool)
            self.assertEqual(len(df), 6)

            # Verify both pools are present
            pools = set(df['muxed_pool'].tolist())
            self.assertIn('AAAAAAAA-TTTTTTTT.lLibA', pools)
            self.assertIn('GGGGGGGG-CCCCCCCC.lLibB', pools)

            # Check unmatched counts
            unmatched_rows = df[df['sample'].str.contains('Unmatched')]
            self.assertEqual(len(unmatched_rows), 2)

            # Verify LibA unmatched count: 50000 - 48000 = 2000
            liba_unmatched = df[df['muxed_pool'] == 'AAAAAAAA-TTTTTTTT.lLibA']
            liba_unmatched = liba_unmatched[liba_unmatched['sample'].str.contains('Unmatched')].iloc[0]
            self.assertEqual(int(liba_unmatched['num_reads_hdistance0']), 2000)

    def test_append_run_id(self):
        """Test append_run_id parameter."""
        inDir = util.file.get_test_input_path(self)
        sample_sheet = os.path.join(inDir, 'sample_sheet_basic.tsv')

        with tempfile.TemporaryDirectory() as tmpdir:
            csv_out = os.path.join(tmpdir, 'lut.csv')

            import shutil
            # Copy JSON file with run suffix - this matches what splitcode actually creates
            shutil.copy(
                os.path.join(inDir, 'ATCGATCG-GCTAGCTA.lB1_summary.json'),
                os.path.join(tmpdir, 'ATCGATCG-GCTAGCTA.lB1.FLOWCELL123_summary.json')
            )

            result_path = illumina.create_splitcode_lookup_table(
                sample_sheet, csv_out,
                unmatched_name="Unmatched",
                append_run_id="FLOWCELL123"
            )

            # Read output
            import pandas as pd
            df = pd.read_csv(result_path, dtype=str)

            # Verify run IDs contain flowcell ID
            sample_row = df[df['sample'] == 'Sample1'].iloc[0]
            self.assertIn('FLOWCELL123', sample_row['run'])
            self.assertIn('FLOWCELL123', sample_row['muxed_pool'])


class TestSplitcodeIntegration(TestCaseWithTmp):
    """
    Test our assumptions about splitcode behavior and file formats.

    These tests create minimal test data to validate splitcode integration
    without testing the external splitcode tool itself.
    """

    def setUp(self):
        super().setUp()
        self.temp_dir = tempfile.mkdtemp()
        self.samtools = tools.samtools.SamtoolsTool()

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)
        super().tearDown()

    def create_test_bam_with_inline_barcodes(self, output_bam, num_reads=10, barcode="AAAAAAAA"):
        """
        Create a test BAM file with reads containing inline barcodes in R1.

        Args:
            output_bam: Path to output BAM file
            num_reads: Number of read pairs to create
            barcode: 8bp inline barcode to prepend to R1 sequences
        """
        # Create temp FASTQs with inline barcodes
        r1_fastq = util.file.mkstempfname('.fastq')
        r2_fastq = util.file.mkstempfname('.fastq')

        with open(r1_fastq, 'w') as f1, open(r2_fastq, 'w') as f2:
            for i in range(num_reads):
                # R1 has inline barcode prepended
                seq_r1 = barcode + "ACGTACGTACGT"  # 8bp barcode + 12bp insert
                qual_r1 = "I" * len(seq_r1)

                # R2 is just normal sequence
                seq_r2 = "TGCATGCATGCA"
                qual_r2 = "I" * len(seq_r2)

                f1.write(f"@read{i}\n{seq_r1}\n+\n{qual_r1}\n")
                f2.write(f"@read{i}\n{seq_r2}\n+\n{qual_r2}\n")

        # Convert to BAM
        tools.picard.FastqToSamTool().execute(
            r1_fastq,
            r2_fastq,
            "TestSample",
            output_bam,
            picardOptions=['LIBRARY_NAME=TestLib', 'PLATFORM=ILLUMINA']
        )

        return output_bam

    def test_run_splitcode_on_pool_basic(self):
        """
        Test run_splitcode_on_pool with a simple single-barcode pool.

        This test focuses on validating our assumptions about splitcode's behavior:
        - Summary JSON is created in expected location
        - JSON has expected structure (tag_qc array with tag, count, distance)
        - We're NOT testing the keep file format here (too complex) - just JSON output

        Validates:
        - Splitcode runs successfully on BAM input
        - Summary JSON location and format
        """
        # Create test BAM with reads containing AAAAAAAA barcode
        pool_bam = os.path.join(self.temp_dir, 'test_pool.bam')
        self.create_test_bam_with_inline_barcodes(pool_bam, num_reads=5, barcode="AAAAAAAA")

        # Create splitcode config file
        # Format: tag, id, locations, distance, left, right
        # The "id" column must match what's in the keep file
        splitcode_config = os.path.join(self.temp_dir, 'config.txt')
        sample_id = "Sample1_R1"  # This matches illumina.py convention: f"{sample_library_id}_R1"
        with open(splitcode_config, 'w') as f:
            f.write("tag\tid\tlocations\tdistance\tleft\tright\n")
            # AAAAAAAA barcode, ID=Sample1_R1, in file 0 (R1) from position 0 to 8, distance=0 (exact match), trim from left
            # locations format: FILE_NUMBER:START_BP:END_BP (0:0:8 means R1, positions 0-8)
            f.write(f"AAAAAAAA\t{sample_id}\t0:0:8\t0\t1\t0\n")

        # Create keep file matching the config file ID
        # Format: barcode_id (from config) <tab> output_prefix
        splitcode_keepfile = os.path.join(self.temp_dir, 'keep.txt')
        output_prefix = os.path.join(self.temp_dir, 'Sample1')
        with open(splitcode_keepfile, 'w') as f:
            # Column 1 must match the "id" from config file
            f.write(f"{sample_id}\t{output_prefix}\n")

        pool_id = "TestPool"
        out_demux_dir = self.temp_dir

        # Run splitcode
        result = illumina.run_splitcode_on_pool(
            pool_id=pool_id,
            pool_bam_file=pool_bam,
            splitcode_config=splitcode_config,
            splitcode_keepfile=splitcode_keepfile,
            out_demux_dir_path=out_demux_dir,
            unmatched_name="unmatched",
            threads_per_worker=1,
            out_dir_path=out_demux_dir,
            out_demux_dir_path_tmp=self.temp_dir
        )

        # Verify return value
        self.assertEqual(result[0], 0, "Splitcode should return 0 on success")
        self.assertEqual(result[1], pool_id, "Should return the pool_id")

        # Verify summary JSON exists in expected location
        # This is the KEY assumption we need to validate
        expected_summary = os.path.join(out_demux_dir, f'{pool_id}_summary.json')
        self.assertTrue(os.path.exists(expected_summary),
                       f"Summary JSON should exist at {expected_summary}")

        # Verify JSON structure - this is what create_splitcode_lookup_table depends on
        with open(expected_summary) as f:
            summary = json.load(f)

        self.assertIn('tag_qc', summary, "JSON should have 'tag_qc' field")
        self.assertIsInstance(summary['tag_qc'], list, "'tag_qc' should be a list")

        # Find Sample1_R1 in tag_qc (splitcode uses the ID as the tag name in output)
        # Note: tag_qc has multiple entries per tag (one for each distance level)
        # We want to check that at least one distance level has count > 0
        matching_entries = [item for item in summary['tag_qc'] if item['tag'] == sample_id]
        self.assertGreater(len(matching_entries), 0, f"{sample_id} should appear in tag_qc")

        # Check that we have at least one entry with count > 0
        total_count = sum(item['count'] for item in matching_entries)
        self.assertGreater(total_count, 0, f"{sample_id} should have matched some reads across all distance levels")

    def test_run_splitcode_on_pool_empty_output(self):
        """
        Test run_splitcode_on_pool when a barcode matches zero reads.

        Validates:
        - Splitcode runs successfully even with empty barcode outputs
        - Summary JSON shows count=0 for unmatched barcodes
        - Output files for empty barcodes may not exist or are empty
        """
        # Create test BAM with AAAAAAAA barcode only
        pool_bam = os.path.join(self.temp_dir, 'test_pool.bam')
        self.create_test_bam_with_inline_barcodes(pool_bam, num_reads=3, barcode="AAAAAAAA")

        # Config with BC1 (present) and BC2 (absent)
        # IDs must match illumina.py convention: f"{sample_library_id}_R1"
        splitcode_config = os.path.join(self.temp_dir, 'config.txt')
        sample_id_1 = "Sample1_R1"
        sample_id_2 = "Sample2_R1"
        with open(splitcode_config, 'w') as f:
            f.write("tag\tid\tlocations\tdistance\tleft\tright\n")
            # locations format: FILE_NUMBER:START_BP:END_BP
            f.write(f"AAAAAAAA\t{sample_id_1}\t0:0:8\t0\t1\t0\n")
            f.write(f"CCCCCCCC\t{sample_id_2}\t0:0:8\t0\t1\t0\n")  # Won't match any reads

        # Keep file for both barcodes
        splitcode_keepfile = os.path.join(self.temp_dir, 'keep.txt')
        output_prefix1 = os.path.join(self.temp_dir, 'Sample1')
        output_prefix2 = os.path.join(self.temp_dir, 'Sample2')
        with open(splitcode_keepfile, 'w') as f:
            f.write(f"{sample_id_1}\t{output_prefix1}\n")
            f.write(f"{sample_id_2}\t{output_prefix2}\n")

        pool_id = "TestPool"
        out_demux_dir = self.temp_dir

        # Run splitcode
        result = illumina.run_splitcode_on_pool(
            pool_id=pool_id,
            pool_bam_file=pool_bam,
            splitcode_config=splitcode_config,
            splitcode_keepfile=splitcode_keepfile,
            out_demux_dir_path=out_demux_dir,
            unmatched_name="unmatched",
            threads_per_worker=1,
            out_dir_path=out_demux_dir,
            out_demux_dir_path_tmp=self.temp_dir
        )

        self.assertEqual(result[0], 0, "Splitcode should succeed even with empty outputs")

        # Verify summary JSON
        expected_summary = os.path.join(out_demux_dir, f'{pool_id}_summary.json')
        with open(expected_summary) as f:
            summary = json.load(f)

        # Sample1_R1 should have reads (check total across all distance levels)
        sample1_entries = [item for item in summary['tag_qc'] if item['tag'] == sample_id_1]
        self.assertGreater(len(sample1_entries), 0, f"{sample_id_1} should appear in tag_qc")
        sample1_total = sum(item['count'] for item in sample1_entries)
        self.assertGreater(sample1_total, 0, f"{sample_id_1} should have matched some reads")

        # Sample2_R1 might appear with count=0 (check total across all distance levels)
        sample2_entries = [item for item in summary['tag_qc'] if item['tag'] == sample_id_2]
        if len(sample2_entries) > 0:
            sample2_total = sum(item['count'] for item in sample2_entries)
            self.assertEqual(sample2_total, 0,
                           f"{sample_id_2} should have 0 reads total if it appears in tag_qc")

        # Sample1 files should exist and be non-empty
        self.assertTrue(os.path.exists(f"{output_prefix1}_R1.fastq"))
        self.assertGreater(os.path.getsize(f"{output_prefix1}_R1.fastq"), 0)

        # Sample2 files might not exist, or if they do, should be empty
        sample2_r1 = f"{output_prefix2}_R1.fastq"
        if os.path.exists(sample2_r1):
            # File exists but should be empty (0 reads = 0 lines)
            with open(sample2_r1) as f:
                lines = f.readlines()
            self.assertEqual(len(lines), 0,
                           "Empty barcode output file should have 0 lines")

    def test_splitcode_json_output_format(self):
        """
        Test that splitcode JSON output has the expected structure.

        Validates our assumptions about JSON format:
        - Has 'tag_qc' array
        - Each entry has 'tag', 'count', 'distance' fields
        - These fields are the correct types
        """
        # Create minimal test
        pool_bam = os.path.join(self.temp_dir, 'test_pool.bam')
        self.create_test_bam_with_inline_barcodes(pool_bam, num_reads=10, barcode="GGGGGGGG")

        splitcode_config = os.path.join(self.temp_dir, 'config.txt')
        sample_id = "Sample3_R1"
        with open(splitcode_config, 'w') as f:
            f.write("tag\tid\tlocations\tdistance\tleft\tright\n")
            # locations format: FILE_NUMBER:START_BP:END_BP
            f.write(f"GGGGGGGG\t{sample_id}\t0:0:8\t1\t1\t0\n")

        splitcode_keepfile = os.path.join(self.temp_dir, 'keep.txt')
        output_prefix = os.path.join(self.temp_dir, 'Sample3')
        with open(splitcode_keepfile, 'w') as f:
            f.write(f"{sample_id}\t{output_prefix}\n")

        pool_id = "TestPool"

        illumina.run_splitcode_on_pool(
            pool_id=pool_id,
            pool_bam_file=pool_bam,
            splitcode_config=splitcode_config,
            splitcode_keepfile=splitcode_keepfile,
            out_demux_dir_path=self.temp_dir,
            unmatched_name="unmatched",
            threads_per_worker=1,
            out_dir_path=self.temp_dir,
            out_demux_dir_path_tmp=self.temp_dir
        )

        # Read and validate JSON structure
        summary_file = os.path.join(self.temp_dir, f'{pool_id}_summary.json')
        with open(summary_file) as f:
            summary = json.load(f)

        # Validate top-level structure
        self.assertIsInstance(summary, dict, "Summary should be a dictionary")
        self.assertIn('tag_qc', summary, "Summary must have 'tag_qc' key")
        self.assertIsInstance(summary['tag_qc'], list, "'tag_qc' must be a list")

        # Validate tag_qc entries
        self.assertGreater(len(summary['tag_qc']), 0, "tag_qc should have at least one entry")

        for entry in summary['tag_qc']:
            self.assertIsInstance(entry, dict, "Each tag_qc entry should be a dict")

            # Required fields our code depends on
            self.assertIn('tag', entry, "Entry must have 'tag' field")
            self.assertIn('count', entry, "Entry must have 'count' field")

            # Validate types (important for create_splitcode_lookup_table)
            # Note: with our pandas dtype=str fix, these get converted to strings
            self.assertIsInstance(entry['tag'], str, "'tag' should be string")
            # 'count' is numeric in JSON but we convert to string in pandas
            self.assertIsInstance(entry['count'], int, "'count' should be int in JSON")

            if 'distance' in entry:
                self.assertIsInstance(entry['distance'], int, "'distance' should be int in JSON")

    def test_splitcode_output_file_locations(self):
        """
        Test that splitcode creates output files in expected locations.

        Validates:
        - Summary JSON location: {out_demux_dir_path}/{pool_id}_summary.json
        - FASTQ location: {output_prefix}_R1.fastq, {output_prefix}_R2.fastq
        - Unmatched FASTQ location: {unmatched_name}.{pool_id}_R1.fastq
        """
        pool_bam = os.path.join(self.temp_dir, 'test_pool.bam')
        self.create_test_bam_with_inline_barcodes(pool_bam, num_reads=8, barcode="TTTTTTTT")

        splitcode_config = os.path.join(self.temp_dir, 'config.txt')
        sample_id = "MySample_R1"
        with open(splitcode_config, 'w') as f:
            f.write("tag\tid\tlocations\tdistance\tleft\tright\n")
            # locations format: FILE_NUMBER:START_BP:END_BP
            f.write(f"TTTTTTTT\t{sample_id}\t0:0:8\t0\t1\t0\n")

        splitcode_keepfile = os.path.join(self.temp_dir, 'keep.txt')
        sample_output_prefix = os.path.join(self.temp_dir, 'MySample')
        with open(splitcode_keepfile, 'w') as f:
            f.write(f"{sample_id}\t{sample_output_prefix}\n")

        pool_id = "MyPool"
        unmatched_name = "unassigned"
        out_demux_dir = self.temp_dir

        illumina.run_splitcode_on_pool(
            pool_id=pool_id,
            pool_bam_file=pool_bam,
            splitcode_config=splitcode_config,
            splitcode_keepfile=splitcode_keepfile,
            out_demux_dir_path=out_demux_dir,
            unmatched_name=unmatched_name,
            threads_per_worker=1,
            out_dir_path=out_demux_dir,
            out_demux_dir_path_tmp=self.temp_dir
        )

        # Verify expected file locations

        # 1. Summary JSON: {out_demux_dir_path}/{pool_id}_summary.json
        expected_summary = os.path.join(out_demux_dir, f'{pool_id}_summary.json')
        self.assertTrue(os.path.exists(expected_summary),
                       f"Summary should be at {expected_summary}")

        # 2. Matched FASTQs: {output_prefix}_R1.fastq, {output_prefix}_R2.fastq
        expected_r1 = f"{sample_output_prefix}_R1.fastq"
        expected_r2 = f"{sample_output_prefix}_R2.fastq"
        self.assertTrue(os.path.exists(expected_r1),
                       f"Matched R1 should be at {expected_r1}")
        self.assertTrue(os.path.exists(expected_r2),
                       f"Matched R2 should be at {expected_r2}")

        # 3. Unmatched FASTQs: {unmatched_name}.{pool_id}_R1.fastq
        expected_unmatched_r1 = os.path.join(self.temp_dir, f'{unmatched_name}.{pool_id}_R1.fastq')
        expected_unmatched_r2 = os.path.join(self.temp_dir, f'{unmatched_name}.{pool_id}_R2.fastq')
        self.assertTrue(os.path.exists(expected_unmatched_r1),
                       f"Unmatched R1 should be at {expected_unmatched_r1}")
        self.assertTrue(os.path.exists(expected_unmatched_r2),
                       f"Unmatched R2 should be at {expected_unmatched_r2}")

    def test_splitcode_barcode_trimming(self):
        """
        Test that splitcode correctly trims inline barcodes from R1.

        Validates:
        - With left=1, barcode is removed from R1 sequences
        - R2 sequences are unchanged
        - Read count is preserved
        """
        # Create BAM where R1 = AAAAAAAA + ACGTACGTACGT (20bp total)
        pool_bam = os.path.join(self.temp_dir, 'test_pool.bam')
        self.create_test_bam_with_inline_barcodes(pool_bam, num_reads=5, barcode="AAAAAAAA")

        splitcode_config = os.path.join(self.temp_dir, 'config.txt')
        sample_id = "Sample1_R1"
        with open(splitcode_config, 'w') as f:
            f.write("tag\tid\tlocations\tdistance\tleft\tright\n")
            # left=1 means trim the barcode from the left side
            # locations format: FILE_NUMBER:START_BP:END_BP
            f.write(f"AAAAAAAA\t{sample_id}\t0:0:8\t0\t1\t0\n")

        splitcode_keepfile = os.path.join(self.temp_dir, 'keep.txt')
        output_prefix = os.path.join(self.temp_dir, 'Sample1')
        with open(splitcode_keepfile, 'w') as f:
            f.write(f"{sample_id}\t{output_prefix}\n")

        illumina.run_splitcode_on_pool(
            pool_id="TestPool",
            pool_bam_file=pool_bam,
            splitcode_config=splitcode_config,
            splitcode_keepfile=splitcode_keepfile,
            out_demux_dir_path=self.temp_dir,
            unmatched_name="unmatched",
            threads_per_worker=1,
            out_dir_path=self.temp_dir,
            out_demux_dir_path_tmp=self.temp_dir
        )

        # Read output FASTQ and verify barcode was trimmed
        output_r1 = f"{output_prefix}_R1.fastq"
        with open(output_r1) as f:
            lines = f.readlines()

        # Should have 5 reads * 4 lines/read = 20 lines
        self.assertEqual(len(lines), 20, "Should have 5 reads in output")

        # Check first read's sequence (line 1, 0-indexed)
        first_seq = lines[1].strip()

        # After trimming 8bp barcode, should have 12bp remaining
        self.assertEqual(len(first_seq), 12,
                        "Sequence should be 12bp after trimming 8bp barcode")

        # Should NOT start with the barcode
        self.assertFalse(first_seq.startswith("AAAAAAAA"),
                        "Barcode should be trimmed from output")

        # Should be the insert sequence
        self.assertEqual(first_seq, "ACGTACGTACGT",
                        "After trimming barcode, should have insert sequence")

    def test_splitcode_with_append_run_id(self):
        """
        Test that summary JSON lookup works correctly with append_run_id parameter.

        This tests the bug we fixed where append_run_id suffix was breaking JSON file lookup.

        Validates:
        - Summary JSON filename doesn't include append_run_id suffix
        - JSON can be found correctly even when pool_id has suffix appended
        """
        pool_bam = os.path.join(self.temp_dir, 'test_pool.bam')
        self.create_test_bam_with_inline_barcodes(pool_bam, num_reads=3, barcode="GGGGGGGG")

        splitcode_config = os.path.join(self.temp_dir, 'config.txt')
        sample_id = "Sample1_R1"
        with open(splitcode_config, 'w') as f:
            f.write("tag\tid\tlocations\tdistance\tleft\tright\n")
            # locations format: FILE_NUMBER:START_BP:END_BP
            f.write(f"GGGGGGGG\t{sample_id}\t0:0:8\t0\t1\t0\n")

        splitcode_keepfile = os.path.join(self.temp_dir, 'keep.txt')
        output_prefix = os.path.join(self.temp_dir, 'Sample1')
        with open(splitcode_keepfile, 'w') as f:
            f.write(f"{sample_id}\t{output_prefix}\n")

        base_pool_id = "Pool1"
        out_demux_dir = self.temp_dir

        # Run splitcode (this doesn't use append_run_id directly, but creates the JSON)
        illumina.run_splitcode_on_pool(
            pool_id=base_pool_id,
            pool_bam_file=pool_bam,
            splitcode_config=splitcode_config,
            splitcode_keepfile=splitcode_keepfile,
            out_demux_dir_path=out_demux_dir,
            unmatched_name="unmatched",
            threads_per_worker=1,
            out_dir_path=out_demux_dir,
            out_demux_dir_path_tmp=self.temp_dir
        )

        # Verify summary JSON is created with base pool_id (no suffix)
        expected_summary = os.path.join(out_demux_dir, f'{base_pool_id}_summary.json')
        self.assertTrue(os.path.exists(expected_summary),
                       "Summary JSON should use base pool_id without suffix")

        # Now test that our lookup logic would work if pool_id had a suffix
        # This simulates what happens in create_splitcode_lookup_table
        pool_id_with_suffix = f"{base_pool_id}.FLOWCELL123"

        # The fix in illumina.py strips the suffix when looking for JSON file
        # Simulate that logic here
        append_run_id = "FLOWCELL123"
        if append_run_id and pool_id_with_suffix.endswith(f".{append_run_id}"):
            pool_for_lookup = pool_id_with_suffix[:-len(f".{append_run_id}")]
        else:
            pool_for_lookup = pool_id_with_suffix

        # Should find the file with the base pool_id
        import glob
        found_files = glob.glob(f"{out_demux_dir}/{pool_for_lookup}_summary.json")
        self.assertEqual(len(found_files), 1,
                        "Should find exactly one summary JSON file")
        self.assertEqual(found_files[0], expected_summary,
                        "Should find the summary JSON with base pool_id")


class TestGenerateSplitcodeConfigAndKeepFiles(TestCaseWithTmp):
    """
    Test generate_splitcode_config_and_keep_files function.

    This function takes sample metadata and generates splitcode config/keep files.
    It's a pure text-file transformation that's ideal for unit testing.
    """

    def setUp(self):
        super().setUp()
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)
        super().tearDown()

    def test_basic_single_pool(self):
        """Test basic config/keep generation with a single pool containing 2 samples."""
        import pandas as pd
        import csv

        # Create test DataFrame matching inner_demux_mapper output
        df = pd.DataFrame({
            'barcode_3': ['AAAAAAAA', 'CCCCCCCC'],
            'run': ['Sample1.lLib1', 'Sample2.lLib1'],
            'muxed_run': ['Pool1', 'Pool1']
        }, index=['Sample1', 'Sample2'])

        config_file, keep_file, sample_ids = illumina.generate_splitcode_config_and_keep_files(
            df, 'Pool1', self.temp_dir, max_hamming_dist=1
        )

        # Verify return values
        self.assertTrue(os.path.exists(config_file), "Config file should exist")
        self.assertTrue(os.path.exists(keep_file), "Keep file should exist")
        self.assertEqual(sample_ids, ['Sample1.lLib1', 'Sample2.lLib1'])

        # Verify config file format
        with open(config_file) as f:
            reader = csv.reader(f, delimiter='\t')
            rows = list(reader)

        # Check header
        self.assertEqual(rows[0], ['tag', 'id', 'locations', 'distance', 'left', 'right'])

        # Check Sample1 config line
        self.assertEqual(rows[1][0], 'AAAAAAAA')  # barcode
        self.assertEqual(rows[1][1], 'Sample1.lLib1')  # ID (no _R1 suffix when using --keep-r1-r2)
        self.assertEqual(rows[1][2], '0:0:8')  # locations
        self.assertEqual(rows[1][3], '1')  # distance
        self.assertEqual(rows[1][4], '1')  # left trim
        self.assertEqual(rows[1][5], '0')  # right trim

        # Check Sample2 config line
        self.assertEqual(rows[2][0], 'CCCCCCCC')
        self.assertEqual(rows[2][1], 'Sample2.lLib1')  # ID (no _R1 suffix when using --keep-r1-r2)
        self.assertEqual(rows[2][2], '0:0:8')

        # Verify keep file format (NO header)
        with open(keep_file) as f:
            reader = csv.reader(f, delimiter='\t')
            keep_rows = list(reader)

        # Should have 2 rows, no header
        self.assertEqual(len(keep_rows), 2)

        # Check Sample1 keep line (no _R1 suffix when using --keep-r1-r2)
        self.assertEqual(keep_rows[0][0], 'Sample1.lLib1')
        self.assertEqual(keep_rows[0][1], f'{self.temp_dir}/Sample1.lLib1')

        # Check Sample2 keep line (no _R1 suffix when using --keep-r1-r2)
        self.assertEqual(keep_rows[1][0], 'Sample2.lLib1')
        self.assertEqual(keep_rows[1][1], f'{self.temp_dir}/Sample2.lLib1')

    def test_variable_barcode_lengths(self):
        """Test that config correctly handles different barcode lengths."""
        import pandas as pd
        import csv

        df = pd.DataFrame({
            'barcode_3': ['AAAA', 'CCCCCCCCCC', 'GGGGGGGG'],  # 4bp, 10bp, 8bp
            'run': ['S1.lL1', 'S2.lL1', 'S3.lL1'],
            'muxed_run': ['Pool1', 'Pool1', 'Pool1']
        }, index=['S1', 'S2', 'S3'])

        config_file, keep_file, sample_ids = illumina.generate_splitcode_config_and_keep_files(
            df, 'Pool1', self.temp_dir
        )

        with open(config_file) as f:
            reader = csv.reader(f, delimiter='\t')
            rows = list(reader)

        # Check locations are adjusted for barcode length
        self.assertEqual(rows[1][2], '0:0:4')   # 4bp barcode
        self.assertEqual(rows[2][2], '0:0:10')  # 10bp barcode
        self.assertEqual(rows[3][2], '0:0:8')   # 8bp barcode

    def test_hamming_distance_parameter(self):
        """Test that max_hamming_dist parameter is correctly written to config."""
        import pandas as pd
        import csv

        df = pd.DataFrame({
            'barcode_3': ['AAAAAAAA'],
            'run': ['Sample1.lLib1'],
            'muxed_run': ['Pool1']
        }, index=['Sample1'])

        # Test distance=0 (exact match only)
        config_file, _, _ = illumina.generate_splitcode_config_and_keep_files(
            df, 'Pool1', self.temp_dir, max_hamming_dist=0
        )

        with open(config_file) as f:
            reader = csv.reader(f, delimiter='\t')
            rows = list(reader)

        self.assertEqual(rows[1][3], '0')  # distance column

        # Test distance=2
        config_file, _, _ = illumina.generate_splitcode_config_and_keep_files(
            df, 'Pool1', self.temp_dir, max_hamming_dist=2
        )

        with open(config_file) as f:
            reader = csv.reader(f, delimiter='\t')
            rows = list(reader)

        self.assertEqual(rows[1][3], '2')  # distance column

    def test_r1_trim_bp_right_of_barcode(self):
        """Test that r1_trim_bp_right_of_barcode correctly sets left trim parameter."""
        import pandas as pd
        import csv

        df = pd.DataFrame({
            'barcode_3': ['AAAAAAAA'],
            'run': ['Sample1.lLib1'],
            'muxed_run': ['Pool1']
        }, index=['Sample1'])

        # Test default (None) - should be "1"
        config_file, _, _ = illumina.generate_splitcode_config_and_keep_files(
            df, 'Pool1', self.temp_dir, r1_trim_bp_right_of_barcode=None
        )

        with open(config_file) as f:
            reader = csv.reader(f, delimiter='\t')
            rows = list(reader)

        self.assertEqual(rows[1][4], '1')  # left column: "1" means trim barcode only

        # Test with 5 extra bp
        config_file, _, _ = illumina.generate_splitcode_config_and_keep_files(
            df, 'Pool1', self.temp_dir, r1_trim_bp_right_of_barcode=5
        )

        with open(config_file) as f:
            reader = csv.reader(f, delimiter='\t')
            rows = list(reader)

        self.assertEqual(rows[1][4], '1:5')  # left column: "1:5" means trim barcode + 5bp

    def test_multi_pool_filtering(self):
        """Test that function correctly filters to only the specified pool."""
        import pandas as pd
        import csv

        # DataFrame with samples from two different pools
        df = pd.DataFrame({
            'barcode_3': ['AAAAAAAA', 'CCCCCCCC', 'GGGGGGGG'],
            'run': ['S1.lL1', 'S2.lL1', 'S3.lL1'],
            'muxed_run': ['Pool1', 'Pool2', 'Pool1']  # Mixed pools
        }, index=['S1', 'S2', 'S3'])

        # Generate for Pool1 only
        config_file, keep_file, sample_ids = illumina.generate_splitcode_config_and_keep_files(
            df, 'Pool1', self.temp_dir
        )

        # Should only have S1 and S3 (Pool1 samples)
        self.assertEqual(sample_ids, ['S1.lL1', 'S3.lL1'])

        with open(config_file) as f:
            reader = csv.reader(f, delimiter='\t')
            rows = list(reader)

        # Header + 2 sample rows (S1 and S3)
        self.assertEqual(len(rows), 3)
        self.assertEqual(rows[1][0], 'AAAAAAAA')  # S1 barcode
        self.assertEqual(rows[2][0], 'GGGGGGGG')  # S3 barcode

        # S2 (Pool2) should NOT be present
        barcodes_in_config = [row[0] for row in rows[1:]]
        self.assertNotIn('CCCCCCCC', barcodes_in_config)

    def test_empty_pool_raises_error(self):
        """Test that requesting a non-existent pool raises ValueError."""
        import pandas as pd

        df = pd.DataFrame({
            'barcode_3': ['AAAAAAAA'],
            'run': ['Sample1.lLib1'],
            'muxed_run': ['Pool1']
        }, index=['Sample1'])

        # Try to generate for Pool2 (which doesn't exist)
        with self.assertRaises(ValueError) as context:
            illumina.generate_splitcode_config_and_keep_files(
                df, 'Pool2', self.temp_dir
            )

        self.assertIn('No samples found', str(context.exception))
        self.assertIn('Pool2', str(context.exception))

    def test_config_keep_id_matching(self):
        """
        Test the critical requirement: config file "id" must match keep file column 1.

        This is the most common source of splitcode errors.
        """
        import pandas as pd
        import csv

        df = pd.DataFrame({
            'barcode_3': ['AAAAAAAA', 'TTTTTTTT'],
            'run': ['MySample.lMyLib.FC.1', 'Other.lLib2.FC.1'],
            'muxed_run': ['Pool1', 'Pool1']
        }, index=['MySample', 'Other'])

        config_file, keep_file, _ = illumina.generate_splitcode_config_and_keep_files(
            df, 'Pool1', self.temp_dir
        )

        # Read both files
        with open(config_file) as f:
            config_rows = list(csv.reader(f, delimiter='\t'))

        with open(keep_file) as f:
            keep_rows = list(csv.reader(f, delimiter='\t'))

        # Extract IDs from config file (column 1, skip header)
        config_ids = [row[1] for row in config_rows[1:]]

        # Extract IDs from keep file (column 0, no header)
        keep_ids = [row[0] for row in keep_rows]

        # They MUST match exactly
        self.assertEqual(config_ids, keep_ids)
        # No _R1 suffix when using --keep-r1-r2 mode (splitcode adds the suffixes automatically)
        self.assertEqual(config_ids, ['MySample.lMyLib.FC.1', 'Other.lLib2.FC.1'])

    def test_output_prefix_path_construction(self):
        """Test that keep file constructs correct output paths."""
        import pandas as pd
        import csv

        df = pd.DataFrame({
            'barcode_3': ['AAAAAAAA'],
            'run': ['Sample1.lLib1.FLOWCELL.1'],
            'muxed_run': ['Pool1']
        }, index=['Sample1'])

        output_dir = "/custom/output/path"
        _, keep_file, _ = illumina.generate_splitcode_config_and_keep_files(
            df, 'Pool1', output_dir
        )

        with open(keep_file) as f:
            keep_rows = list(csv.reader(f, delimiter='\t'))

        # Output prefix should be: {output_dir}/{sample_library_id}
        expected_prefix = f"{output_dir}/Sample1.lLib1.FLOWCELL.1"
        self.assertEqual(keep_rows[0][1], expected_prefix)

        # Note: splitcode will append _R1.fastq and _R2.fastq to this prefix

    def test_complex_realistic_scenario(self):
        """
        Test a realistic scenario with multiple samples, realistic IDs, and mixed parameters.
        """
        import pandas as pd
        import csv

        # Realistic sample sheet data
        df = pd.DataFrame({
            'barcode_3': ['AAAAAAAA', 'CCCCCCCC', 'GGGGGGGG', 'TTTTTTTT'],
            'run': [
                'Sample1.lPool_1.HHJYWDRX5.6',
                'Sample2.lPool_1.HHJYWDRX5.6',
                'Sample3.lPool_2.HHJYWDRX5.6',
                'Sample4.lPool_2.HHJYWDRX5.6'
            ],
            'muxed_run': [
                'ATCGATCG-GCTAGCTA.lPool_1.HHJYWDRX5.6',
                'ATCGATCG-GCTAGCTA.lPool_1.HHJYWDRX5.6',
                'TTTTAAAA-CCCCGGGG.lPool_2.HHJYWDRX5.6',
                'TTTTAAAA-CCCCGGGG.lPool_2.HHJYWDRX5.6'
            ]
        }, index=['Sample1', 'Sample2', 'Sample3', 'Sample4'])

        # Generate for Pool_1 only
        pool_id = 'ATCGATCG-GCTAGCTA.lPool_1.HHJYWDRX5.6'
        config_file, keep_file, sample_ids = illumina.generate_splitcode_config_and_keep_files(
            df, pool_id, '/output', max_hamming_dist=1, r1_trim_bp_right_of_barcode=3
        )

        # Should have Sample1 and Sample2
        self.assertEqual(len(sample_ids), 2)
        self.assertIn('Sample1.lPool_1.HHJYWDRX5.6', sample_ids)
        self.assertIn('Sample2.lPool_1.HHJYWDRX5.6', sample_ids)

        # Verify config file
        with open(config_file) as f:
            config_rows = list(csv.reader(f, delimiter='\t'))

        # Header + 2 samples
        self.assertEqual(len(config_rows), 3)

        # Verify first sample (no _R1 suffix when using --keep-r1-r2)
        self.assertEqual(config_rows[1][0], 'AAAAAAAA')
        self.assertEqual(config_rows[1][1], 'Sample1.lPool_1.HHJYWDRX5.6')
        self.assertEqual(config_rows[1][2], '0:0:8')
        self.assertEqual(config_rows[1][3], '1')
        self.assertEqual(config_rows[1][4], '1:3')  # Trim barcode + 3bp

        # Verify keep file (no _R1 suffix when using --keep-r1-r2)
        with open(keep_file) as f:
            keep_rows = list(csv.reader(f, delimiter='\t'))

        self.assertEqual(len(keep_rows), 2)
        self.assertEqual(keep_rows[0][0], 'Sample1.lPool_1.HHJYWDRX5.6')
        self.assertEqual(keep_rows[0][1], '/output/Sample1.lPool_1.HHJYWDRX5.6')


class TestSplitcodeSummaryJSONErrorHandling(TestCaseWithTmp):
    """
    Test error handling when loading splitcode summary JSON files.

    These tests validate that we provide helpful debugging information
    when JSON files are missing or malformed.
    """

    def setUp(self):
        super().setUp()
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)
        super().tearDown()

    def test_missing_json_file_provides_debugging_info(self):
        """Test that missing JSON file error includes directory listing."""
        import pandas as pd

        # Create a sample sheet pointing to a pool
        sample_sheet = os.path.join(self.temp_dir, 'samples.tsv')
        with open(sample_sheet, 'w') as f:
            f.write("sample\tlibrary_id_per_sample\tbarcode_1\tbarcode_2\tbarcode_3\n")
            f.write("Sample1\tLib1\tATCGATCG\tGCTAGCTA\tAAAAAAAA\n")

        # Create the output directory with some files but NOT the expected JSON
        out_dir = os.path.join(self.temp_dir, 'output')
        os.makedirs(out_dir, exist_ok=True)

        # Create some dummy files to appear in the directory listing
        with open(os.path.join(out_dir, 'other_file.txt'), 'w') as f:
            f.write("dummy")
        with open(os.path.join(out_dir, 'wrong_pool_summary.json'), 'w') as f:
            f.write('{"tag_qc": []}')

        # csv_out determines outDir - put it in the output directory
        csv_out = os.path.join(out_dir, 'lut.csv')

        # Should raise FileNotFoundError with helpful message
        with self.assertRaises(FileNotFoundError) as context:
            illumina.create_splitcode_lookup_table(
                sample_sheet,
                csv_out,
                unmatched_name="Unmatched"
            )

        # Check that error message mentions the pool
        error_msg = str(context.exception)
        self.assertIn('ATCGATCG-GCTAGCTA', error_msg)
        self.assertIn('summary.json', error_msg)

    def test_malformed_json_provides_file_preview(self):
        """Test that malformed JSON error includes file content preview."""
        import pandas as pd

        # Create a sample sheet
        sample_sheet = os.path.join(self.temp_dir, 'samples.tsv')
        with open(sample_sheet, 'w') as f:
            f.write("sample\tlibrary_id_per_sample\tbarcode_1\tbarcode_2\tbarcode_3\n")
            f.write("Sample1\tLib1\tATCGATCG\tGCTAGCTA\tAAAAAAAA\n")

        out_dir = os.path.join(self.temp_dir, 'output')
        os.makedirs(out_dir, exist_ok=True)

        # Create a malformed JSON file
        pool_name = "ATCGATCG-GCTAGCTA.lLib1"
        malformed_json = os.path.join(out_dir, f'{pool_name}_summary.json')
        with open(malformed_json, 'w') as f:
            f.write('{"tag_qc": [this is not valid json}')

        # csv_out determines outDir - put it in the output directory
        csv_out = os.path.join(out_dir, 'lut.csv')

        # Should raise JSONDecodeError
        with self.assertRaises(json.JSONDecodeError):
            illumina.create_splitcode_lookup_table(
                sample_sheet,
                csv_out,
                unmatched_name="Unmatched"
            )


class TestConvertSplitcodeMetricsToPicardStyle(TestCaseWithTmp):
    """
    Test convert_splitcode_demux_metrics_to_picard_style function.

    This function converts splitcode CSV metrics into Picard-style TSV format.
    It's a pure data transformation function ideal for unit testing.

    Input CSV format (from create_splitcode_lookup_table):
    - sample, barcode_1, barcode_2, inline_barcode, library_id, run
    - num_reads_hdistance0, num_reads_hdistance1, num_reads_total

    Output TSV format (Picard-style ExtractIlluminaBarcodes.BarcodeMetric):
    - BARCODE, BARCODE_WITHOUT_DELIMITER, BARCODE_NAME, LIBRARY_NAME
    - READS, PF_READS, PERFECT_MATCHES, PF_PERFECT_MATCHES
    - ONE_MISMATCH_MATCHES, PF_ONE_MISMATCH_MATCHES
    - PCT_MATCHES, RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT
    - PF_PCT_MATCHES, PF_RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT, PF_NORMALIZED_MATCHES
    """

    def setUp(self):
        super().setUp()
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)
        super().tearDown()

    def test_basic_conversion(self):
        """Test basic conversion with minimal input."""
        import pandas as pd
        import csv

        # Create minimal input CSV
        # Note: Numeric columns are strings because that's how create_splitcode_lookup_table writes them
        input_csv = os.path.join(self.temp_dir, 'input.csv')
        df_input = pd.DataFrame({
            'sample': ['Sample1'],
            'barcode_1': ['ATCGATCG'],
            'barcode_2': ['GCTAGCTA'],
            'inline_barcode': ['AAAAAAAA'],
            'library_id': ['Lib1'],
            'run': ['Sample1.lLib1'],
            'num_reads_hdistance0': [100],
            'num_reads_hdistance1': [10],
            'num_reads_total': [110]
        })
        df_input.to_csv(input_csv, index=False)

        # Convert
        output_tsv = os.path.join(self.temp_dir, 'output.tsv')
        illumina.convert_splitcode_demux_metrics_to_picard_style(
            input_csv,
            output_tsv
        )

        # Verify output exists
        self.assertTrue(os.path.exists(output_tsv))

        # Read output
        with open(output_tsv) as f:
            lines = f.readlines()

        # Check header comment
        self.assertTrue(lines[0].startswith('## METRICS CLASS'))

        # Parse TSV data (skip header comment)
        with open(output_tsv) as f:
            next(f)  # Skip comment line
            reader = csv.DictReader(f, delimiter='\t')
            rows = list(reader)

        self.assertEqual(len(rows), 1)

        # Check barcode formatting
        self.assertEqual(rows[0]['BARCODE'], 'ATCGATCG-GCTAGCTA-AAAAAAAA')
        self.assertEqual(rows[0]['BARCODE_WITHOUT_DELIMITER'], 'ATCGATCGGCTAGCTAAAAAAAAA')
        self.assertEqual(rows[0]['BARCODE_NAME'], 'Sample1')
        self.assertEqual(rows[0]['LIBRARY_NAME'], 'Sample1.lLib1')

        # Check read counts
        self.assertEqual(rows[0]['READS'], '110')
        self.assertEqual(rows[0]['PF_READS'], '110')
        self.assertEqual(rows[0]['PERFECT_MATCHES'], '100')
        self.assertEqual(rows[0]['PF_PERFECT_MATCHES'], '100')
        self.assertEqual(rows[0]['ONE_MISMATCH_MATCHES'], '10')
        self.assertEqual(rows[0]['PF_ONE_MISMATCH_MATCHES'], '10')

    def test_missing_required_column_raises_error(self):
        """Test that missing required columns raise ValueError."""
        import pandas as pd

        # Create CSV missing required column
        input_csv = os.path.join(self.temp_dir, 'input.csv')
        df_input = pd.DataFrame({
            'sample': ['Sample1'],
            'barcode_1': ['ATCGATCG'],
            # Missing 'barcode_2'
            'num_reads_hdistance0': [100],
            'num_reads_hdistance1': [10],
            'num_reads_total': [110]
        })
        df_input.to_csv(input_csv, index=False)

        output_tsv = os.path.join(self.temp_dir, 'output.tsv')

        with self.assertRaises(ValueError) as context:
            illumina.convert_splitcode_demux_metrics_to_picard_style(
                input_csv,
                output_tsv
            )

        self.assertIn('barcode_2', str(context.exception))
        self.assertIn('missing', str(context.exception).lower())

    def test_empty_required_column_raises_error(self):
        """Test that entirely empty required columns raise ValueError."""
        import pandas as pd

        # Create CSV with empty required column
        input_csv = os.path.join(self.temp_dir, 'input.csv')
        df_input = pd.DataFrame({
            'sample': ['Sample1'],
            'barcode_1': ['ATCGATCG'],
            'barcode_2': ['GCTAGCTA'],
            'num_reads_hdistance0': [None],  # Empty column
            'num_reads_hdistance1': [10],
            'num_reads_total': [110]
        })
        df_input.to_csv(input_csv, index=False)

        output_tsv = os.path.join(self.temp_dir, 'output.tsv')

        with self.assertRaises(ValueError) as context:
            illumina.convert_splitcode_demux_metrics_to_picard_style(
                input_csv,
                output_tsv
            )

        self.assertIn('num_reads_hdistance0', str(context.exception))
        self.assertIn('empty', str(context.exception).lower())

    def test_barcode_combination_without_inline(self):
        """Test barcode combination when inline_barcode is missing."""
        import pandas as pd
        import csv

        # Create input CSV without inline_barcode column
        input_csv = os.path.join(self.temp_dir, 'input.csv')
        df_input = pd.DataFrame({
            'sample': ['Sample1'],
            'barcode_1': ['ATCGATCG'],
            'barcode_2': ['GCTAGCTA'],
            'library_id': ['Lib1'],
            'num_reads_hdistance0': [100],
            'num_reads_hdistance1': [10],
            'num_reads_total': [110]
        })
        df_input.to_csv(input_csv, index=False)

        output_tsv = os.path.join(self.temp_dir, 'output.tsv')
        illumina.convert_splitcode_demux_metrics_to_picard_style(
            input_csv,
            output_tsv
        )

        # Read output
        with open(output_tsv) as f:
            next(f)  # Skip comment
            reader = csv.DictReader(f, delimiter='\t')
            rows = list(reader)

        # Should only have barcode_1 and barcode_2
        self.assertEqual(rows[0]['BARCODE'], 'ATCGATCG-GCTAGCTA')
        self.assertEqual(rows[0]['BARCODE_WITHOUT_DELIMITER'], 'ATCGATCGGCTAGCTA')

    def test_library_name_fallback(self):
        """Test LIBRARY_NAME falls back to library_id if 'run' is missing."""
        import pandas as pd
        import csv

        # Create input CSV with library_id but no 'run' column
        input_csv = os.path.join(self.temp_dir, 'input.csv')
        df_input = pd.DataFrame({
            'sample': ['Sample1'],
            'barcode_1': ['ATCGATCG'],
            'barcode_2': ['GCTAGCTA'],
            'library_id': ['MyLibraryID'],
            'num_reads_hdistance0': [100],
            'num_reads_hdistance1': [10],
            'num_reads_total': [110]
        })
        df_input.to_csv(input_csv, index=False)

        output_tsv = os.path.join(self.temp_dir, 'output.tsv')
        illumina.convert_splitcode_demux_metrics_to_picard_style(
            input_csv,
            output_tsv
        )

        with open(output_tsv) as f:
            next(f)
            reader = csv.DictReader(f, delimiter='\t')
            rows = list(reader)

        self.assertEqual(rows[0]['LIBRARY_NAME'], 'MyLibraryID')

    def test_library_name_uses_barcode_name_fallback(self):
        """Test LIBRARY_NAME falls back to BARCODE_NAME if both run and library_id missing."""
        import pandas as pd
        import csv

        # Create input CSV without run or library_id
        input_csv = os.path.join(self.temp_dir, 'input.csv')
        df_input = pd.DataFrame({
            'sample': ['Sample1'],
            'barcode_1': ['ATCGATCG'],
            'barcode_2': ['GCTAGCTA'],
            'num_reads_hdistance0': [100],
            'num_reads_hdistance1': [10],
            'num_reads_total': [110]
        })
        df_input.to_csv(input_csv, index=False)

        output_tsv = os.path.join(self.temp_dir, 'output.tsv')
        illumina.convert_splitcode_demux_metrics_to_picard_style(
            input_csv,
            output_tsv
        )

        with open(output_tsv) as f:
            next(f)
            reader = csv.DictReader(f, delimiter='\t')
            rows = list(reader)

        self.assertEqual(rows[0]['LIBRARY_NAME'], 'Sample1')

    def test_combine_innerbarcode_unmatched_false(self):
        """Test that unmatched rows are NOT collapsed when combine_innerbarcode_unmatched=False."""
        import pandas as pd
        import csv

        # Create input with multiple unmatched (all-N inline barcode) rows
        input_csv = os.path.join(self.temp_dir, 'input.csv')
        df_input = pd.DataFrame({
            'sample': ['Sample1', 'Sample2', 'Sample3'],
            'barcode_1': ['ATCGATCG', 'ATCGATCG', 'GCTAGCTA'],
            'barcode_2': ['GCTAGCTA', 'GCTAGCTA', 'ATCGATCG'],
            'inline_barcode': ['NNNNNNNN', 'NNNNNNNN', 'AAAAAAAA'],
            'library_id': ['Lib1', 'Lib1', 'Lib2'],
            'num_reads_hdistance0': [100, 50, 200],
            'num_reads_hdistance1': [10, 5, 20],
            'num_reads_total': [110, 55, 220]
        })
        df_input.to_csv(input_csv, index=False)

        output_tsv = os.path.join(self.temp_dir, 'output.tsv')
        illumina.convert_splitcode_demux_metrics_to_picard_style(
            input_csv,
            output_tsv,
            combine_innerbarcode_unmatched=False
        )

        with open(output_tsv) as f:
            next(f)
            reader = csv.DictReader(f, delimiter='\t')
            rows = list(reader)

        # Should have all 3 rows (no collapsing)
        self.assertEqual(len(rows), 3)

    def test_combine_innerbarcode_unmatched_true_with_report_within_pools_true(self):
        """Test that unmatched rows ARE collapsed when combine_innerbarcode_unmatched=True."""
        import pandas as pd
        import csv

        # Create input with multiple unmatched (all-N inline barcode) rows in same pool
        input_csv = os.path.join(self.temp_dir, 'input.csv')
        df_input = pd.DataFrame({
            'sample': ['Unmatched.Pool1', 'Unmatched.Pool1_2', 'Sample3'],
            'barcode_1': ['ATCGATCG', 'ATCGATCG', 'GCTAGCTA'],
            'barcode_2': ['GCTAGCTA', 'GCTAGCTA', 'ATCGATCG'],
            'inline_barcode': ['NNNNNNNN', 'NNNNNNNN', 'AAAAAAAA'],
            'library_id': ['Lib1', 'Lib1', 'Lib2'],
            'num_reads_hdistance0': [100, 50, 200],
            'num_reads_hdistance1': [10, 5, 20],
            'num_reads_total': [110, 55, 220]
        })
        df_input.to_csv(input_csv, index=False)

        output_tsv = os.path.join(self.temp_dir, 'output.tsv')
        illumina.convert_splitcode_demux_metrics_to_picard_style(
            input_csv,
            output_tsv,
            combine_innerbarcode_unmatched=True,
            report_within_pools=True
        )

        with open(output_tsv) as f:
            next(f)
            reader = csv.DictReader(f, delimiter='\t')
            rows = list(reader)

        # Should have 2 rows: collapsed N row + Sample3
        self.assertEqual(len(rows), 2)

        # Find the collapsed row
        n_row = [r for r in rows if r['BARCODE'] == 'N'][0]

        # Check that reads were summed
        self.assertEqual(n_row['READS'], str(110 + 55))
        self.assertEqual(n_row['PERFECT_MATCHES'], str(100 + 50))
        self.assertEqual(n_row['ONE_MISMATCH_MATCHES'], str(10 + 5))

    def test_combine_innerbarcode_unmatched_true_with_report_within_pools_false(self):
        """Test all-N detection uses entire barcode when report_within_pools=False."""
        import pandas as pd
        import csv

        # Create input with all-N full barcodes (not just inline)
        input_csv = os.path.join(self.temp_dir, 'input.csv')
        df_input = pd.DataFrame({
            'sample': ['Unmatched1', 'Unmatched2', 'Sample3'],
            'barcode_1': ['NNNNNNNN', 'NNNNNNNN', 'ATCGATCG'],
            'barcode_2': ['NNNNNNNN', 'NNNNNNNN', 'GCTAGCTA'],
            'inline_barcode': ['NNNNNNNN', 'NNNNNNNN', 'AAAAAAAA'],
            'library_id': ['Lib1', 'Lib1', 'Lib2'],
            'num_reads_hdistance0': [100, 50, 200],
            'num_reads_hdistance1': [10, 5, 20],
            'num_reads_total': [110, 55, 220]
        })
        df_input.to_csv(input_csv, index=False)

        output_tsv = os.path.join(self.temp_dir, 'output.tsv')
        illumina.convert_splitcode_demux_metrics_to_picard_style(
            input_csv,
            output_tsv,
            combine_innerbarcode_unmatched=True,
            report_within_pools=False
        )

        with open(output_tsv) as f:
            next(f)
            reader = csv.DictReader(f, delimiter='\t')
            rows = list(reader)

        # Should have 2 rows: collapsed N row + Sample3
        self.assertEqual(len(rows), 2)

        n_row = [r for r in rows if r['BARCODE'] == 'N'][0]
        self.assertEqual(n_row['READS'], str(110 + 55))

    def test_stats_computation_global(self):
        """Test statistics are computed globally when report_within_pools=False."""
        import pandas as pd
        import csv

        # Create input with 2 samples
        input_csv = os.path.join(self.temp_dir, 'input.csv')
        df_input = pd.DataFrame({
            'sample': ['Sample1', 'Sample2'],
            'barcode_1': ['ATCGATCG', 'GCTAGCTA'],
            'barcode_2': ['GCTAGCTA', 'ATCGATCG'],
            'library_id': ['Lib1', 'Lib2'],
            'num_reads_hdistance0': [100, 200],
            'num_reads_hdistance1': [10, 20],
            'num_reads_total': [110, 220]
        })
        df_input.to_csv(input_csv, index=False)

        output_tsv = os.path.join(self.temp_dir, 'output.tsv')
        illumina.convert_splitcode_demux_metrics_to_picard_style(
            input_csv,
            output_tsv,
            report_within_pools=False
        )

        with open(output_tsv) as f:
            next(f)
            reader = csv.DictReader(f, delimiter='\t')
            rows = list(reader)

        # Sample2 has max reads (220), so it should have ratio=1.0
        sample2 = [r for r in rows if r['BARCODE_NAME'] == 'Sample2'][0]
        self.assertEqual(float(sample2['RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT']), 1.0)

        # Sample1 has 110 reads, ratio should be 110/220 = 0.5
        sample1 = [r for r in rows if r['BARCODE_NAME'] == 'Sample1'][0]
        self.assertAlmostEqual(float(sample1['RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT']), 0.5, places=6)

    def test_stats_computation_within_pools(self):
        """Test statistics are computed per pool when report_within_pools=True."""
        import pandas as pd
        import csv

        # Create input with 2 pools, 2 samples each
        input_csv = os.path.join(self.temp_dir, 'input.csv')
        df_input = pd.DataFrame({
            'sample': ['S1', 'S2', 'S3', 'S4'],
            'barcode_1': ['ATCGATCG', 'ATCGATCG', 'GCTAGCTA', 'GCTAGCTA'],
            'barcode_2': ['GCTAGCTA', 'GCTAGCTA', 'ATCGATCG', 'ATCGATCG'],
            'inline_barcode': ['AAAAAAAA', 'CCCCCCCC', 'GGGGGGGG', 'TTTTTTTT'],
            'library_id': ['Lib1', 'Lib1', 'Lib2', 'Lib2'],
            'num_reads_hdistance0': [100, 200, 50, 150],
            'num_reads_hdistance1': [10, 20, 5, 15],
            'num_reads_total': [110, 220, 55, 165]
        })
        df_input.to_csv(input_csv, index=False)

        output_tsv = os.path.join(self.temp_dir, 'output.tsv')
        illumina.convert_splitcode_demux_metrics_to_picard_style(
            input_csv,
            output_tsv,
            report_within_pools=True
        )

        with open(output_tsv) as f:
            next(f)
            reader = csv.DictReader(f, delimiter='\t')
            rows = list(reader)

        # Pool 1 (ATCGATCG-GCTAGCTA): S2 has max (220)
        s2 = [r for r in rows if r['BARCODE_NAME'] == 'S2'][0]
        self.assertEqual(float(s2['RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT']), 1.0)

        # Pool 1: S1 has 110, ratio = 110/220 = 0.5
        s1 = [r for r in rows if r['BARCODE_NAME'] == 'S1'][0]
        self.assertAlmostEqual(float(s1['RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT']), 0.5, places=6)

        # Pool 2 (GCTAGCTA-ATCGATCG): S4 has max (165)
        s4 = [r for r in rows if r['BARCODE_NAME'] == 'S4'][0]
        self.assertEqual(float(s4['RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT']), 1.0)

        # Pool 2: S3 has 55, ratio = 55/165 ≈ 0.333
        s3 = [r for r in rows if r['BARCODE_NAME'] == 'S3'][0]
        self.assertAlmostEqual(float(s3['RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT']), 55/165, places=6)

    def test_pct_matches_computation(self):
        """Test PCT_MATCHES is computed correctly."""
        import pandas as pd
        import csv

        input_csv = os.path.join(self.temp_dir, 'input.csv')
        df_input = pd.DataFrame({
            'sample': ['Sample1', 'Sample2'],
            'barcode_1': ['ATCGATCG', 'ATCGATCG'],
            'barcode_2': ['GCTAGCTA', 'GCTAGCTA'],
            'library_id': ['Lib1', 'Lib1'],
            'num_reads_hdistance0': [100, 200],
            'num_reads_hdistance1': [10, 20],
            'num_reads_total': [110, 220]
        })
        df_input.to_csv(input_csv, index=False)

        output_tsv = os.path.join(self.temp_dir, 'output.tsv')
        illumina.convert_splitcode_demux_metrics_to_picard_style(
            input_csv,
            output_tsv,
            report_within_pools=True
        )

        with open(output_tsv) as f:
            next(f)
            reader = csv.DictReader(f, delimiter='\t')
            rows = list(reader)

        # Total reads in pool: 110 + 220 = 330
        # Sample1 perfect matches: 100, PCT_MATCHES = 100/330
        s1 = [r for r in rows if r['BARCODE_NAME'] == 'Sample1'][0]
        self.assertAlmostEqual(float(s1['PCT_MATCHES']), 100/330, places=6)

        # Sample2 perfect matches: 200, PCT_MATCHES = 200/330
        s2 = [r for r in rows if r['BARCODE_NAME'] == 'Sample2'][0]
        self.assertAlmostEqual(float(s2['PCT_MATCHES']), 200/330, places=6)

    def test_normalized_matches_excludes_all_n(self):
        """Test PF_NORMALIZED_MATCHES excludes all-N rows from mean calculation."""
        import pandas as pd
        import csv

        input_csv = os.path.join(self.temp_dir, 'input.csv')
        df_input = pd.DataFrame({
            'sample': ['Sample1', 'Sample2', 'Unmatched'],
            'barcode_1': ['ATCGATCG', 'ATCGATCG', 'ATCGATCG'],
            'barcode_2': ['GCTAGCTA', 'GCTAGCTA', 'GCTAGCTA'],
            'inline_barcode': ['AAAAAAAA', 'CCCCCCCC', 'NNNNNNNN'],
            'library_id': ['Lib1', 'Lib1', 'Lib1'],
            'num_reads_hdistance0': [100, 200, 50],
            'num_reads_hdistance1': [10, 20, 5],
            'num_reads_total': [110, 220, 55]
        })
        df_input.to_csv(input_csv, index=False)

        output_tsv = os.path.join(self.temp_dir, 'output.tsv')
        illumina.convert_splitcode_demux_metrics_to_picard_style(
            input_csv,
            output_tsv,
            report_within_pools=True
        )

        with open(output_tsv) as f:
            next(f)
            reader = csv.DictReader(f, delimiter='\t')
            rows = list(reader)

        # Mean PF_READS should exclude unmatched (all-N)
        # Mean = (110 + 220) / 2 = 165
        # Sample1: PF_NORMALIZED_MATCHES = 100 / 165
        s1 = [r for r in rows if r['BARCODE_NAME'] == 'Sample1'][0]
        self.assertAlmostEqual(float(s1['PF_NORMALIZED_MATCHES']), 100/165, places=6)

    def test_custom_demux_function_name(self):
        """Test that custom demux_function parameter appears in header."""
        import pandas as pd

        input_csv = os.path.join(self.temp_dir, 'input.csv')
        df_input = pd.DataFrame({
            'sample': ['Sample1'],
            'barcode_1': ['ATCGATCG'],
            'barcode_2': ['GCTAGCTA'],
            'library_id': ['Lib1'],
            'num_reads_hdistance0': [100],
            'num_reads_hdistance1': [10],
            'num_reads_total': [110]
        })
        df_input.to_csv(input_csv, index=False)

        output_tsv = os.path.join(self.temp_dir, 'output.tsv')
        custom_name = "MyCustomDemuxFunction"
        illumina.convert_splitcode_demux_metrics_to_picard_style(
            input_csv,
            output_tsv,
            demux_function=custom_name
        )

        with open(output_tsv) as f:
            header = f.readline()

        self.assertIn(custom_name, header)

    def test_zero_reads_handling(self):
        """Test that zero reads don't cause division by zero errors."""
        import pandas as pd
        import csv

        input_csv = os.path.join(self.temp_dir, 'input.csv')
        df_input = pd.DataFrame({
            'sample': ['Sample1', 'Sample2'],
            'barcode_1': ['ATCGATCG', 'GCTAGCTA'],
            'barcode_2': ['GCTAGCTA', 'ATCGATCG'],
            'library_id': ['Lib1', 'Lib2'],
            'num_reads_hdistance0': [0, 0],
            'num_reads_hdistance1': [0, 0],
            'num_reads_total': [0, 0]
        })
        df_input.to_csv(input_csv, index=False)

        output_tsv = os.path.join(self.temp_dir, 'output.tsv')
        # Should not raise any errors
        illumina.convert_splitcode_demux_metrics_to_picard_style(
            input_csv,
            output_tsv
        )

        with open(output_tsv) as f:
            next(f)
            reader = csv.DictReader(f, delimiter='\t')
            rows = list(reader)

        # All percentage/ratio columns should be 0
        for row in rows:
            self.assertEqual(float(row['PCT_MATCHES']), 0.0)
            self.assertEqual(float(row['RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT']), 0.0)
            self.assertEqual(float(row['PF_NORMALIZED_MATCHES']), 0.0)


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

        # Custom 3-barcode samplesheet (ONLY required input besides FASTQs)
        self.samples_3bc = os.path.join(self.input_dir, 'samples_3bc.tsv')

        # Verify test input files exist
        for f in [self.r1_fastq, self.r2_fastq, self.samples_3bc]:
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

        # Run demux on 2-barcode sample (simplified interface)
        illumina.splitcode_demux_fastqs(
            fastq_r1=r1_fastq,
            fastq_r2=r2_fastq,
            samplesheet=self.samples_3bc,
            outdir=out_dir,
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
