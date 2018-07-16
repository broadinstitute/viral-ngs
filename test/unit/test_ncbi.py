# Unit tests for ncbi.py

__author__ = "tomkinsc@broadinstitute.org"

# built-ins
import os
import tempfile
import unittest
import argparse

# module-specific
import ncbi
import util.file
from test import assert_equal_bam_reads, TestCaseWithTmp, assert_equal_contents, assert_md5_equal_to_line_in_file


class TestCommandHelp(unittest.TestCase):

    def test_help_parser_for_each_command(self):
        for cmd_name, parser_fun in ncbi.__commands__:
            parser = parser_fun(argparse.ArgumentParser())
            helpstring = parser.format_help()

class TestFeatureTransfer(TestCaseWithTmp):
    def setUp(self):
        super(TestFeatureTransfer, self).setUp()
        self.input_dir = util.file.get_test_input_path(self)

    def test_synthetic_feature_table(self):
        input_dir = os.path.join(self.input_dir, "synthetic", "input")
        expected_dir = os.path.join(self.input_dir, "synthetic", "expected")
        temp_dir = tempfile.gettempdir()

        in_tbl = os.path.join(input_dir,"ref.tbl")
        out_tbl = os.path.join(temp_dir,"sample.tbl")
        expected = os.path.join(expected_dir, "mapped.tbl")

        ncbi.tbl_transfer_prealigned( 
                                    os.path.join(input_dir,"aligned_1.fasta"),
                                    os.path.join(input_dir,"ref.fasta"), 
                                    [in_tbl], 
                                    temp_dir)

        assert_equal_contents(self, out_tbl, expected)

    def test_synthetic_feature_table_oob_clip(self):
        input_dir    = os.path.join(self.input_dir, "synthetic_oob_clip", "input")
        expected_dir = os.path.join(self.input_dir, "synthetic_oob_clip", "expected")
        temp_dir     = tempfile.gettempdir()

        in_tbl = os.path.join(input_dir,"ref.tbl")
        out_tbl = os.path.join(temp_dir,"sample.tbl")
        expected = os.path.join(expected_dir, "mapped.tbl")

        ncbi.tbl_transfer_prealigned( 
                                    os.path.join(input_dir,"aligned_1.fasta"),
                                    os.path.join(input_dir,"ref.fasta"), 
                                    [in_tbl], 
                                    temp_dir,
                                    oob_clip=True)

        assert_equal_contents(self, out_tbl, expected)

    def test_synthetic_feature_table_ignore_ambiguous_edges(self):
        input_dir    = os.path.join(self.input_dir, "synthetic_ignore_ambig_edges", "input")
        expected_dir = os.path.join(self.input_dir, "synthetic_ignore_ambig_edges", "expected")
        temp_dir     = tempfile.gettempdir()

        in_tbl = os.path.join(input_dir,"ref.tbl")
        out_tbl = os.path.join(temp_dir,"sample.tbl")
        expected = os.path.join(expected_dir, "mapped.tbl")

        ncbi.tbl_transfer_prealigned( 
                                    os.path.join(input_dir,"aligned_1.fasta"),
                                    os.path.join(input_dir,"ref.fasta"), 
                                    [in_tbl], 
                                    temp_dir,
                                    ignore_ambig_feature_edge=True)

        assert_equal_contents(self, out_tbl, expected)

    def test_lasv_oob_clip(self):
        input_dir    = os.path.join(self.input_dir, "lasv", "input")
        expected_dir = os.path.join(self.input_dir, "lasv", "expected")
        temp_dir     = tempfile.gettempdir()

        infastas = [os.path.join(input_dir, f) for f in [
                "align_mafft-ref-lasv-ISTH2376_1.fasta",
                "align_mafft-ref-lasv-ISTH2376_2.fasta"
            ]
        ]

        intables =[os.path.join(input_dir, f) for f in [
                "KM821997.1.tbl",
                "KM821998.1.tbl"
            ]
        ]

        out_table_names = [
            "LASV_NGA_2018_0026-1.tbl",
            "LASV_NGA_2018_0026-2.tbl",
            "LASV_NGA_2018_0097-1.tbl",
            "LASV_NGA_2018_0097-2.tbl",
            "LASV_NGA_2018_0541-1.tbl",
            "LASV_NGA_2018_0541-2.tbl",
            "LASV_NGA_2018_0611-1.tbl",
            "LASV_NGA_2018_0611-2.tbl",
            "LASV_NGA_2018_0664-1.tbl",
            "LASV_NGA_2018_0664-2.tbl",
            "LASV_NGA_2018_0959-1.tbl",
            "LASV_NGA_2018_0959-2.tbl",
            "LASV_NGA_2018_0998-1.tbl",
            "LASV_NGA_2018_0998-2.tbl",
            "LASV_NGA_2018_1024-1.tbl",
            "LASV_NGA_2018_1024-2.tbl",
            "LASV_NGA_2018_1079-1.tbl",
            "LASV_NGA_2018_1079-2.tbl",
            "LASV_NGA_2018_1177-1.tbl",
            "LASV_NGA_2018_1177-2.tbl",
            "LASV_NGA_2018_1375-1.tbl",
            "LASV_NGA_2018_1375-2.tbl",
            "LASV_NGA_2018_1381-1.tbl",
            "LASV_NGA_2018_1381-2.tbl",
            "LASV_NGA_2018_1392-1.tbl",
            "LASV_NGA_2018_1392-2.tbl",
            "LASV_NGA_2018_1643-1.tbl",
            "LASV_NGA_2018_1643-2.tbl"
        ]
        out_tbls =[os.path.join(temp_dir, f) for f in out_table_names]
        expected_tbls = [os.path.join(expected_dir, f) for f in out_table_names]


        for i in range(0, len(infastas)):
            ncbi.tbl_transfer_prealigned( 
                                        infastas[i],
                                        os.path.join(input_dir,"ref-lasv-ISTH2376.fasta"), 
                                        intables, 
                                        temp_dir,
                                        oob_clip=True)

        for i in range(0,len(out_table_names)):
            out_tbl = out_tbls[i]
            expected_tbl = expected_tbls[i]
            assert_equal_contents(self, out_tbl, expected_tbl)
