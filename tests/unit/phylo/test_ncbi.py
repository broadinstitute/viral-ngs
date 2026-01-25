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
import phylo.genbank
from test import assert_equal_bam_reads, TestCaseWithTmp, assert_equal_contents, assert_md5_equal_to_line_in_file


class TestCommandHelp(unittest.TestCase):

    def test_help_parser_for_each_command(self):
        for cmd_name, parser_fun in ncbi.__commands__:
            parser = parser_fun(argparse.ArgumentParser())
            helpstring = parser.format_help()

class TestFeatureReader(TestCaseWithTmp):
    def setUp(self):
        super(TestFeatureReader, self).setUp()
        self.input_dir = util.file.get_test_input_path(self)

    def test_read_seq_id_simple(self):
        accessions = ('GU481072.1', 'GU481073.1',
            'KM821772.1', 'KM821773.1')
        for acc in accessions:
            self.assertEqual(acc, phylo.genbank.get_feature_table_id(os.path.join(self.input_dir, acc+'.tbl')))

    def test_read_seq_id_different_fnames(self):
        self.assertEqual('KM821998.1', phylo.genbank.get_feature_table_id(os.path.join(self.input_dir,
            'test1-S.tbl')))
        self.assertEqual('KM821997.1', phylo.genbank.get_feature_table_id(os.path.join(self.input_dir,
            'test2-L.tbl')))

    def test_read_seq_id_refseq(self):
        self.assertEqual('NC_026438.1', phylo.genbank.get_feature_table_id(os.path.join(self.input_dir,
            'NC_026438.1.tbl')))

    def test_read_seq_id_ddbj(self):
        self.assertEqual('LC889323.1', phylo.genbank.get_feature_table_id(os.path.join(self.input_dir,
            'LC889323.1.tbl')))

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

    def test_severely_truncated_assembly_oob_clip(self):
        """Test annotation transfer with severely truncated assembly (adenovirus-like scenario).

        This reproduces issue #67 where a ~50% truncated assembly creates invalid
        feature table syntax with coordinates but no feature keys.
        """
        input_dir = os.path.join(self.input_dir, "adenovirus_truncated", "input")
        expected_dir = os.path.join(self.input_dir, "adenovirus_truncated", "expected")
        temp_dir = tempfile.gettempdir()

        in_tbl = os.path.join(input_dir, "ref.tbl")
        out_tbl = os.path.join(temp_dir, "sample.tbl")
        expected = os.path.join(expected_dir, "mapped.tbl")

        ncbi.tbl_transfer_prealigned(
            os.path.join(input_dir, "aligned_1.fasta"),
            os.path.join(input_dir, "ref.fasta"),
            [in_tbl],
            temp_dir,
            oob_clip=True
        )

        # Verify output is valid feature table syntax (no blank lines with coordinates)
        with open(out_tbl, 'r') as f:
            for line in f:
                line = line.rstrip('\n')
                if not line or line.startswith('>') or line.startswith('\t'):
                    continue
                # Lines with coordinates must have a feature key (third column)
                parts = line.split('\t')
                if len(parts) >= 2 and parts[0] and parts[1]:
                    self.assertTrue(len(parts) >= 3 and parts[2],
                        f"Invalid line (coordinates without feature key): {line}")

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
