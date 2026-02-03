# Unit tests for ncbi.py

__author__ = "tomkinsc@broadinstitute.org"

# built-ins
import os
import tempfile
import unittest
import argparse

# module-specific
import viral_ngs.ncbi
import viral_ngs.core.file
import viral_ngs.phylo.genbank
from tests import (assert_equal_bam_reads, TestCaseWithTmp, assert_equal_contents,
                   assert_md5_equal_to_line_in_file, assert_valid_feature_table)


class TestCommandHelp(unittest.TestCase):

    def test_help_parser_for_each_command(self):
        for cmd_name, parser_fun in ncbi.__commands__:
            parser = parser_fun(argparse.ArgumentParser())
            helpstring = parser.format_help()

class TestFeatureReader(TestCaseWithTmp):
    def setUp(self):
        super(TestFeatureReader, self).setUp()
        self.input_dir = viral_ngs.core.file.get_test_input_path(self)

    def test_read_seq_id_simple(self):
        accessions = ('GU481072.1', 'GU481073.1',
            'KM821772.1', 'KM821773.1')
        for acc in accessions:
            self.assertEqual(acc, viral_ngs.phylo.genbank.get_feature_table_id(os.path.join(self.input_dir, acc+'.tbl')))

    def test_read_seq_id_different_fnames(self):
        self.assertEqual('KM821998.1', viral_ngs.phylo.genbank.get_feature_table_id(os.path.join(self.input_dir,
            'test1-S.tbl')))
        self.assertEqual('KM821997.1', viral_ngs.phylo.genbank.get_feature_table_id(os.path.join(self.input_dir,
            'test2-L.tbl')))

    def test_read_seq_id_refseq(self):
        self.assertEqual('NC_026438.1', viral_ngs.phylo.genbank.get_feature_table_id(os.path.join(self.input_dir,
            'NC_026438.1.tbl')))

    def test_read_seq_id_ddbj(self):
        self.assertEqual('LC889323.1', viral_ngs.phylo.genbank.get_feature_table_id(os.path.join(self.input_dir,
            'LC889323.1.tbl')))

    def test_seq_location_str_format(self):
        """Test SeqLocation string format - continuation lines should not have trailing tabs."""
        from phylo.feature_table import SeqLocation, SeqPosition

        # First line with feature type - should have 3 columns
        first = SeqLocation(SeqPosition(100), SeqPosition(200), feature_type='CDS')
        self.assertEqual(str(first), "100\t200\tCDS")

        # Continuation line without feature type - should have 2 columns, NO trailing tab
        continuation = SeqLocation(SeqPosition(300), SeqPosition(400), feature_type=None)
        self.assertEqual(str(continuation), "300\t400")  # NOT "300\t400\t"

class TestFeatureTransfer(TestCaseWithTmp):
    def setUp(self):
        super(TestFeatureTransfer, self).setUp()
        self.input_dir = viral_ngs.core.file.get_test_input_path(self)

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

        # Validate with table2asn
        assert_valid_feature_table(
            self,
            out_tbl,
            os.path.join(input_dir, "aligned_1.fasta"),
            temp_dir
        )

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

        # Verify output is valid feature table syntax
        # First interval line has 3 fields (start, end, feature_type)
        # Continuation lines have exactly 2 fields (start, end) - NO trailing tab
        with open(out_tbl, 'r') as f:
            for line in f:
                line = line.rstrip('\n')
                if not line or line.startswith('>') or line.startswith('\t'):
                    continue
                # Lines with coordinates must be properly tab-delimited
                parts = line.split('\t')
                if len(parts) >= 2 and parts[0] and parts[1]:
                    # Must have 2 or 3 columns
                    self.assertTrue(len(parts) in (2, 3),
                        f"Invalid line (must have 2 or 3 tab-delimited fields): {line}")
                    # If 3 columns, the 3rd must not be empty (that's a trailing tab)
                    if len(parts) == 3:
                        self.assertNotEqual(parts[2], '',
                            f"Line has trailing tab (continuation lines should have 2 columns): {line}")

        assert_equal_contents(self, out_tbl, expected)

        # Validate with table2asn
        assert_valid_feature_table(
            self,
            out_tbl,
            os.path.join(input_dir, "aligned_1.fasta"),
            temp_dir
        )

    def test_partial_symbols_column_placement(self):
        """Test that < only appears in column 1 and > only appears in column 2.

        Issue #72: Partial stop coordinates incorrectly used '<' instead of '>'.
        Per NCBI spec, '<' always goes in column 1, '>' always goes in column 2,
        regardless of strand.
        """
        input_dir = os.path.join(self.input_dir, "negative_strand_partial", "input")
        expected_dir = os.path.join(self.input_dir, "negative_strand_partial", "expected")
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

        # Verify partial symbols are in correct columns per NCBI spec
        with open(out_tbl, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.rstrip('\n')
                if not line or line.startswith('>Feature') or line.startswith('\t'):
                    continue
                parts = line.split('\t')
                if len(parts) >= 2:
                    col1, col2 = parts[0], parts[1]
                    # < must only appear in column 1 (start position)
                    self.assertFalse(col2.startswith('<'),
                        f"Line {line_num}: '<' in column 2 is invalid per NCBI spec (must be '>'): {line}")
                    # > must only appear in column 2 (end position)
                    self.assertFalse(col1.startswith('>'),
                        f"Line {line_num}: '>' in column 1 is invalid per NCBI spec (must be '<'): {line}")

        assert_equal_contents(self, out_tbl, expected)

        # Validate with table2asn
        assert_valid_feature_table(
            self,
            out_tbl,
            os.path.join(input_dir, "aligned_1.fasta"),
            temp_dir
        )

    def test_internal_partials_dropped(self):
        """Test that multi-interval features with internal partials are dropped.

        Issue #74: table2asn rejects features where partial symbols appear
        on internal intervals rather than the true 5'/3' ends. When a
        multi-interval feature is truncated such that a partial symbol would
        appear on an internal interval, the entire feature should be dropped.
        """
        input_dir = os.path.join(self.input_dir, "internal_partials", "input")
        expected_dir = os.path.join(self.input_dir, "internal_partials", "expected")
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

        # Feature with internal partials should be dropped entirely
        assert_equal_contents(self, out_tbl, expected)

        # Validate with table2asn - should pass since feature was dropped
        assert_valid_feature_table(
            self,
            out_tbl,
            os.path.join(input_dir, "aligned_1.fasta"),
            temp_dir
        )

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
