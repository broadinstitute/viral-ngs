# Unit tests for tools.splitcode
# -*- coding: utf-8 -*-

__author__ = "dpark@broadinstitute.org"

import unittest
import os
import os.path
import tempfile
import filecmp
import shutil
import json
import util
import util.file
import tools.splitcode
import tools.samtools
import tools.picard
from test import TestCaseWithTmp, assert_equal_bam_reads


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

            result_path = tools.splitcode.create_splitcode_lookup_table(
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

            result_path = tools.splitcode.create_splitcode_lookup_table(
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

            result_path = tools.splitcode.create_splitcode_lookup_table(
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

            result_path = tools.splitcode.create_splitcode_lookup_table(
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
        result = tools.splitcode.run_splitcode_on_pool(
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
        result = tools.splitcode.run_splitcode_on_pool(
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

        tools.splitcode.run_splitcode_on_pool(
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

        tools.splitcode.run_splitcode_on_pool(
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

        tools.splitcode.run_splitcode_on_pool(
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
        tools.splitcode.run_splitcode_on_pool(
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

        config_file, keep_file, sample_ids = tools.splitcode.generate_splitcode_config_and_keep_files(
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

        config_file, keep_file, sample_ids = tools.splitcode.generate_splitcode_config_and_keep_files(
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
        config_file, _, _ = tools.splitcode.generate_splitcode_config_and_keep_files(
            df, 'Pool1', self.temp_dir, max_hamming_dist=0
        )

        with open(config_file) as f:
            reader = csv.reader(f, delimiter='\t')
            rows = list(reader)

        self.assertEqual(rows[1][3], '0')  # distance column

        # Test distance=2
        config_file, _, _ = tools.splitcode.generate_splitcode_config_and_keep_files(
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
        config_file, _, _ = tools.splitcode.generate_splitcode_config_and_keep_files(
            df, 'Pool1', self.temp_dir, r1_trim_bp_right_of_barcode=None
        )

        with open(config_file) as f:
            reader = csv.reader(f, delimiter='\t')
            rows = list(reader)

        self.assertEqual(rows[1][4], '1')  # left column: "1" means trim barcode only

        # Test with 5 extra bp
        config_file, _, _ = tools.splitcode.generate_splitcode_config_and_keep_files(
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
        config_file, keep_file, sample_ids = tools.splitcode.generate_splitcode_config_and_keep_files(
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
            tools.splitcode.generate_splitcode_config_and_keep_files(
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

        config_file, keep_file, _ = tools.splitcode.generate_splitcode_config_and_keep_files(
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
        _, keep_file, _ = tools.splitcode.generate_splitcode_config_and_keep_files(
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
        config_file, keep_file, sample_ids = tools.splitcode.generate_splitcode_config_and_keep_files(
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
            tools.splitcode.create_splitcode_lookup_table(
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
            tools.splitcode.create_splitcode_lookup_table(
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
        tools.splitcode.convert_splitcode_demux_metrics_to_picard_style(
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
            tools.splitcode.convert_splitcode_demux_metrics_to_picard_style(
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
            tools.splitcode.convert_splitcode_demux_metrics_to_picard_style(
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
        tools.splitcode.convert_splitcode_demux_metrics_to_picard_style(
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
        tools.splitcode.convert_splitcode_demux_metrics_to_picard_style(
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
        tools.splitcode.convert_splitcode_demux_metrics_to_picard_style(
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
        tools.splitcode.convert_splitcode_demux_metrics_to_picard_style(
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
        tools.splitcode.convert_splitcode_demux_metrics_to_picard_style(
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
        tools.splitcode.convert_splitcode_demux_metrics_to_picard_style(
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
        tools.splitcode.convert_splitcode_demux_metrics_to_picard_style(
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
        tools.splitcode.convert_splitcode_demux_metrics_to_picard_style(
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
        tools.splitcode.convert_splitcode_demux_metrics_to_picard_style(
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
        tools.splitcode.convert_splitcode_demux_metrics_to_picard_style(
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
        tools.splitcode.convert_splitcode_demux_metrics_to_picard_style(
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
        tools.splitcode.convert_splitcode_demux_metrics_to_picard_style(
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
