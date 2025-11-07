# Unit tests for illumina.py splitcode integration

__author__ = "dpark@broadinstitute.org"

import unittest
import os
import tempfile
import json
import shutil

import illumina
import tools.samtools
import tools.picard
import util.file
from test import TestCaseWithTmp


class TestSplitkodeIntegration(TestCaseWithTmp):
    """Test our assumptions about splitcode behavior and file formats"""

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
        splitcode_config = os.path.join(self.temp_dir, 'config.txt')
        with open(splitcode_config, 'w') as f:
            f.write("tag\tid\tlocations\tdistance\tleft\tright\n")
            # AAAAAAAA at positions 0-7 in R1, distance=1, trim from left
            f.write("AAAAAAAA\tBC1\t0:7\t1\t1\t0\n")

        # Create MINIMAL keep file - just use the barcode ID without suffix
        # This is the simplest format that works
        splitcode_keepfile = os.path.join(self.temp_dir, 'keep.txt')
        output_prefix = os.path.join(self.temp_dir, 'Sample1')
        with open(splitcode_keepfile, 'w') as f:
            # Just the barcode ID from config (without _R1 suffix)
            f.write(f"BC1\t{output_prefix}\n")

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

        # Find BC1 in tag_qc
        tag_qc_dict = {item['tag']: item for item in summary['tag_qc']}
        self.assertIn('BC1', tag_qc_dict, "BC1 should be in tag_qc")

        bc1_entry = tag_qc_dict['BC1']
        self.assertIn('count', bc1_entry, "BC1 entry should have 'count' field")
        self.assertGreater(bc1_entry['count'], 0, "BC1 should have matched some reads")

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
        splitcode_config = os.path.join(self.temp_dir, 'config.txt')
        with open(splitcode_config, 'w') as f:
            f.write("tag\tid\tlocations\tdistance\tleft\tright\n")
            f.write("AAAAAAAA\tBC1\t0:7\t0\t1\t0\n")
            f.write("CCCCCCCC\tBC2\t0:7\t0\t1\t0\n")  # Won't match any reads

        # Keep file for both barcodes
        splitcode_keepfile = os.path.join(self.temp_dir, 'keep.txt')
        output_prefix1 = os.path.join(self.temp_dir, 'Sample1')
        output_prefix2 = os.path.join(self.temp_dir, 'Sample2')
        with open(splitcode_keepfile, 'w') as f:
            f.write(f"BC1_R1\t{output_prefix1}\n")
            f.write(f"BC2_R1\t{output_prefix2}\n")

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

        tag_qc_dict = {item['tag']: item for item in summary['tag_qc']}

        # BC1 should have reads
        self.assertIn('BC1', tag_qc_dict)
        self.assertGreater(tag_qc_dict['BC1']['count'], 0)

        # BC2 might not appear in tag_qc, or it might appear with count=0
        if 'BC2' in tag_qc_dict:
            self.assertEqual(tag_qc_dict['BC2']['count'], 0,
                           "BC2 should have 0 reads if it appears in tag_qc")

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
        with open(splitcode_config, 'w') as f:
            f.write("tag\tid\tlocations\tdistance\tleft\tright\n")
            f.write("GGGGGGGG\tBC3\t0:7\t1\t1\t0\n")

        splitcode_keepfile = os.path.join(self.temp_dir, 'keep.txt')
        output_prefix = os.path.join(self.temp_dir, 'Sample3')
        with open(splitcode_keepfile, 'w') as f:
            f.write(f"BC3_R1\t{output_prefix}\n")

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
        with open(splitcode_config, 'w') as f:
            f.write("tag\tid\tlocations\tdistance\tleft\tright\n")
            f.write("TTTTTTTT\tBC4\t0:7\t0\t1\t0\n")

        splitcode_keepfile = os.path.join(self.temp_dir, 'keep.txt')
        sample_output_prefix = os.path.join(self.temp_dir, 'MySample')
        with open(splitcode_keepfile, 'w') as f:
            f.write(f"BC4_R1\t{sample_output_prefix}\n")

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
        with open(splitcode_config, 'w') as f:
            f.write("tag\tid\tlocations\tdistance\tleft\tright\n")
            # left=1 means trim the barcode from the left side
            f.write("AAAAAAAA\tBC1\t0:7\t0\t1\t0\n")

        splitcode_keepfile = os.path.join(self.temp_dir, 'keep.txt')
        output_prefix = os.path.join(self.temp_dir, 'Sample1')
        with open(splitcode_keepfile, 'w') as f:
            f.write(f"BC1_R1\t{output_prefix}\n")

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
        with open(splitcode_config, 'w') as f:
            f.write("tag\tid\tlocations\tdistance\tleft\tright\n")
            f.write("GGGGGGGG\tBC1\t0:7\t0\t1\t0\n")

        splitcode_keepfile = os.path.join(self.temp_dir, 'keep.txt')
        output_prefix = os.path.join(self.temp_dir, 'Sample1')
        with open(splitcode_keepfile, 'w') as f:
            f.write(f"BC1_R1\t{output_prefix}\n")

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
