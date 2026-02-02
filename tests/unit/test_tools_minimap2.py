# Unit tests for minimap2 aligner

__author__ = "dpark@broadinstitute.org"

import os
import os.path
import shutil
import tempfile
import subprocess
import util.file
import tools.minimap2
import tools.samtools
from test import TestCaseWithTmp

class TestToolMinimap2(TestCaseWithTmp):

    def setUp(self):
        super(TestToolMinimap2, self).setUp()
        self.mm2 = tools.minimap2.Minimap2()
        self.samtools = tools.samtools.SamtoolsTool()
        self.tempDir = tempfile.mkdtemp()
        root_input_dir = util.file.get_test_input_path()
        self.ref_fasta = os.path.join(root_input_dir, '5kb_human_from_chr6.fasta')

    def test_human_bam(self):
        in_bam = os.path.join(util.file.get_test_input_path(), "TestDepleteHuman", 'test-reads-human.bam')
        with util.file.tempfname('.bam') as outfile:
            self.mm2.align_bam(in_bam, self.ref_fasta, outfile, options=['-x', 'sr'])
            self.assertEqual(self.samtools.count(outfile), 20)

    def test_ebola_bam(self):
        in_bam = os.path.join(util.file.get_test_input_path(), 'G5012.3.testreads.bam')
        ref_fasta = os.path.join(util.file.get_test_input_path(), 'ebov-makona.fasta')
        with util.file.tempfname('.bam') as outfile:
            self.mm2.align_bam(in_bam, ref_fasta, outfile)
            self.assertEqual(self.samtools.count(outfile), 18843)

    def test_corrupt_bam(self):
        # pipe.poll() should raise an exception
        in_bam = os.path.join(util.file.get_test_input_path(), 'broken.bam')
        with util.file.tempfname('.bam') as outfile:
            self.assertRaises(subprocess.CalledProcessError, self.mm2.align_bam, in_bam, self.ref_fasta, outfile)


class TestMinimap2Idxstats(TestCaseWithTmp):
    """Tests for the Minimap2.idxstats() method - PAF-based read counting."""

    def setUp(self):
        super(TestMinimap2Idxstats, self).setUp()
        self.mm2 = tools.minimap2.Minimap2()
        self.samtools = tools.samtools.SamtoolsTool()
        self.tempDir = tempfile.mkdtemp()
        self.root_input_dir = util.file.get_test_input_path()

    def test_idxstats_basic(self):
        """Basic functionality: single reference, verify idxstats format and counts."""
        in_bam = os.path.join(self.root_input_dir, 'G5012.3.testreads.bam')
        ref_fasta = os.path.join(self.root_input_dir, 'ebov-makona.fasta')

        with util.file.tempfname('.idxstats.txt') as out_stats:
            self.mm2.idxstats(in_bam, ref_fasta, out_stats)

            # Verify file exists and has content
            self.assertTrue(os.path.isfile(out_stats))
            self.assertGreater(os.path.getsize(out_stats), 0)

            # Parse and validate idxstats format
            with open(out_stats, 'rt') as f:
                lines = [line.strip().split('\t') for line in f if line.strip()]

            # Should have at least one line (the reference sequence)
            self.assertGreater(len(lines), 0)

            # Each line should have 4 tab-separated columns
            for line in lines:
                self.assertEqual(len(line), 4, f"Expected 4 columns, got {len(line)}: {line}")

            # First line should be the reference: name, length, mapped, unmapped
            ref_name, ref_len, mapped_count, unmapped_count = lines[0]
            self.assertEqual(ref_name, 'KJ660346.2')  # ebov-makona reference name
            self.assertGreater(int(ref_len), 0)
            self.assertGreater(int(mapped_count), 18000)  # Should be close to 18843
            self.assertEqual(int(unmapped_count), 0)  # Always 0 in our implementation

    def test_idxstats_with_readlist(self):
        """Test optional read list output."""
        in_bam = os.path.join(self.root_input_dir, 'G5012.3.testreads.bam')
        ref_fasta = os.path.join(self.root_input_dir, 'ebov-makona.fasta')

        with util.file.tempfname('.idxstats.txt') as out_stats:
            with util.file.tempfname('.readlist.txt') as out_readlist:
                self.mm2.idxstats(in_bam, ref_fasta, out_stats, outReadlist=out_readlist)

                # Verify read list file is created
                self.assertTrue(os.path.isfile(out_readlist))
                self.assertGreater(os.path.getsize(out_readlist), 0)

                # Parse read list
                with open(out_readlist, 'rt') as f:
                    read_ids = [line.strip() for line in f if line.strip()]

                # Should have reads
                self.assertGreater(len(read_ids), 0)

                # All read IDs should be unique
                self.assertEqual(len(read_ids), len(set(read_ids)),
                                 "Read list should contain unique read IDs only")

                # Verify stats count (counts R1 + R2 separately, like samtools idxstats)
                with open(out_stats, 'rt') as f:
                    first_line = f.readline().strip().split('\t')
                    mapped_count = int(first_line[2])
                self.assertGreater(mapped_count, 18000)

                # Readlist contains unique read pairs (~half of mapped count for paired-end data)
                self.assertGreater(len(read_ids), 9000)

    def test_idxstats_no_readlist(self):
        """Test that readlist output can be skipped (outReadlist=None)."""
        in_bam = os.path.join(self.root_input_dir, 'G5012.3.testreads.bam')
        ref_fasta = os.path.join(self.root_input_dir, 'ebov-makona.fasta')

        with util.file.tempfname('.idxstats.txt') as out_stats:
            # Call without readlist
            self.mm2.idxstats(in_bam, ref_fasta, out_stats, outReadlist=None)

            # idxstats should still work
            self.assertTrue(os.path.isfile(out_stats))
            with open(out_stats, 'rt') as f:
                lines = [line.strip().split('\t') for line in f if line.strip()]
            self.assertGreater(len(lines), 0)

    def test_idxstats_empty_input(self):
        """Edge case: empty BAM input should produce zero counts."""
        in_bam = os.path.join(self.root_input_dir, 'empty.bam')
        ref_fasta = os.path.join(self.root_input_dir, 'ebov-makona.fasta')

        with util.file.tempfname('.idxstats.txt') as out_stats:
            self.mm2.idxstats(in_bam, ref_fasta, out_stats)

            # Should have output with zero counts
            self.assertTrue(os.path.isfile(out_stats))
            with open(out_stats, 'rt') as f:
                lines = [line.strip().split('\t') for line in f if line.strip()]

            # Reference should be listed but with 0 mapped count
            self.assertEqual(len(lines), 1)
            ref_name, ref_len, mapped_count, unmapped_count = lines[0]
            self.assertEqual(int(mapped_count), 0)

    def test_idxstats_empty_reference(self):
        """Edge case: empty FASTA reference should handle gracefully."""
        in_bam = os.path.join(self.root_input_dir, 'G5012.3.testreads.bam')
        ref_fasta = os.path.join(self.root_input_dir, 'empty.fasta')

        with util.file.tempfname('.idxstats.txt') as out_stats:
            self.mm2.idxstats(in_bam, ref_fasta, out_stats)

            # Should produce empty output (no references = no lines)
            self.assertTrue(os.path.isfile(out_stats))
            with open(out_stats, 'rt') as f:
                lines = [line.strip() for line in f if line.strip()]
            self.assertEqual(len(lines), 0)

    def test_idxstats_multi_reference(self):
        """Multi-reference: counts should be correctly separated per reference sequence."""
        in_bam = os.path.join(self.root_input_dir, 'TestMinimap2Idxstats', 'multi-viral-reads.bam')
        ref_fasta = os.path.join(self.root_input_dir, 'TestMinimap2Idxstats', 'multi-viral-refs.fasta')

        with util.file.tempfname('.idxstats.txt') as out_stats:
            self.mm2.idxstats(in_bam, ref_fasta, out_stats)

            # Parse output
            with open(out_stats, 'rt') as f:
                lines = [line.strip().split('\t') for line in f if line.strip()]

            # Should have 4 references (HCV, Tomato mosaic, Poliovirus, Coxsackievirus)
            self.assertEqual(len(lines), 4)

            # Validate format for each line
            ref_names = []
            total_mapped = 0
            for line in lines:
                self.assertEqual(len(line), 4)
                ref_name, ref_len, mapped_count, unmapped_count = line
                ref_names.append(ref_name)
                self.assertGreater(int(ref_len), 0)
                self.assertGreaterEqual(int(mapped_count), 0)
                self.assertEqual(int(unmapped_count), 0)
                total_mapped += int(mapped_count)

            # Reference names should include our expected sequences
            expected_refs = ['EU255973.1', 'KJ207374.1', 'AY184220.1', 'KM609478.1']
            for expected in expected_refs:
                self.assertIn(expected, ref_names)
