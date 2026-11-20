# Unit tests for bbmap aligner and related tools

__author__ = "ilya@broadinstitute.org"

import unittest
import os.path
import shutil
import viral_ngs.core
import viral_ngs.core.bbmap
import viral_ngs.core.samtools
import pysam
from tests import TestCaseWithTmp, assert_md5_equal_to_line_in_file


class TestToolBBMap(TestCaseWithTmp):

    def setUp(self):
        super(TestToolBBMap, self).setUp()
        self.bbmap = viral_ngs.core.bbmap.BBMapTool()
        self.bbmap.install()
        self.samtools = viral_ngs.core.samtools.SamtoolsTool()

    def test_align(self):
        orig_ref = os.path.join(viral_ngs.core.file.get_test_input_path(), 'ebola.fasta')
        inRef = viral_ngs.core.file.mkstempfname('.fasta')
        shutil.copyfile(orig_ref, inRef)
        reads = os.path.join(viral_ngs.core.file.get_test_input_path(self), 'ebov_reads.bam')
        outBam = viral_ngs.core.file.mkstempfname('.bam')
        # Use conservative memory and single thread for CI compatibility
        # BBMap alignment needs more memory than bbnorm for index building
        self.bbmap.align(inBam=reads, refFasta=inRef, outBam=outBam,
                         threads=1, Xmx='1g')
        self.assertTrue(os.path.isfile(outBam))
        self.assertTrue(os.path.getsize(outBam))

    def test_bbnorm_paired_interleaved(self):
        """Test bbnorm on paired-end interleaved FASTQ file."""
        # Create interleaved FASTQ from TestRmdupUnaligned BAM
        inBam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'TestRmdupUnaligned', 'input.bam')
        inFastq = viral_ngs.core.file.mkstempfname('.fastq')
        outFastq = viral_ngs.core.file.mkstempfname('.fastq')

        # Convert BAM to interleaved FASTQ
        self.samtools.bam2fq(inBam, inFastq)
        self.assertGreater(os.path.getsize(inFastq), 0)

        # Count input reads
        with open(inFastq, 'rt') as f:
            input_lines = sum(1 for _ in f)
        input_reads = input_lines // 4

        # Run bbnorm (auto-detects interleaved format)
        # Use low memory and single thread for CI compatibility
        with viral_ngs.core.file.tmp_dir('_bbnorm_test') as tmpdir:
            self.bbmap.bbnorm(inFastq, outFastq, tmpdir=tmpdir, target=50,
                              threads=1, memory='250m')

        # Verify output exists and has reads
        self.assertTrue(os.path.isfile(outFastq))
        self.assertGreater(os.path.getsize(outFastq), 0)

        # Count output reads
        with open(outFastq, 'rt') as f:
            output_lines = sum(1 for _ in f)
        output_reads = output_lines // 4

        # Output should have reads and be <= input
        self.assertGreater(output_reads, 0)
        self.assertLessEqual(output_reads, input_reads)

    def test_bbnorm_single_end(self):
        """Test bbnorm on single-end FASTQ file."""
        # Create single-end FASTQ from TestPerSample/in.2libs3rgs.bam
        inBam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'TestPerSample', 'in.2libs3rgs.bam')
        inFastq = viral_ngs.core.file.mkstempfname('.fastq')
        outFastq = viral_ngs.core.file.mkstempfname('.fastq')

        # Convert BAM to FASTQ (single-end)
        self.samtools.bam2fq(inBam, inFastq)
        self.assertGreater(os.path.getsize(inFastq), 0)

        # Count input reads
        with open(inFastq, 'rt') as f:
            input_lines = sum(1 for _ in f)
        input_reads = input_lines // 4

        # Run bbnorm (auto-detects single-end format)
        # Use low memory and single thread for CI compatibility
        with viral_ngs.core.file.tmp_dir('_bbnorm_test') as tmpdir:
            self.bbmap.bbnorm(inFastq, outFastq, tmpdir=tmpdir, target=50,
                              threads=1, memory='250m')

        # Verify output exists and has reads
        self.assertTrue(os.path.isfile(outFastq))
        self.assertGreater(os.path.getsize(outFastq), 0)

        # Count output reads
        with open(outFastq, 'rt') as f:
            output_lines = sum(1 for _ in f)
        output_reads = output_lines // 4

        # Output should have reads and be <= input
        self.assertGreater(output_reads, 0)
        self.assertLessEqual(output_reads, input_reads)

    # --- clumpify ---

    def _write_dup_fastqs(self, out1, out2=None, interleaved=None):
        """Write a tiny fixture: 4 read pairs where pair B duplicates pair A exactly.

        Returns the number of read pairs written (4).
        """
        recs = [
            ('readA', 'ACGTACGTACGTACGTACGTACGTACGTACGTAAAA', 'ACGTACGTACGTACGTACGTACGTACGTACGTTTTT'),
            ('readB', 'ACGTACGTACGTACGTACGTACGTACGTACGTAAAA', 'ACGTACGTACGTACGTACGTACGTACGTACGTTTTT'),
            ('readC', 'GGGGCCCCGGGGCCCCGGGGCCCCGGGGCCCCACGT', 'TTTTAAAATTTTAAAATTTTAAAATTTTAAAACGTA'),
            ('readD', 'TACGTACGATCGATCGATCGATCGATCGATCGGGGG', 'CATCATCATCATCATCATCATCATCATCATCAAAAA'),
        ]

        def rec(name, mate, seq):
            return '@{}/{}\n{}\n+\n{}\n'.format(name, mate, seq, 'I' * len(seq))

        if interleaved:
            with open(interleaved, 'wt') as inter:
                for name, seq1, seq2 in recs:
                    inter.write(rec(name, 1, seq1))
                    inter.write(rec(name, 2, seq2))
        if out2 is not None:
            with open(out1, 'wt') as f1, open(out2, 'wt') as f2:
                for name, seq1, seq2 in recs:
                    f1.write(rec(name, 1, seq1))
                    f2.write(rec(name, 2, seq2))
        elif out1 is not None:
            with open(out1, 'wt') as f1:
                for name, seq1, _ in recs:
                    f1.write(rec(name, 1, seq1))
        return len(recs)

    @staticmethod
    def _fastq_read_names(fastq):
        names = []
        with open(fastq, 'rt') as f:
            for i, line in enumerate(f):
                if i % 4 == 0:
                    names.append(line.rstrip('\n')[1:])
        return names

    def test_clumpify_dedupe_removes_exact_duplicates(self):
        """clumpify dedupe collapses an exactly-duplicated pair and leaves the rest."""
        in1 = viral_ngs.core.file.mkstempfname('.1.fastq')
        in2 = viral_ngs.core.file.mkstempfname('.2.fastq')
        self._write_dup_fastqs(in1, in2)

        with viral_ngs.core.file.tmp_dir('_clumpify_test') as tmpdir:
            out1 = os.path.join(tmpdir, 'out.1.fastq')
            out2 = os.path.join(tmpdir, 'out.2.fastq')
            self.bbmap.clumpify(in1, out1, inFastq2=in2, outFastq2=out2,
                                threads=1, memory='250m')

            names1 = self._fastq_read_names(out1)
            names2 = self._fastq_read_names(out2)

        # 4 pairs in, exactly one duplicate pair collapsed -> 3 pairs out
        self.assertEqual(len(names1), 3)
        self.assertEqual(len(names2), 3)
        # mates stay in step
        self.assertEqual([n[:-2] for n in names1], [n[:-2] for n in names2])
        # the unique pairs all survive; exactly one of the A/B twins is dropped
        bases = set(n[:-2] for n in names1)
        self.assertIn('readC', bases)
        self.assertIn('readD', bases)
        self.assertEqual(len({'readA', 'readB'} & bases), 1)

    def test_clumpify_single_end(self):
        """clumpify dedupes single-end input (the case M-Vicuna silently no-op'd)."""
        inFastq = viral_ngs.core.file.mkstempfname('.fastq')
        self._write_dup_fastqs(inFastq)

        with viral_ngs.core.file.tmp_dir('_clumpify_test') as tmpdir:
            outFastq = os.path.join(tmpdir, 'out.fastq')
            self.bbmap.clumpify(inFastq, outFastq, threads=1, memory='250m')
            names = self._fastq_read_names(outFastq)

        self.assertEqual(len(names), 3)

    def test_clumpify_paired_interleaved(self):
        """clumpify on an interleaved FASTQ derived from a real BAM."""
        inBam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'TestRmdupUnaligned', 'input.bam')
        inFastq = viral_ngs.core.file.mkstempfname('.fastq')
        self.samtools.bam2fq(inBam, inFastq)
        self.assertGreater(os.path.getsize(inFastq), 0)

        with open(inFastq, 'rt') as f:
            input_reads = sum(1 for _ in f) // 4

        with viral_ngs.core.file.tmp_dir('_clumpify_test') as tmpdir:
            outFastq = os.path.join(tmpdir, 'out.fastq')
            self.bbmap.clumpify(inFastq, outFastq, interleaved=True,
                                threads=1, memory='250m')
            self.assertTrue(os.path.isfile(outFastq))
            with open(outFastq, 'rt') as f:
                output_reads = sum(1 for _ in f) // 4

        self.assertGreater(output_reads, 0)
        self.assertLess(output_reads, input_reads)

    def test_clumpify_paired_split(self):
        """clumpify with separate in1/in2 keeps mate files the same length."""
        inBam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'TestRmdupUnaligned', 'input.bam')
        in1 = viral_ngs.core.file.mkstempfname('.1.fastq')
        in2 = viral_ngs.core.file.mkstempfname('.2.fastq')
        self.samtools.bam2fq(inBam, in1, in2)

        with viral_ngs.core.file.tmp_dir('_clumpify_test') as tmpdir:
            out1 = os.path.join(tmpdir, 'out.1.fastq')
            out2 = os.path.join(tmpdir, 'out.2.fastq')
            self.bbmap.clumpify(in1, out1, inFastq2=in2, outFastq2=out2,
                                threads=1, memory='250m')
            n1 = len(self._fastq_read_names(out1))
            n2 = len(self._fastq_read_names(out2))

        self.assertGreater(n1, 0)
        self.assertEqual(n1, n2)

    def test_clumpify_subs_and_containment(self):
        """The sensitivity knobs reach the command line and still produce valid output."""
        inBam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'TestRmdupUnaligned', 'input.bam')
        inFastq = viral_ngs.core.file.mkstempfname('.fastq')
        self.samtools.bam2fq(inBam, inFastq)

        with viral_ngs.core.file.tmp_dir('_clumpify_test') as tmpdir:
            strict = os.path.join(tmpdir, 'strict.fastq')
            loose = os.path.join(tmpdir, 'loose.fastq')
            self.bbmap.clumpify(inFastq, strict, interleaved=True, subs=0,
                                threads=1, memory='250m')
            self.bbmap.clumpify(inFastq, loose, interleaved=True, subs=5,
                                containment=True, threads=1, memory='250m')
            with open(strict, 'rt') as f:
                strict_reads = sum(1 for _ in f) // 4
            with open(loose, 'rt') as f:
                loose_reads = sum(1 for _ in f) // 4

        # a looser substitution budget can only remove at least as many reads
        self.assertGreater(strict_reads, 0)
        self.assertGreater(loose_reads, 0)
        self.assertLessEqual(loose_reads, strict_reads)
