# Unit tests for Novoalign aligner

__author__ = "dpark@broadinstitute.org"

import unittest
import os.path
import shutil
import viral_ngs.core
import viral_ngs.core.novoalign
import viral_ngs.core.samtools
import pysam
from test import TestCaseWithTmp, assert_md5_equal_to_line_in_file


class TestToolNovoalign(TestCaseWithTmp):

    def setUp(self):
        super(TestToolNovoalign, self).setUp()
        self.novoalign = viral_ngs.core.novoalign.NovoalignTool()
        self.novoalign.install()
        self.samtools = viral_ngs.core.samtools.SamtoolsTool()

    def test_index(self):
        orig_ref = os.path.join(viral_ngs.core.file.get_test_input_path(), 'ebola.fasta')
        inRef = viral_ngs.core.file.mkstempfname('.fasta')
        shutil.copyfile(orig_ref, inRef)
        self.novoalign.index_fasta(inRef)
        outfile = inRef[:-6] + '.nix'

        self.assertTrue(os.path.isfile(outfile))
        self.assertTrue(os.path.getsize(outfile))

    def test_align(self):
        orig_ref = os.path.join(viral_ngs.core.file.get_test_input_path(), 'ebola.fasta')
        inRef = viral_ngs.core.file.mkstempfname('.fasta')
        shutil.copyfile(orig_ref, inRef)
        self.novoalign.index_fasta(inRef)
        reads = os.path.join(viral_ngs.core.file.get_test_input_path(self), 'ebov_reads.bam')
        outBam = viral_ngs.core.file.mkstempfname('.bam')
        self.novoalign.execute(reads, inRef, outBam)
        self.assertTrue(os.path.isfile(outBam))
        self.assertTrue(os.path.getsize(outBam))
        self.assertTrue(os.path.isfile(outBam[:-1] + 'i'))

    def test_align_filter(self):
        orig_ref = os.path.join(viral_ngs.core.file.get_test_input_path(), 'ebola.fasta')
        inRef = viral_ngs.core.file.mkstempfname('.fasta')
        shutil.copyfile(orig_ref, inRef)
        self.novoalign.index_fasta(inRef)
        reads = os.path.join(viral_ngs.core.file.get_test_input_path(self), 'ebov_reads.bam')
        outBam = viral_ngs.core.file.mkstempfname('.bam')
        self.novoalign.execute(reads, inRef, outBam, min_qual=1)
        self.assertTrue(os.path.isfile(outBam))
        self.assertTrue(os.path.getsize(outBam))
        self.assertTrue(os.path.isfile(outBam[:-1] + 'i'))

    def test_multi_read_groups(self):
        orig_ref = os.path.join(viral_ngs.core.file.get_test_input_path(), 'G5012.3.fasta')
        inRef = viral_ngs.core.file.mkstempfname('.fasta')
        shutil.copyfile(orig_ref, inRef)
        self.novoalign.index_fasta(inRef)

        # align with Novoalign (BAM input, BAM output)
        reads = os.path.join(viral_ngs.core.file.get_test_input_path(), 'G5012.3.subset.bam')
        outBam = viral_ngs.core.file.mkstempfname('.bam')
        self.novoalign.execute(reads, inRef, outBam)
        self.assertTrue(os.path.isfile(outBam))
        self.assertTrue(os.path.getsize(outBam))
        self.assertTrue(os.path.isfile(outBam[:-1] + 'i'))
        sam_in = viral_ngs.core.file.mkstempfname('.in.sam')
        sam_out = viral_ngs.core.file.mkstempfname('.out.sam')
        self.samtools.view([], reads, sam_in)
        self.samtools.view([], outBam, sam_out)

        # assert that all reads are present in output
        self.assertEqual(self.samtools.count(reads), self.samtools.count(outBam))

        # assert that all read groups are described exactly the same
        # in the output header (not necessarily same order, but same content)
        orig_rgs = self.samtools.getReadGroups(reads)
        new_rgs = self.samtools.getReadGroups(outBam)
        self.assertTrue(len(orig_rgs) > 1)
        self.assertTrue(len(new_rgs) > 1)
        for rgid in new_rgs.keys():
            self.assertIn(rgid, orig_rgs)
            orig_rg = orig_rgs[rgid]
            new_rg = new_rgs[rgid]
            if 'DT' in orig_rg and 'DT' in new_rg:
                # don't bother comparing the date because Novoalign seems
                # to enjoy doing senseless time zone conversions on the
                # ISO date string
                del orig_rg['DT']
                del new_rg['DT']
            self.assertEqual(orig_rgs[rgid], new_rgs[rgid])

        # if any RGs are missing, it should be because they were never really there
        for rgid in orig_rgs.keys():
            if rgid not in new_rgs:
                self.assertEqual(0, self.samtools.count(reads, ['-r', rgid]))

        # assert that all reads retained their original RG assignments
        read_to_rg = {}
        read_seen = set()
        with open(sam_in, 'rt') as inf:
            for read in inf:
                read = read.rstrip('\n').split('\t')
                tags = [t[5:] for t in read[11:] if t.startswith('RG:Z:')]
                self.assertTrue(len(tags) == 1)
                read_to_rg[read[0]] = tags[0]
        with open(sam_out, 'rt') as inf:
            for read in inf:
                read = read.rstrip('\n').split('\t')
                tags = [t[5:] for t in read[11:] if t.startswith('RG:Z:')]
                self.assertTrue(len(tags) == 1)
                self.assertIn(read[0], read_to_rg)
                self.assertEqual(tags[0], read_to_rg[read[0]])
                read_seen.add(read[0])
        self.assertEqual(len(read_seen), len(read_to_rg))

        # clean up
        for fn in (sam_in, sam_out, outBam, inRef):
            os.unlink(fn)

    def test_multi_read_groups_filter(self):
        orig_ref = os.path.join(viral_ngs.core.file.get_test_input_path(), 'G5012.3.fasta')
        inRef = viral_ngs.core.file.mkstempfname('.fasta')
        shutil.copyfile(orig_ref, inRef)
        self.novoalign.index_fasta(inRef)

        # align with Novoalign (BAM input, BAM output)
        reads = os.path.join(viral_ngs.core.file.get_test_input_path(), 'G5012.3.subset.bam')
        outBam = viral_ngs.core.file.mkstempfname('.bam')
        self.novoalign.execute(reads, inRef, outBam, min_qual=1)
        self.assertTrue(os.path.isfile(outBam))
        self.assertTrue(os.path.getsize(outBam))
        self.assertTrue(os.path.isfile(outBam[:-1] + 'i'))
