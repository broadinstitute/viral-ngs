# Unit tests for intrahost.py

__author__ = "dpark@broadinstitute.org"

# built-ins
from collections import OrderedDict
import os
import os.path
import shutil
import tempfile
import itertools
import argparse
import unittest

# third-party
import Bio
import Bio.SeqRecord
import Bio.Seq

# module-specific
import intrahost
import util.file
import util.vcf
import test
from intrahost import AlleleFieldParser
import interhost
import tools.mafft

class TestPerSample(test.TestCaseWithTmp):
    ''' This tests step 1 of the iSNV calling process
        (intrahost.vphaser_one_sample), which runs V-Phaser2 on
        a single sample, reformats the output slightly, and performs
        strand-bias filtering and adds library-bias statistics.
    '''

    def test_vphaser_one_sample(self):
        # Files here were created as follows:
        # - in.bam was copied from input directory for TestVPhaser2; see notes
        #   there on how it was created.
        # - ref.fasta was created by making two identical chromosomes, chr1
        #   and chr2, with the sequence from West Nile virus isolate
        #   WNV-1/US/BID-V7821/2011. That genome is a guess for the reference
        #   of the V-Phaser 2 test file because BLAST matches the reads to
        #   West Nile virus and that isolate has the size reported in the bam file.
        #   Note that ref.fasta is not exactly the consensus of the reads in in.bam;
        #   for example, pos 660 is C in ref.fasta, but more reads have T there
        #   than C in in.bam. So we are actually testing the case that
        #   V-Phaser 2 consensus != our consensus.
        myInputDir = util.file.get_test_input_path(self)
        inBam = os.path.join(myInputDir, 'in.bam')
        refFasta = os.path.join(myInputDir, 'ref.fasta')
        outTab = util.file.mkstempfname('.txt')
        intrahost.vphaser_one_sample(inBam, refFasta, outTab, vphaserNumThreads=4, minReadsEach=6, maxBias=3)
        expected = os.path.join(myInputDir, 'vphaser_one_sample_expected.txt')
        self.assertEqualContents(outTab, expected)

    def test_vphaser_one_sample_indels(self):
        # Files here were created as follows:
        # ref.indels.fasta is Seed-stock-137_S2_L001_001.fasta
        # in.indels.bam was created from Seed-stock-137_S2_L001_001.mapped.bam
        #     as follows:
        # Created two .sam files using samtools view, restricting to ranges
        # 6811-7011 and 13081-13281, respectively. Paired reads were removed
        # from those files by throwing out the second occurence of any read name
        # and anding the flag fields with 16. Then, a random 90% of reads were
        # removed, except that any reads containing the indel variants at
        # 6911 and 13181 were kept. Then the resulting 2 files were combined.
        myInputDir = util.file.get_test_input_path(self)
        inBam = os.path.join(myInputDir, 'in.indels.bam')
        refFasta = os.path.join(myInputDir, 'ref.indels.fasta')
        outTab = util.file.mkstempfname('.txt')
        intrahost.vphaser_one_sample(inBam, refFasta, outTab, vphaserNumThreads=4, minReadsEach=0)
        expected = os.path.join(myInputDir, 'vphaser_one_sample_indels_expected.txt')
        self.assertEqualContents(outTab, expected)

    def test_vphaser_one_sample_2libs(self):
        # in.2libs.bam was created by "manually" editing in.bam and moving about
        # 1/3 of the reads to ReadGroup2.
        myInputDir = util.file.get_test_input_path(self)
        inBam = os.path.join(myInputDir, 'in.2libs.bam')
        refFasta = os.path.join(myInputDir, 'ref.fasta')
        outTab = util.file.mkstempfname('.txt')
        intrahost.vphaser_one_sample(inBam, refFasta, outTab, vphaserNumThreads=4, minReadsEach=6, maxBias=3)
        expected = os.path.join(myInputDir, 'vphaser_one_sample_2libs_expected.txt')
        self.assertEqualContents(outTab, expected)

    def test_vphaser_one_sample_3libs_and_chi2(self):
        # In addition to testing that we can handle 3 libraries, this is testing
        #    the chi2_contingency approximation to fisher_exact. The 4th, 5th,
        #    and 6th rows have large enough minor allele count that their
        #    p-values are calculated using the chi2 approximation. The other
        #    rows are testing the 2 x 3 case of fisher_exact.
        # in.3libs.bam was created by "manually" editing in.2libs.bam and moving
        # about 1/2 of the reads in ReadGroup2 to ReadGroup3.
        myInputDir = util.file.get_test_input_path(self)
        inBam = os.path.join(myInputDir, 'in.3libs.bam')
        refFasta = os.path.join(myInputDir, 'ref.fasta')
        outTab = util.file.mkstempfname('.txt')
        intrahost.vphaser_one_sample(inBam, refFasta, outTab, vphaserNumThreads=4, minReadsEach=6, maxBias=3)
        expected = os.path.join(myInputDir, 'vphaser_one_sample_3libs_expected.txt')
        self.assertEqualContents(outTab, expected)
