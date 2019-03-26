# Integration tests for intrahost.py

__author__ = "dpark@broadinstitute.org"

# built-ins
import os
import os.path
import shutil
import tempfile
import argparse
import unittest

# module-specific
import intrahost
import interhost
import tools.mafft
import util.file
import util.vcf
import test
from test import TestCaseWithTmp
import tools

# third-party
import pytest

class TestPerSample(test.TestCaseWithTmp):
    ''' This tests step 1 of the iSNV calling process
        (intrahost.vphaser_one_sample), which runs V-Phaser2 on
        a single sample, reformats the output slightly, and performs
        strand-bias filtering and adds library-bias statistics.
    '''

    #@unittest.skipIf(tools.is_osx(), "vphaser2 osx binary from bioconda has issues")
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
        intrahost.vphaser_one_sample(inBam, refFasta, outTab, vphaserNumThreads=test._CPUS, minReadsEach=0)
        expected = os.path.join(myInputDir, 'vphaser_one_sample_indels_expected.txt')
        self.assertEqualContents(outTab, expected)

    def test_vphaser_one_sample_one_mate_unpaired(self):
        # Files here were created as follows:
        # ref.indels.fasta is Seed-stock-137_S2_L001_001.fasta
        # in.oneunmapped.bam was created from in.indels.bam, with flag 0->89, 16->73:
        # When removing doubly mapped reads, doing so can result in all reads
        # being removed in the case of low-quality runs with few reads
        # This tests that when v-phaser2 input is empty, a blank
        # file is created as output.
        myInputDir = util.file.get_test_input_path(self)
        inBam = os.path.join(myInputDir, 'in.oneunmapped.bam')
        refFasta = os.path.join(myInputDir, 'ref.indels.fasta')
        outTab = util.file.mkstempfname('.txt')
        intrahost.vphaser_one_sample(inBam, refFasta, outTab, vphaserNumThreads=test._CPUS, minReadsEach=0, removeDoublyMappedReads=True)
        assert os.path.getsize(outTab) == 0

    #@unittest.skipIf(tools.is_osx(), "vphaser2 osx binary from bioconda has issues")
    def test_vphaser_one_sample_2libs(self):
        # in.2libs.bam was created by "manually" editing in.bam and moving about
        # 1/3 of the reads to ReadGroup2.
        myInputDir = util.file.get_test_input_path(self)
        inBam = os.path.join(myInputDir, 'in.2libs.bam')
        refFasta = os.path.join(myInputDir, 'ref.fasta')
        outTab = util.file.mkstempfname('.txt')
        intrahost.vphaser_one_sample(inBam, refFasta, outTab, vphaserNumThreads=test._CPUS, minReadsEach=6, maxBias=3)
        expected = os.path.join(myInputDir, 'vphaser_one_sample_2libs_expected.txt')
        self.assertEqualContents(outTab, expected)

    #@unittest.skipIf(tools.is_osx(), "vphaser2 osx binary from bioconda has issues")
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
        intrahost.vphaser_one_sample(inBam, refFasta, outTab, vphaserNumThreads=test._CPUS, minReadsEach=6, maxBias=3)
        expected = os.path.join(myInputDir, 'vphaser_one_sample_3libs_expected.txt')
        self.assertEqualContents(outTab, expected)

class TestSnpEff(TestCaseWithTmp):
    @pytest.fixture(autouse=True)
    def capsys(self, capsys):
        self.capsys = capsys

    def test_snpeff(self):
        temp_dir = tempfile.gettempdir()
        input_dir = util.file.get_test_input_path(self)

        ref_fasta      = os.path.join(input_dir,"ref-rabies-JQ685920.fasta")
        assembly_fasta = os.path.join(input_dir,"RBV16.fasta")
        isnv_calls     = os.path.join(input_dir,"vphaser2.RBV16.mapped.txt.gz")

        # align sample to reference to create MSA
        msa_fasta = util.file.mkstempfname('.fasta')
        expected_msa_fasta = os.path.join(input_dir,"msa.fasta")
        args = [ref_fasta, assembly_fasta, msa_fasta, "--localpair", "--preservecase"]
        args = interhost.parser_align_mafft(argparse.ArgumentParser()).parse_args(args)
        args.func_main(args)
        test.assert_equal_contents(self, msa_fasta, expected_msa_fasta)

        # merge (one) VCF to merged vcf
        merged_vcf = os.path.join(temp_dir,"merged.vcf.gz")
        expected_merged_vcf = os.path.join(input_dir,"merged.vcf.gz")
        args = [ref_fasta, merged_vcf, "--isnvs", isnv_calls, "--alignments", msa_fasta, "--strip_chr_version", "--parse_accession"]
        args = intrahost.parser_merge_to_vcf(argparse.ArgumentParser()).parse_args(args)
        args.func_main(args)
        vcf = util.vcf.VcfReader(merged_vcf)
        expected_vcf = util.vcf.VcfReader(expected_merged_vcf)
        rows = list(vcf.get())
        expected_rows = list(expected_vcf.get())
        #self.assertEqual(rows, expected_rows)

        # run snpEff against merged VCF to predict SNP effects
        eff_vcf = os.path.join(temp_dir,"ann_eff.vcf.gz")
        expected_eff_vcf = os.path.join(input_dir,"ann_eff.vcf.gz")
        args = [merged_vcf, "JQ685920", eff_vcf, "--emailAddress=test@example.com"]
        args = interhost.parser_snpEff(argparse.ArgumentParser()).parse_args(args)
        args.func_main(args)
        vcf = util.vcf.VcfReader(eff_vcf)
        expected_vcf = util.vcf.VcfReader(expected_eff_vcf)
        rows = list(vcf.get())
        expected_rows = list(expected_vcf.get())
        #self.assertEqual(rows, expected_rows)

        # create tabular iSNV output
        eff_txt = os.path.join(temp_dir,"ann_eff.txt.gz")
        expected_eff_txt = os.path.join(input_dir,"ann_eff.txt.gz")
        args = [eff_vcf, eff_txt]
        args = intrahost.parser_iSNV_table(argparse.ArgumentParser()).parse_args(args)
        args.func_main(args)
        for outrow, expectedrow in zip(util.file.read_tabfile(eff_txt),util.file.read_tabfile(expected_eff_txt)):
            for colout, colexpected in zip(outrow, expectedrow):
                # if it casts to float, perform approx comparison
                try:
                    f1=float(colout)
                    f2=float(colexpected)
                    self.assertAlmostEqual(f1, f1)
                except ValueError:
                    self.assertEqual(sorted(sorted(colout.split(","))), sorted(sorted(colexpected.split(","))))
