# Unit tests for taxon_filter.py

__author__ = "dpark@broadinstitute.org, irwin@broadinstitute.org," \
    + "hlevitin@broadinstitute.org"

import unittest
import os
import tempfile
import shutil
import filecmp
import subprocess

import argparse

import read_utils
import taxon_filter
import util.file
import util.misc
import tools.last
import tools.bmtagger
import tools.blast
from test import assert_equal_contents, assert_equal_bam_reads, assert_md5_equal_to_line_in_file, TestCaseWithTmp


class TestCommandHelp(unittest.TestCase):

    def test_help_parser_for_each_command(self):
        for cmd_name, parser_fun in taxon_filter.__commands__:
            parser = parser_fun(argparse.ArgumentParser())
            helpstring = parser.format_help()



class TestFilterLastal(TestCaseWithTmp):

    def setUp(self):
        TestCaseWithTmp.setUp(self)
        polio_fasta = os.path.join(
            util.file.get_test_input_path(),
            'TestMetagenomicsViralMix', 'db', 'library', 'Viruses', 'Poliovirus_uid15288', 'NC_002058.ffn'
        )
        dbDir = tempfile.mkdtemp()
        self.lastdb_path = tools.last.Lastdb().build_database(polio_fasta, os.path.join(dbDir, 'NC_002058'))

    def test_filter_lastal_bam_polio(self):
        inBam = os.path.join(util.file.get_test_input_path(), 'TestDepleteHuman', 'expected', 'test-reads.blastn.bam')
        outBam = util.file.mkstempfname('-out-taxfilt.bam')
        args = taxon_filter.parser_filter_lastal_bam(argparse.ArgumentParser()).parse_args([
            inBam, self.lastdb_path, outBam])
        args.func_main(args)
        expectedOut = os.path.join(util.file.get_test_input_path(), 'TestDepleteHuman', 'expected', 'test-reads.taxfilt.imperfect.bam')
        assert_equal_bam_reads(self, outBam, expectedOut)

    def test_lastal_empty_input(self):
        empty_bam = os.path.join(util.file.get_test_input_path(), 'empty.bam')
        outBam = util.file.mkstempfname('-out-taxfilt.bam')
        taxon_filter.filter_lastal_bam(
            empty_bam,
            self.lastdb_path,
            outBam
        )
        assert_equal_bam_reads(self, outBam, empty_bam)

    def test_lastal_empty_output(self):
        empty_bam = os.path.join(util.file.get_test_input_path(), 'empty.bam')
        in_bam = os.path.join(util.file.get_test_input_path(), 'TestDepleteHuman', 'test-reads-human.bam')
        outBam = util.file.mkstempfname('-out-taxfilt.bam')
        taxon_filter.filter_lastal_bam(
            in_bam,
            self.lastdb_path,
            outBam
        )
        assert_equal_bam_reads(self, outBam, empty_bam)


class TestBmtagger(TestCaseWithTmp):
    '''
        How test data was created:
          exported 5kb region of chr6
          created pan-viral fasta file from all NCBI viral accessions
          used wgsim to create simulated reads
    '''

    def setUp(self):
        TestCaseWithTmp.setUp(self)
        self.tempDir = tempfile.mkdtemp()
        myInputDir = os.path.join(util.file.get_test_input_path(), 'TestDepleteHuman')
        ref_fasta = os.path.join(myInputDir, '5kb_human_from_chr6.fasta')
        self.database_prefix_path = os.path.join(self.tempDir, "5kb_human_from_chr6")
        self.bmtooldb_path = tools.bmtagger.BmtoolTool().build_database(ref_fasta, self.database_prefix_path + ".bitmask")
        self.srprismdb_path = tools.bmtagger.SrprismTool().build_database(ref_fasta, self.database_prefix_path + ".srprism")

    def test_deplete_bmtagger_bam(self):
        inBam = os.path.join(util.file.get_test_input_path(), 'TestDepleteHuman', 'test-reads.bam')
        outBam = util.file.mkstempfname('-out.bam')
        args = taxon_filter.parser_deplete_bam_bmtagger(argparse.ArgumentParser()).parse_args([
            inBam, self.database_prefix_path, outBam])
        args.func_main(args)
        expectedOut = os.path.join(util.file.get_test_input_path(), 'TestDepleteHuman', 'expected', 'test-reads.bmtagger.bam')
        assert_equal_bam_reads(self, outBam, expectedOut)

    def test_bmtagger_empty_input(self):
        empty_bam = os.path.join(util.file.get_test_input_path(), 'empty.bam')
        out_bam = util.file.mkstempfname('-out.bam')
        args = taxon_filter.parser_deplete_bam_bmtagger(argparse.ArgumentParser()).parse_args([
            empty_bam, self.database_prefix_path, out_bam])
        args.func_main(args)
        assert_equal_bam_reads(self, out_bam, empty_bam)

    def test_bmtagger_empty_output(self):
        empty_bam = os.path.join(util.file.get_test_input_path(), 'empty.bam')
        in_bam = os.path.join(util.file.get_test_input_path(), 'TestDepleteHuman', 'test-reads-human.bam')
        out_bam = util.file.mkstempfname('-out.bam')
        args = taxon_filter.parser_deplete_bam_bmtagger(argparse.ArgumentParser()).parse_args([
            in_bam, self.database_prefix_path, out_bam])
        args.func_main(args)
        assert_equal_bam_reads(self, out_bam, empty_bam)


class TestBlastnDbBuild(TestCaseWithTmp):

    def test_blastn_db_build(self):
        commonInputDir = util.file.get_test_input_path()
        refFasta = os.path.join(commonInputDir, 'ebola.fasta')

        myInputDir = util.file.get_test_input_path(self)
        tempDir = tempfile.mkdtemp()

        output_prefix = self.__class__.__name__

        args = taxon_filter.parser_blastn_build_db(argparse.ArgumentParser()).parse_args(
            [
                # input fasta
                refFasta,
                # output directory
                tempDir,
                "--outputFilePrefix",
                output_prefix
            ]
        )
        args.func_main(args)

        # nhr=header. nin=index, nsq=sequence
        for ext in [".nhr", ".nsq"]: # ".nin" can change
            assert_equal_contents(
                self, os.path.join(tempDir, output_prefix + ext),
                os.path.join(myInputDir, "expected", output_prefix + ext)
            )


class TestBmtaggerDbBuild(TestCaseWithTmp):

    def test_bmtagger_db_build(self):
        commonInputDir = util.file.get_test_input_path()
        refFasta = os.path.join(commonInputDir, 'ebola.fasta')

        myInputDir = util.file.get_test_input_path(self)
        tempDir = tempfile.mkdtemp()

        output_prefix = self.__class__.__name__

        args = taxon_filter.parser_bmtagger_build_db(argparse.ArgumentParser()).parse_args(
            [
                # input fasta
                refFasta,
                # output directory
                tempDir,
                "--outputFilePrefix",
                output_prefix
            ]
        )
        args.func_main(args)

        for ext in [
            ".bitmask", ".srprism.amp", ".srprism.idx", ".srprism.imp", ".srprism.pmp", ".srprism.rmp",
            ".srprism.ss", ".srprism.ssa", ".srprism.ssd"
        ]:
            assert_equal_contents(
                self, os.path.join(tempDir, output_prefix + ext),
                os.path.join(myInputDir, "expected", output_prefix + ext)
            )

        for ext in [".srprism.map"]:
            assert_md5_equal_to_line_in_file(self, os.path.join(tempDir, output_prefix + ext), os.path.join(myInputDir, "expected", output_prefix + ext+".md5"))


class TestLastalDbBuild(TestCaseWithTmp):

    def test_lastal_db_build(self):
        commonInputDir = util.file.get_test_input_path()
        refFasta = os.path.join(commonInputDir, 'ebola.fasta')

        myInputDir = util.file.get_test_input_path(self)
        tempDir = tempfile.mkdtemp()

        output_prefix = self.__class__.__name__

        args = taxon_filter.parser_lastal_build_db(argparse.ArgumentParser()).parse_args(
            [
                # input fasta
                refFasta,
                # output directory
                tempDir,
                "--outputFilePrefix",
                output_prefix
            ]
        )
        args.func_main(args)

        for ext in [".bck", ".des", ".prj", ".sds", ".ssp", ".suf", ".tis"]:
            assert_equal_contents(
                self, os.path.join(tempDir, output_prefix + ext),
                os.path.join(myInputDir, "expected", output_prefix + ext)
            )


class TestDepleteBlastnBam(TestCaseWithTmp):
    '''
        How test data was created:
        humanChr1Subset.fa has 200 bases from human chr1
        humanChr9Subset.fa has 200 bases from human chr9
        in.fastq "reads" are from humanChr[19]Subset.fa and ebola genome,
        with arbitrary quality scores.
    '''

    def setUp(self):
        TestCaseWithTmp.setUp(self)
        self.tempDir = tempfile.mkdtemp()
        commonInputDir = util.file.get_test_input_path()
        ref_fasta = os.path.join(commonInputDir, 'TestDepleteHuman', '5kb_human_from_chr6.fasta')
        self.database_prefix_path = os.path.join(self.tempDir, "5kb_human_from_chr6")

        # create blast db
        self.blastdb_path = tools.blast.MakeblastdbTool().build_database(ref_fasta, self.database_prefix_path)

    def test_deplete_blastn_bam(self):
        tempDir = tempfile.mkdtemp()
        myInputDir = util.file.get_test_input_path(self)

        # Make blast databases
        makeblastdbPath = tools.blast.MakeblastdbTool().install_and_get_path()
        dbnames = ['humanChr1Subset.fa', 'humanChr9Subset.fa']
        refDbs = []
        for dbname in dbnames:
            refDb = os.path.join(tempDir, dbname)
            os.symlink(os.path.join(myInputDir, dbname), refDb)
            refDbs.append(refDb)
            util.misc.run_and_print([makeblastdbPath, '-dbtype', 'nucl', '-in', refDb], check=True)

        # convert the input fastq's to a bam
        inFastq1 = os.path.join(myInputDir, "in1.fastq")
        inFastq2 = os.path.join(myInputDir, "in2.fastq")
        inBam = os.path.join(tempDir, 'in.bam')
        parser = read_utils.parser_fastq_to_bam(argparse.ArgumentParser())
        args = parser.parse_args(
            [
                inFastq1,
                inFastq2,
                inBam,
                '--sampleName',
                'FreeSample',
                '--JVMmemory',
                '1g',
                '--picardOptions',
                'LIBRARY_NAME=Alexandria',
                'PLATFORM=9.75',
                'SEQUENCING_CENTER=KareemAbdul-Jabbar',
            ]
        )
        args.func_main(args)

        # Run deplete_blastn_bam
        outBam = os.path.join(tempDir, 'out.bam')
        args = taxon_filter.parser_deplete_blastn_bam(argparse.ArgumentParser()).parse_args(
            [inBam, refDbs[0], refDbs[1], outBam, "--chunkSize", "1"]
        )
        args.func_main(args)

        # samtools view for out.sam and compare to expected
        outSam = os.path.join(tempDir, 'out.sam')
        samtools = tools.samtools.SamtoolsTool()
        samtools.view(['-h'], outBam, outSam)

        with open(outSam, "r") as outSamFile:
            for line in outSamFile.readlines():
                print(line)

        # the header field ordering may be different with Java 1.8
        self.assertTrue(
            filecmp.cmp(
                outSam,
                os.path.join(myInputDir, 'expected.sam'),
                shallow=False
            ) or filecmp.cmp(
                outSam,
                os.path.join(myInputDir, 'expected_1_8.sam'),
                shallow=False
            ) or filecmp.cmp(
                outSam,
                os.path.join(myInputDir, 'expected_alt_v1.5.sam'),
                shallow=False
            ) or filecmp.cmp(
                outSam,
                os.path.join(myInputDir, 'expected_1_8_v1.5.sam'),
                shallow=False
            )
        )

    def test_blastn_empty_input(self):
        empty_bam = os.path.join(util.file.get_test_input_path(), 'empty.bam')
        out_bam = util.file.mkstempfname('-out.bam')
        taxon_filter.multi_db_deplete_bam(
            empty_bam,
            [self.blastdb_path],
            taxon_filter.deplete_blastn_bam,
            out_bam
        )
        assert_equal_bam_reads(self, out_bam, empty_bam)

    def test_blastn_empty_output(self):
        in_bam = os.path.join(util.file.get_test_input_path(), 'TestDepleteHuman', 'test-reads-human.bam')
        out_bam = util.file.mkstempfname('-out.bam')
        empty_bam = os.path.join(util.file.get_test_input_path(), 'empty.bam')
        taxon_filter.multi_db_deplete_bam(
            in_bam,
            [self.blastdb_path],
            taxon_filter.deplete_blastn_bam,
            out_bam
        )
        assert_equal_bam_reads(self, out_bam, empty_bam)

