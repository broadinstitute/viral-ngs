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
import classify.last
import classify.bmtagger
import classify.blast
from test import assert_equal_contents, assert_equal_bam_reads, assert_md5_equal_to_line_in_file, TestCaseWithTmp


class TestCommandHelp(unittest.TestCase):

    def test_help_parser_for_each_command(self):
        for cmd_name, parser_fun in taxon_filter.__commands__:
            parser = parser_fun(argparse.ArgumentParser())
            helpstring = parser.format_help()



class TestFilterLastal(TestCaseWithTmp):

    def setUp(self):
        TestCaseWithTmp.setUp(self)
        self.polio_fasta = os.path.join(
            util.file.get_test_input_path(),
            'TestMetagenomicsViralMix', 'db', 'library', 'Viruses', 'Enterovirus_C', 'GCF_000861165.1_ViralProj15288_genomic.fna'
        )
        dbDir = tempfile.mkdtemp()
        self.lastdb_path = classify.last.Lastdb().build_database(self.polio_fasta, os.path.join(dbDir, 'NC_002058'))

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

    def test_lastal_unbuilt_db(self):
        empty_bam = os.path.join(util.file.get_test_input_path(), 'empty.bam')
        in_bam = os.path.join(util.file.get_test_input_path(), 'TestDepleteHuman', 'test-reads-human.bam')
        outBam = util.file.mkstempfname('-out-taxfilt.bam')
        taxon_filter.filter_lastal_bam(
            in_bam,
            self.polio_fasta,
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
        ref_fasta = os.path.join(util.file.get_test_input_path(), '5kb_human_from_chr6.fasta')
        self.database_prefix_path = os.path.join(self.tempDir, "5kb_human_from_chr6")
        taxon_filter.bmtagger_build_db(ref_fasta, self.tempDir, "5kb_human_from_chr6", word_size=8)

    def test_deplete_bmtagger_bam(self):
        inBam = os.path.join(util.file.get_test_input_path(), 'TestDepleteHuman', 'test-reads.bam')
        outBam = util.file.mkstempfname('-out.bam')
        args = taxon_filter.parser_deplete_bam_bmtagger(argparse.ArgumentParser()).parse_args([
            inBam, self.database_prefix_path, outBam, '--srprismMemory', '1500'])
        args.func_main(args)
        expectedOut = os.path.join(util.file.get_test_input_path(), 'TestDepleteHuman', 'expected', 'test-reads.bmtagger.bam')
        assert_equal_bam_reads(self, outBam, expectedOut)

    @unittest.skip("too slow for real word size of 18bp")
    def test_deplete_bmtagger_fasta_db(self):
        inBam = os.path.join(util.file.get_test_input_path(), 'TestDepleteHuman', 'test-reads.bam')
        ref_fasta = os.path.join(util.file.get_test_input_path(), '5kb_human_from_chr6.fasta')
        outBam = util.file.mkstempfname('-out.bam')
        args = taxon_filter.parser_deplete_bam_bmtagger(argparse.ArgumentParser()).parse_args([
            inBam, ref_fasta, outBam, '--srprismMemory', '1500'])
        args.func_main(args)
        expectedOut = os.path.join(util.file.get_test_input_path(), 'TestDepleteHuman', 'expected', 'test-reads.bmtagger.bam')
        assert_equal_bam_reads(self, outBam, expectedOut)

    def test_deplete_bmtagger_tar_db(self):
        inBam = os.path.join(util.file.get_test_input_path(), 'TestDepleteHuman', 'test-reads.bam')
        outBam = util.file.mkstempfname('-out.bam')
        tar_db_tgz = util.file.mkstempfname('.db.tar.gz')
        cmd = ['tar', '-C', os.path.dirname(self.database_prefix_path), '-cvzf', tar_db_tgz, '.']
        subprocess.check_call(cmd)
        args = taxon_filter.parser_deplete_bam_bmtagger(argparse.ArgumentParser()).parse_args([
            inBam, tar_db_tgz, outBam, '--srprismMemory', '1500'])
        args.func_main(args)
        expectedOut = os.path.join(util.file.get_test_input_path(), 'TestDepleteHuman', 'expected', 'test-reads.bmtagger.bam')
        assert_equal_bam_reads(self, outBam, expectedOut)
        os.unlink(tar_db_tgz)

    def test_bmtagger_empty_input(self):
        empty_bam = os.path.join(util.file.get_test_input_path(), 'empty.bam')
        out_bam = util.file.mkstempfname('-out.bam')
        args = taxon_filter.parser_deplete_bam_bmtagger(argparse.ArgumentParser()).parse_args([
            empty_bam, self.database_prefix_path, out_bam, '--srprismMemory', '1500'])
        args.func_main(args)
        assert_equal_bam_reads(self, out_bam, empty_bam)

    def test_bmtagger_empty_output(self):
        empty_bam = os.path.join(util.file.get_test_input_path(), 'empty.bam')
        in_bam = os.path.join(util.file.get_test_input_path(), 'TestDepleteHuman', 'test-reads-human.bam')
        out_bam = util.file.mkstempfname('-out.bam')
        args = taxon_filter.parser_deplete_bam_bmtagger(argparse.ArgumentParser()).parse_args([
            in_bam, self.database_prefix_path, out_bam, '--srprismMemory', '1500'])
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

    def test_blastn_db_build_gz(self):
        commonInputDir = util.file.get_test_input_path()
        refFasta = os.path.join(commonInputDir, 'ebola.fasta.gz')

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

        commonInputDir = util.file.get_test_input_path()

        refFasta = os.path.join(commonInputDir, 'ebola.fasta.lz4')
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


class TestBmtaggerDbBuild(TestCaseWithTmp):

    def test_bmtagger_db_build(self):
        os.environ.pop('TMPDIR', None)
        util.file.set_tmp_dir(None)
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
                output_prefix,
                "--word_size",
                "8",
            ]
        )
        args.func_main(args)

        for ext in [
            ".bitmask", ".srprism.amp", ".srprism.imp", ".srprism.pmp", ".srprism.rmp",
            ".srprism.ssa", ".srprism.ssd"
        ]:
            assert_equal_contents(
                self, os.path.join(tempDir, output_prefix + ext),
                os.path.join(myInputDir, "expected", output_prefix + ext)
            )

        for ext in [".srprism.map", ".srprism.idx", ".srprism.ss"]:
            assert_md5_equal_to_line_in_file(self, os.path.join(tempDir, output_prefix + ext), os.path.join(myInputDir, "expected", output_prefix + ext+".md5"))

    def test_bmtagger_db_build_gz(self):
        commonInputDir = util.file.get_test_input_path()
        refFasta = os.path.join(commonInputDir, 'ebola.fasta.gz')
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
                output_prefix,
                "--word_size",
                "8",
            ]
        )
        args.func_main(args)
        refFasta = os.path.join(commonInputDir, 'ebola.fasta.lz4')
        args = taxon_filter.parser_bmtagger_build_db(argparse.ArgumentParser()).parse_args(
            [
                # input fasta
                refFasta,
                # output directory
                tempDir,
                "--outputFilePrefix",
                output_prefix,
                "--word_size",
                "8",
            ]
        )
        args.func_main(args)


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

        for ext in [".bck", ".des", ".prj", ".sds", ".ssp", ".tis"]:
            assert_equal_contents(
                self, os.path.join(tempDir, output_prefix + ext),
                os.path.join(myInputDir, "expected", output_prefix + ext)
            )

        for ext in [".suf"]:
            assert_md5_equal_to_line_in_file(self, os.path.join(tempDir, output_prefix + ext), os.path.join(myInputDir, "expected", output_prefix + ext+".md5"))

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
        ref_fasta = os.path.join(commonInputDir, '5kb_human_from_chr6.fasta')
        self.database_prefix_path = os.path.join(self.tempDir, "5kb_human_from_chr6")

        # create blast db
        self.blastdb_path = classify.blast.MakeblastdbTool().build_database(ref_fasta, self.database_prefix_path)

        # create multiple dbs
        self.blastdbs_multi = []
        for db in ['humanChr1Subset.fa', 'humanChr9Subset.fa']:
            dbPath = classify.blast.MakeblastdbTool().build_database(
                os.path.join(util.file.get_test_input_path(self), db),
                os.path.join(self.tempDir, db[:-3]))
            self.blastdbs_multi.append(dbPath)

        # tar one db, but not the other
        tar_db_tgz = util.file.mkstempfname('-humanChr9Subset.blastn.db.tar.gz')
        cmd = ['tar', '-C', self.tempDir, '-cvzf', tar_db_tgz]
        for ext in ('nhr', 'nin', 'nsq'):
            cmd.append('humanChr9Subset.'+ext)
        subprocess.check_call(cmd)
        self.blastdbs_multi[1] = tar_db_tgz
        for ext in ('nhr', 'nin', 'nsq'):
            os.unlink(os.path.join(self.tempDir, 'humanChr9Subset.'+ext))

    def test_deplete_blastn_bam(self):
        tempDir = tempfile.mkdtemp()
        myInputDir = util.file.get_test_input_path(self)

        # Run deplete_blastn_bam
        inBam = os.path.join(myInputDir, 'in.bam')
        outBam = os.path.join(tempDir, 'out.bam')
        args = taxon_filter.parser_deplete_blastn_bam(argparse.ArgumentParser()).parse_args(
            [inBam] + self.blastdbs_multi + [outBam, "--chunkSize", "0"]
        )
        args.func_main(args)

        # samtools view for out.sam and compare to expected
        outSam = os.path.join(tempDir, 'out.sam')
        samtools = tools.samtools.SamtoolsTool()
        samtools.view(['-h'], outBam, outSam)

        #with open(outSam, "r") as outSamFile:
        #    for line in outSamFile.readlines():
        #        print(line)

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

    def test_deplete_blastn_bam_chunked(self):
        tempDir = tempfile.mkdtemp()
        myInputDir = util.file.get_test_input_path(self)

        # Run deplete_blastn_bam
        inBam = os.path.join(myInputDir, 'in.bam')
        outBam = os.path.join(tempDir, 'out.bam')
        args = taxon_filter.parser_deplete_blastn_bam(argparse.ArgumentParser()).parse_args(
            [inBam] + self.blastdbs_multi + [outBam, "--chunkSize", "1"]
        )
        args.func_main(args)

        # samtools view for out.sam and compare to expected
        outSam = os.path.join(tempDir, 'out.sam')
        samtools = tools.samtools.SamtoolsTool()
        samtools.view(['-h'], outBam, outSam)

        #with open(outSam, "r") as outSamFile:
        #    for line in outSamFile.readlines():
        #        print(line)

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
        self.assertEqual(0, tools.samtools.SamtoolsTool().count(out_bam))

    def test_blastn_empty_output(self):
        in_bam = os.path.join(util.file.get_test_input_path(), 'TestDepleteHuman', 'test-reads-human.bam')
        out_bam = util.file.mkstempfname('-out.bam')
        taxon_filter.multi_db_deplete_bam(
            in_bam,
            [self.blastdb_path],
            taxon_filter.deplete_blastn_bam,
            out_bam
        )
        self.assertEqual(0, tools.samtools.SamtoolsTool().count(out_bam))
