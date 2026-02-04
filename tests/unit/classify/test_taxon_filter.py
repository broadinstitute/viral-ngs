# Unit tests for taxon_filter.py

__author__ = "dpark@broadinstitute.org, irwin@broadinstitute.org," \
    + "hlevitin@broadinstitute.org"

import unittest
import glob
import os, os.path
import platform
import tempfile
import shutil
import filecmp
import subprocess

import argparse

from viral_ngs import read_utils
from viral_ngs import taxon_filter
from viral_ngs.core import file as util_file
from viral_ngs.core import misc as util_misc
from viral_ngs.core import samtools

from viral_ngs.classify import last
from viral_ngs.classify import bmtagger
from viral_ngs.classify import blast

from tests import assert_equal_contents, assert_equal_bam_reads, assert_md5_equal_to_line_in_file, TestCaseWithTmp

# Skip bmtagger tests on ARM platforms - bmtagger is x86-only (no ARM64 builds in bioconda)
IS_ARM = platform.machine() in ('arm64', 'aarch64')
SKIP_BMTAGGER_REASON = "bmtagger requires x86-only bioconda package (not available on ARM)"


class TestCommandHelp(unittest.TestCase):

    def test_help_parser_for_each_command(self):
        for cmd_name, parser_fun in taxon_filter.__commands__:
            parser = parser_fun(argparse.ArgumentParser())
            helpstring = parser.format_help()



class TestFilterLastal(TestCaseWithTmp):

    def setUp(self):
        TestCaseWithTmp.setUp(self)
        self.polio_fasta = os.path.join(
            util_file.get_test_input_path(),
            'TestMetagenomicsViralMix', 'db', 'library', 'Viruses', 'Enterovirus_C', 'GCF_000861165.1_ViralProj15288_genomic.fna'
        )
        dbDir = tempfile.mkdtemp()
        self.lastdb_path = last.Lastdb().build_database(self.polio_fasta, os.path.join(dbDir, 'NC_002058'))

    def test_filter_lastal_bam_polio(self):
        inBam = os.path.join(util_file.get_test_input_path(), 'TestDepleteHuman', 'expected', 'test-reads.blastn.bam')
        outBam = util_file.mkstempfname('-out-taxfilt.bam')
        args = taxon_filter.parser_filter_lastal_bam(argparse.ArgumentParser()).parse_args([
            inBam, self.lastdb_path, outBam])
        args.func_main(args)
        expectedOut = os.path.join(util_file.get_test_input_path(), 'TestDepleteHuman', 'expected', 'test-reads.taxfilt.imperfect-2.bam')
        assert_equal_bam_reads(self, outBam, expectedOut)

    def test_lastal_empty_input(self):
        empty_bam = os.path.join(util_file.get_test_input_path(), 'empty.bam')
        outBam = util_file.mkstempfname('-out-taxfilt.bam')
        taxon_filter.filter_lastal_bam(
            empty_bam,
            self.lastdb_path,
            outBam
        )
        assert_equal_bam_reads(self, outBam, empty_bam)

    def test_lastal_empty_output(self):
        empty_bam = os.path.join(util_file.get_test_input_path(), 'empty.bam')
        in_bam = os.path.join(util_file.get_test_input_path(), 'TestDepleteHuman', 'test-reads-human.bam')
        outBam = util_file.mkstempfname('-out-taxfilt.bam')
        taxon_filter.filter_lastal_bam(
            in_bam,
            self.lastdb_path,
            outBam
        )
        assert_equal_bam_reads(self, outBam, empty_bam)

    def test_lastal_unbuilt_db(self):
        empty_bam = os.path.join(util_file.get_test_input_path(), 'empty.bam')
        in_bam = os.path.join(util_file.get_test_input_path(), 'TestDepleteHuman', 'test-reads-human.bam')
        outBam = util_file.mkstempfname('-out-taxfilt.bam')
        taxon_filter.filter_lastal_bam(
            in_bam,
            self.polio_fasta,
            outBam
        )
        assert_equal_bam_reads(self, outBam, empty_bam)


@unittest.skipIf(IS_ARM, SKIP_BMTAGGER_REASON)
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
        ref_fasta = os.path.join(util_file.get_test_input_path(), '5kb_human_from_chr6.fasta')
        self.database_prefix_path = os.path.join(self.tempDir, "5kb_human_from_chr6")
        taxon_filter.bmtagger_build_db(ref_fasta, self.tempDir, "5kb_human_from_chr6", word_size=8)

    def test_deplete_bmtagger_bam(self):
        inBam = os.path.join(util_file.get_test_input_path(), 'TestDepleteHuman', 'test-reads.bam')
        outBam = util_file.mkstempfname('-out.bam')
        args = taxon_filter.parser_deplete_bam_bmtagger(argparse.ArgumentParser()).parse_args([
            inBam, self.database_prefix_path, outBam, '--srprismMemory', '1500'])
        args.func_main(args)
        expectedOut = os.path.join(util_file.get_test_input_path(), 'TestDepleteHuman', 'expected', 'test-reads.bmtagger.bam')
        assert_equal_bam_reads(self, outBam, expectedOut)

    @unittest.skip("too slow for real word size of 18bp")
    def test_deplete_bmtagger_fasta_db(self):
        inBam = os.path.join(util_file.get_test_input_path(), 'TestDepleteHuman', 'test-reads.bam')
        ref_fasta = os.path.join(util_file.get_test_input_path(), '5kb_human_from_chr6.fasta')
        outBam = util_file.mkstempfname('-out.bam')
        args = taxon_filter.parser_deplete_bam_bmtagger(argparse.ArgumentParser()).parse_args([
            inBam, ref_fasta, outBam, '--srprismMemory', '1500'])
        args.func_main(args)
        expectedOut = os.path.join(util_file.get_test_input_path(), 'TestDepleteHuman', 'expected', 'test-reads.bmtagger.bam')
        assert_equal_bam_reads(self, outBam, expectedOut)

    def test_deplete_bmtagger_tar_db(self):
        inBam = os.path.join(util_file.get_test_input_path(), 'TestDepleteHuman', 'test-reads.bam')
        outBam = util_file.mkstempfname('-out.bam')
        tar_db_tgz = util_file.mkstempfname('.db.tar.gz')
        cmd = ['tar', '-C', os.path.dirname(self.database_prefix_path), '-cvzf', tar_db_tgz, '.']
        subprocess.check_call(cmd)
        args = taxon_filter.parser_deplete_bam_bmtagger(argparse.ArgumentParser()).parse_args([
            inBam, tar_db_tgz, outBam, '--srprismMemory', '1500'])
        args.func_main(args)
        expectedOut = os.path.join(util_file.get_test_input_path(), 'TestDepleteHuman', 'expected', 'test-reads.bmtagger.bam')
        assert_equal_bam_reads(self, outBam, expectedOut)
        os.unlink(tar_db_tgz)

    def test_bmtagger_empty_input(self):
        empty_bam = os.path.join(util_file.get_test_input_path(), 'empty.bam')
        out_bam = util_file.mkstempfname('-out.bam')
        args = taxon_filter.parser_deplete_bam_bmtagger(argparse.ArgumentParser()).parse_args([
            empty_bam, self.database_prefix_path, out_bam, '--srprismMemory', '1500'])
        args.func_main(args)
        assert_equal_bam_reads(self, out_bam, empty_bam)

    def test_bmtagger_empty_output(self):
        empty_bam = os.path.join(util_file.get_test_input_path(), 'empty.bam')
        in_bam = os.path.join(util_file.get_test_input_path(), 'TestDepleteHuman', 'test-reads-human.bam')
        out_bam = util_file.mkstempfname('-out.bam')
        args = taxon_filter.parser_deplete_bam_bmtagger(argparse.ArgumentParser()).parse_args([
            in_bam, self.database_prefix_path, out_bam, '--srprismMemory', '1500'])
        args.func_main(args)
        assert_equal_bam_reads(self, out_bam, empty_bam)


class TestBlastnDbBuild(TestCaseWithTmp):

    def test_blastn_db_build(self):
        commonInputDir = util_file.get_test_input_path()
        refFasta = os.path.join(commonInputDir, 'ebola.fasta')

        myInputDir = util_file.get_test_input_path(self)
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
        commonInputDir = util_file.get_test_input_path()
        refFasta = os.path.join(commonInputDir, 'ebola.fasta.gz')

        myInputDir = util_file.get_test_input_path(self)
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

        commonInputDir = util_file.get_test_input_path()

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


@unittest.skipIf(IS_ARM, SKIP_BMTAGGER_REASON)
class TestBmtaggerDbBuild(TestCaseWithTmp):

    def test_bmtagger_db_build(self):
        os.environ.pop('TMPDIR', None)
        util_file.set_tmp_dir(None)
        commonInputDir = util_file.get_test_input_path()
        refFasta = os.path.join(commonInputDir, 'ebola.fasta')

        myInputDir = util_file.get_test_input_path(self)
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
        commonInputDir = util_file.get_test_input_path()
        refFasta = os.path.join(commonInputDir, 'ebola.fasta.gz')
        myInputDir = util_file.get_test_input_path(self)
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
        commonInputDir = util_file.get_test_input_path()
        refFasta = os.path.join(commonInputDir, 'ebola.fasta')

        myInputDir = util_file.get_test_input_path(self)
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
            self.assertGreater(os.path.getsize(os.path.join(tempDir, output_prefix + ext)), 0)
            #assert_equal_contents(
            #    self, os.path.join(tempDir, output_prefix + ext),
            #    os.path.join(myInputDir, "expected", output_prefix + ext)
            #)

        for ext in [".suf"]:
            self.assertGreater(os.path.getsize(os.path.join(tempDir, output_prefix + ext)), 0)
        #    assert_md5_equal_to_line_in_file(self, os.path.join(tempDir, output_prefix + ext), os.path.join(myInputDir, "expected", output_prefix + ext+".md5"))

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
        commonInputDir = util_file.get_test_input_path()
        ref_fasta = os.path.join(commonInputDir, '5kb_human_from_chr6.fasta')
        self.database_prefix_path = os.path.join(self.tempDir, "5kb_human_from_chr6")

        # create blast db
        self.blastdb_path = blast.MakeblastdbTool().build_database(ref_fasta, self.database_prefix_path)

        # create multiple dbs
        self.blastdbs_multi = []
        for db in ['humanChr1Subset.fa', 'humanChr9Subset.fa']:
            dbPath = blast.MakeblastdbTool().build_database(
                os.path.join(util_file.get_test_input_path(self), db),
                os.path.join(self.tempDir, db[:-3]))
            self.blastdbs_multi.append(dbPath)

        # tar one db, but not the other
        tar_db_tgz = util_file.mkstempfname('-humanChr9Subset.blastn.db.tar.gz')
        cmd = ['tar', '-C', self.tempDir, '-cvzf', tar_db_tgz] + list(os.path.basename(f) for f in glob.glob(os.path.join(self.tempDir, "humanChr9Subset.n*")))
        subprocess.check_call(cmd)
        self.blastdbs_multi[1] = tar_db_tgz
        for idx in glob.glob(os.path.join(self.tempDir, "humanChr9Subset.n*")):
            os.unlink(idx)

    def test_deplete_blastn_bam(self):
        tempDir = tempfile.mkdtemp()
        myInputDir = util_file.get_test_input_path(self)

        # Run deplete_blastn_bam
        inBam = os.path.join(myInputDir, 'in.bam')
        outBam = os.path.join(tempDir, 'out.bam')
        args = taxon_filter.parser_deplete_blastn_bam(argparse.ArgumentParser()).parse_args(
            [inBam] + self.blastdbs_multi + [outBam, "--chunkSize", "0"]
        )
        args.func_main(args)

        # samtools view for out.sam and compare to expected
        outSam = os.path.join(tempDir, 'out.sam')
        samtools = samtools.SamtoolsTool()
        samtools.view(['-h'], outBam, outSam)

        assert_equal_bam_reads(self,
            outSam,
            os.path.join(myInputDir, 'expected.sam'))

    def test_deplete_blastn_bam_chunked(self):
        tempDir = tempfile.mkdtemp()
        myInputDir = util_file.get_test_input_path(self)

        # Run deplete_blastn_bam
        inBam = os.path.join(myInputDir, 'in.bam')
        outBam = os.path.join(tempDir, 'out.bam')
        args = taxon_filter.parser_deplete_blastn_bam(argparse.ArgumentParser()).parse_args(
            [inBam] + self.blastdbs_multi + [outBam, "--chunkSize", "1"]
        )
        args.func_main(args)

        # samtools view for out.sam and compare to expected
        outSam = os.path.join(tempDir, 'out.sam')
        samtools = samtools.SamtoolsTool()
        samtools.view(['-h'], outBam, outSam)

        assert_equal_bam_reads(self,
            outSam,
            os.path.join(myInputDir, 'expected.sam'))

    def test_blastn_empty_input(self):
        empty_bam = os.path.join(util_file.get_test_input_path(), 'empty.bam')
        out_bam = util_file.mkstempfname('-out.bam')
        taxon_filter.multi_db_deplete_bam(
            empty_bam,
            [self.blastdb_path],
            taxon_filter.deplete_blastn_bam,
            out_bam
        )
        self.assertEqual(0, samtools.SamtoolsTool().count(out_bam))

    def test_blastn_empty_output(self):
        in_bam = os.path.join(util_file.get_test_input_path(), 'TestDepleteHuman', 'test-reads-human.bam')
        out_bam = util_file.mkstempfname('-out.bam')
        taxon_filter.multi_db_deplete_bam(
            in_bam,
            [self.blastdb_path],
            taxon_filter.deplete_blastn_bam,
            out_bam
        )
        self.assertEqual(0, samtools.SamtoolsTool().count(out_bam))


class TestDepleteMinimap2Bam(TestCaseWithTmp):
    '''
    Tests for minimap2-based read depletion.
    Reuses test data from TestDepleteBlastnBam/TestBmtagger.
    '''

    def setUp(self):
        TestCaseWithTmp.setUp(self)
        self.tempDir = tempfile.mkdtemp()
        commonInputDir = util_file.get_test_input_path()
        # Reuse existing test reference fasta
        self.ref_fasta = os.path.join(commonInputDir, '5kb_human_from_chr6.fasta')

    def test_deplete_minimap2_bam(self):
        '''Basic depletion test - remove human reads from mixed input'''
        inBam = os.path.join(util_file.get_test_input_path(), 'TestDepleteHuman', 'test-reads.bam')
        outBam = util_file.mkstempfname('-out.bam')
        args = taxon_filter.parser_deplete_minimap2_bam(argparse.ArgumentParser()).parse_args([
            inBam, self.ref_fasta, outBam])
        args.func_main(args)
        # Should have fewer reads than input after depletion
        samtools = samtools.SamtoolsTool()
        self.assertLess(samtools.count(outBam), samtools.count(inBam))

    def test_minimap2_empty_input(self):
        '''Empty input BAM should produce empty output BAM'''
        empty_bam = os.path.join(util_file.get_test_input_path(), 'empty.bam')
        out_bam = util_file.mkstempfname('-out.bam')
        taxon_filter.multi_db_deplete_bam(
            empty_bam,
            [self.ref_fasta],
            taxon_filter.deplete_minimap2_bam,
            out_bam
        )
        self.assertEqual(0, samtools.SamtoolsTool().count(out_bam))

    def test_minimap2_empty_output(self):
        '''All-human input should result in empty output after depletion'''
        in_bam = os.path.join(util_file.get_test_input_path(), 'TestDepleteHuman', 'test-reads-human.bam')
        out_bam = util_file.mkstempfname('-out.bam')
        taxon_filter.multi_db_deplete_bam(
            in_bam,
            [self.ref_fasta],
            taxon_filter.deplete_minimap2_bam,
            out_bam
        )
        self.assertEqual(0, samtools.SamtoolsTool().count(out_bam))


class TestDepletePipeline(TestCaseWithTmp):
    '''
    Tests for the full deplete pipeline (main_deplete) including minimap2 depletion.
    '''

    def setUp(self):
        TestCaseWithTmp.setUp(self)
        self.tempDir = tempfile.mkdtemp()
        commonInputDir = util_file.get_test_input_path()
        self.ref_fasta = os.path.join(commonInputDir, '5kb_human_from_chr6.fasta')

    def test_deplete_pipeline_with_minimap(self):
        '''Test full deplete pipeline with minimap2 depletion enabled'''
        inBam = os.path.join(util_file.get_test_input_path(), 'TestDepleteHuman', 'test-reads.bam')
        revertBam = util_file.mkstempfname('-revert.bam')
        minimapBam = util_file.mkstempfname('-minimap.bam')
        bwaBam = util_file.mkstempfname('-bwa.bam')
        bmtaggerBam = util_file.mkstempfname('-bmtagger.bam')
        blastnBam = util_file.mkstempfname('-blastn.bam')

        args = taxon_filter.parser_deplete(argparse.ArgumentParser()).parse_args([
            inBam,
            revertBam,
            minimapBam,
            bwaBam,
            bmtaggerBam,
            blastnBam,
            '--minimapDbs', self.ref_fasta
        ])
        args.func_main(args)

        samtools = samtools.SamtoolsTool()
        # Verify output files exist and have fewer reads than input
        self.assertTrue(os.path.exists(minimapBam))
        self.assertLess(samtools.count(minimapBam), samtools.count(inBam))

    def test_deplete_pipeline_empty_minimap_dbs(self):
        '''Test deplete pipeline with empty minimapDbs (default behavior)'''
        inBam = os.path.join(util_file.get_test_input_path(), 'TestDepleteHuman', 'test-reads.bam')
        revertBam = util_file.mkstempfname('-revert.bam')
        minimapBam = util_file.mkstempfname('-minimap.bam')
        bwaBam = util_file.mkstempfname('-bwa.bam')
        bmtaggerBam = util_file.mkstempfname('-bmtagger.bam')
        blastnBam = util_file.mkstempfname('-blastn.bam')

        args = taxon_filter.parser_deplete(argparse.ArgumentParser()).parse_args([
            inBam,
            revertBam,
            minimapBam,
            bwaBam,
            bmtaggerBam,
            blastnBam,
            '--bwaDbs', self.ref_fasta  # Use bwa instead of minimap
        ])
        args.func_main(args)

        samtools = samtools.SamtoolsTool()
        # Verify minimap output exists (should be copy of input since no minimap dbs)
        self.assertTrue(os.path.exists(minimapBam))
        # With no minimap dbs, minimapBam should have same count as reverted input
        # BWA depletion happens after, so bwaBam should have fewer reads
        self.assertLess(samtools.count(bwaBam), samtools.count(inBam))
