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


class TestTrimmomatic(TestCaseWithTmp):

    def test_trimmomatic(self):
        myInputDir = util.file.get_test_input_path(self)
        inFastq1 = os.path.join(myInputDir, 'in1.fastq')
        inFastq2 = os.path.join(myInputDir, 'in2.fastq')
        pairedOutFastq1 = util.file.mkstempfname()
        pairedOutFastq2 = util.file.mkstempfname()
        clipFasta = os.path.join(myInputDir, 'clip.fasta')
        parser = taxon_filter.parser_trim_trimmomatic(argparse.ArgumentParser())
        args = parser.parse_args([inFastq1, inFastq2, pairedOutFastq1, pairedOutFastq2, clipFasta])
        args.func_main(args)

        # Check that results match expected
        expected1Fastq = os.path.join(myInputDir, 'expected1.fastq')
        expected2Fastq = os.path.join(myInputDir, 'expected2.fastq')
        assert_equal_contents(self, pairedOutFastq1, expected1Fastq)
        assert_equal_contents(self, pairedOutFastq2, expected2Fastq)


class TestFilterLastal(TestCaseWithTmp):

    def test_filter_lastal_fastq_ebola(self):
        # Create refDbs
        commonInputDir = util.file.get_test_input_path()
        myInputDir = util.file.get_test_input_path(self)
        refFasta = os.path.join(commonInputDir, 'ebola.fasta')
        dbsDir = tempfile.mkdtemp()
        refDbs = os.path.join(dbsDir, 'ebola')
        lastdbPath = tools.last.Lastdb().install_and_get_path()
        subprocess.check_call([lastdbPath, refDbs, refFasta])

        # Call main_filter_lastal
        inFastq = os.path.join(myInputDir, 'in.fastq')
        outFastq = util.file.mkstempfname('.fastq')
        args = taxon_filter.parser_filter_lastal(argparse.ArgumentParser()).parse_args([inFastq, refDbs, outFastq])
        args.func_main(args)

        # Check that results match expected
        expectedFastq = os.path.join(myInputDir, 'expected.fastq')
        assert_equal_contents(self, outFastq, expectedFastq)

    def test_filter_lastal_bam_polio(self):

        # Build polio database
        polio_fasta = os.path.join(
            util.file.get_test_input_path(),
            'TestMetagenomicsViralMix', 'db', 'library', 'Viruses', 'Poliovirus_uid15288', 'NC_002058.ffn'
        )
        dbDir = tempfile.mkdtemp()
        lastdb_path = tools.last.Lastdb().build_database(polio_fasta, os.path.join(dbDir, 'NC_002058'))

        # Call main_filter_lastal_bam
        inBam = os.path.join(util.file.get_test_input_path(), 'TestDepleteHuman', 'expected', 'test-reads.blastn.bam')
        outBam = util.file.mkstempfname('-out-taxfilt.bam', directory=dbDir)
        args = taxon_filter.parser_filter_lastal_bam(argparse.ArgumentParser()).parse_args([
            inBam, lastdb_path, outBam])
        args.func_main(args)

        # Check that results match expected
        expectedOut = os.path.join(util.file.get_test_input_path(), 'TestDepleteHuman', 'expected', 'test-reads.taxfilt.imperfect.bam')
        assert_equal_bam_reads(self, outBam, expectedOut)


class TestBmtagger(TestCaseWithTmp):
    """
    How test data was created:
      humanChr1Subset.fa has 200 bases from human chr1
      humanChr9Subset.fa has 200 bases from human chr9
      bmtool -d humanChr1Subset.fa -o humanChr1Subset.bitmask -w 8
      bmtool -d humanChr9Subset.fa -o humanChr9Subset.bitmask -w 8
      in[12].fastq "reads" are from humanChr[19]Subset.fa and ebola genome,
          with arbitrary quality scores.
    """

    def setUp(self):
        TestCaseWithTmp.setUp(self)
        self.tempDir = tempfile.mkdtemp()
        myInputDir = util.file.get_test_input_path(self)

        for db in ['humanChr1Subset', 'humanChr9Subset']:
            # .map file is > 100M, so recreate instead of copying
            dbfa = os.path.join(myInputDir, db + '.fa')
            dbsrprism = os.path.join(self.tempDir, db + '.srprism')
            srprismdb_path = tools.bmtagger.SrprismTool().build_database(dbfa, dbsrprism)
            # .bitmask and .srprism.* files must be in same dir, so copy
            shutil.copy(os.path.join(myInputDir, db + '.bitmask'), self.tempDir)

    def test_partition_bmtagger(self):
        outMatch = [os.path.join(self.tempDir, 'outMatch.{}.fastq'.format(n)) for n in '12']
        outNoMatch = [os.path.join(self.tempDir, 'outNoMatch.{}.fastq'.format(n)) for n in '12']
        myInputDir = util.file.get_test_input_path(self)
        args = taxon_filter.parser_partition_bmtagger(argparse.ArgumentParser()).parse_args(
            [
                os.path.join(myInputDir, 'in1.fastq'), os.path.join(myInputDir, 'in2.fastq'), os.path.join(
                    self.tempDir, 'humanChr1Subset'
                ), os.path.join(self.tempDir, 'humanChr9Subset'), '--outMatch', outMatch[0], outMatch[
                    1], '--outNoMatch', outNoMatch[0], outNoMatch[1]
            ]
        )
        args.func_main(args)

        # Compare to expected
        for case in ['Match.1', 'Match.2', 'NoMatch.1', 'NoMatch.2']:
            assert_equal_contents(
                self, os.path.join(self.tempDir, 'out' + case + '.fastq'),
                os.path.join(myInputDir, 'expected.' + case + '.fastq')
            )

    def test_deplete_bmtagger(self):
        myInputDir = util.file.get_test_input_path(self)
        args = taxon_filter.parser_partition_bmtagger(argparse.ArgumentParser()).parse_args(
            [
                os.path.join(myInputDir, 'in1.fastq'), os.path.join(myInputDir, 'in2.fastq'), os.path.join(
                    self.tempDir, 'humanChr1Subset'
                ), os.path.join(self.tempDir, 'humanChr9Subset'), '--outNoMatch', os.path.join(
                    self.tempDir, 'deplete.1.fastq'), os.path.join(self.tempDir, 'deplete.2.fastq')
            ]
        )
        args.func_main(args)

        # Compare to expected
        for case in ['1', '2']:
            assert_equal_contents(
                self, os.path.join(self.tempDir, 'deplete.' + case + '.fastq'),
                os.path.join(myInputDir, 'expected.NoMatch.' + case + '.fastq')
            )


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


class TestDepleteHuman(TestCaseWithTmp):
    '''
        How test data was created:
          exported 5kb region of chr6
          created pan-viral fasta file from all NCBI viral accessions
          used wgsim to create simulated reads
    '''

    def setUp(self):
        TestCaseWithTmp.setUp(self)
        self.tempDir = tempfile.mkdtemp()
        myInputDir = util.file.get_test_input_path(self)
        ref_fasta = os.path.join(myInputDir, '5kb_human_from_chr6.fasta')
        self.database_prefix_path = os.path.join(self.tempDir, "5kb_human_from_chr6")
        polio_fasta = os.path.join(
            util.file.get_test_input_path(),
            'TestMetagenomicsViralMix', 'db', 'library', 'Viruses', 'Poliovirus_uid15288', 'NC_002058.ffn'
        )

        # create blast db
        self.blastdb_path = tools.blast.MakeblastdbTool().build_database(ref_fasta, self.database_prefix_path)

        # create bmtagger db
        self.bmtooldb_path = tools.bmtagger.BmtoolTool().build_database(ref_fasta, self.database_prefix_path + ".bitmask")
        self.srprismdb_path = tools.bmtagger.SrprismTool().build_database(ref_fasta, self.database_prefix_path + ".srprism")

        # create last db
        self.lastdb_path = tools.last.Lastdb().build_database(polio_fasta, os.path.join(self.tempDir, 'polio'))

    def test_deplete_human(self):
        myInputDir = util.file.get_test_input_path(self)

        # Run deplete_human
        args = taxon_filter.parser_deplete_human(argparse.ArgumentParser()).parse_args(
            [
                os.path.join(myInputDir, 'test-reads.bam'),
                # output files
                os.path.join(self.tempDir, 'test-reads.bmtagger.bam'),
                os.path.join(self.tempDir, 'test-reads.rmdup.bam'),
                os.path.join(self.tempDir, 'test-reads.blastn.bam'),
                "--taxfiltBam", os.path.join(self.tempDir, 'test-reads.taxfilt.imperfect.bam'),
                # DBs
                "--blastDbs", self.blastdb_path,
                "--bmtaggerDbs", self.database_prefix_path,
                "--lastDb", self.lastdb_path,
                "--threads", "4"
            ]
        )
        args.func_main(args)

        # Compare to expected
        for fname in [
            'test-reads.bmtagger.bam',
            'test-reads.rmdup.bam', 'test-reads.blastn.bam',
            'test-reads.taxfilt.imperfect.bam'
        ]:
            assert_equal_bam_reads(self, os.path.join(self.tempDir, fname), os.path.join(myInputDir, 'expected', fname))

    def test_deplete_human_aligned_input(self):
        myInputDir = util.file.get_test_input_path(self)

        # Run deplete_human
        args = taxon_filter.parser_deplete_human(argparse.ArgumentParser()).parse_args(
            [
                os.path.join(myInputDir, 'test-reads-aligned.bam'),
                # output files
                os.path.join(self.tempDir, 'test-reads.revert.bam'),
                os.path.join(self.tempDir, 'test-reads.bmtagger.bam'),
                os.path.join(self.tempDir, 'test-reads.rmdup.bam'),
                os.path.join(self.tempDir, 'test-reads.blastn.bam'),
                "--taxfiltBam", os.path.join(self.tempDir, 'test-reads.taxfilt.imperfect.bam'),
                # DBs
                "--blastDbs", self.blastdb_path,
                "--bmtaggerDbs", self.database_prefix_path,
                "--lastDb", self.lastdb_path,
                "--threads", "4"
            ]
        )
        args.func_main(args)

        # Compare to expected
        for fname in [
            'test-reads.revert.bam', 'test-reads.bmtagger.bam',
            'test-reads.rmdup.bam', 'test-reads.blastn.bam',
            'test-reads.taxfilt.imperfect.bam'
        ]:
            assert_equal_bam_reads(self, os.path.join(self.tempDir, fname), os.path.join(myInputDir, 'aligned-expected', fname))

    def test_deplete_empty(self):
        myInputDir = util.file.get_test_input_path(self)
        empty_bam = os.path.join(util.file.get_test_input_path(), 'empty.bam')

        # Run deplete_human
        args = taxon_filter.parser_deplete_human(argparse.ArgumentParser()).parse_args(
            [
                empty_bam,
                # output files
                os.path.join(self.tempDir, 'deplete-empty.bmtagger.bam'),
                os.path.join(self.tempDir, 'deplete-empty.rmdup.bam'),
                os.path.join(self.tempDir, 'deplete-empty.blastn.bam'),
                "--taxfiltBam", os.path.join(self.tempDir, 'deplete-empty.taxfilt.bam'),
                # DBs
                "--blastDbs", self.blastdb_path,
                "--bmtaggerDbs", self.database_prefix_path,
                "--lastDb", self.lastdb_path,
                "--threads", "4"
            ]
        )
        args.func_main(args)

        # Compare to expected
        for fname in [
            'deplete-empty.bmtagger.bam',
            'deplete-empty.rmdup.bam', 'deplete-empty.blastn.bam',
            'deplete-empty.taxfilt.bam'
        ]:
            assert_equal_bam_reads(self, os.path.join(self.tempDir, fname), empty_bam)

    def test_revert_empty_input(self):
        empty_bam = os.path.join(util.file.get_test_input_path(), 'empty.bam')
        tools.picard.RevertSamTool().execute(
            empty_bam,
            util.file.mkstempfname(),
            picardOptions=['SORT_ORDER=queryname', 'SANITIZE=true']
        )

    def test_bmtagger_empty_input(self):
        empty_bam = os.path.join(util.file.get_test_input_path(), 'empty.bam')
        taxon_filter.multi_db_deplete_bam(
            empty_bam,
            [self.database_prefix_path],
            taxon_filter.deplete_bmtagger_bam,
            util.file.mkstempfname()
        )

    def test_bmtagger_empty_output(self):
        in_bam = os.path.join(util.file.get_test_input_path(self), 'test-reads-human.bam')
        taxon_filter.multi_db_deplete_bam(
            in_bam,
            [self.database_prefix_path],
            taxon_filter.deplete_bmtagger_bam,
            util.file.mkstempfname()
        )

    def test_mvicuna_empty_input(self):
        empty_bam = os.path.join(util.file.get_test_input_path(), 'empty.bam')
        read_utils.rmdup_mvicuna_bam(
            empty_bam,
            util.file.mkstempfname()
        )

    def test_blastn_empty_input(self):
        empty_bam = os.path.join(util.file.get_test_input_path(), 'empty.bam')
        taxon_filter.multi_db_deplete_bam(
            empty_bam,
            [self.blastdb_path],
            taxon_filter.deplete_blastn_bam,
            util.file.mkstempfname()
        )

    def test_blastn_empty_output(self):
        in_bam = os.path.join(util.file.get_test_input_path(self), 'test-reads-human.bam')
        taxon_filter.multi_db_deplete_bam(
            in_bam,
            [self.blastdb_path],
            taxon_filter.deplete_blastn_bam,
            util.file.mkstempfname()
        )

    def test_lastal_empty_input(self):
        empty_bam = os.path.join(util.file.get_test_input_path(), 'empty.bam')
        taxon_filter.filter_lastal_bam(
            empty_bam,
            self.lastdb_path,
            util.file.mkstempfname()
        )

    def test_lastal_empty_output(self):
        in_bam = os.path.join(util.file.get_test_input_path(self), 'test-reads-human.bam')
        taxon_filter.filter_lastal_bam(
            in_bam,
            self.lastdb_path,
            util.file.mkstempfname()
        )


class TestDepleteBlastnBam(TestCaseWithTmp):
    '''
        How test data was created:
        humanChr1Subset.fa has 200 bases from human chr1
        humanChr9Subset.fa has 200 bases from human chr9
        in.fastq "reads" are from humanChr[19]Subset.fa and ebola genome,
        with arbitrary quality scores.
    '''

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
