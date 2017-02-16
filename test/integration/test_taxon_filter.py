# Integration tests for taxon_filter.py

__author__ = "dpark@broadinstitute.org"

import unittest
import os
import tempfile
import subprocess

import argparse

import taxon_filter
import util.file
import tools.last
import tools.bmtagger
import tools.blast
from test import assert_equal_bam_reads, TestCaseWithTmp


class TestDepleteHuman(TestCaseWithTmp):
    '''
        This class should move to test/integration. 

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
