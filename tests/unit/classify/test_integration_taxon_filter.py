# Integration tests for taxon_filter.py

__author__ = "dpark@broadinstitute.org"

import unittest
import os
import tempfile
import subprocess
import shutil

import argparse

from viral_ngs import taxon_filter
from viral_ngs.core import file as util_file
from viral_ngs.classify import last
from viral_ngs.classify import bmtagger
from viral_ngs.classify import blast
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
        myInputDir = util_file.get_test_input_path()
        ref_fasta = os.path.join(myInputDir, '5kb_human_from_chr6.fasta')
        self.database_prefix_path = os.path.join(self.tempDir, "5kb_human_from_chr6")

        # create blast db
        self.blastdb_path = blast.MakeblastdbTool().build_database(ref_fasta, self.database_prefix_path)

        # create bmtagger db
        taxon_filter.bmtagger_build_db(ref_fasta, self.tempDir, "5kb_human_from_chr6", word_size=8)

    def test_deplete_human(self):
        myInputDir = util_file.get_test_input_path(self)

        # Run deplete_human
        args = taxon_filter.parser_deplete(argparse.ArgumentParser()).parse_args(
            [
                os.path.join(myInputDir, 'test-reads.bam'),
                # output files
                os.path.join(self.tempDir, 'test-reads.revert.bam'),
                os.path.join(self.tempDir, 'test-reads.minimap.bam'),
                os.path.join(self.tempDir, 'test-reads.bwa.bam'),
                os.path.join(self.tempDir, 'test-reads.bmtagger.bam'),
                os.path.join(self.tempDir, 'test-reads.blastn.bam'),
                # DBs
                "--blastDbs", self.blastdb_path,
                "--bmtaggerDbs", self.database_prefix_path,
                "--chunkSize", "0",
                "--srprismMemory", '1500',
            ]
        )
        args.func_main(args)

        # Compare to expected
        for fname in [
            'test-reads.revert.bam',
            'test-reads.bwa.bam',
            'test-reads.bmtagger.bam',
            'test-reads.blastn.bam'
        ]:
            assert_equal_bam_reads(self, os.path.join(self.tempDir, fname), os.path.join(myInputDir, 'expected', fname))

    def test_deplete_human_aligned_input(self):
        myInputDir = util_file.get_test_input_path(self)

        # Run deplete_human
        args = taxon_filter.parser_deplete(argparse.ArgumentParser()).parse_args(
            [
                os.path.join(myInputDir, 'test-reads-aligned.bam'),
                # output files
                os.path.join(self.tempDir, 'test-reads.revert.bam'),
                os.path.join(self.tempDir, 'test-reads.minimap.bam'),
                os.path.join(self.tempDir, 'test-reads.bwa.bam'),
                os.path.join(self.tempDir, 'test-reads.bmtagger.bam'),
                os.path.join(self.tempDir, 'test-reads.blastn.bam'),
                # DBs
                "--blastDbs", self.blastdb_path,
                "--bmtaggerDbs", self.database_prefix_path,
                "--srprismMemory", '1500',
            ]
        )
        args.func_main(args)

        # Compare to expected
        for fname in [
            'test-reads.revert.bam',
            'test-reads.bwa.bam',
            'test-reads.bmtagger.bam',
            'test-reads.blastn.bam'
        ]:
            assert_equal_bam_reads(self, os.path.join(self.tempDir, fname), os.path.join(myInputDir, 'aligned-expected', fname))

    def test_deplete_empty(self):
        empty_bam = os.path.join(util_file.get_test_input_path(), 'empty.bam')

        # Run deplete_human
        args = taxon_filter.parser_deplete(argparse.ArgumentParser()).parse_args(
            [
                empty_bam,
                # output files
                os.path.join(self.tempDir, 'deplete-empty.revert.bam'),
                os.path.join(self.tempDir, 'deplete-empty.minimap.bam'),
                os.path.join(self.tempDir, 'deplete-empty.bwa.bam'),
                os.path.join(self.tempDir, 'deplete-empty.bmtagger.bam'),
                os.path.join(self.tempDir, 'deplete-empty.blastn.bam'),
                # DBs
                "--blastDbs", self.blastdb_path,
                "--bmtaggerDbs", self.database_prefix_path,
                "--srprismMemory", '1500',
            ]
        )
        args.func_main(args)

        # Compare to expected
        for fname in [
            'deplete-empty.bmtagger.bam', 'deplete-empty.blastn.bam',
        ]:
            assert_equal_bam_reads(self, os.path.join(self.tempDir, fname), empty_bam)
