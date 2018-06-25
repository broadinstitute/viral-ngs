"""Unit tests for kmers.py"""

__author__ = "ilya@broadinstitute.org"

import unittest
import os
import tempfile
import shutil
import filecmp
import subprocess

import argparse

import kmers
import util.cmd
import util.file
import util.misc
from test import assert_equal_contents, assert_equal_bam_reads, assert_md5_equal_to_line_in_file, TestCaseWithTmp

class TestKmc(TestCaseWithTmp):

    """Test the tool wrapper for KMC kmer counter"""

    def setUp(self):
        super(TestKmc, self).setUp()

    def test_simple_kmer_extraction(self):
        with util.file.tmp_dir(suffix='kmctest') as t_dir:
            simple_fasta = os.path.join(util.file.get_test_input_path(self), 'simple.fasta')
            kmer_db = os.path.join(t_dir, 'kmer_db')
            util.cmd.run_cmd(kmers, 'build_kmer_db', [simple_fasta, kmer_db, '-k', 4])
            kmers_fasta = os.path.join(t_dir, 'kmers.txt')
            kmers_fasta_exp = os.path.join(util.file.get_test_input_path(self), 'simple.fasta.kmers.k4.txt')
            util.cmd.run_cmd(kmers, 'dump_kmers', [kmer_db, kmers_fasta])
            assert_equal_contents(self, kmers_fasta, kmers_fasta_exp)

    def test_read_filtering(self):
        with util.file.tmp_dir(suffix='kmctest') as t_dir:
            simple_fasta = os.path.join(util.file.get_test_input_path(self), 'simple.fasta')
            kmer_db = os.path.join(t_dir, 'kmer_db')
            util.cmd.run_cmd(kmers, 'build_kmer_db', [simple_fasta, kmer_db, '-k', 4])


            filt_fasta = os.path.join(util.file.get_test_input_path(self), 'filt.fasta')
            filt_fasta_out = os.path.join(t_dir, 'filtered.fasta')
            util.cmd.run_cmd(kmers, 'filter_by_kmers', [kmer_db, filt_fasta, filt_fasta_out, '--readMinOccs', 1])
