"""Unit tests for kmers.py"""

__author__ = "ilya@broadinstitute.org"

import unittest
import os
import tempfile
import shutil
import filecmp
import subprocess
import collections
import operator
import functools


import argparse

import kmers
import util.cmd
import util.file
import util.misc
import tools.kmc
from test import assert_equal_contents, assert_equal_bam_reads, assert_md5_equal_to_line_in_file, TestCaseWithTmp

import Bio.SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

class TestKmc(TestCaseWithTmp):

    """Test the tool wrapper for KMC kmer counter"""

    def setUp(self):
        super(TestKmc, self).setUp()

    def _get_seq(self, s):
        if isinstance(s, Seq): return str(s)
        if isinstance(s, SeqRecord): return str(s.seq)
        return s

    def _revcomp(self, kmer):
        return str(Seq(kmer, IUPAC.unambiguous_dna).reverse_complement())
    
    def _canonize(self, kmer):
        return min(kmer, self._revcomp(kmer))

    def _get_seq_kmers(self, seqs, k, single_strand):
        """Get kmers of seq(s)"""
        for seq in util.misc.make_seq(seqs):
            seq = self._get_seq(seq)
            for i in range(len(seq)-k+1):
                kmer = seq[i:i+k]
                yield kmer if single_strand else self._canonize(kmer)

    def _get_seq_kmer_counts(self, seqs, args):
        """Get kmer counts of seq(s)"""
        counts = collections.Counter(self._get_seq_kmers(seqs, args.kmer_size, args.single_strand))
        if any((args.min_occs, args.max_occs, args.counter_cap)):
            counts = dict((kmer, min(count, args.counter_cap or count)) \
                          for kmer, count in counts.items() \
                          if (count >= (args.min_occs or count)) and \
                          (count <= (args.max_occs or count)))
        return counts

    def _make_seq_recs(self, seqs):
        return [SeqRecord(Seq(seq, IUPAC.unambiguous_dna),
                          id='seq_%d'.format(i), name='seq_%d'.format(i), 
                          description='sequence number %d'.format(i)) 
                for i, seq in enumerate(util.misc.make_seq(seqs))]

    def _write_seqs_to_fasta(self, seqs, seqs_fasta):
        Bio.SeqIO.write(self._make_seq_recs(seqs),
                        seqs_fasta, 'fasta')

    def test_kmer_extraction(self):

        test_data = (
            ('A'*15, '-k 4'),
            ('T'*15, '-k 4' ),
            ([], '-k 1'),
            (['TCGA'*3, 'ATTT'*5], '-k 7'),
            (['TCGA'*3, 'ATTT'*5], '-k 31'),
        )

        for seqs, opts in test_data:
            with util.file.tmp_dir(suffix='kmctest') as t_dir:
                seqs_fasta = os.path.join(t_dir, 'seqs.fasta')
                self._write_seqs_to_fasta(seqs, seqs_fasta)
                kmer_db = os.path.join(t_dir, 'kmer_db')
                util.cmd.run_cmd(kmers, 'build_kmer_db', opts.split() + [seqs_fasta, kmer_db])
                args = util.cmd.parse_cmd(kmers, 'build_kmer_db', opts.split() + [seqs_fasta, kmer_db])

                kmers_txt = os.path.join(t_dir, 'kmers.txt')
                util.cmd.run_cmd(kmers, 'dump_kmer_counts', [kmer_db, kmers_txt])
                assert tools.kmc.KmcTool().read_kmer_counts(kmers_txt) == \
                    self._get_seq_kmer_counts(seqs, args)

    # to test:
    #   empty bam, empty fasta, empty both
    #   getting kmers from bam (here just convert the fasta to it?), from fastq,
    #   from multiple files.

    def _filter_seqs(self, kmer_counts, seqs, args):
        seqs_out = []
        for seq in util.misc.make_seq(seqs):
            kmerCounts = self._get_seq_kmer_counts(seqs)

    def test_read_filtering(self):
        with util.file.tmp_dir(suffix='kmctest') as t_dir:
            simple_fasta = self.input('simple.fasta')
            kmer_db = os.path.join(t_dir, 'kmer_db')
            util.cmd.run_cmd(kmers, 'build_kmer_db', [simple_fasta, kmer_db, '-k', 4])


            filt_fasta = os.path.join(util.file.get_test_input_path(self), 'filt.fasta')
            filt_fasta_out = os.path.join(t_dir, 'filtered.fasta')
            util.cmd.run_cmd(kmers, 'filter_by_kmers', [kmer_db, filt_fasta, filt_fasta_out, '--readMinOccs', 1])
