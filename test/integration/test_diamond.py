# Integration tests for diamond

from builtins import super
import fnmatch
import os
from os.path import join
import sys
import tempfile
import unittest
from Bio import SeqIO
import util.file
import tools.diamond
import tools.picard
from test import TestCaseWithTmp
from test.integration.snake import SnakemakeRunner


def find_files(root_dir, filt):
    matches = []
    for root, dirnames, filenames in os.walk(root_dir):
        for filename in fnmatch.filter(filenames, '*.ffn'):
            yield join(root, filename)


class TestDiamondBase(TestCaseWithTmp):

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.diamond = tools.diamond.Diamond()
        cls.diamond.install()

    @classmethod
    def input_bam(cls, sample_name):
        bam = util.file.mkstempfname('.bam')
        picard = tools.picard.FastqToSamTool()
        picard.execute(cls.fastqs[0], cls.fastqs[1], 'zaire_ebola', bam)
        return bam

    @classmethod
    def input_fastqs(cls):
        fastq1 = util.file.mkstempfname('.1.fastq')
        fastq2 = util.file.mkstempfname('.2.fastq')
        picard = tools.picard.SamToFastqTool()
        picard.execute(cls.bam, fastq1, fastq2)
        return fastq1, fastq2

    @classmethod
    def build_db(cls):
        db = util.file.mkstempfname('.dmnd')
        translated = util.file.mkstempfname('.fa')
        lib_dir = join(cls.db_dir, 'library')
        with open(translated, "w") as f_out:
            for fname in find_files(cls.db_dir, '*.ffn'):
                with open(fname) as f:
                    for seq_record in SeqIO.parse(f, 'fasta'):
                        seq_record.seq = seq_record.seq.translate()
                        SeqIO.write(seq_record, f_out, 'fasta')
        assert cls.diamond.build(db, [translated]).returncode == 0
        return db


class CommonTests(object):

    @unittest.skipIf(sys.version_info < (3,2), "Python version is too old for snakemake.")
    def test_pipes(self):
        """Build diamond db and execute Snakemake pipeline."""
        runner = SnakemakeRunner()
        runner.set_override_config({
            'diamond_db': self.db,
        })
        runner.setup()
        runner.link_samples([self.bam], destination='source')
        runner.create_sample_files(sample_files=['samples_metagenomics'])

        diamond_output = join(runner.workdir, runner.data_dir, runner.config['subdirs']['metagenomics'],
                             '.'.join([os.path.splitext(os.path.basename(self.bam))[0], 'diamond.m8.gz']))
        runner.run([diamond_output])
        self.assertGreater(os.path.getsize(diamond_output), 0)

    def test_diamond(self):
        bin = join(util.file.get_project_path(), 'metagenomics.py')
        out_m8 = util.file.mkstempfname('.m8.gz')
        cmd = [bin, 'diamond', self.bam, self.db, out_m8]
        util.misc.run_and_print(cmd, check=True)
        self.assertGreater(os.path.getsize(out_m8), 0)


class TestDiamondTiny(TestDiamondBase, CommonTests):

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.data_dir = join(util.file.get_test_input_path(), 'TestToolDiamond')
        cls.db_dir = os.path.join(cls.data_dir, 'db')
        cls.db = cls.build_db()
        cls.fastqs = [os.path.join(cls.data_dir, f)
                      for f in ['zaire_ebola.1.fastq', 'zaire_ebola.2.fastq']]
        cls.bam = cls.input_bam('zaire_ebola')


class TestDiamondVirusMix(TestDiamondBase, CommonTests):

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.data_dir = join(util.file.get_test_input_path(), 'TestKrakenViralMix')
        cls.db_dir = os.path.join(cls.data_dir, 'db')
        cls.db = cls.build_db()
        cls.bam = join(cls.data_dir, 'test-reads.bam')
        cls.fastqs = cls.input_fastqs()
