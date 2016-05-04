# Integration tests for kraken

from builtins import super
import argparse
import os.path
from os.path import join
import sys
import tempfile
import unittest
import metagenomics
import util.file
import tools.kraken
import tools.krona
import tools.picard
from test import TestCaseWithTmp
from test.integration.snake import SnakemakeRunner


class TestKrakenBase(TestCaseWithTmp):

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.kraken = tools.kraken.Kraken()
        cls.kraken.install()

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
    def build_kraken_db(cls):
        db = tempfile.mkdtemp('kraken_db')
        for d in ['library', 'taxonomy']:
            realpath = os.path.join(cls.db_dir, d)
            name = os.path.join(db, d)
            os.symlink(realpath, name)

        # Minimizer len corresponds to memory/disk usage of index.
        assert cls.kraken.build(db, options={
            '--minimizer-len': 10,
            '--build': None,
        }).returncode == 0
        return db


class TestKronaBase(object):
    TAXONOMY_FILES = (
        'gi_taxid_nucl.dmp',
        'gi_taxid_prot.dmp',
        'names.dmp',
        'nodes.dmp',
        )

    @classmethod
    def setUpClass(cls):
        cls.krona = tools.krona.Krona()
        cls.krona.install()

    @classmethod
    def build_krona_db(cls):
        cls.data_dir = join(util.file.get_test_input_path(), 'TestKrakenViralMix')
        cls.db_dir = os.path.join(cls.data_dir, 'db')

        db = tempfile.mkdtemp('krona_db')
        for d in TestKronaBase.TAXONOMY_FILES:
            src = join(cls.db_dir, 'taxonomy', d)
            dest = os.path.join(db, d)
            os.symlink(src, dest)
        cls.krona.create_db(db)
        return db


class CommonTests(object):

    @unittest.skipIf(sys.version_info < (3,2), "Python version is too old for snakemake.")
    def test_pipes(self):
        runner = SnakemakeRunner()
        runner.set_override_config({
            'kraken_db': self.db,
            'diamond_db': None,
        })
        runner.setup()
        runner.link_samples([self.bam], destination='source')
        runner.create_sample_files(sample_files=['samples_metagenomics'])

        kraken_out = join(runner.workdir, runner.data_dir, runner.config['subdirs']['metagenomics'],
                             '.'.join([os.path.splitext(os.path.basename(self.bam))[0], 'kraken.report']))

        krona_out = join(runner.workdir, runner.data_dir, runner.config['subdirs']['metagenomics'],
                             '.'.join([os.path.splitext(os.path.basename(self.bam))[0], 'kraken.krona.html']))

        # runner.run(['all_metagenomics'])
        runner.run([kraken_out])
        self.assertGreater(os.path.getsize(kraken_out), 0)
        # self.assertGreater(os.path.getsize(krona_out), 0)

    def test_kraken(self):
        bin = join(util.file.get_project_path(), 'metagenomics.py')
        out_report = util.file.mkstempfname('.report')
        out_reads = util.file.mkstempfname('.reads.gz')
        cmd = [self.bam, self.db, '--outReport', out_report, '--outReads', out_reads]
        parser = metagenomics.parser_kraken(argparse.ArgumentParser())
        args = parser.parse_args(cmd)
        args.func_main(args)

        self.assertGreater(os.path.getsize(out_report), 0)
        self.assertGreater(os.path.getsize(out_reads), 0)


class TestKrakenTiny(TestKrakenBase, CommonTests):

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.data_dir = join(util.file.get_test_input_path(), 'TestToolKraken')
        cls.db_dir = os.path.join(cls.data_dir, 'db')
        cls.db = cls.build_kraken_db()
        cls.fastqs = [os.path.join(cls.data_dir, f)
                      for f in ['zaire_ebola.1.fastq', 'zaire_ebola.2.fastq']]
        cls.bam = cls.input_bam('zaire_ebola')

    def test_kraken_tool(self):
        outdir = tempfile.mkdtemp('-kraken')
        out = join(outdir, 'zaire_ebola.kraken')
        out_filtered = join(outdir, 'zaire_ebola.filtered-kraken')
        out_report = join(outdir, 'zaire_ebola.kraken-report')
        self.assertEqual(0, self.kraken.classify(self.db, self.fastqs, out).returncode)
        result = self.kraken.execute(
            'kraken-filter', self.db, out_filtered, [out],
            options={'--threshold': 0.05})
        self.assertEqual(0, result.returncode)
        result = self.kraken.execute(
            'kraken-report', self.db, out_report, [out_filtered])
        self.assertEqual(0, result.returncode)
        self.assertGreater(os.path.getsize(out_report), 0)
        self.assertGreater(os.path.getsize(out_filtered), 0)

class TestKrakenVirusMix(TestKrakenBase, CommonTests):

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.data_dir = join(util.file.get_test_input_path(), 'TestKrakenViralMix')
        cls.db_dir = os.path.join(cls.data_dir, 'db')
        cls.db = cls.build_kraken_db()
        cls.bam = join(cls.data_dir, 'test-reads.bam')
        cls.fastqs = cls.input_fastqs()


class TestKrakenKrona(TestKrakenBase, TestKronaBase):

    @classmethod
    def setUpClass(cls):
        super(TestKrakenKrona, cls).setUpClass()
        TestKronaBase.setUpClass()
        cls.data_dir = join(util.file.get_test_input_path(), 'TestKrakenViralMix')
        cls.db_dir = os.path.join(cls.data_dir, 'db')
        cls.db = cls.build_kraken_db()
        cls.krona_db = cls.build_krona_db()
        cls.bam = join(cls.data_dir, 'test-reads.bam')
        cls.fastqs = cls.input_fastqs()

    def test_kraken_krona(self):
        bin = join(util.file.get_project_path(), 'metagenomics.py')
        out_report = util.file.mkstempfname('.report')
        out_reads = util.file.mkstempfname('.reads.gz')

        cmd = [self.bam, self.db, '--outReport', out_report, '--outReads', out_reads]
        parser = metagenomics.parser_kraken(argparse.ArgumentParser())
        args = parser.parse_args(cmd)
        args.func_main(args)

        out_html = util.file.mkstempfname('.krona.html')
        parser = metagenomics.parser_krona(argparse.ArgumentParser())
        args = parser.parse_args([out_reads, self.krona_db, out_html])
        args.func_main(args)
