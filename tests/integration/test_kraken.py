# Integration tests for kraken

from builtins import super
import os.path
from os.path import join
import sys
import tempfile
import unittest
import util.file
import tools.kraken
import tools.picard
from test import TestCaseWithTmp
from test.integration.snake import SnakemakeRunner


class TestKrakenTiny(TestCaseWithTmp):

    def setUp(self):
        super().setUp()

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.data_dir = join(util.file.get_test_input_path(), 'TestToolKraken')
        cls.db_dir = os.path.join(cls.data_dir, 'db')
        cls.kraken = tools.kraken.Kraken()
        cls.kraken.install()
        cls.db = cls.build_kraken_db()
        cls.bam = cls.input_bam()

    @classmethod
    def input_bam(cls):
        inputs = [os.path.join(cls.data_dir, f)
                  for f in ['zaire_ebola.1.fastq', 'zaire_ebola.2.fastq']]
        bam = util.file.mkstempfname('.bam')

        picard = tools.picard.FastqToSamTool()
        picard.execute(inputs[0], inputs[1], 'zaire_ebola', bam)
        return bam

    @classmethod
    def build_kraken_db(cls):
        db = os.path.join(tempfile.tempdir, 'db')
        os.mkdir(db)
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

    @unittest.skipIf(sys.version_info < (3,2), "Python version is too old for snakemake.")
    def test_pipes(self):
        """Build kraken db and execute Snakemake pipeline."""
        runner = SnakemakeRunner()
        runner.set_override_config({
            'kraken_db': self.db,
            'diamond_db': None,
        })
        runner.setup()
        runner.link_samples([self.bam], destination='source')
        runner.create_sample_files(sample_files=['samples_metagenomics'])

        kraken_output = join(runner.workdir, runner.data_dir, runner.config['subdirs']['metagenomics'],
                             '.'.join([os.path.splitext(os.path.basename(self.bam))[0], 'kraken.report']))
        runner.run([kraken_output])

    def test_kraken(self):
        bin = join(util.file.get_project_path(), 'metagenomics.py')
        out_report = util.file.mkstempfname('.report')
        out_reads = util.file.mkstempfname('.reads.gz')
        cmd = [bin, 'kraken', self.bam, self.db, '--outReport', out_report, '--outReads', out_reads]
        util.misc.run_and_print(cmd, check=True)


if __name__ == '__main__':
    unittest.main()
