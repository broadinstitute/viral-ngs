# Integration tests for kraken

from builtins import super
import os.path
from os.path import join
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
        self.kraken = tools.kraken.Kraken()
        self.kraken.install()

    @property
    def data_dir(self):
        return join(util.file.get_test_input_path(), 'TestToolKraken')

    @property
    def db_dir(self):
        return os.path.join(self.data_dir, 'db')

    def build_kraken_db(self):
        db = os.path.join(tempfile.tempdir, 'db')
        os.mkdir(db)
        for d in ['library', 'taxonomy']:
            realpath = os.path.join(self.db_dir, d)
            name = os.path.join(db, d)
            os.symlink(realpath, name)

        # Minimizer len corresponds to memory/disk usage of index.
        self.assertEqual(0, self.kraken.build(db, options={
            '--minimizer-len': 10,
            '--build': None,
        }).returncode)
        return db

    def test_pipes(self):
        """Build kraken db and execute Snakemake pipeline."""
        db = self.build_kraken_db()

        inputs = [os.path.join(self.data_dir, f)
                  for f in ['zaire_ebola.1.fastq', 'zaire_ebola.2.fastq']]
        bam = util.file.mkstempfname('.bam')

        picard = tools.picard.FastqToSamTool()
        picard.execute(inputs[0], inputs[1], 'zaire_ebola', bam)

        runner = SnakemakeRunner()
        runner.set_override_config({
            'kraken_db': db,
            'diamond_db': None,
        })
        runner.setup()
        runner.link_samples([bam], destination='source')
        runner.create_sample_files(sample_files=['samples_metagenomics'])

        kraken_output = join(runner.workdir, runner.data_dir, runner.config['subdirs']['metagenomics'],
                             '.'.join([os.path.splitext(os.path.basename(bam))[0], 'kraken.report']))
        runner.run([kraken_output])


if __name__ == '__main__':
    unittest.main()
