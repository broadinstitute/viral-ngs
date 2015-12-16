# Unit tests for kraken

#from builtins import super
import os.path
import tempfile
import unittest
from util import file
import tools.kraken
from test import TestCaseWithTmp


class TestToolKraken(TestCaseWithTmp):

    def setUp(self):
        #super().setUp()
        super(TestToolKraken, self).setUp()
        self.kraken = tools.kraken.Kraken()
        self.kraken.install()

    @property
    def data_dir(self):
        return file.get_test_input_path(self)

    @property
    def db_dir(self):
        return os.path.join(self.data_dir, 'db')

    def test_create_and_classify(self):
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
        expected_files = [
            'database.idx',
            'database.jdb',
            'database.kdb',
            'gi2seqid.map',
            'lca.complete',
            'seqid2taxid.map'
        ]
        for f in expected_files:
            path = os.path.join(db, f)
            self.assertTrue(os.path.isfile(path))
        inputs = [os.path.join(self.data_dir, f)
                  for f in ['zaire_ebola.1.fastq', 'zaire_ebola.2.fastq']]
        output = os.path.join(tempfile.tempdir, 'zaire_ebola.kraken')
        self.assertEqual(0, self.kraken.classify(db, args=inputs, options={
            '--output': output
        }).returncode)
        result = self.kraken.execute(
            'kraken-filter', db, args=[output],
            options={'--threshold': 0.05})
        self.assertEqual(0, result.returncode)
        result = self.kraken.execute(
            'kraken-report', db, args=[output])
        self.assertEqual(0, result.returncode)


if __name__ == '__main__':
    unittest.main()
