# Unit tests for diamond

#from builtins import super
import glob
import os.path
import tempfile
import unittest
from util import file
import tools.diamond
from test import TestCaseWithTmp


class TestToolDiamond(TestCaseWithTmp):

    def setUp(self):
        #super().setUp()
        super(TestToolDiamond, self).setUp()
        self.diamond = tools.diamond.Diamond()
        self.diamond.install()

    @property
    def data_dir(self):
        return file.get_test_input_path(self)

    @property
    def db_dir(self):
        return os.path.join(self.data_dir, 'db')

    def test_create_and_blastx(self):
        db = os.path.join(tempfile.tempdir, 'db.dmnd')

        # To be replaced with recursive glob in Python 3.5.
        protein_fastas = glob.glob('{}/library/*/*/*.ffn'.format(self.db_dir))
        self.assertEqual(0, self.diamond.build(db, protein_fastas).returncode)
        self.assertTrue(os.path.isfile(db))

        inputs = [os.path.join(self.data_dir, f)
                  for f in ['zaire_ebola.1.fastq', 'zaire_ebola.2.fastq']]
        output = os.path.join(tempfile.tempdir, 'zaire_ebola.daa')
        output_tab = os.path.join(tempfile.tempdir, 'zaire_ebola.m8')

        tmpdir = tempfile.gettempdir()
        self.assertEqual(
            0, self.diamond.blastx(db, inputs, output,
                                   options={'--tmpdir': tmpdir}).returncode)
        self.assertEqual(0, self.diamond.view(output, output_tab).returncode)


if __name__ == '__main__':
    unittest.main()
