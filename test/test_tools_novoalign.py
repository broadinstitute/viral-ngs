# Unit tests for Novoalign aligner

__author__ = "dpark@broadinstitute.org"

import unittest, os.path, shutil
import util.file, tools.novoalign
from test import TestCaseWithTmp

class TestToolNovoalign(TestCaseWithTmp) :

    def setUp(self):
        super(TestToolNovoalign, self).setUp()
        self.novoalign = tools.novoalign.NovoalignTool()
        self.novoalign.install()

    def test_index(self) :
        orig_ref = os.path.join(util.file.get_test_input_path(),
            'ebola.fasta')
        inRef = util.file.mkstempfname('.fasta')
        shutil.copyfile(orig_ref, inRef)
        self.novoalign.index_fasta(inRef)
        outfile = inRef[:-6] + '.nix'
        self.assertTrue(os.path.isfile(outfile))
        self.assertTrue(os.path.getsize(outfile))
