# Unit tests for vphaser tool

__author__ = "irwin@broadinstitute.org"

import intrahost
import unittest, argparse, os, tempfile
import util.file
from tools.vphaser2 import Vphaser2Tool
from test import TestCaseWithTmp

class TestVPhaser2(TestCaseWithTmp) :
    def test_vphaser2(self) :
        inBam = os.path.join(util.file.get_binaries_path(),
                              'V-Phaser-2.0', 'TestData',
                              '4528.454.indelRealigned.bam')
        outDir = tempfile.mkdtemp()
        Vphaser2Tool().execute(inBam, outDir, numThreads = 8)
        filesToCompare = [
                          'V4528_assembly.5000.5499.region',
                          'V4528_assembly.covplot.R',
                          'V4528_assembly.eb',
                          # p-values in following vary from run to run; exclude.
                          #'V4528_assembly.fdr.var.txt',
                          #'V4528_assembly.nofdr.var.txt',
                          #'V4528_assembly.var.raw.txt',
                         ]
        myComparisonDir = util.file.get_test_input_path(self)
        for fileName in filesToCompare :
            self.assertEqualContents(os.path.join(myComparisonDir, fileName),
                                     os.path.join(outDir, fileName))
