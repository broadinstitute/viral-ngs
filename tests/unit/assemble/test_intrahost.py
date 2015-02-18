# Unit tests for intrahost.py

__author__ = "irwin@broadinstitute.org"

import intrahost
import unittest, argparse, os, tempfile
import util.file, tools.vphaser
from test import TestCaseWithTmp

class TestCommandHelp(unittest.TestCase):
    def test_help_parser_for_each_command(self):
        for cmd_name, parser_fun in intrahost.__commands__:
            parser = parser_fun(argparse.ArgumentParser())
            helpstring = parser.format_help()

class TestVPhaser(TestCaseWithTmp) :
    def test_vphaser(self) :
        inBam = os.path.join(util.file.get_binaries_path(),
                              'V-Phaser-2.0', 'TestData', '4528.454.indelRealigned.bam')
        outDir = tempfile.mkdtemp()
        tools.vphaser.Vphaser2Tool().execute(inBam, outDir, numThreads = 8)
        filesToCompare = [
                          'V4528_assembly.5000.5499.region',
                          'V4528_assembly.covplot.R',
                          'V4528_assembly.eb',
                          # Some p-values in the following vary from run to run, so exclude for now.
                          #'V4528_assembly.fdr.var.txt',
                          #'V4528_assembly.nofdr.var.txt',
                          #'V4528_assembly.var.raw.txt',
                         ]
        myComparisonDir = util.file.get_test_input_path(self)
        for fileName in filesToCompare :
            self.assertEqualContents(os.path.join(myComparisonDir, fileName),
                                     os.path.join(outDir, fileName))
