# Integration tests for assembly.py

__author__ = "dpark@broadinstitute.org"

import assembly, util.cmd, util.file, tools.novoalign
import unittest
import os, shutil, tempfile, argparse, itertools
from test import TestCaseWithTmp

class TestAssemble(TestCaseWithTmp):
    ''' Test the de novo assembly pipeline '''

    def test_ref_assisted_assembly(self):
        novoalign = tools.novoalign.NovoalignTool()
        novoalign.install()
        
        orig_ref = os.path.join(util.file.get_test_input_path(), 'ebov-makona.fasta')
        refGenome = util.file.mkstempfname('.ref.fasta')
        shutil.copyfile(orig_ref, refGenome)
        novoalign.index_fasta(inRef)
        inBam = os.path.join(util.file.get_test_input_path(), 'G5012.3.testreads.bam')
        outFasta = util.file.mkstempfname('.refined.fasta')
        
        args = [refGenome, inBam, outFasta,
                "--chr_names", 'G5012.3',
                "--min_coverage", '3',
                "--novo_params", "-r Random -l 30 -g 40 -x 20 -t 502"]
        args = assembly.parser_refine_assembly().parse_args(args)
        args.func_main(args)
        self.assertTrue(os.path.isfile(outFasta))
        self.assertTrue(os.path.getsize(outFasta) > 1000)

# in order to test the actual de novo pipeline, we need to add a clip db for trimmomatic
# then we should test from G5012.3.testreads.bam all the way through the assembly pipe

