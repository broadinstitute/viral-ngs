# Integration tests for assembly.py

__author__ = "dpark@broadinstitute.org"

import assembly
import util.cmd
import util.file
import tools.novoalign
import Bio.SeqIO
import unittest
import os
import shutil
import tempfile
import argparse
import itertools
from test import TestCaseWithTmp


class TestAssemble(TestCaseWithTmp):
    ''' Test the de novo assembly pipeline '''

    def test_ref_assisted_assembly(self):
        novoalign = tools.novoalign.NovoalignTool()
        novoalign.install()

        # prep inputs
        orig_ref = os.path.join(util.file.get_test_input_path(), 'ebov-makona.fasta')
        refGenome = util.file.mkstempfname('.ref.fasta')
        shutil.copyfile(orig_ref, refGenome)
        novoalign.index_fasta(refGenome)
        inBam = os.path.join(util.file.get_test_input_path(), 'G5012.3.testreads.bam')
        outFasta = util.file.mkstempfname('.refined.fasta')

        # run refine_assembly
        args = [refGenome, inBam, outFasta, "--chr_names", 'G5012.3', "--min_coverage", '3', "--novo_params",
                "-r Random -l 30 -g 40 -x 20 -t 502"]
        args = assembly.parser_refine_assembly().parse_args(args)
        args.func_main(args)
        self.assertTrue(os.path.isfile(outFasta))
        self.assertTrue(os.path.getsize(outFasta) > 1000)

        # check assembly quality
        with open(outFasta, 'rt') as inf:
            seq = Bio.SeqIO.read(inf, 'fasta')
            self.assertGreater(len(seq), 17000)
            self.assertGreater(assembly.unambig_count(seq.seq), len(seq) * 0.95)

# in order to test the actual de novo pipeline, we need to add a clip db for trimmomatic
# then we should test from G5012.3.testreads.bam all the way through the assembly pipe

class TestRefineAssembly(TestCaseWithTmp):
    def test_ebov_refine1(self):
        inDir = util.file.get_test_input_path(self)
        inFasta = os.path.join(inDir, 'impute.ebov.fasta')
        imputeFasta = util.file.mkstempfname('.imputed.fasta')
        refine1Fasta = util.file.mkstempfname('.refine1.fasta')
        shutil.copy(inFasta, imputeFasta)
        tools.picard.CreateSequenceDictionaryTool().execute(imputeFasta, overwrite=True)
        tools.novoalign.NovoalignTool().index_fasta(imputeFasta)
        assembly.refine_assembly(
            imputeFasta,
            os.path.join(util.file.get_test_input_path(), 'G5012.3.testreads.bam'),
            refine1Fasta,
            # normally -r Random, but for unit tests, we want deterministic behavior
            novo_params='-r None -l 30 -x 20 -t 502',
            min_coverage=2,
            threads=4)
        actual = str(Bio.SeqIO.read(refine1Fasta, 'fasta').seq)
        expected = str(Bio.SeqIO.read(os.path.join(inDir, 'expected.ebov.refine1.fasta'), 'fasta').seq)
        self.assertEqual(actual, expected)

    def test_ebov_refine2(self):
        inDir = util.file.get_test_input_path(self)
        inFasta = os.path.join(inDir, 'expected.ebov.refine1.fasta')
        refine1Fasta = util.file.mkstempfname('.refine1.fasta')
        refine2Fasta = util.file.mkstempfname('.refine2.fasta')
        shutil.copy(inFasta, refine1Fasta)
        tools.picard.CreateSequenceDictionaryTool().execute(refine1Fasta, overwrite=True)
        tools.novoalign.NovoalignTool().index_fasta(refine1Fasta)
        assembly.refine_assembly(
            refine1Fasta,
            os.path.join(util.file.get_test_input_path(), 'G5012.3.testreads.bam'),
            refine2Fasta,
            # normally -r Random, but for unit tests, we want deterministic behavior
            novo_params='-r None -l 40 -x 20 -t 100',
            min_coverage=3,
            threads=4)
        actual = str(Bio.SeqIO.read(refine2Fasta, 'fasta').seq)
        expected = str(Bio.SeqIO.read(os.path.join(inDir, 'expected.ebov.refine2.fasta'), 'fasta').seq)
        self.assertEqual(actual, expected)

