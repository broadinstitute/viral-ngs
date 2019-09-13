# Integration tests for assembly.py

__author__ = "dpark@broadinstitute.org"

import os
import shutil
import tempfile
import itertools
import unittest
import pytest

import Bio.SeqIO

import assembly
import util.cmd
import util.file
import tools.novoalign
from test import TestCaseWithTmp, _CPUS


@unittest.skip("redundant, and takes 1 minute")
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
        inBam = os.path.join(util.file.get_test_input_path(), 'G5012.3.mini.bam')
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
            os.path.join(util.file.get_test_input_path(), 'G5012.3.mini.bam'),
            refine1Fasta,
            # normally -r Random, but for unit tests, we want deterministic behavior
            novo_params='-r None -l 30 -x 20 -t 502',
            min_coverage=2)
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
            os.path.join(util.file.get_test_input_path(), 'G5012.3.mini.bam'),
            refine2Fasta,
            # normally -r Random, but for unit tests, we want deterministic behavior
            novo_params='-r None -l 40 -x 20 -t 100',
            min_coverage=3)
        actual = str(Bio.SeqIO.read(refine2Fasta, 'fasta').seq)
        expected = str(Bio.SeqIO.read(os.path.join(inDir, 'expected.ebov.refine2.fasta'), 'fasta').seq)
        self.assertEqual(actual, expected)


class TestOrderOrientAndImputeFromReference(TestCaseWithTmp):
    # common setup
    def setUp(self):
        super(TestOrderOrientAndImputeFromReference, self).setUp()
        self.inDir = util.file.get_test_input_path(self)

        self.refFasta = os.path.join(self.inDir, 'ref.influenza_partial.fasta')
        self.outOrientFasta = util.file.mkstempfname('.fasta')
        assembly.order_and_orient(
            os.path.join(self.inDir, 'contigs.influenza.fasta'),
            self.refFasta,
            self.outOrientFasta)

    # common teardown
    def tearDown(self):
        os.unlink(self.outOrientFasta)
        super(TestOrderOrientAndImputeFromReference, self).tearDown()

    def test_impute_from_oriented_muscle(self):
        self.influenza_impute("muscle")

    def test_impute_from_oriented_mafft(self):
        self.influenza_impute("mafft")

    def test_impute_from_oriented_mummer(self):
        self.influenza_impute("mummer")

    # common impute function using the specified aligner
    def influenza_impute(self, aligner):
        outImputeFasta = util.file.mkstempfname('.fasta')
        expected = os.path.join(self.inDir, 'expected.influenza.impute.'+aligner+'.fasta')
        # ensure we can run impute_from_reference from the output of order_and_orient
        # without errors, but don't actually check the output
        assembly.impute_from_reference(
            self.outOrientFasta,
            self.refFasta,
            outImputeFasta,
            minLengthFraction=0.8,
            minUnambig=0.2,
            replaceLength=5,
            newName='test_influenza.genome',
            aligner=aligner)

        # if we were interested in checking the output...
        # self.assertEqualContents(
        #     outImputeFasta,
        #     expected
        # )
