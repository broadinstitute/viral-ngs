# Integration tests for viral_ngs.assembly.py

__author__ = "dpark@broadinstitute.org"

import os
import shutil
import tempfile
import itertools
import unittest
import pytest

import Bio.SeqIO

import viral_ngs.assembly
import viral_ngs.core.cmd
import viral_ngs.core.file
import viral_ngs.core.novoalign
from tests import TestCaseWithTmp, _CPUS, IS_ARM

SKIP_X86_ONLY_REASON = "novoalign requires x86-only bioconda package (not available on ARM)"


@unittest.skip("redundant, and takes 1 minute")
class TestAssemble(TestCaseWithTmp):
    ''' Test the de novo assembly pipeline '''

    def test_ref_assisted_assembly(self):
        novoalign = viral_ngs.core.novoalign.NovoalignTool()
        novoalign.install()

        # prep inputs
        orig_ref = os.path.join(viral_ngs.core.file.get_test_input_path(), 'ebov-makona.fasta')
        refGenome = viral_ngs.core.file.mkstempfname('.ref.fasta')
        shutil.copyfile(orig_ref, refGenome)
        novoalign.index_fasta(refGenome)
        inBam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'G5012.3.mini.bam')
        outFasta = viral_ngs.core.file.mkstempfname('.refined.fasta')

        # run refine_assembly
        args = [refGenome, inBam, outFasta, "--chr_names", 'G5012.3', "--min_coverage", '3', "--novo_params",
                "-r Random -l 30 -g 40 -x 20 -t 502"]
        args = viral_ngs.assembly.parser_refine_assembly().parse_args(args)
        args.func_main(args)
        self.assertTrue(os.path.isfile(outFasta))
        self.assertTrue(os.path.getsize(outFasta) > 1000)

        # check assembly quality
        with open(outFasta, 'rt') as inf:
            seq = Bio.SeqIO.read(inf, 'fasta')
            self.assertGreater(len(seq), 17000)
            self.assertGreater(viral_ngs.assembly.unambig_count(seq.seq), len(seq) * 0.95)

# in order to test the actual de novo pipeline, we need to add a clip db for trimmomatic
# then we should test from G5012.3.testreads.bam all the way through the assembly pipe

@unittest.skipIf(IS_ARM, SKIP_X86_ONLY_REASON)
class TestRefineAssembly(TestCaseWithTmp):
    def test_ebov_refine1(self):
        inDir = viral_ngs.core.file.get_test_input_path(self)
        inFasta = os.path.join(inDir, 'impute.ebov.fasta')
        imputeFasta = viral_ngs.core.file.mkstempfname('.imputed.fasta')
        refine1Fasta = viral_ngs.core.file.mkstempfname('.refine1.fasta')
        shutil.copy(inFasta, imputeFasta)
        viral_ngs.core.picard.CreateSequenceDictionaryTool().execute(imputeFasta, overwrite=True)
        viral_ngs.core.novoalign.NovoalignTool().index_fasta(imputeFasta)
        viral_ngs.assembly.refine_assembly(
            imputeFasta,
            os.path.join(viral_ngs.core.file.get_test_input_path(), 'G5012.3.mini.bam'),
            refine1Fasta,
            # normally -r Random, but for unit tests, we want deterministic behavior
            novo_params='-r None -l 30 -x 20 -t 502',
            min_coverage=2)
        actual = str(Bio.SeqIO.read(refine1Fasta, 'fasta').seq)
        expected = str(Bio.SeqIO.read(os.path.join(inDir, 'expected.ebov.refine1.new.fasta'), 'fasta').seq)
        self.assertEqual(actual, expected)

    def test_ebov_refine2(self):
        inDir = viral_ngs.core.file.get_test_input_path(self)
        inFasta = os.path.join(inDir, 'expected.ebov.refine1.fasta')
        refine1Fasta = viral_ngs.core.file.mkstempfname('.refine1.fasta')
        refine2Fasta = viral_ngs.core.file.mkstempfname('.refine2.fasta')
        shutil.copy(inFasta, refine1Fasta)
        viral_ngs.core.picard.CreateSequenceDictionaryTool().execute(refine1Fasta, overwrite=True)
        viral_ngs.core.novoalign.NovoalignTool().index_fasta(refine1Fasta)
        viral_ngs.assembly.refine_assembly(
            refine1Fasta,
            os.path.join(viral_ngs.core.file.get_test_input_path(), 'G5012.3.mini.bam'),
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
        self.inDir = viral_ngs.core.file.get_test_input_path(self)

        self.refFasta = os.path.join(self.inDir, 'ref.influenza_partial.fasta')
        self.outOrientFasta = viral_ngs.core.file.mkstempfname('.fasta')
        viral_ngs.assembly.order_and_orient(
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
        outImputeFasta = viral_ngs.core.file.mkstempfname('.fasta')
        expected = os.path.join(self.inDir, 'expected.influenza.impute.'+aligner+'.fasta')
        # ensure we can run impute_from_reference from the output of order_and_orient
        # without errors, but don't actually check the output
        viral_ngs.assembly.impute_from_reference(
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
