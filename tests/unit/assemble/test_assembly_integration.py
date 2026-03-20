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
import viral_ngs.core.minimap2
import viral_ngs.core.picard
import viral_ngs.core.samtools
from tests import TestCaseWithTmp, _CPUS, IS_ARM


def _align_with_minimap2(refFasta, inBam):
    """Align reads to reference with minimap2, return indexed BAM path."""
    mm2 = viral_ngs.core.minimap2.Minimap2()
    samtools = viral_ngs.core.samtools.SamtoolsTool()
    outBam = viral_ngs.core.file.mkstempfname('.aligned.bam')
    mm2.align_bam(inBam, refFasta, outBam, options=['-x', 'sr'])
    samtools.index(outBam)
    return outBam


class TestRefineAssembly(TestCaseWithTmp):
    def test_ebov_refine1(self):
        inDir = viral_ngs.core.file.get_test_input_path(self)
        inFasta = os.path.join(inDir, 'impute.ebov.fasta')
        imputeFasta = viral_ngs.core.file.mkstempfname('.imputed.fasta')
        refine1Fasta = viral_ngs.core.file.mkstempfname('.refine1.fasta')
        shutil.copy(inFasta, imputeFasta)
        viral_ngs.core.picard.CreateSequenceDictionaryTool().execute(imputeFasta, overwrite=True)
        viral_ngs.core.samtools.SamtoolsTool().faidx(imputeFasta, overwrite=True)
        inBam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'G5012.3.mini.bam')
        alignedBam = _align_with_minimap2(imputeFasta, inBam)
        viral_ngs.assembly.refine_assembly(
            imputeFasta, inBam, refine1Fasta,
            already_realigned_bam=alignedBam,
            min_coverage=2)
        actual = str(Bio.SeqIO.read(refine1Fasta, 'fasta').seq)
        expected = str(Bio.SeqIO.read(os.path.join(inDir, 'expected.ebov.refine1.freebayes.fasta'), 'fasta').seq)
        self.assertEqual(actual, expected)

    def test_ebov_refine2(self):
        inDir = viral_ngs.core.file.get_test_input_path(self)
        inFasta = os.path.join(inDir, 'expected.ebov.refine1.fasta')
        refine1Fasta = viral_ngs.core.file.mkstempfname('.refine1.fasta')
        refine2Fasta = viral_ngs.core.file.mkstempfname('.refine2.fasta')
        shutil.copy(inFasta, refine1Fasta)
        viral_ngs.core.picard.CreateSequenceDictionaryTool().execute(refine1Fasta, overwrite=True)
        viral_ngs.core.samtools.SamtoolsTool().faidx(refine1Fasta, overwrite=True)
        inBam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'G5012.3.mini.bam')
        alignedBam = _align_with_minimap2(refine1Fasta, inBam)
        viral_ngs.assembly.refine_assembly(
            refine1Fasta, inBam, refine2Fasta,
            already_realigned_bam=alignedBam,
            min_coverage=3)
        actual = str(Bio.SeqIO.read(refine2Fasta, 'fasta').seq)
        expected = str(Bio.SeqIO.read(os.path.join(inDir, 'expected.ebov.refine2.freebayes.fasta'), 'fasta').seq)
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
