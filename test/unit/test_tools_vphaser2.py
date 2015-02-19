# Unit tests for vphaser tool

__author__ = "irwin@broadinstitute.org"

import unittest, os, cPickle
import util.file
from tools.vphaser2 import Vphaser2Tool
from test import TestCaseWithTmp

class TestVPhaser2(TestCaseWithTmp) :
    def test_vphaser2(self) :
        myInputDir = util.file.get_test_input_path(self)
        inBam = os.path.join(myInputDir, 'in.bam')
        recs = list(Vphaser2Tool().iterate(inBam, numThreads = 8))
        expectedRecs = cPickle.load(open(os.path.join(myInputDir, 'expected.cp')))
        # Vphaser2 p-val calculation is unstable and sometimes varies from
        # run to run, so exclude it from comparison.
        self.assertEqual([rec[:4] + rec[5:] for rec in recs],
                         [rec[:4] + rec[5:] for rec in expectedRecs])
        """
        Creation of in.bam:
        Start with test file that ships with V-Phaser 2.
        cp 4528.454.indelRealigned.bam orig.bam
        samtools index orig.bam
        samtools view -h orig.bam V4528_assembly:1-100 V4528_assembly:950-1050 >c1.sam
        samtools view -h orig.bam V4528_assembly:8950-9050 V4528_assembly:9900-9950 >c2.sam
        Change all occurences of V4528_assembly in c1.sam to chr1,
        Change all occurences of V4528_assembly in c2.sam to chr2,
        Move the @SQ line from c2.sam to c1.sam and delete header of c2.sam.
        cat c1.sam c2.sam >new.sam
        samtools view -bh new.sam >new.bam
        
        Creation of expected.cp:
        cPickle.dump(list(Vphaser2Tool().iterate(inBam, numThreads = 8)),
                     open('expected.cp', 'w'))
        """
