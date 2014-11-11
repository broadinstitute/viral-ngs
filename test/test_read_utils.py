# Unit tests for read_utils.py

__author__ = "irwin@broadinstitute.org"

import unittest, os, tempfile
import util, util.file, read_utils, tools, tools.samtools
from test import assert_equal_contents, TestCaseWithTmp

class TestPurgeUnmated(TestCaseWithTmp) :
    def test_purge_unmated(self) :
        tempDir = tempfile.mkdtemp()
        myInputDir = util.file.get_test_input_path(self)
        inFastq1 = os.path.join(myInputDir, 'in1.fastq')
        inFastq2 = os.path.join(myInputDir, 'in2.fastq')
        outFastq1 = util.file.mkstempfname()
        outFastq2 = util.file.mkstempfname()
        parser = read_utils.parser_purge_unmated()
        args = parser.parse_args([inFastq1, inFastq2, outFastq1, outFastq2])
        read_utils.main_purge_unmated(args)

        # Check that results match expected
        expected1Fastq = os.path.join(myInputDir, 'expected1.fastq')
        expected2Fastq = os.path.join(myInputDir, 'expected2.fastq')
        assert_equal_contents(self, outFastq1, expected1Fastq)
        assert_equal_contents(self, outFastq2, expected2Fastq)

class TestFastqToFasta(TestCaseWithTmp) :
    def test_fastq_to_fasta(self) :
        tempDir = tempfile.mkdtemp()
        myInputDir = util.file.get_test_input_path(self)
        inFastq = os.path.join(myInputDir, 'in.fastq')
        outFasta = util.file.mkstempfname()
        parser = read_utils.parser_fastq_to_fasta()
        args = parser.parse_args([inFastq, outFasta])
        read_utils.main_fastq_to_fasta(args)

        # Check that results match expected
        expectedFasta = os.path.join(myInputDir, 'expected.fasta')
        assert_equal_contents(self, outFasta, expectedFasta)

class TestFastqBam(TestCaseWithTmp) :
    'Class for testing fastq <-> bam conversions'
    def test_fastq_bam(self) :
        tempDir = tempfile.mkdtemp()
        myInputDir = util.file.get_test_input_path(self)
        
        # Define file names
        inFastq1 = os.path.join(myInputDir, 'in1.fastq')
        inFastq2 = os.path.join(myInputDir, 'in2.fastq')
        inHeader = os.path.join(myInputDir, 'inHeader.txt')
        expectedSam = os.path.join(myInputDir, 'expected.sam')
        expectedFastq1 = os.path.join(myInputDir, 'expected.fastq1')
        outBamCmd = util.file.mkstempfname() + '.bam'
        outBamTxt = util.file.mkstempfname() + '.bam'
        outSam = util.file.mkstempfname() + '.sam'
        outFastq1 = util.file.mkstempfname() + '.fastq'
        outFastq2 = util.file.mkstempfname() + '.fastq'
        outHeader = util.file.mkstempfname() + '.txt'
        
        # in1.fastq, in2.fastq -> out.bam; header params from command-line
        parser = read_utils.parser_fastq_to_bam()
        args = parser.parse_args([inFastq1, inFastq2, outBamCmd,
            '--sampleName', 'FreeSample',
            '--JVMmemory', '1g',
            '--picardOptions',
            'LIBRARY_NAME=Alexandria',
            'PLATFORM=9.75',
            'SEQUENCING_CENTER=KareemAbdul-Jabbar',
            ])
        read_utils.main_fastq_to_bam(args)

        # samtools view for out.sam and compare to expected
        samtoolsPath = tools.samtools.SamtoolsTool().install_and_get_path()
        viewCmd = '{samtoolsPath} view -h {outBamCmd} > {outSam}'.format(
            **locals())
        assert not os.system(viewCmd)
        assert_equal_contents(self, outSam, expectedSam)

        # in1.fastq, in2.fastq, inHeader.txt -> out.bam; header from txt
        parser = read_utils.parser_fastq_to_bam()
        args = parser.parse_args([inFastq1, inFastq2, outBamTxt,
            '--header', inHeader])
        read_utils.main_fastq_to_bam(args)

        # out.bam -> out1.fastq, out2.fastq, outHeader.txt; trim 1 base from 1
        parser = read_utils.parser_bam_to_fastq()
        args = parser.parse_args([outBamTxt, outFastq1, outFastq2,
            '--outHeader', outHeader,
            '--JVMmemory', '1g',
            '--picardOptions', 'READ1_TRIM=1',
            ])
        read_utils.main_bam_to_fastq(args)

        # compare to out1.fastq, out2.fastq, outHeader.txt to in and expected
        assert_equal_contents(self, outFastq1, expectedFastq1) # 1 base trimmed
        assert_equal_contents(self, outFastq2, inFastq2)
        assert_equal_contents(self, outHeader, inHeader)

if __name__ == '__main__':
    unittest.main()
