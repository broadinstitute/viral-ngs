# Unit tests for read_utils.py

__author__ = "irwin@broadinstitute.org"

import unittest, os, tempfile, argparse, filecmp
import util, util.file, read_utils, tools, tools.samtools
from test import TestCaseWithTmp

class TestCommandHelp(unittest.TestCase):
    def test_help_parser_for_each_command(self):
        for cmd_name, parser_fun in read_utils.__commands__:
            parser = parser_fun(argparse.ArgumentParser())
            helpstring = parser.format_help()

class TestPurgeUnmated(TestCaseWithTmp) :
    def test_purge_unmated(self) :
        myInputDir = util.file.get_test_input_path(self)
        inFastq1 = os.path.join(myInputDir, 'in1.fastq')
        inFastq2 = os.path.join(myInputDir, 'in2.fastq')
        outFastq1 = util.file.mkstempfname('.fastq')
        outFastq2 = util.file.mkstempfname('.fastq')
        parser = read_utils.parser_purge_unmated(argparse.ArgumentParser())
        args = parser.parse_args([inFastq1, inFastq2, outFastq1, outFastq2])
        args.func_main(args)

        # Check that results match expected
        expected1Fastq = os.path.join(myInputDir, 'expected1.fastq')
        expected2Fastq = os.path.join(myInputDir, 'expected2.fastq')
        self.assertEqualContents(outFastq1, expected1Fastq)
        self.assertEqualContents(outFastq2, expected2Fastq)

    # test on FASTQs with read IDs in the style of SRA fastq-dump
    def test_purge_unmated_sra(self):
        myInputDir = util.file.get_test_input_path(self)
        inFastq1 = os.path.join(myInputDir, 'in_sra1.fastq')
        inFastq2 = os.path.join(myInputDir, 'in_sra2.fastq')
        outFastq1 = util.file.mkstempfname('.fastq')
        outFastq2 = util.file.mkstempfname('.fastq')
        parser = read_utils.parser_purge_unmated(argparse.ArgumentParser())
        args = parser.parse_args(['--regex', '^@(\S+).[1|2] .*', inFastq1, inFastq2, outFastq1, outFastq2])
        args.func_main(args)

        # The expected outputs are identical to the previous case.
        expected1Fastq = os.path.join(myInputDir, 'expected1.fastq')
        expected2Fastq = os.path.join(myInputDir, 'expected2.fastq')
        self.assertEqualContents(outFastq1, expected1Fastq)
        self.assertEqualContents(outFastq2, expected2Fastq)

class TestFastqToFasta(TestCaseWithTmp) :
    def test_fastq_to_fasta(self) :
        myInputDir = util.file.get_test_input_path(self)
        inFastq = os.path.join(myInputDir, 'in.fastq')
        outFasta = util.file.mkstempfname('.fasta')
        parser = read_utils.parser_fastq_to_fasta(argparse.ArgumentParser())
        args = parser.parse_args([inFastq, outFasta])
        args.func_main(args)

        # Check that results match expected
        expectedFasta = os.path.join(myInputDir, 'expected.fasta')
        self.assertEqualContents(outFasta, expectedFasta)

class TestFastqBam(TestCaseWithTmp) :
    'Class for testing fastq <-> bam conversions'
    def test_fastq_bam(self) :
        myInputDir = util.file.get_test_input_path(self)
        
        # Define file names
        inFastq1 = os.path.join(myInputDir, 'in1.fastq')
        inFastq2 = os.path.join(myInputDir, 'in2.fastq')
        inHeader = os.path.join(myInputDir, 'inHeader.txt')
        expected1_7Sam = os.path.join(myInputDir, 'expected.java1_7.sam')
        expected1_8Sam = os.path.join(myInputDir, 'expected.java1_8.sam')
        expectedFastq1 = os.path.join(myInputDir, 'expected.fastq1')
        outBamCmd = util.file.mkstempfname('.bam')
        outBamTxt = util.file.mkstempfname('.bam')
        outSam = util.file.mkstempfname('.sam')
        outFastq1 = util.file.mkstempfname('.fastq')
        outFastq2 = util.file.mkstempfname('.fastq')
        outHeader = util.file.mkstempfname('.txt')
        
        # in1.fastq, in2.fastq -> out.bam; header params from command-line
        parser = read_utils.parser_fastq_to_bam(argparse.ArgumentParser())
        args = parser.parse_args([inFastq1, inFastq2, outBamCmd,
            '--sampleName', 'FreeSample',
            '--JVMmemory', '1g',
            '--picardOptions',
            'LIBRARY_NAME=Alexandria',
            'PLATFORM=9.75',
            'SEQUENCING_CENTER=KareemAbdul-Jabbar',
            ])
        args.func_main(args)

        # samtools view for out.sam and compare to expected
        samtools = tools.samtools.SamtoolsTool()
        samtools.view(['-h'], outBamCmd, outSam)
        # picard.sam.FastqToSam outputs header fields in different order for
        #    java version 1.8 vs 1.7/1.6, so compare both
        self.assertTrue(filecmp.cmp(outSam, expected1_7Sam, shallow=False) or
                        filecmp.cmp(outSam, expected1_8Sam, shallow=False))
 
        # in1.fastq, in2.fastq, inHeader.txt -> out.bam; header from txt
        parser = read_utils.parser_fastq_to_bam(argparse.ArgumentParser())
        args = parser.parse_args([inFastq1, inFastq2, outBamTxt,
            '--header', inHeader])
        args.func_main(args)

        # out.bam -> out1.fastq, out2.fastq, outHeader.txt; trim 1 base from 1
        parser = read_utils.parser_bam_to_fastq(argparse.ArgumentParser())
        args = parser.parse_args([outBamTxt, outFastq1, outFastq2,
            '--outHeader', outHeader,
            '--JVMmemory', '1g',
            '--picardOptions', 'READ1_TRIM=1',
            ])
        args.func_main(args)

        # compare to out1.fastq, out2.fastq, outHeader.txt to in and expected
        self.assertEqualContents(outFastq1, expectedFastq1) # 1 base trimmed
        self.assertEqualContents(outFastq2, inFastq2)
        self.assertEqualContents(outHeader, inHeader)
        

class TestSplitReads(TestCaseWithTmp) :
    'Test various options of split_reads command.'
    def test_max_reads(self) :
        'Test splitting fastq using --maxReads option, with indexLen 1.'
        myInputDir = util.file.get_test_input_path(self)
        inFastq = os.path.join(myInputDir, 'in.fastq')
        outPrefix = util.file.mkstempfname()
        
        # Split
        parser = read_utils.parser_split_reads(argparse.ArgumentParser())
        args = parser.parse_args([inFastq, outPrefix, '--maxReads', '4',
                                  '--indexLen', '1'])
        args.func_main(args)
        
        # Check that results match expected
        expectedFastq1 = os.path.join(myInputDir, 'expected.fastq.1')
        expectedFastq2 = os.path.join(myInputDir, 'expected.fastq.2')
        self.assertEqualContents(outPrefix + '1', expectedFastq1)
        self.assertEqualContents(outPrefix + '2', expectedFastq2)

    def test_num_chunks(self) :
        'Test spliting fastq.gz using --numChunks option, with default indexLen.'
        myInputDir = util.file.get_test_input_path(self)
        inFastq = os.path.join(myInputDir, 'in.fastq.gz')
        outPrefix = util.file.mkstempfname()

        # Split
        parser = read_utils.parser_split_reads(argparse.ArgumentParser())
        args = parser.parse_args([inFastq, outPrefix, '--numChunks', '3'])
        args.func_main(args)
        
        # Check that results match expected
        expectedFastq1 = os.path.join(myInputDir, 'expected.fastq.01')
        expectedFastq2 = os.path.join(myInputDir, 'expected.fastq.02')
        expectedFastq3 = os.path.join(myInputDir, 'expected.fastq.03')
        self.assertEqualContents(outPrefix + '01', expectedFastq1)
        self.assertEqualContents(outPrefix + '02', expectedFastq2)
        self.assertEqualContents(outPrefix + '03', expectedFastq3)

    def test_fasta(self) :
        'Test splitting fasta file.'
        myInputDir = util.file.get_test_input_path(self)
        inFasta = os.path.join(myInputDir, 'in.fasta')
        outPrefix = util.file.mkstempfname()

        # Split
        parser = read_utils.parser_split_reads(argparse.ArgumentParser())
        args = parser.parse_args([inFasta, outPrefix, '--numChunks', '2',
                                  '--format', 'fasta'])
        args.func_main(args)
        
        # Check that results match expected
        expectedFasta1 = os.path.join(myInputDir, 'expected.fasta.01')
        expectedFasta2 = os.path.join(myInputDir, 'expected.fasta.02')
        self.assertEqualContents(outPrefix + '01', expectedFasta1)
        self.assertEqualContents(outPrefix + '02', expectedFasta2)
        

class TestMvicuna(TestCaseWithTmp) :
    """
    Input consists of 3 read pairs.
    Second read pair is identical to first.
    Third read pair has same 5' read as first, but different 3' read.
    What Mvicuna did was create paired output files in which the 2nd read
        was deleted. It created an empty unpaired output file. Although
        it initially created the postDupRm pair, it renamed them to the output
        pair.
    [IJ:]I have no idea if this is the correct behavior, but test checks that it
        doesn't change.
    """
    def test_mvicuna(self) :
        tempDir = tempfile.mkdtemp()
        myInputDir = util.file.get_test_input_path(self)
        
        # Run mvicuna
        inFastq1 = os.path.join(myInputDir, 'in.1.fastq')
        inFastq2 = os.path.join(myInputDir, 'in.2.fastq')
        pairedOutFastq1 = os.path.join(tempDir, 'pairedOut.1.fastq')
        pairedOutFastq2 = os.path.join(tempDir, 'pairedOut.2.fastq')
        unpairedOutFastq = os.path.join(tempDir, 'unpairedOut.fastq')
        args = read_utils.parser_dup_remove_mvicuna(argparse.ArgumentParser()).parse_args(
            [inFastq1, inFastq2,
             pairedOutFastq1, pairedOutFastq2,
             '--unpairedOutFastq', unpairedOutFastq])
        args.func_main(args)
        
        # Compare to expected
        for filename in ['pairedOut.1.fastq', 'pairedOut.2.fastq',
                         'unpairedOut.fastq'] :
            self.assertEqualContents(
                os.path.join(tempDir, filename),
                os.path.join(myInputDir, 'expected_' + filename))

if __name__ == '__main__':
    unittest.main()
