#!/usr/bin/python

__author__ = "tomkinsc@broadinstitute.org"

# built-ins
import unittest, os, argparse, pickle, shutil, tempfile
from collections import OrderedDict
import logging

# module-specific
from test import TestCaseWithTmp
import ncbi
import util.file

log = logging.getLogger(__name__)

class TestNcbiFetch(TestCaseWithTmp) :
    def setUp(self):
        super(TestNcbiFetch, self).setUp()

        # these are Orungo accessions
        self.accessions = [ "JQ610675.1",
                              "JQ610676.1",
                              "JQ610677.1",
                              "JQ610678.1",
                              "JQ610679.1",
                              "JQ610680.1",
                              "JQ610681.1",
                              "JQ610682.1",
                              "JQ610683.1",
                              "JQ610684.1"] 

        self.myInputDir = util.file.get_test_input_path(self)

    def perform_download_and_check(self, parser_func, additional_args, expected_files, null_files):
        
        tempDir = tempfile.gettempdir()

        args = [ "viral-ngs-test@example.com",
                tempDir]
        args.extend(self.accessions)
        args.extend(additional_args)           

        args = parser_func(argparse.ArgumentParser()).parse_args(args)
        args.func_main(args)

        # check that each file that each expected file was downloaded
        # and that the contents match what they should be
        for fileName in expected_files:
            createdFilePath = os.path.join(tempDir, fileName)
            log.info("createdFilePath: {}".format(createdFilePath))
            assert os.path.exists(createdFilePath), "File that should have been created does not exist: %s" % createdFilePath
            self.assertEqualContents(createdFilePath, os.path.join(self.myInputDir, fileName))

        for fileName in null_files:
            shouldNotExistFilePath = os.path.join(tempDir, fileName)
            assert not os.path.exists(shouldNotExistFilePath), "File exists but it should not: %s" % shouldNotExistFilePath

class TestFastaFetch(TestNcbiFetch):
    def setUp(self):
        super(TestFastaFetch, self).setUp()

    def test_download(self):
        args = []
        expectedFiles = [ a+".fasta" for a in self.accessions ]
        null_files = []      

        self.perform_download_and_check(ncbi.parser_fetch_fastas, 
            additional_args=args, expected_files=expectedFiles, null_files=null_files)

    def test_concat(self):
        args = ["--combinedFilePrefix", "orungo"]
        expectedFiles = [ "orungo.fasta" ]
        null_files = []
        
        self.perform_download_and_check(ncbi.parser_fetch_fastas, 
            additional_args=args, expected_files=expectedFiles, null_files=null_files)

    def test_removal_of_intermediates(self):
        args = ["--combinedFilePrefix", "orungo", "--removeSeparateFiles"]
        expectedFiles = [ "orungo.fasta" ]
        null_files = [ a+".fasta" for a in self.accessions ]
        
        self.perform_download_and_check(ncbi.parser_fetch_fastas, 
            additional_args=args, expected_files=expectedFiles, null_files=null_files)

    def test_individual_preexistance(self):
        # since the arguments are positional, including an accession here makes a duplicate that should
        # raise an Error
        args = [self.accessions[0]]
        args.extend(["--combinedFilePrefix", "orungo"])
        expectedFiles = [ "orungo.fasta" ]
        null_files = []
        
        with self.assertRaises(AssertionError):
            self.perform_download_and_check(ncbi.parser_fetch_fastas, 
                additional_args=args, expected_files=expectedFiles, null_files=null_files)

    def test_combined_preexistance(self):
        args = ["--combinedFilePrefix", "orungo"]
        expectedFiles = [ "orungo.fasta" ]
        null_files = []
        
        # call once to create the combined file
        self.perform_download_and_check(ncbi.parser_fetch_fastas, 
            additional_args=args, expected_files=expectedFiles, null_files=null_files)

        # an error should be raised the second time the call is made
        with self.assertRaises(AssertionError):
            self.perform_download_and_check(ncbi.parser_fetch_fastas, 
                additional_args=args, expected_files=expectedFiles, null_files=null_files)

    def test_overwrite(self):
        args = ["--combinedFilePrefix", "orungo", "--forceOverwrite"]
        expectedFiles = [ "orungo.fasta" ]
        null_files = []
        
        # call once to create the combined file
        self.perform_download_and_check(ncbi.parser_fetch_fastas, 
            additional_args=args, expected_files=expectedFiles, null_files=null_files)

        # no error should be raised the second time the call is made
        self.perform_download_and_check(ncbi.parser_fetch_fastas, 
            additional_args=args, expected_files=expectedFiles, null_files=null_files)

    def test_different_file_extension(self):
        args = ["--fileExt", "fa", "--combinedFilePrefix", "orungo"]
        expectedFiles = [ a+".fa" for a in self.accessions ]
        expectedFiles.append("orungo.fa")
        null_files = []      

        self.perform_download_and_check(ncbi.parser_fetch_fastas, 
            additional_args=args, expected_files=expectedFiles, null_files=null_files)


class TestFeatureTableFetch(TestNcbiFetch):
    def setUp(self):
        super(TestFeatureTableFetch, self).setUp()

    def test_download(self):
        args = []
        expectedFiles = [ a+".tbl" for a in self.accessions ]
        null_files = []      

        self.perform_download_and_check(ncbi.parser_fetch_feature_tables, 
            additional_args=args, expected_files=expectedFiles, null_files=null_files)

    def test_concat(self):
        args = ["--combinedFilePrefix", "orungo"]
        expectedFiles = [ "orungo.tbl" ]
        null_files = []
        
        self.perform_download_and_check(ncbi.parser_fetch_feature_tables, 
            additional_args=args, expected_files=expectedFiles, null_files=null_files)

    def test_removal_of_intermediates(self):
        args = ["--combinedFilePrefix", "orungo", "--removeSeparateFiles"]
        expectedFiles = [ "orungo.tbl" ]
        null_files = [ a+".tbl" for a in self.accessions ]
        
        self.perform_download_and_check(ncbi.parser_fetch_feature_tables, 
            additional_args=args, expected_files=expectedFiles, null_files=null_files)

    def test_individual_preexistance(self):
        # since the arguments are positional, including an accession here makes a duplicate that should
        # raise an Error
        args = [self.accessions[0]]
        args.extend(["--combinedFilePrefix", "orungo"])
        expectedFiles = [ "orungo.tbl" ]
        null_files = []
        
        with self.assertRaises(AssertionError):
            self.perform_download_and_check(ncbi.parser_fetch_feature_tables, 
                additional_args=args, expected_files=expectedFiles, null_files=null_files)

    def test_combined_preexistance(self):
        args = ["--combinedFilePrefix", "orungo"]
        expectedFiles = [ "orungo.tbl" ]
        null_files = []
        
        # call once to create the combined file
        self.perform_download_and_check(ncbi.parser_fetch_feature_tables, 
            additional_args=args, expected_files=expectedFiles, null_files=null_files)

        # an error should be raised the second time the call is made
        with self.assertRaises(AssertionError):
            self.perform_download_and_check(ncbi.parser_fetch_feature_tables, 
                additional_args=args, expected_files=expectedFiles, null_files=null_files)

    def test_overwrite(self):
        args = ["--combinedFilePrefix", "orungo", "--forceOverwrite"]
        expectedFiles = [ "orungo.tbl" ]
        null_files = []
        
        # call once to create the combined file
        self.perform_download_and_check(ncbi.parser_fetch_feature_tables, 
            additional_args=args, expected_files=expectedFiles, null_files=null_files)

        # no error should be raised the second time the call is made
        self.perform_download_and_check(ncbi.parser_fetch_feature_tables, 
            additional_args=args, expected_files=expectedFiles, null_files=null_files)

    def test_different_file_extension(self):
        args = ["--fileExt", "table", "--combinedFilePrefix", "orungo"]
        expectedFiles = [ a+".table" for a in self.accessions ]
        expectedFiles.append("orungo.table")
        null_files = []      

        self.perform_download_and_check(ncbi.parser_fetch_feature_tables, 
            additional_args=args, expected_files=expectedFiles, null_files=null_files)

class TestGenbankRecordFetch(TestNcbiFetch):
    def setUp(self):
        super(TestGenbankRecordFetch, self).setUp()

    def test_download(self):
        args = []
        expectedFiles = [ a+".gbk" for a in self.accessions ]
        null_files = []      

        self.perform_download_and_check(ncbi.parser_fetch_genbank_records, 
            additional_args=args, expected_files=expectedFiles, null_files=null_files)

    def test_concat(self):
        args = ["--combinedFilePrefix", "orungo"]
        expectedFiles = [ "orungo.gbk" ]
        null_files = []
        
        self.perform_download_and_check(ncbi.parser_fetch_genbank_records, 
            additional_args=args, expected_files=expectedFiles, null_files=null_files)

    def test_removal_of_intermediates(self):
        args = ["--combinedFilePrefix", "orungo", "--removeSeparateFiles"]
        expectedFiles = [ "orungo.gbk" ]
        null_files = [ a+".gbk" for a in self.accessions ]
        
        self.perform_download_and_check(ncbi.parser_fetch_genbank_records, 
            additional_args=args, expected_files=expectedFiles, null_files=null_files)

    def test_individual_preexistance(self):
        # since the arguments are positional, including an accession here makes a duplicate that should
        # raise an Error
        args = [self.accessions[0]]
        args.extend(["--combinedFilePrefix", "orungo"])
        expectedFiles = [ "orungo.gbk" ]
        null_files = []
        
        with self.assertRaises(AssertionError):
            self.perform_download_and_check(ncbi.parser_fetch_genbank_records, 
                additional_args=args, expected_files=expectedFiles, null_files=null_files)

    def test_combined_preexistance(self):
        args = ["--combinedFilePrefix", "orungo"]
        expectedFiles = [ "orungo.gbk" ]
        null_files = []
        
        # call once to create the combined file
        self.perform_download_and_check(ncbi.parser_fetch_genbank_records, 
            additional_args=args, expected_files=expectedFiles, null_files=null_files)

        # an error should be raised the second time the call is made
        with self.assertRaises(AssertionError):
            self.perform_download_and_check(ncbi.parser_fetch_genbank_records, 
                additional_args=args, expected_files=expectedFiles, null_files=null_files)

    def test_overwrite(self):
        args = ["--combinedFilePrefix", "orungo", "--forceOverwrite"]
        expectedFiles = [ "orungo.gbk" ]
        null_files = []
        
        # call once to create the combined file
        self.perform_download_and_check(ncbi.parser_fetch_genbank_records, 
            additional_args=args, expected_files=expectedFiles, null_files=null_files)

        # no error should be raised the second time the call is made
        self.perform_download_and_check(ncbi.parser_fetch_genbank_records, 
            additional_args=args, expected_files=expectedFiles, null_files=null_files)

    def test_different_file_extension(self):
        args = ["--fileExt", "gb", "--combinedFilePrefix", "orungo"]
        expectedFiles = [ a+".gb" for a in self.accessions ]
        expectedFiles.append("orungo.gb")
        null_files = []      

        self.perform_download_and_check(ncbi.parser_fetch_genbank_records, 
            additional_args=args, expected_files=expectedFiles, null_files=null_files)
