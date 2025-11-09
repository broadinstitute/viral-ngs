# Unit tests for tools.fastqc

__author__ = "dpark@broadinstitute.org"

import unittest
import os
import zipfile
from html.parser import HTMLParser
import util.file
import tools.fastqc
from test import TestCaseWithTmp


class HTMLValidator(HTMLParser):
    """Simple HTML validator to check if HTML can be parsed without errors"""
    def __init__(self):
        super().__init__()
        self.valid = True
        self.error_msg = None

    def error(self, message):
        self.valid = False
        self.error_msg = message


class TestToolFastQC(TestCaseWithTmp):

    def _validate_html(self, html_path):
        """Validate that a file contains parseable HTML"""
        with open(html_path, 'rt') as f:
            html_content = f.read()

        parser = HTMLValidator()
        try:
            parser.feed(html_content)
            self.assertTrue(parser.valid, f"HTML parsing failed: {parser.error_msg}")
        except Exception as e:
            self.fail(f"HTML parsing raised exception: {e}")

    def _validate_zip(self, zip_path):
        """Validate that a file is a valid zip archive"""
        try:
            with zipfile.ZipFile(zip_path, 'r') as zf:
                # testzip() returns None if the archive is valid
                result = zf.testzip()
                self.assertIsNone(result, f"ZIP validation failed: {result}")
        except zipfile.BadZipFile as e:
            self.fail(f"Invalid ZIP file: {e}")

    def test_fastqc_nonempty_bam(self):
        """Test FastQC on a non-empty BAM file"""
        in_bam = os.path.join(util.file.get_test_input_path(), 'G5012.3.subset.bam')
        out_html = util.file.mkstempfname('.html')
        out_zip = util.file.mkstempfname('.zip')

        fastqc = tools.fastqc.FastQC()
        fastqc.execute(in_bam, out_html, out_zip=out_zip)

        # Verify HTML file is created and non-empty
        self.assertTrue(os.path.exists(out_html), "HTML output file should exist")
        self.assertGreater(os.path.getsize(out_html), 0, "HTML file should not be empty")

        # Validate HTML is parseable
        self._validate_html(out_html)

        # Verify ZIP file is created and is valid
        self.assertTrue(os.path.exists(out_zip), "ZIP output file should exist")
        self.assertGreater(os.path.getsize(out_zip), 0, "ZIP file should not be empty")

        # Validate ZIP file
        self._validate_zip(out_zip)

        # Verify ZIP contains FastQC output files
        with zipfile.ZipFile(out_zip, 'r') as zf:
            namelist = zf.namelist()
            self.assertGreater(len(namelist), 0, "ZIP should contain FastQC output files")

    def test_fastqc_empty_bam(self):
        """Test FastQC on an empty BAM file (zero reads)"""
        in_bam = os.path.join(util.file.get_test_input_path(), 'empty.bam')
        out_html = util.file.mkstempfname('.html')
        out_zip = util.file.mkstempfname('.zip')

        fastqc = tools.fastqc.FastQC()
        fastqc.execute(in_bam, out_html, out_zip=out_zip)

        # Verify HTML file is created and non-empty
        self.assertTrue(os.path.exists(out_html), "HTML output file should exist")
        self.assertGreater(os.path.getsize(out_html), 0, "HTML file should not be empty")

        # Validate HTML is parseable
        self._validate_html(out_html)

        # Verify HTML contains expected message about empty input
        with open(out_html, 'rt') as f:
            html_content = f.read()
            self.assertIn('zero reads', html_content.lower(), "HTML should mention zero reads")

        # Verify ZIP file is created and is a valid (empty) zip
        self.assertTrue(os.path.exists(out_zip), "ZIP output file should exist")

        # A valid empty zip file should be 22 bytes (zip structure overhead)
        self.assertEqual(os.path.getsize(out_zip), 22, "Empty ZIP should be 22 bytes (valid zip structure)")

        # Validate ZIP file
        self._validate_zip(out_zip)

        # Verify ZIP has no files (empty)
        with zipfile.ZipFile(out_zip, 'r') as zf:
            namelist = zf.namelist()
            self.assertEqual(len(namelist), 0, "Empty BAM should produce empty ZIP (0 files)")

    def test_fastqc_without_zip(self):
        """Test FastQC when --out_zip is not specified"""
        in_bam = os.path.join(util.file.get_test_input_path(), 'G5012.3.subset.bam')
        out_html = util.file.mkstempfname('.html')

        fastqc = tools.fastqc.FastQC()
        fastqc.execute(in_bam, out_html, out_zip=None)

        # Verify HTML file is created
        self.assertTrue(os.path.exists(out_html), "HTML output file should exist")
        self.assertGreater(os.path.getsize(out_html), 0, "HTML file should not be empty")

        # Validate HTML is parseable
        self._validate_html(out_html)
