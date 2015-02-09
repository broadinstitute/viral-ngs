# Unit tests for tools.PicardTools

__author__ = "dpark@broadinstitute.org"

import unittest, os, tempfile, shutil
import util, util.file, tools, tools.picard
from test import TestCaseWithTmp

class TestToolPicard(TestCaseWithTmp):
    def test_fasta_index(self):
        orig_ref = os.path.join(util.file.get_test_input_path(self), 'in.fasta')
        expected_dict = os.path.join(util.file.get_test_input_path(self), 'in.dict')
        picard_index = tools.picard.CreateSequenceDictionaryTool()
        with open(expected_dict, 'rt') as inf:
            expected_first3 = [x.strip().split('\t')[:3] for x in inf.readlines()]
        for ext in ('.fasta', '.fa'):
            inRef = util.file.mkstempfname(ext)
            shutil.copyfile(orig_ref, inRef)
            outDict = inRef[:-len(ext)] + '.dict'
            
            picard_index.execute(inRef)
            
            # the dict files will not be exactly the same, just the first 3 cols
            with open(outDict, 'rt') as inf:
                actual_first3 = [x.strip().split('\t')[:3] for x in inf.readlines()]
            self.assertEqual(actual_first3, expected_first3)


if __name__ == '__main__':
    unittest.main()
