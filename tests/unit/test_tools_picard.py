# Unit tests for tools.PicardTools

__author__ = "dpark@broadinstitute.org"

import unittest
import re
import os
import tempfile
import shutil
import viral_ngs.core
import viral_ngs.core
import viral_ngs.core
import viral_ngs.core.picard
import viral_ngs.core.samtools
from test import TestCaseWithTmp


class TestToolPicard(TestCaseWithTmp):

    def test_fasta_index(self):
        orig_ref = os.path.join(viral_ngs.core.file.get_test_input_path(self), 'in.fasta')
        expected_dict = os.path.join(viral_ngs.core.file.get_test_input_path(self), 'in.dict')
        picard_index = viral_ngs.core.picard.CreateSequenceDictionaryTool()
        with open(expected_dict, 'rt') as inf:
            expected_first3 = [x.strip().split('\t')[:3] for x in inf.readlines()]
        for ext in ('.fasta', '.fa'):
            inRef = viral_ngs.core.file.mkstempfname(ext)
            shutil.copyfile(orig_ref, inRef)
            outDict = inRef[:-len(ext)] + '.dict'

            picard_index.execute(inRef, overwrite=True)

            # the dict files will not be exactly the same, just the first 3 cols
            with open(outDict, 'rt') as inf:
                # .replace("VN:1.5","VN:1.4") ==> because the header version may be 1.[4|5|6]
                actual_first3 = [re.sub(r'VN:1\.[4-6]','VN:1.4',x.strip()).split('\t')[:3] for x in inf.readlines()]
            self.assertEqual(actual_first3, expected_first3)

    def test_messy_fasta_index(self):
        orig_ref = os.path.join(viral_ngs.core.file.get_test_input_path(self), 'messy-headers.fasta')
        picard_index = viral_ngs.core.picard.CreateSequenceDictionaryTool()
        with viral_ngs.core.file.tempfname('.fasta') as inRef:
            shutil.copyfile(orig_ref, inRef)
            picard_index.execute(inRef, overwrite=True)
        with open(inRef[:-6] + '.dict', 'rt') as inf:
            seqnames = set()
            for line in inf:
                if line.startswith('@SQ'):
                    seq_name = dict(x.split(':', maxsplit=1) for x in line.strip().split('\t')[1:])['SN']
                    # old versions of code cut this off at "Influenza"
                    self.assertGreater(len(seq_name), 50)
                    seqnames.add(seq_name)
            # require that all sequence names are unique
            self.assertEqual(len(seqnames), 8)

    def test_sam_downsample(self):
        desired_count = 100
        tolerance = 0.02

        in_sam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'G5012.3.subset.bam')
        out_bam = viral_ngs.core.file.mkstempfname('.bam')

        downsamplesam = viral_ngs.core.picard.DownsampleSamTool()
        samtools = viral_ngs.core.samtools.SamtoolsTool()

        downsamplesam.downsample_to_approx_count(in_sam, out_bam, desired_count)

        assert samtools.count(out_bam) in range(
            int(desired_count - (desired_count * tolerance)), int(desired_count + (desired_count * tolerance))+1
        ), "Downsampled bam file does not contain the expected number of reads within tolerance: %s" % tolerance

    def test_revert_bam_empty_input(self):
        empty_bam = os.path.join(viral_ngs.core.file.get_test_input_path(), 'empty.bam')
        out_bam = viral_ngs.core.file.mkstempfname()
        viral_ngs.core.picard.RevertSamTool().execute(
            empty_bam,
            out_bam,
            picardOptions=['SORT_ORDER=queryname', 'SANITIZE=true']
        )
        samtools = viral_ngs.core.samtools.SamtoolsTool()
        assert samtools.count(out_bam) == 0
