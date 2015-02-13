# Unit tests for Novoalign aligner

__author__ = "mlin@dnanexus.com"

import unittest, os.path, shutil
import util.file, tools.mosaik
from test import TestCaseWithTmp

class TestToolMosaik(TestCaseWithTmp) :

    def setUp(self):
        super(TestToolMosaik, self).setUp()
        self.mosaik = tools.mosaik.MosaikTool()
        self.mosaik.install()

    def test_get_networkFile(self):
        nndir = self.mosaik.get_networkFile()
        assert os.path.exists(os.path.join(nndir, "2.1.26.pe.100.0065.ann"))

    # TO DO: further testing of Mosaik invocations
    # system($mosaikpath."MosaikBuild -q ".$option{fq}." -q2 ".$option{fq2}." -out $readdat -st ".$option{st}." -mfl ".$option{mfl});
    # system($mosaikpath."MosaikAligner -in $readdat -out $output -ia $refdat -hs ".$option{hs}." -act ".$option{act}." -mm 500 -mmp ".$option{mmp}." -minp ".$option{minp}." -ms ".$option{ms}." -mms ".$option{mms}." -gop ".$option{gop}." -hgop ".$option{hgop}." -gep ".$option{gep}."$bw -m ".$option{m}." -annpe ".$option{annpe}." -annse ".$option{annse});
