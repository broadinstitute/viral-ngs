# Unit tests for snakemake pipelines

__author__ = "dpark@broadinstitute.org"

import util.cmd
import util.file
import unittest
import argparse
import sys
import os
import subprocess
import shutil
import tempfile
import argparse
import itertools
from test import TestCaseWithTmp

if sys.version_info >= (3, 2):
    import snakemake


def setup_dummy_simple(sample_names=('G1234', 'G5678', 'G3671.1_r1', 'G3680-1_4', '9876', 'x.y-7b')):
    ''' Set up a very simple project directory with empty input files. '''

    workdir = tempfile.mkdtemp()
    os.mkdir(os.path.join(workdir, 'data'))
    os.mkdir(os.path.join(workdir, 'ref_genome_dir'))
    os.mkdir(os.path.join(workdir, 'data', '00_raw'))
    os.mkdir(os.path.join(workdir, 'data', '01_per_sample'))
    os.mkdir(os.path.join(workdir, 'log'))
    os.mkdir(os.path.join(workdir, 'reports'))
    os.mkdir(os.path.join(workdir, 'tmp'))

    for s in sample_names:
        with open(os.path.join(workdir, 'data', '01_per_sample', s + '.raw.bam'), 'wt') as outf:
            pass
        with open(os.path.join(workdir, 'data', '00_raw', s + '.bam'), 'wt') as outf:
            pass
    for fn in ('samples-assembly.txt', 'samples-depletion.txt', 'samples-runs.txt', 'samples-assembly-failures.txt',
               'samples-metagenomics.txt'):
        with open(os.path.join(workdir, fn), 'wt') as outf:
            for s in sample_names:
                outf.write(s + '\n')

    shutil.copy(os.path.join(util.file.get_project_path(), 'pipes', 'Snakefile'), workdir)
    shutil.copy(os.path.join(util.file.get_project_path(), 'pipes', 'config.yaml'), workdir)

    os.symlink(util.file.get_project_path(), os.path.join(workdir, 'bin'))

    return workdir


@unittest.skipIf(sys.version_info < (3, 2), "python version is too old for snakemake")
class TestSimpleDryRuns(TestCaseWithTmp):

    def setUp(self):
        super(TestSimpleDryRuns, self).setUp()
        self.workdir = setup_dummy_simple()
        self.env = {'GATK_PATH': os.environ.get('GATK_PATH'), 'NOVOALIGN_PATH': os.environ.get('NOVOALIGN_PATH')}

    def tearDown(self):
        for k, v in self.env.items():
            if v:
                os.environ[k] = v
        super(TestSimpleDryRuns, self).tearDown()

    def test_dryrun_all(self):
        ''' Test that the "all" rule dryruns properly '''
        self.assertTrue(snakemake.snakemake(
            os.path.join(self.workdir, 'Snakefile'),
            #configfile=os.path.join(self.workdir, 'config.yaml'),
            workdir=self.workdir,
            dryrun=True))
        self.assertTrue(snakemake.snakemake(
            os.path.join(self.workdir, 'Snakefile'),
            #configfile=os.path.join(self.workdir, 'config.yaml'),
            workdir=self.workdir,
            dryrun=True,
            targets=['all']))

    def test_dryrun_all_assemble(self):
        ''' Test that the "all_assemble" rule dryruns properly '''
        self.assertTrue(snakemake.snakemake(
            os.path.join(self.workdir, 'Snakefile'),
            #configfile=os.path.join(self.workdir, 'config.yaml'),
            workdir=self.workdir,
            dryrun=True,
            targets=['all_assemble']))

    def test_dryrun_all_deplete(self):
        ''' Test that the "all_deplete" rule dryruns properly '''
        self.assertTrue(snakemake.snakemake(
            os.path.join(self.workdir, 'Snakefile'),
            #configfile=os.path.join(self.workdir, 'config.yaml'),
            workdir=self.workdir,
            dryrun=True,
            targets=['all_deplete']))

    def test_dryrun_all_metagenomics(self):
        ''' Test that the "all_metagenomics" rule dryruns properly '''
        self.assertTrue(snakemake.snakemake(
            os.path.join(self.workdir, 'Snakefile'),
            #configfile=os.path.join(self.workdir, 'config.yaml'),
            workdir=self.workdir,
            dryrun=True,
            targets=['all_metagenomics']))
