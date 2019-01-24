# Unit tests for snakemake pipelines

__author__ = "dpark@broadinstitute.org"

import sys
import os
import subprocess
import shutil
import itertools

import pytest
import yaml

import util.cmd
import util.file


needs_snakemake = pytest.mark.skipif(
    sys.version_info < (3, 5),
    reason='python version is too old for snakemake')


if sys.version_info >= (3, 5):
    import snakemake


def add_to_sample_list(workdir, sample_list_name, sample_names):
    with open(os.path.join(workdir, 'samples-{}.txt'.format(sample_list_name)), 'a') as outf:
        for sample_name in sample_names:
            outf.write(sample_name + '\n')


def setup_dummy_simple(workdir, sample_names=('G1234', 'G5678', 'G3671.1_r1', 'G3680-1_4', '9876', 'x.y-7b')):
    ''' Set up a very simple project directory with empty input files. '''

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
    for name in ('assembly', 'depletion', 'runs', 'assembly-failures', 'metagenomics'):
        add_to_sample_list(workdir, name, sample_names)

    shutil.copy(os.path.join(util.file.get_project_path(), 'pipes', 'Snakefile'), workdir)

    with open(os.path.join(util.file.get_project_path(), 'pipes', 'config.yaml')) as f:
        config = yaml.load(f)

    def translate_remote_s3(uri):
      remote_path = uri[5:]
      fake_s3_root = os.path.join(util.file.get_project_path(), 'test', 'input', 's3')
      local_path = os.path.join(fake_s3_root, remote_path)
      return local_path

    for k, v in config.items():
        if isinstance(v, str):
            if v.startswith('s3://'):
                config[k] = translate_remote_s3(v)

        if util.misc.is_nonstr_iterable(v):
            for i, vv in enumerate(v):
                if isinstance(vv, str):
                  if vv.startswith('s3://'):
                      v[i] = translate_remote_s3(vv)
    with open(os.path.join(workdir, 'config.yaml'), 'w') as f:
        yaml.dump(config, f)

    os.symlink(util.file.get_project_path(), os.path.join(workdir, 'bin'))
    return workdir


@pytest.fixture
def workdir(request, tmpdir_function):
    env = {'GATK_PATH': os.environ.get('GATK_PATH'), 'NOVOALIGN_PATH': os.environ.get('NOVOALIGN_PATH')}
    setup_dummy_simple(tmpdir_function)
    yield tmpdir_function
    for k, v in env.items():
        if v:
            os.environ[k] = v


def call_snakemake(workdir, targets=None):
    return snakemake.snakemake(
        os.path.join(workdir, 'Snakefile'),
        configfile=os.path.join(workdir, 'config.yaml'),
        workdir=workdir,
        dryrun=True,
        targets=targets)


@needs_snakemake
def test_dryrun_all(workdir):
    ''' Test that the "all" rule dryruns properly '''
    assert call_snakemake(workdir)
    assert call_snakemake(workdir, ['all'])


@needs_snakemake
def test_dryrun_all_assemble(workdir):
    ''' Test that the "all_assemble" rule dryruns properly '''
    assert call_snakemake(workdir, ['all_assemble'])

@needs_snakemake
def test_dryrun_all_deplete(workdir):
    ''' Test that the "all_deplete" rule dryruns properly '''
    assert call_snakemake(workdir, ['all_deplete'])

@needs_snakemake
def test_dryrun_all_metagenomics(workdir):
    ''' Test that the "all_metagenomics" rule dryruns properly '''
    assert call_snakemake(workdir, ['all_metagenomics'])

@needs_snakemake
def test_missing_merge_inputs(workdir):
    add_to_sample_list(workdir, 'assembly', 'G_missing')
    assert call_snakemake(workdir, ['all_assemble']) == False
