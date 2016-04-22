"""Snakemake integration test runner."""
from __future__ import print_function
import copy
import os
from os.path import join
import shutil
import tempfile
import subprocess
import util.file
import util.misc
import yaml


def merge_yaml_dicts(*dicts):
    """Merge yaml dicts.

    Only handle top-level yaml dicts. Simply uses dict.update and ignores None
    dicts.
    """
    out = copy.deepcopy(dicts[0])
    for d in dicts[1:]:
        if d:
            out.update(d)
    return out


def touch(f):
    open(f, 'a').close()


class SnakemakeRunner(object):
    """Generates a working directory for snakemake integration tests.
    """

    # Keys in config for sample file lists
    SAMPLE_FILE_KEYS = {
        'samples_depletion',
        'samples_assembly',
        'samples_metagenomics',
        'samples_assembly_failures',
        'samples_per_run',
    }

    def __init__(self, workdir=None):
        self.workdir = workdir
        self.samples = set()
        self.config = None

    @property
    def data_dir(self):
        return join(self.workdir, self.config['data_dir'])

    def set_override_config(self, config):
        """Sets a dict of keys to override in base config."""
        self.override_config = config

    def setup(self):
        """Create working directory, subdirs, and config."""
        if not self.workdir:
            self.workdir = tempfile.mkdtemp('-snakemake')
        elif not os.path.isdir(self.workdir):
            os.makedirs(self.workdir)
        self.bindir = join(self.workdir, 'bin')
        self.root = util.file.get_project_path()
        os.symlink(self.root, self.bindir)
        os.symlink(join(self.root, 'pipes', 'Snakefile'),
                   join(self.workdir, 'Snakefile'))
        with open(join(self.root, 'pipes', 'config.yaml')) as f:
            config = yaml.load(f)
        if self.override_config:
            self.config = merge_yaml_dicts(config, self.override_config)
        else:
            self.config = config
        with open(join(self.workdir, 'config.yaml'), 'w') as f:
            yaml.dump(self.config, f)
        self.create_subdirs()

    def create_subdirs(self):
        """Create the data subdirs."""
        for _, subdir in self.config['subdirs'].items():
            os.makedirs(join(self.data_dir, subdir))


    def link_samples(self, samples, destination='source'):
        """Links samples files in data destination dir."""
        for sample in samples:
            link = join(self.data_dir, self.config['subdirs'][destination],
                        os.path.basename(sample))
            os.symlink(sample, link)
            self.samples.add(sample)


    def create_sample_files(self, samples=None, sample_files=None):
        """Creates files for sample lists.

        If samples is None, add all samples to all selected sample_files. If
        sample_files is None, write to all sample files. Afterwards, touch all
        sample files that weren't written to.
        """
        samples = samples or self.samples
        all_sample_files = [self.config[key] for key in
                            SnakemakeRunner.SAMPLE_FILE_KEYS]

        if not sample_files:
            sample_files = all_sample_files
        for sample_key in sample_files:
            sample_file = self.config[sample_key]
            with open(join(self.workdir, sample_file), 'w') as f:
                for sample in samples:
                    print(os.path.splitext(os.path.basename(sample))[0], file=f)

        for sample_file in set(all_sample_files) - set(sample_files):
            touch(join(self.workdir, sample_file))

    def run(self, rules=None):
        """Run snakemake with extra verbosity. """
        cmd = ['snakemake', '--verbose', '--reason', '--printshellcmds']
        if rules:
            cmd.extend(rules)
        try:
            util.misc.run_and_print(cmd, check=True, cwd=self.workdir)
        except subprocess.CalledProcessError as e:
            print(e.output.decode('utf-8'))
            raise
